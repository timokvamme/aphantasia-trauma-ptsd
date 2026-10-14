"""
Step 9: compare every result table of this run (output/*.csv) with the published results (expected/*.csv).

A number matches when it agrees to 6 significant digits; text cells must be identical. The participant-level file
output/scored.csv is not compared. Prints "REPRODUCTION OK" and exits 0 when every file matches.
Usage: python scripts/09_check_reproduction.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
expected = sorted((ROOT / "expected").glob("*.csv"))
problems = []
for e in expected:
    o = ROOT / "output" / e.name
    if not o.exists():
        problems.append(f"{e.name}: not produced by this run")
        continue
    a, b = pd.read_csv(e), pd.read_csv(o)
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        problems.append(f"{e.name}: shape or columns differ ({a.shape} vs {b.shape})")
        continue
    for c in a.columns:
        x, y = a[c], b[c]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            ok = np.allclose(x.astype(float), y.astype(float), rtol=1e-6, atol=1e-9, equal_nan=True)
        else:
            ok = (x.astype(object).where(x.notna(), "") == y.astype(object).where(y.notna(), "")).all()
        if not ok:
            problems.append(f"{e.name}: column '{c}' differs")
if not expected:
    problems.append("expected/ holds no result files to compare with")
print(f"{len(expected)} result files compared with expected/")
print("\n".join(problems) if problems else "REPRODUCTION OK")
sys.exit(1 if problems else 0)
