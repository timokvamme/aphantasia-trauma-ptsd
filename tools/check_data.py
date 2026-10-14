"""
Check that data/trauma_survey_anonymised.csv contains only the allowlisted columns and no identifying content:
  - every column is in the allowlist of tools/anonymise.py and described in data/codebook.csv
  - pid matches pr0000 / an0000 and src is prolific / an; every other column is numeric; nobody is under 18
  - no cell looks like an IP address, e-mail address, date or time, or a long ID string
Prints "DATA OK" and exits 0, or lists the problems and exits 1.
Usage: python tools/check_data.py
"""
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from anonymise import ITEM, KEEP  # noqa: E402

df = pd.read_csv(ROOT / "data" / "trauma_survey_anonymised.csv", dtype=str)
book = pd.read_csv(ROOT / "data" / "codebook.csv")
problems = []
allowed = {"pid", "src", *KEEP.values()}
extra = [c for c in df.columns if c not in allowed and not ITEM.match(c)]
if extra:
    problems.append(f"columns outside the allowlist: {extra}")
if set(df.columns) != set(book.column):
    problems.append("codebook and data columns differ")
if not df.pid.str.fullmatch(r"(pr|an)\d{4}").all():
    problems.append("pid values outside the pr0000/an0000 pattern")
if (pd.to_numeric(df.age, errors="coerce") < 18).any():
    problems.append("participants under 18 in the file")
if not df.src.isin(["prolific", "an"]).all():
    problems.append("src values outside prolific/an")
for c in df.columns.drop(["pid", "src"]):
    bad = pd.to_numeric(df[c], errors="coerce").isna() & df[c].notna()
    if bad.any():
        problems.append(f"{c}: {int(bad.sum())} non-numeric cells")
text = df.to_csv(index=False)
for name, pat in [("IP address", r"\b(?:\d{1,3}\.){3}\d{1,3}\b"), ("e-mail", r"[\w.+-]+@[\w-]+\.[\w.]+"),
                  ("date", r"\b\d{1,4}[/-]\d{1,2}[/-]\d{1,4}\b"), ("time", r"\b\d{1,2}:\d{2}\b"),
                  ("long ID", r"\b[0-9a-fA-F]{16,}\b")]:
    n = len(re.findall(pat, text))
    if n:
        problems.append(f"{n} cells look like a {name}")
print(f"{len(df)} rows, {df.shape[1]} columns checked")
print("\n".join(problems) if problems else "DATA OK")
sys.exit(1 if problems else 0)
