"""
Step 7: Supplementary Figure S2, change in scaled CFI at each invariance step, per questionnaire and comparison
(script 03 outputs).
Usage: python scripts/07_figure_invariance.py
"""

import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1].as_posix()
CAT = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a"}          # same palette as the other figures
TXT1, TXT2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": TXT2,
                     "xtick.color": TXT2, "ytick.color": TXT2, "axes.titlecolor": TXT1, "figure.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False})
NAME = {"DES": "DES (dissociation)", "MAIA": "MAIA-2 (interoception)", "IAMI": "IAMI (involuntary memory)",
        "CTQ": "CTQ (childhood trauma)", "PHQ": "PHQ-9 (depression)", "GAD": "GAD-7 (anxiety)", "PCL": "PCL-5 (PTSD)"}
CAT["violet"] = "#4a3aa7"
STEP = {"thresholds": ("Equal thresholds", CAT["violet"], "^"), "metric": ("+ equal loadings", CAT["blue"], "o"),
        "scalar": ("+ equal intercepts (scalar)", CAT["orange"], "s"),
        "partial_scalar": ("Partial scalar (items freed)", CAT["aqua"], "D")}
CON = {"control": "Congenital aphantasia vs Prolific controls",
       "an_with_imagery": "Congenital aphantasia vs Aphantasia Network members with imagery"}

REFS = ["control", "an_with_imagery"]
F = pd.concat(pd.read_csv(f"{ROOT}/output/03_invariance_fit_{r}.csv") for r in REFS)
fr = pd.concat(pd.read_csv(f"{ROOT}/output/03_invariance_freed_{r}.csv") for r in REFS).set_index(["contrast", "instrument"])
order = [k for k in NAME if k in set(F.instrument)][::-1]
lo = min(-.03, F.d_cfi.min() - .005)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True, gridspec_kw={"wspace": .15})
for ax, (c, title) in zip(axes, CON.items()):
    ax.axvspan(lo, -.01, color="#fbeceb", lw=0)
    ax.axvline(-.01, color=TXT2, lw=1, ls=(0, (3, 2)))
    for step, (lab, col, mk) in STEP.items():
        s = F[(F.contrast == c) & (F.model == step)].set_index("instrument").reindex(order)
        ax.scatter(s.d_cfi, range(len(order)), s=60, color=col, marker=mk, edgecolor="white", linewidth=1.5, zorder=3, label=lab)
    for j, k in enumerate(order):
        n = fr.n_freed.get((c, k), 0)
        if n:
            cfi_ok = F[(F.contrast == c) & (F.instrument == k) & (F.model == "scalar")].d_cfi.iloc[0] >= -.01
            why = ",\nfails on RMSEA only" if cfi_ok else ",\nstill not invariant"
            ax.text(.021, j, f"{n} of {fr.n_items[(c, k)]} items freed" + (why if fr.still_fails[(c, k)] else ""),
                    va="center", fontsize=7.5, color=TXT2)
    ax.set_xlim(lo, .05); ax.set_yticks(range(len(order)), [NAME[k] for k in order]); ax.tick_params(length=0)
    ax.set_xlabel("Change in scaled CFI at this step"); ax.grid(axis="x", color=GRID, lw=.6); ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=9.5, fontweight="bold")
axes[0].text(lo + .004, (len(order) - 1) / 2, "not invariant", fontsize=7.5, color=TXT2, va="center", ha="center", rotation=90)
axes[1].legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0, -0.16), ncol=4)
fig.savefig(f"{ROOT}/figures/Supplementary_Figure_S2.png", dpi=200, bbox_inches="tight"); plt.close(fig)
print("wrote figures/Supplementary_Figure_S2.png")
