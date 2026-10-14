"""
Step 6: Supplementary Figure S1, factor structure of the questionnaire items (script 02 outputs).
Usage: python scripts/06_figure_factor_structure.py
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from pathlib import Path
ROOT = Path(__file__).resolve().parents[1].as_posix()
# --- palette ---
CAT = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100", "magenta": "#e87ba4",
       "green": "#008300", "violet": "#4a3aa7", "red": "#e34948"}
INST = {"MAIA": CAT["blue"], "DES": CAT["orange"], "CTQ": CAT["aqua"], "PHQ": CAT["yellow"], "IAMI": CAT["magenta"],
        "GAD": CAT["green"], "PCL": CAT["red"], "PSIQ": CAT["violet"]}
INST_NAME = {"PSIQ": "Psi-Q (imagery)", "DES": "DES (dissociation)", "MAIA": "MAIA-2 (interoception)",
             "IAMI": "IAMI (involuntary memory)", "CTQ": "CTQ (childhood trauma)", "PHQ": "PHQ-9 (depression)",
             "GAD": "GAD-7 (anxiety)", "PCL": "PCL-5 (PTSD)"}
DIV = LinearSegmentedColormap.from_list("div", ["#e34948", "#f0efec", "#2a78d6"])
TXT1, TXT2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": TXT2,
                     "xtick.color": TXT2, "ytick.color": TXT2, "axes.titlecolor": TXT1, "figure.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False})
# plain names for the pooled factors (named from the content of the items loading on each)
PLAIN = {"F1": "Imagery", "F2": "Dissociation", "F3": "Depression and anxiety", "F4": "Involuntary future thoughts",
         "F5": "Body awareness", "F6": "Emotional abuse and neglect", "F7": "PTSD re-experiencing",
         "F8": "Not distracting from the body", "F9": "Sexual abuse", "F10": "Physical abuse",
         "F11": "Involuntary past memories", "F12": "Not worrying about the body", "F13": "Self-regulation"}
plain = lambda s: PLAIN[s.split(":")[0]]
# guard: each name assumes the questionnaire that supplies most of the factor's salient items (script 02 labels such as
# "F3:PHQ9+GAD7+PCL7"); a different factor order, e.g. after a package update, stops the script instead of mislabelling
MAIN = {"F1": "PSIQ", "F2": "DES", "F3": "PHQ", "F4": "IAMI", "F5": "MAIA", "F6": "CTQ", "F7": "PCL", "F8": "MAIA",
        "F9": "CTQ", "F10": "CTQ", "F11": "IAMI", "F12": "MAIA", "F13": "MAIA"}

# ---------------- Supplementary Figure S1: factor structure ----------------
ld = pd.read_csv(f"{ROOT}/output/02_loadings_pooled.csv", index_col=0)
bad_names = [c for c in ld.columns.drop("instrument") if not c.split(":")[1].startswith(MAIN[c.split(":")[0]])]
assert not bad_names, f"factor order changed, check PLAIN and MAIN: {bad_names}"
order = ["PSIQ", "DES", "MAIA", "IAMI", "CTQ", "PHQ", "GAD", "PCL"]
ld["o"] = ld.instrument.map(order.index)
ld = ld.sort_values("o", kind="stable")
M = ld.drop(columns=["instrument", "o"])
col_o = sorted(M.columns, key=lambda c: np.average(np.arange(len(M)), weights=M[c].abs() ** 4))
cg = pd.read_csv(f"{ROOT}/output/02_congruence.csv", index_col=0).tucker_phi
null = pd.read_csv(f"{ROOT}/output/02_congruence_null.csv").mean_phi

fig, (a0, a1) = plt.subplots(1, 2, figsize=(11, 6.4), gridspec_kw={"width_ratios": [1.35, 1], "wspace": .8})
im = a0.imshow(M[col_o].values, aspect="auto", cmap=DIV, vmin=-1, vmax=1, interpolation="nearest")
bounds = ld.groupby("instrument", sort=False).size().cumsum()
start = 0
for k, end in bounds.items():
    a0.add_patch(plt.Rectangle((-1.4, start - .5), .7, end - start, color=INST[k], clip_on=False, lw=0))
    a0.text(-1.8, (start + end) / 2 - .5, INST_NAME[k], ha="right", va="center", fontsize=8, color=TXT2)
    a0.axhline(end - .5, color="white", lw=1.5)
    start = end
a0.set_xticks(range(len(col_o)), [plain(c) for c in col_o], rotation=55, ha="right", fontsize=8)
a0.set_yticks([]); a0.spines[:].set_visible(False)
a0.set_title("a", loc="left", fontsize=10, fontweight="bold")
cb = fig.colorbar(im, ax=a0, fraction=.035, pad=.02); cb.set_label("Loading", fontsize=8); cb.outline.set_visible(False)

cg = cg.sort_values()
a1.axvspan(np.percentile(null, 5), np.percentile(null, 95), color=GRID, lw=0)
a1.axvline(null.mean(), color=TXT2, lw=1, ls=(0, (3, 2)))
a1.scatter(cg.values, range(len(cg)), s=55, color=CAT["blue"], edgecolor="white", linewidth=1.5, zorder=3)
a1.set_yticks(range(len(cg)), [f"Factor {i.split(':')[0][1:]} ({i.split(':')[1]})" for i in cg.index], fontsize=8)
a1.set_xlim(.55, 1); a1.set_xlabel("Tucker congruence,\ncontrols vs congenital aphantasia")
a1.text(null.mean() - .012, .5, "two random halves\nof the controls", ha="right", va="center", fontsize=7.5, color=TXT2)
a1.grid(axis="x", color=GRID, lw=.6); a1.set_axisbelow(True); a1.tick_params(length=0)
a1.set_title("b", loc="left", fontsize=10, fontweight="bold")
fig.savefig(f"{ROOT}/figures/Supplementary_Figure_S1.png", dpi=200, bbox_inches="tight"); plt.close(fig)

print("wrote figures/Supplementary_Figure_S1.png")
