"""
Step 8: Figures 1 to 4 (script 04 outputs).
  Figure 1  exposure, mood, post-traumatic stress, dissociation and involuntary memory in congenital aphantasia
  Figure 2  each PCL-5 item at equal exposure, against both reference groups
  Figure 3  depression and post-traumatic stress against childhood trauma, by group
  Figure 4  share above the PCL-5 cut-off, and odds ratios with increasing adjustment
Usage: python scripts/08_figures_trauma.py   (run script 04 first)
"""
import warnings

import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

from pathlib import Path
warnings.simplefilter("ignore", pd.errors.PerformanceWarning)   # many derived columns are added one at a time
ROOT = Path(__file__).resolve().parents[1].as_posix()
CAT = {"blue": "#2a78d6", "orange": "#eb6834", "aqua": "#1baf7a", "yellow": "#eda100", "violet": "#4a3aa7"}
TXT1, TXT2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": TXT2,
                     "xtick.color": TXT2, "ytick.color": TXT2, "axes.titlecolor": TXT1, "figure.facecolor": "white",
                     "axes.spines.top": False, "axes.spines.right": False})
REF = {"control": ("vs Prolific controls", CAT["blue"], "o"), "an_imagers": ("vs Aphantasia Network members with imagery", CAT["aqua"], "s")}
GRP = {"control": ("Prolific controls", CAT["blue"]), "an_imagers": ("Aphantasia Network members with imagery", CAT["aqua"]),
       "congenital": ("Congenital aphantasia", CAT["orange"])}
A = pd.read_csv(f"{ROOT}/output/04_all.csv", low_memory=False)
B = pd.read_csv(f"{ROOT}/output/04_symptom_form.csv")
B = B[B.group == "congenital"]

# ---------------- t1: profile ----------------
SECTIONS = [("Trauma exposure", [("Childhood trauma (CTQ)", "A", "Childhood trauma, total"),
                                 ("Traumatic events (LEC-5)", "A", "Traumatic events experienced or witnessed")]),
            ("Mood", [("Depression (PHQ-9)", "B", "Depression"), ("Anxiety (GAD-7)", "B", "Anxiety")]),
            ("Post-traumatic stress (PCL-5)", [("Intrusion", "B", "Intrusion (B)"), ("Avoidance", "B", "Avoidance (C)"),
                                               ("Negative cognition and mood", "B", "Negative cognition and mood (D)"),
                                               ("Arousal", "B", "Arousal (E)")]),
            ("Dissociation (DES)", [("Depersonalisation and derealisation", "B", "Depersonalisation/derealisation"),
                                    ("Amnesia", "B", "Amnesia"), ("Absorption, unbiased items", "B", "Absorption, unbiased items"),
                                    ("No memory of important events (item 9)", "B", "DES 9: no memory of important events")]),
            ("Involuntary memory (IAMI)", [("Past memories", "B", "Involuntary past memories (IAMI)"),
                                           ("Future thoughts", "B", "Involuntary future thoughts (IAMI)")])]
rows, y, ticks, heads = [], 0, [], []
for sec, items in SECTIONS:
    heads.append((y, sec)); y += 1
    for lab, part, key in items:
        for r_ in REF:
            if part == "A":
                s = A[(A.part == "A_exposure") & (A.measure == key) & (A.vs == r_)].iloc[0]
                rows.append((y, r_, s.d, s.lo, s.hi))
            else:
                s = B[(B.measure == key) & (B.vs == r_)].iloc[0]
                rows.append((y, r_, s.est, s.lo, s.hi))
        ticks.append((y, lab)); y += 1
    y += .4
fig, ax = plt.subplots(figsize=(7.6, 8.2))
ax.axvline(0, color=TXT2, lw=1)
for r_, (lab, col, mk) in REF.items():
    R = [r for r in rows if r[1] == r_]
    off = -.14 if r_ == "control" else .14
    yy = np.array([r[0] for r in R]) + off
    ax.hlines(yy, [r[3] for r in R], [r[4] for r in R], color=col, lw=2, alpha=.55)
    ax.scatter([r[2] for r in R], yy, s=46, color=col, marker=mk, edgecolor="white", linewidth=1.2, zorder=3, label=lab)
ax.set_yticks([t[0] for t in ticks], [t[1] for t in ticks]); ax.tick_params(length=0)
for yh, sec in heads:
    ax.text(-1.6, yh, sec, fontweight="bold", color=TXT1, va="center", ha="left", fontsize=9)
ax.set_ylim(y - .2, -.8); ax.set_xlim(-1.6, 1.25)
ax.set_xlabel("Congenital aphantasia minus reference group (SD)")
ax.grid(axis="x", color=GRID, lw=.6); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=8, loc="upper right")
# no figure title: the caption describes the figure, and a title that states the result reads as interpretation
fig.savefig(f"{ROOT}/figures/Figure_1.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------------- t2: PCL-5 items ----------------
WORD = {1: "Unwanted memories", 2: "Disturbing dreams", 3: "Reliving it (flashback)", 4: "Upset by reminders",
        5: "Physical reactions to reminders", 6: "Avoiding memories, thoughts", 7: "Avoiding external reminders",
        8: "Trouble remembering the event", 9: "Negative beliefs", 10: "Blame", 11: "Strong negative feelings",
        12: "Loss of interest", 13: "Feeling distant from others", 14: "Trouble feeling positive",
        15: "Irritable, angry", 16: "Taking risks", 17: "Superalert", 18: "Jumpy, easily startled",
        19: "Difficulty concentrating", 20: "Trouble sleeping"}
CL = [("Intrusion", 1, 5), ("Avoidance", 6, 7), ("Negative cognition and mood", 8, 14), ("Arousal", 15, 20)]
fig, ax = plt.subplots(figsize=(7.6, 7.4))
ax.axvline(0, color=TXT2, lw=1)
for k, (name, a, b) in enumerate(CL):
    if k % 2 == 0:
        ax.axhspan(a - 1.5, b - .5, color="#f6f5f2", lw=0, zorder=0)
    ax.text(1.02, (a + b) / 2 - 1, name, transform=ax.get_yaxis_transform(), va="center", fontsize=8, color=TXT2)
for r_, (lab, col, mk) in REF.items():
    s = B[(B.vs == r_) & B["var"].str.match(r"PCL_Q\d+$")].copy()
    s["i"] = s["var"].str[5:].astype(int)
    yy = s.i - 1 + (-.15 if r_ == "control" else .15)
    ax.hlines(yy, s.lo, s.hi, color=col, lw=2, alpha=.55)
    ax.scatter(s.est, yy, s=44, color=col, marker=mk, edgecolor="white", linewidth=1.2, zorder=3, label=lab)
ax.set_yticks(range(20), [f"{i}. {WORD[i]}" for i in range(1, 21)]); ax.set_ylim(19.6, -.6); ax.tick_params(length=0)
ax.set_xlabel("Congenital aphantasia minus reference group, at equal exposure (SD)")
ax.grid(axis="x", color=GRID, lw=.6); ax.set_axisbelow(True)
ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.savefig(f"{ROOT}/figures/Figure_2.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------------- t3: slopes ----------------
# rebuilds only the columns this panel needs from scored.csv (the same definitions as script 04)
d = pd.read_csv(f"{ROOT}/output/scored.csv", low_memory=False)
cong = d.grp_ext.isin(["cong_floor", "cong_above"])
d["g"] = np.select([d.grp_ext == "control", cong, d.grp_ext == "an_with_imagery"], ["control", "congenital", "an_imagers"], "other")
d["CTQ_total"] = d[["CTQ_EmoAb", "CTQ_PhyAb", "CTQ_SexAb", "CTQ_EmoNeg", "CTQ_PhyNeg"]].sum(axis=1, min_count=5)
C = A[A.part == "C_slopes"]
fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"wspace": .28})
for ax, out, lab in zip(axes, ["PHQ9", "PCL5"], ["Depression (PHQ-9)", "Post-traumatic stress (PCL-5)"]):
    for g in ["control", "an_imagers", "congenital"]:
        s = d[d.g == g].dropna(subset=["CTQ_total", out])
        s = s.assign(bin=pd.qcut(s.CTQ_total.rank(method="first"), 5, labels=False))
        m = s.groupby("bin").agg(x=("CTQ_total", "mean"), y=(out, "mean"), se=(out, lambda v: v.std() / np.sqrt(len(v))))
        f = smf.ols(f"{out} ~ CTQ_total", data=s).fit()
        xs = np.linspace(s.CTQ_total.quantile(.02), s.CTQ_total.quantile(.98), 50)
        ax.plot(xs, f.params.Intercept + f.params.CTQ_total * xs, color=GRP[g][1], lw=2)
        ax.errorbar(m.x, m.y, yerr=1.96 * m.se, fmt="o", ms=6, color=GRP[g][1], mec="white", mew=1.2, elinewidth=1.2,
                    label=f"{GRP[g][0]} (n = {len(s)})")
    ax.set_title(f"{'ab'[list(axes).index(ax)]}", loc="left", fontsize=10, fontweight="bold")
    ax.set_xlabel("Childhood trauma (CTQ total)"); ax.set_ylabel(lab)
    ax.grid(color=GRID, lw=.5); ax.set_axisbelow(True)
axes[0].legend(frameon=False, fontsize=8, loc="upper left")
fig.savefig(f"{ROOT}/figures/Figure_3.png", dpi=200, bbox_inches="tight"); plt.close(fig)

# ---------------- t4: detection ----------------
F = pd.read_csv(f"{ROOT}/output/04_detection.csv")
F = F[F.cutoff_on == "PCL-5 total"]
ADJ = {"age_z + female": "Age and gender", "CTQ_total_z + LEC_direct_z + age_z + female": "+ trauma exposure",
       "CTQ_total_z + LEC_direct_z + PHQ9_z + GAD7_z + age_z + female": "+ depression and anxiety"}
fig, (a0, a1) = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [.8, 1.2], "wspace": .45})
rates = {"control": F[F.vs == "control"].rate_ref.iloc[0], "an_imagers": F[F.vs == "an_imagers"].rate_ref.iloc[0],
         "congenital": F.rate_aph.iloc[0]}
for j, g in enumerate(["control", "an_imagers", "congenital"]):
    a0.barh(j, rates[g] * 100, color=GRP[g][1], height=.6)
    a0.text(rates[g] * 100 + 1, j, f"{rates[g] * 100:.0f}%", va="center", fontsize=9, color=TXT1)
a0.set_yticks(range(3), [GRP[g][0] for g in ["control", "an_imagers", "congenital"]]); a0.invert_yaxis()
a0.set_xlabel("Above the PCL-5 cut-off of 33 (%)"); a0.set_xlim(0, 42); a0.tick_params(length=0)
a0.set_title("a", loc="left", fontsize=10, fontweight="bold")
a1.axvline(1, color=TXT2, lw=1)
for k, (adj, lab) in enumerate(ADJ.items()):
    for r_, (rl, col, mk) in REF.items():
        s = F[(F.vs == r_) & (F.adjusted_for == adj)].iloc[0]
        yy = k + (-.13 if r_ == "control" else .13)
        a1.hlines(yy, s.lo, s.hi, color=col, lw=2, alpha=.55)
        a1.scatter(s.OR, yy, s=46, color=col, marker=mk, edgecolor="white", linewidth=1.2, zorder=3, label=rl if k == 0 else None)
a1.set_xscale("log"); a1.xaxis.set_minor_formatter(plt.NullFormatter()); a1.set_xticks([.25, .5, 1, 2], ["0.25", "0.5", "1", "2"]); a1.set_xlim(.22, 2.2)
a1.set_yticks(range(3), list(ADJ.values())); a1.invert_yaxis(); a1.tick_params(length=0)
a1.text(-.02, 1.02, "Adjusted for", transform=a1.transAxes, ha="right", va="bottom", fontsize=8.5, color=TXT2)
a1.set_xlabel("Odds ratio, congenital aphantasia (log scale)"); a1.grid(axis="x", color=GRID, lw=.6); a1.set_axisbelow(True)
a1.legend(frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(.5, -.22), ncol=2)
a1.set_title("b", loc="left", fontsize=10, fontweight="bold")
fig.savefig(f"{ROOT}/figures/Figure_4.png", dpi=200, bbox_inches="tight"); plt.close(fig)
print("wrote figures/Figure_1.png to Figure_4.png")
