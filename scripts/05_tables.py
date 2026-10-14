"""
Step 5: the paper's tables and the mediation estimates, as .csv (for checking) and as one readable .md file.

  table_participants   participant characteristics by group (M (SD) or n (%), and n with complete data)
  table_dsm_criteria   DSM-5 symptom criteria on the PCL-5: % meeting each, odds ratios at equal exposure
  table_mediation      indirect associations of childhood trauma with depression, by group
  table_alpha          Cronbach's alpha of each scale in each group
  table_biased_items   items answered differently in aphantasia at the same level of their factor

Input : output/scored.csv, output/02_item_bias.csv, output/04_dsm_criteria.csv, output/04_routes.csv
Output: output/05_table_*.csv, output/05_tables.md
"""
from pathlib import Path

import numpy as np
import warnings

import pandas as pd

warnings.simplefilter("ignore", pd.errors.PerformanceWarning)   # many derived columns are added one at a time
ROOT = Path(__file__).resolve().parents[1].as_posix()
d = pd.read_csv(f"{ROOT}/output/scored.csv", low_memory=False)
cong = d.grp_ext.isin(["cong_floor", "cong_above"])
d["g"] = np.select([cong, d.grp_ext == "control", d.grp_ext == "an_with_imagery"],
                   ["congenital", "control", "an_imagers"], "other")
G = ["congenital", "control", "an_imagers"]
COLS = {"congenital": "Congenital aphantasia", "control": "Prolific controls", "an_imagers": "Network members with imagery"}
d = d[d.g != "other"]

# ---------------- scores as reported in the paper ----------------
item = lambda pre, idx: [f"{pre}{i}" for i in idx]
D_ = lambda idx: d[item("DES_Q", idx)].mean(axis=1, skipna=False)
PCL = lambda idx: d[item("PCL_Q", idx)].sum(axis=1, min_count=len(idx))
d["CTQ_total"] = d[["CTQ_EmoAb", "CTQ_PhyAb", "CTQ_SexAb", "CTQ_EmoNeg", "CTQ_PhyNeg"]].sum(axis=1, min_count=5)
d["PSIQ_vis"] = d[item("PSIQ_A_Q", range(1, 6))].mean(axis=1, skipna=False)
d["PCL_B"], d["PCL_C"], d["PCL_D"], d["PCL_E"] = PCL(range(1, 6)), PCL([6, 7]), PCL(range(8, 15)), PCL(range(15, 21))
d["DES_amn6"], d["DES_dep6"], d["DES_abs6"] = D_([3, 4, 5, 8, 25, 26]), D_([7, 11, 12, 13, 27, 28]), D_([2, 15, 16, 20, 22, 23])
d["any_event"] = (d.LEC_direct > 0).astype(float).where(d.LEC_direct.notna())
d["phq10"] = (d.PHQ9 >= 10).astype(float).where(d.PHQ9.notna())
d["gad10"] = (d.GAD7 >= 10).astype(float).where(d.GAD7.notna())
d["pcl33"] = (d.PCL5 >= 33).astype(float).where(d.PCL5.notna())
for gname in ["female", "male", "other"]:
    d[f"gender_{gname}"] = (d.gender == gname).astype(float)          # share of the whole group, missing as "not"


def msd(v, digits=1):
    return {g: f"{d.loc[d.g == g, v].mean():.{digits}f} ({d.loc[d.g == g, v].std():.{digits}f})" for g in G}


def npct(v, whole_group=False):
    out = {}
    for g in G:
        s = d.loc[d.g == g, v]
        s = s.fillna(0) if whole_group else s.dropna()
        out[g] = f"{int(s.sum())} ({100 * s.mean():.1f})"
    return out


def nn(v, whole_group=False):
    return " / ".join(str(int((d.g == g).sum()) if whole_group else int(d.loc[d.g == g, v].notna().sum())) for g in G)


ROWS = [("Demographics", None, None, None),
        ("Age, years", "age", msd, False),
        ("Women, n (%)", "gender_female", npct, True), ("Men, n (%)", "gender_male", npct, True),
        ("Another gender, n (%)", "gender_other", npct, True),
        ("Imagery (Psi-Q)", None, None, None),
        ("All modalities", "PSIQ", lambda v: msd(v, 2), False), ("Visual subscale", "PSIQ_vis", lambda v: msd(v, 2), False),
        ("Trauma exposure", None, None, None),
        ("CTQ total (25-125)", "CTQ_total", msd, False), ("Emotional abuse", "CTQ_EmoAb", msd, False),
        ("Physical abuse", "CTQ_PhyAb", msd, False), ("Sexual abuse", "CTQ_SexAb", msd, False),
        ("Emotional neglect", "CTQ_EmoNeg", msd, False), ("Physical neglect", "CTQ_PhyNeg", msd, False),
        ("LEC-5 events (0-17)", "LEC_direct", msd, False), ("At least one event, n (%)", "any_event", npct, False),
        ("Depression and anxiety", None, None, None),
        ("PHQ-9 (0-27)", "PHQ9", msd, False), ("PHQ-9 >= 10, n (%)", "phq10", npct, False),
        ("GAD-7 (0-21)", "GAD7", msd, False), ("GAD-7 >= 10, n (%)", "gad10", npct, False),
        ("Post-traumatic stress (PCL-5)", None, None, None),
        ("Total (0-80)", "PCL5", msd, False), ("Intrusion (0-20)", "PCL_B", msd, False), ("Avoidance (0-8)", "PCL_C", msd, False),
        ("Negative cognition and mood (0-28)", "PCL_D", msd, False), ("Arousal and reactivity (0-24)", "PCL_E", msd, False),
        ("Total >= 33, n (%)", "pcl33", npct, False),
        ("Dissociation (DES, 0-100)", None, None, None),
        ("Amnesia", "DES_amn6", msd, False), ("Depersonalisation and derealisation", "DES_dep6", msd, False),
        ("Absorption", "DES_abs6", msd, False), ("Item 9, no memory of life events", "DES_Q9", msd, False),
        ("Involuntary autobiographical memory (IAMI, 1-5)", None, None, None),
        ("Past memories", "IAMI_past", lambda v: msd(v, 2), False), ("Future thoughts", "IAMI_future", lambda v: msd(v, 2), False)]
rows = []
for lab, v, f, whole in ROWS:
    if v is None:
        rows.append({"Measure": lab, **{COLS[g]: "" for g in G}, "n": ""})
        continue
    cells = f(v, whole) if f is npct else f(v)
    rows.append({"Measure": lab, **{COLS[g]: cells[g] for g in G}, "n": nn(v, whole)})
T1 = pd.DataFrame(rows)
T1.columns = ["Measure"] + [f"{COLS[g]} (n = {int((d.g == g).sum())})" for g in G] + ["n"]

# ---------------- DSM-5 criteria ----------------
H = pd.read_csv(f"{ROOT}/output/04_dsm_criteria.csv")
LAB = {"B intrusion": "B, intrusion (>= 1 of 1-5)", "C avoidance": "C, avoidance (>= 1 of 6-7)",
       "D cognition and mood": "D, cognition/mood (>= 2 of 8-14)", "E arousal": "E, arousal (>= 2 of 15-20)",
       "D and E only": "D and E", "all four (provisional PTSD)": "All (provisional PTSD)"}
ci = lambda r: f"{r.OR_at_equal_exposure:.2f} [{r.lo:.2f}, {r.hi:.2f}]"
rows = []
for crit, lab in LAB.items():
    c, a = H[(H.criterion == crit) & (H.vs == "control")].iloc[0], H[(H.criterion == crit) & (H.vs == "an_imagers")].iloc[0]
    rows.append({"Criterion": lab, "Congenital aphantasia %": f"{100 * c.rate_aph:.1f}", "Prolific controls %": f"{100 * c.rate_ref:.1f}",
                 "Network members %": f"{100 * a.rate_ref:.1f}", "OR vs controls": ci(c), "OR vs Network with imagery": ci(a)})
T3 = pd.DataFrame(rows)

# ---------------- mediation ----------------
Dm = pd.read_csv(f"{ROOT}/output/04_routes.csv")
PATH = {"indirect via intrusions": "Through intrusions", "direct": "Direct"}
rows = []
for spec, title in [("with numbing", "Second mediator: trauma amnesia, numbing and concentration (items 8, 12, 13, 14, 19)"),
                    ("without numbing", "Second mediator: trauma amnesia and concentration (items 8 and 19)")]:
    rows.append({"Path": title, **{COLS[g]: "" for g in G}})
    s = Dm[Dm.spec == spec]
    for path in s.path.unique():
        lab = PATH.get(path, "Through second mediator")
        cell = {g: (lambda r: f"{r.est:.2f} [{r.lo:.2f}, {r.hi:.2f}]")(s[(s.path == path) & (s.group == g)].iloc[0]) for g in G}
        rows.append({"Path": lab, **{COLS[g]: cell[g] for g in G}})
T2 = pd.DataFrame(rows)
nmed = Dm.dropna(subset=["n"]).drop_duplicates("group").set_index("group").n.astype(int)
T2.columns = ["Path"] + [f"{COLS[g]} (n = {nmed[g]})" for g in G]


# ---------------- Cronbach's alpha ----------------
def alpha(cols, g):
    x = d.loc[d.g == g, cols].dropna()
    k = len(cols)
    return k / (k - 1) * (1 - x.var(ddof=1).sum() / x.sum(axis=1).var(ddof=1))


ctq_clin = [c for c in d.columns if c.startswith(("CTQ_A_", "CTQ_B_", "CTQ_C_", "CTQ_D_", "CTQ_E_"))]
SCALES = [("Psi-Q, all items", [c for c in d.columns if c.startswith("PSIQ_") and "_Q" in c]),
          ("Psi-Q, visual subscale", item("PSIQ_A_Q", range(1, 6))),
          ("CTQ total", ctq_clin), *[(n_, [c for c in ctq_clin if c.startswith(p)]) for n_, p in
                                     [("Emotional abuse", "CTQ_A_"), ("Physical abuse", "CTQ_B_"), ("Sexual abuse", "CTQ_C_"),
                                      ("Emotional neglect", "CTQ_D_"), ("Physical neglect", "CTQ_E_")]],
          ("PHQ-9", item("PHQ_Q", range(1, 10))), ("GAD-7", item("GAD_Q", range(1, 8))),
          ("PCL-5 total", item("PCL_Q", range(1, 21))), ("Intrusion", item("PCL_Q", range(1, 6))),
          ("Avoidance", item("PCL_Q", [6, 7])), ("Negative cognition and mood", item("PCL_Q", range(8, 15))),
          ("Arousal and reactivity", item("PCL_Q", range(15, 21))),
          ("DES amnesia", item("DES_Q", [3, 4, 5, 8, 25, 26])),
          ("DES depersonalisation and derealisation", item("DES_Q", [7, 11, 12, 13, 27, 28])),
          ("DES absorption", item("DES_Q", [2, 15, 16, 20, 22, 23])),
          ("IAMI past memories", item("IAMI_Q", [2, 4, 8, 9, 11, 13, 14, 17, 19, 20])),
          ("IAMI future thoughts", item("IAMI_Q", [1, 3, 5, 6, 7, 10, 12, 15, 16, 18]))]
TA = pd.DataFrame([{"Scale (items)": f"{n_} ({len(c)})", **{COLS[g]: f"{alpha(c, g):.2f}" for g in G}} for n_, c in SCALES])

# ---------------- biased items ----------------
F = pd.read_csv(f"{ROOT}/output/02_item_bias.csv", index_col=0)
F = F[F.biased_both].sort_values("dif_vs_prolific")
TB = pd.DataFrame({"Item": F.index, "Against Prolific controls": F.dif_vs_prolific.map("{:+.2f}".format),
                   "Against Network members with imagery": F.dif_vs_an_imagers.map("{:+.2f}".format)})

L = ["# 05 - tables\n"]
for name, t in [("participants", T1), ("dsm_criteria", T3), ("mediation", T2), ("alpha", TA), ("biased_items", TB)]:
    t.to_csv(f"{ROOT}/output/05_table_{name}.csv", index=False)
    L += [f"\n## table_{name}\n", t.to_markdown(index=False)]
Path(f"{ROOT}/output/05_tables.md").write_text("\n".join(L), encoding="utf-8")
print(f"wrote output/05_table_*.csv ({len(L) // 2} tables) and output/05_tables.md")
