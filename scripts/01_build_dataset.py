"""
Step 1: score every questionnaire and define the groups.

Input : data/trauma_survey_anonymised.csv (one row per participant, item responses in the survey's numeric codes)
Output: output/scored.csv          participant-level scores and rescaled items, used by every later script
        output/01_descriptives.md  group sizes and block completeness

Item coding in the data file and the scoring used here:
  Psi-Q  35 items, coded 1-11, rescaled to 0-10; seven modalities A-G (vision, sound, smell, taste, touch, bodily
         sensation, emotional feeling), each the mean of its five items; total = mean of all 35
  DES    28 items, coded 1-11, rescaled to 0-100 (% of the time); subscales are item means (Carlson & Putnam, 1993):
         amnesia 3, 4, 5, 8, 25, 26; depersonalisation/derealisation 7, 11, 12, 13, 27, 28;
         absorption 2, 14, 15, 16, 17, 18, 20, 22, 23
  MAIA-2 37 items, coded 1-6, rescaled to 0-5; reverse-keyed 5-12 and 15 (Mehling et al., 2018)
  CTQ    25 clinical items + 3 minimisation items, coded 1-5; subscale sums; emotional neglect (all items) and
         physical neglect items 2 and 5 reverse-keyed (Bernstein et al., 2003)
  PHQ-9, GAD-7  coded 1-4, rescaled to 0-3, summed
  PCL-5  20 items, coded 1-5, rescaled to 0-4, summed
  IAMI   20 items, coded 1-5, mean; future-thought items 1, 3, 5, 6, 7, 10, 12, 15, 16, 18, past-memory items the rest
  LEC-5  17 events; count of events that happened to the participant (column _1) or that they witnessed (_2)
Out-of-range codes become missing. Scale scores need complete items, except where stated.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1].as_posix()
(Path(ROOT) / "output").mkdir(exist_ok=True)
df = pd.read_csv(f"{ROOT}/data/trauma_survey_anonymised.csv")
num = lambda cols: df[cols].apply(pd.to_numeric, errors="coerce")

# ---------------- items, rescaled to each instrument's published metric ----------------
maia_cols = [f"MAIA_{i}" for i in range(1, 38)]
M = num(maia_cols) - 1                                     # 0-5
M[(M < 0) | (M > 5)] = np.nan
for i in [5, 6, 7, 8, 9, 10, 11, 12, 15]:
    M[f"MAIA_{i}"] = 5 - M[f"MAIA_{i}"]
maia_key = {"MAIA_Not": [1, 2, 3, 4], "MAIA_NDis": [5, 6, 7, 8, 9, 10], "MAIA_NWor": [11, 12, 13, 14, 15],
            "MAIA_Att": list(range(16, 23)), "MAIA_Emo": list(range(23, 28)), "MAIA_Self": [28, 29, 30, 31],
            "MAIA_Body": [32, 33, 34], "MAIA_Trust": [35, 36, 37]}

des_cols = [f"DES_Q{i}" for i in range(1, 29)]
D = (num(des_cols) - 1) * 10                               # 0-100
D[(D < 0) | (D > 100)] = np.nan
des_key = {"DES_Amn": [3, 4, 5, 8, 25, 26], "DES_Dep": [7, 11, 12, 13, 27, 28],
           "DES_Abs": [2, 14, 15, 16, 17, 18, 20, 22, 23]}

ctq_cols = [c for c in df.columns if c.startswith("CTQ_")]
C = num(ctq_cols)
C[(C < 1) | (C > 5)] = np.nan
for c in [f"CTQ_D_Q{i}" for i in range(1, 6)] + ["CTQ_E_Q2", "CTQ_E_Q5"]:
    C[c] = 6 - C[c]
ctq_key = {"CTQ_EmoAb": "CTQ_A", "CTQ_PhyAb": "CTQ_B", "CTQ_SexAb": "CTQ_C", "CTQ_EmoNeg": "CTQ_D",
           "CTQ_PhyNeg": "CTQ_E", "CTQ_Denial": "CTQ_F"}

phq_cols = [f"PHQ_Q{i}" for i in range(1, 10)]
P9 = num(phq_cols) - 1; P9[(P9 < 0) | (P9 > 3)] = np.nan
gad_cols = [f"GAD_Q{i}" for i in range(1, 8)]
G7 = num(gad_cols) - 1; G7[(G7 < 0) | (G7 > 3)] = np.nan
pcl_cols = [f"PCL_Q{i}" for i in range(1, 21)]
PC = num(pcl_cols) - 1; PC[(PC < 0) | (PC > 4)] = np.nan
iami_cols = [f"IAMI_Q{i}" for i in range(1, 21)]
IA = num(iami_cols); IA[(IA < 1) | (IA > 5)] = np.nan
psiq_cols = [c for c in df.columns if c.startswith("PSIQ_")]
PS = num(psiq_cols) - 1; PS[(PS < 0) | (PS > 10)] = np.nan
psiq_mod = {"A": "app", "B": "sound", "C": "smell", "D": "taste", "E": "touch", "F": "body", "G": "emo"}

items = pd.concat([M, D, C, P9, G7, PC, IA, PS], axis=1)

# ---------------- scale scores ----------------
S = pd.DataFrame(index=df.index)
S["pid"] = df["pid"]; S["src"] = df["src"]
S["age"] = df["age"]
S["gender"] = df["gender_code"].map({1: "female", 2: "male", 3: "other"})
S["chg_q1"] = df["imagery_change"]                         # 1 no change, 2 decreased, 3 increased, 4 never had, 5 other
for k, its in maia_key.items():
    S[k] = M[[f"MAIA_{i}" for i in its]].mean(axis=1, skipna=False)
for k, its in des_key.items():
    S[k] = D[[f"DES_Q{i}" for i in its]].mean(axis=1, skipna=False)
S["DES_total"] = D.mean(axis=1, skipna=False)
for k, pre in ctq_key.items():
    S[k] = C[[c for c in ctq_cols if c.startswith(pre)]].sum(axis=1, min_count=5 if pre != "CTQ_F" else 3)
S["PHQ9"] = P9.sum(axis=1, min_count=9)
S["GAD7"] = G7.sum(axis=1, min_count=7)
S["PCL5"] = PC.sum(axis=1, min_count=20)
S["IAMI"] = IA.mean(axis=1, skipna=False)
fut = [f"IAMI_Q{i}" for i in [1, 3, 5, 6, 7, 10, 12, 15, 16, 18]]
S["IAMI_future"] = IA[fut].mean(axis=1, skipna=False)
S["IAMI_past"] = IA[[c for c in iami_cols if c not in fut]].mean(axis=1, skipna=False)
for m, name in psiq_mod.items():
    S[f"PSIQ_{name}"] = PS[[c for c in psiq_cols if c.startswith(f"PSIQ_{m}_")]].mean(axis=1, skipna=False)
S["PSIQ"] = PS.mean(axis=1, skipna=False)                  # 0-10
S["PSIQ_floor"] = (PS.max(axis=1) == 0) & S["PSIQ"].notna()

# LEC-5: events that happened to the participant (_1) or that they witnessed (_2); missing if no LEC answer at all
lec_cols = [c for c in df.columns if c.startswith("LEC_Q")]
hit = lambda c: df[c].eq(1)
S["LEC_direct"] = sum((hit(f"LEC_Q{e}_1") | hit(f"LEC_Q{e}_2")).astype(int) for e in range(1, 18))
S.loc[df[lec_cols].isna().all(axis=1), "LEC_direct"] = np.nan

# ---------------- groups (Psi-Q on the 0-10 metric) ----------------
def group(r):
    if pd.isna(r.PSIQ): return "no_psiq"
    if r.src == "prolific":
        return "control" if r.PSIQ > 4 else "prolific_low"
    if r.PSIQ > 4: return "an_with_imagery"
    if r.chg_q1 == 4: return "cong_floor" if r.PSIQ_floor else "cong_above"
    if r.chg_q1 == 2: return "acquired"
    if r.chg_q1 == 1: return "an_low_nochange"
    return "an_low_other"
S["grp_ext"] = S.apply(group, axis=1)

blocks = {"MAIA": maia_cols, "DES": des_cols, "CTQ": ctq_cols, "PHQ": phq_cols, "GAD": gad_cols,
          "PCL": pcl_cols, "IAMI": iami_cols, "PSIQ": psiq_cols}
for b, cols in blocks.items():
    S[f"ok_{b}"] = items[cols].notna().all(axis=1)

pd.concat([S, items], axis=1).to_csv(f"{ROOT}/output/scored.csv", index=False)

# ---------------- descriptives ----------------
L = ["# 01 - dataset build\n", f"Rows: {S.src.value_counts().to_dict()}, total {len(S)}\n",
     "## Groups (Psi-Q on the 0-10 metric, cut-off 4)\n",
     S.groupby("grp_ext").agg(n=("pid", "size"), psiq=("PSIQ", "mean"), age=("age", "mean"),
                              phq=("PHQ9", "mean"), gad=("GAD7", "mean"), pcl=("PCL5", "mean")).round(2).to_markdown(),
     "\n\n## Complete questionnaire blocks by source\n",
     S.groupby("src")[[f"ok_{b}" for b in blocks]].sum().to_markdown()]
Path(f"{ROOT}/output/01_descriptives.md").write_text("\n".join(L), encoding="utf-8")
print("\n".join(L[:3]))
