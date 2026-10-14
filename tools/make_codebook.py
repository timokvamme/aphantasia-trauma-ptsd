"""
Write data/codebook.csv: one row per column of data/trauma_survey_anonymised.csv.

Item wording is not reproduced, because several questionnaires are copyrighted; the item numbers follow the published
versions cited in the paper.
Usage: python tools/make_codebook.py
"""
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
df = pd.read_csv(ROOT / "data" / "trauma_survey_anonymised.csv")

FIXED = {
    "pid": ("Study code", "", "'pr' (Prolific) or 'an' (Aphantasia Network) plus a running number", ""),
    "src": ("Recruitment source", "", "prolific = general-population sample; an = Aphantasia Network", ""),
    "age": ("Age", "", "years", ""),
    "gender_code": ("Gender", "", "1 female, 2 male, 3 another gender", ""),
    "imagery_change": ("Lifetime change in mental imagery", "",
                       "1 no change, 2 decreased, 3 increased, 4 never had imagery, 5 other", ""),
}
PSIQ_MOD = {"A": "vision", "B": "sound", "C": "smell", "D": "taste", "E": "touch", "F": "bodily sensation", "G": "emotional feeling"}
CTQ_SUB = {"A": "emotional abuse", "B": "physical abuse", "C": "sexual abuse", "D": "emotional neglect",
           "E": "physical neglect", "F": "minimisation/denial"}
LEC_COL = {1: "happened to me", 2: "witnessed it", 3: "learned about it", 4: "part of my job", 5: "not sure",
           6: "does not apply"}
DES_SUB = {**{i: "amnesia" for i in [3, 4, 5, 8, 25, 26]}, **{i: "depersonalisation/derealisation" for i in [7, 11, 12, 13, 27, 28]},
           **{i: "absorption" for i in [2, 14, 15, 16, 17, 18, 20, 22, 23]}}
MAIA_SUB = [("noticing", 1, 4), ("not-distracting", 5, 10), ("not-worrying", 11, 15), ("attention regulation", 16, 22),
            ("emotional awareness", 23, 27), ("self-regulation", 28, 31), ("body listening", 32, 34), ("trusting", 35, 37)]
FUTURE = {1, 3, 5, 6, 7, 10, 12, 15, 16, 18}


def describe(c):
    if c in FIXED:
        return FIXED[c]
    if m := re.match(r"PSIQ_([A-G])_Q(\d+)", c):
        return ("Psi-Q", f"{PSIQ_MOD[m[1]]}, item {m[2]}", "1-11 = 0 (no image at all) to 10 (as vivid as real life)", "")
    if m := re.match(r"DES_Q(\d+)", c):
        i = int(m[1])
        return ("DES", f"item {i}" + (f" ({DES_SUB[i]})" if i in DES_SUB else ""), "1-11 = 0% to 100% of the time, in steps of 10", "")
    if m := re.match(r"MAIA_(\d+)", c):
        i = int(m[1])
        sub = next(n for n, a, b in MAIA_SUB if a <= i <= b)
        return ("MAIA-2", f"item {i} ({sub})", "1-6 = 0 (never) to 5 (always)", "yes" if i in {5, 6, 7, 8, 9, 10, 11, 12, 15} else "")
    if m := re.match(r"IAMI_Q(\d+)", c):
        i = int(m[1])
        return ("IAMI", f"item {i} ({'future thoughts' if i in FUTURE else 'past memories'})", "1 (never) to 5 (every hour)", "")
    if m := re.match(r"CTQ_([A-F])_Q(\d+)", c):
        rev = m[1] == "D" or (m[1] == "E" and m[2] in ("2", "5"))
        return ("CTQ", f"{CTQ_SUB[m[1]]}, item {m[2]}", "1 (never true) to 5 (very often true)", "yes" if rev else "")
    if m := re.match(r"PHQ_Q(\d+)", c):
        return ("PHQ-9", f"item {m[1]}", "1-4 = 0 (not at all) to 3 (nearly every day)", "")
    if m := re.match(r"GAD_Q(\d+)", c):
        return ("GAD-7", f"item {m[1]}", "1-4 = 0 (not at all) to 3 (nearly every day)", "")
    if m := re.match(r"PCL_Q(\d+)", c):
        return ("PCL-5", f"item {m[1]}", "1-5 = 0 (not at all) to 4 (extremely)", "")
    if m := re.match(r"LEC_Q(\d+)_(\d)", c):
        return ("LEC-5", f"event {m[1]}, '{LEC_COL[int(m[2])]}'", "1 = ticked, empty = not ticked (8 cells hold other values from the survey export and are not counted as ticks)", "")
    raise ValueError(f"no codebook entry for column {c}")


rows = [dict(column=c, questionnaire_or_variable=a, item=b, coding=k, reverse_keyed=r) for c in df.columns
        for a, b, k, r in [describe(c)]]
pd.DataFrame(rows).to_csv(ROOT / "data" / "codebook.csv", index=False)
print(f"wrote data/codebook.csv: {len(rows)} columns described (data file has {df.shape[1]})")
