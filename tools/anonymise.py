"""
Build the shared, anonymised data file from the two raw QuestionPro exports.

The raw exports are not part of this repository. They contain IP addresses, e-mail fields, timestamps, country and
region, a response ID, the time taken to complete the survey, free-text answers (occupation, the description of the
worst event, when and why imagery changed) and a dream questionnaire that the paper does not use. This script keeps
an explicit allowlist of columns and nothing else:

  pid              study code: "pr" (Prolific) or "an" (Aphantasia Network) plus a running number
  src              prolific / an
  age              years
  gender_code      1 female, 2 male, 3 another gender (as in the export)
  imagery_change   "Have you noticed any change in your mental imagery abilities over your lifetime?"
                   1 no change, 2 decreased, 3 increased, 4 never had imagery, 5 other
  item responses   PSIQ, DES, MAIA, IAMI, CTQ, PHQ, GAD, PCL (items 1 to 20) and LEC, in the export's numeric codes

Respondents who reported an age under 18 are removed first, so the file holds adults only. Every item value is
converted to a number; the three stray text entries in item columns become missing, exactly as the analysis code
treats them. Row order is the export order, so the shared file reproduces the reported results exactly (bootstrap
and random-split seeds act on row positions).

Usage: python tools/anonymise.py --raw <folder with the two .xlsx exports>
Output: data/trauma_survey_anonymised.csv
"""
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FILES = {"prolific": "Prolific Participants.xlsx", "an": "AN Praticipants.xlsx"}   # original file names
ITEM = re.compile(r"^(PSIQ_[A-G]_Q\d+|DES_Q\d+|MAIA_\d+|IAMI_Q\d+|CTQ_[A-F]_Q\d+|PHQ_Q\d+|GAD_Q\d+|PCL_Q([1-9]|1\d|20)"
                  r"|LEC_Q\d+_\d+)$")
KEEP = {"Age": "age", "Gender": "gender_code", "Changes_Q1": "imagery_change"}


def load(path, src):
    df = pd.read_excel(path)
    df.columns = [c if not str(c).startswith("Unnamed") else str(df.iloc[0][c]) for c in df.columns]
    df = df.iloc[2:].reset_index(drop=True)                     # row 0 item text, row 1 the export's second header
    df.columns = [re.sub(r"^MAIA_Q(\d+)$", r"MAIA_\1", str(c)) for c in df.columns]
    items = [c for c in df.columns if ITEM.match(c)]
    out = pd.DataFrame({"src": [src] * len(df)})
    for old, new in KEEP.items():
        out[new] = pd.to_numeric(df[old], errors="coerce")
    out = pd.concat([out, df[items].apply(pd.to_numeric, errors="coerce")], axis=1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, required=True)
    a = ap.parse_args()
    parts = [load(a.raw / f, src) for src, f in FILES.items()]
    df = pd.concat(parts, ignore_index=True, sort=False)
    df = df[~(df.age < 18)].reset_index(drop=True)            # adults only
    df["pid"] = df.src.str[:2] + df.groupby("src").cumcount().map("{:04d}".format)   # numbered after the removal
    # one column order for both sources: identifiers, demographics, then items by questionnaire
    order = ["PSIQ", "DES", "MAIA", "IAMI", "CTQ", "PHQ", "GAD", "PCL", "LEC"]
    items = sorted([c for c in df.columns if ITEM.match(c)],
                   key=lambda c: (order.index(c.split("_")[0]), [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", c)]))
    df = df[["pid", "src", *KEEP.values(), *items]]
    assert df.pid.is_unique and len(df) == 1533, len(df)
    df.to_csv(ROOT / "data" / "trauma_survey_anonymised.csv", index=False, float_format="%g")
    print(f"wrote data/trauma_survey_anonymised.csv: {len(df)} rows, {df.shape[1]} columns "
          f"({len(items)} item columns)")


if __name__ == "__main__":
    main()
