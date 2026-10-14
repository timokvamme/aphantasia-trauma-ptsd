"""
Step 4: every trauma analysis in the paper.

A. Exposure: childhood trauma (CTQ) and lifetime traumatic events (LEC-5), Cohen d.
B. Symptom form: each PCL-5 cluster and item, DES subscale, depression, anxiety and IAMI score, compared at equal
   exposure (CTQ, LEC-5) and adjusted for age and gender; also for the core aphantasia subgroup.
C. Burden: slopes of depression, anxiety and PCL-5 total on exposure, by group, with an equivalence test (TOST).
D. Indirect associations of childhood trauma with depression: parallel mediation, bootstrap CIs, per group.
E. Trauma amnesia (PCL-5 item 8) as a predictor of depression and numbing.
F. Screening: odds of a PCL-5 total of 33 or more, with increasing adjustment.
H. DSM-5 symptom rule applied to the PCL-5.
I. Involuntary memory (IAMI) and the group difference in intrusions.
K. Screening among participants with clinical-range depression or anxiety.
L. Depersonalisation/derealisation against each distress measure.
M. Imagery vividness and depersonalisation within the groups with imagery.
N. Trauma amnesia with and without general autobiographical gaps (DES item 9).

Groups: Prolific controls (Psi-Q > 4); congenital aphantasia (Aphantasia Network, never had imagery, Psi-Q <= 4);
within it, core = every Psi-Q visual item at 0; Aphantasia Network members with imagery (Psi-Q > 4) as the second
reference group. Group differences are in standard deviations of the reference group of each comparison; slopes and
indirect effects are in standard deviations of the controls.
Input : output/scored.csv.  Output: output/04_*.csv, output/04_trauma_analyses.md
"""
from pathlib import Path

import numpy as np
import warnings

import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm, t as t_dist

warnings.simplefilter("ignore", pd.errors.PerformanceWarning)   # many derived columns are added one at a time
ROOT = Path(__file__).resolve().parents[1].as_posix()
EQ_BOUND = 0.20   # smallest slope difference of interest for the equivalence test (control SD per control SD)
RNG = np.random.default_rng(20261002)
d = pd.read_csv(f"{ROOT}/output/scored.csv", low_memory=False)

# ---------------- groups and scores ----------------
VIS = [f"PSIQ_A_Q{i}" for i in range(1, 6)]
cong = d.grp_ext.isin(["cong_floor", "cong_above"])
d["g"] = np.select([d.grp_ext == "control", cong, d.grp_ext == "an_with_imagery"],
                   ["control", "congenital", "an_imagers"], "other")
d["core"] = cong & d[VIS].max(axis=1).eq(0)
d["CTQ_total"] = d[["CTQ_EmoAb", "CTQ_PhyAb", "CTQ_SexAb", "CTQ_EmoNeg", "CTQ_PhyNeg"]].sum(axis=1, min_count=5)
P = lambda idx: d[[f"PCL_Q{i}" for i in idx]].mean(axis=1, skipna=False)
d["PCL_intrusion"], d["PCL_avoidance"], d["PCL_negmood"], d["PCL_arousal"] = P(range(1, 6)), P([6, 7]), P(range(8, 15)), P(range(15, 21))
d["PCL_imagery_route"] = P(range(1, 6))                       # intrusions, nightmares, flashbacks, cue reactivity
d["PCL_other_route"] = P([8, 12, 13, 14, 19])                 # amnesia, loss of interest, detachment, no positive affect, concentration
d["PCL_numbing"] = P([12, 13, 14])
# total without the two items answered differently in aphantasia (3 and 8), prorated to 20 items
d["PCL_no_biased"] = d[[f"PCL_Q{i}" for i in range(1, 21) if i not in (3, 8)]].sum(axis=1, min_count=18) * 20 / 18
# total without item 8 (trouble remembering the event), which people with aphantasia endorse more
d["PCL_no8"] = d[[f"PCL_Q{i}" for i in range(1, 21) if i != 8]].sum(axis=1, min_count=19) * 20 / 19
D_ = lambda idx: d[[f"DES_Q{i}" for i in idx]].mean(axis=1, skipna=False)
# DES subscales (Carlson & Putnam, 1993); absorption without the three items answered differently in aphantasia
# (14, 17, 18; script 02). Item 9, the fourth such item, belongs to no subscale and is analysed on its own.
d["DES_dep"], d["DES_amn"], d["DES_abs_unbiased"] = D_([7, 11, 12, 13, 27, 28]), D_([3, 4, 5, 8, 25, 26]), D_([2, 15, 16, 20, 22, 23])
d["female"] = (d.gender == "female").astype(float).where(d.gender.notna())   # women vs men and another gender
d["probable_ptsd"] = (d.PCL5 >= 33).astype(float).where(d.PCL5.notna())
d["probable_ptsd_nb"] = (d.PCL_no_biased >= 33).astype(float).where(d.PCL_no_biased.notna())

ref = d.g == "control"
zc = lambda v: (d[v] - d.loc[ref, v].mean()) / d.loc[ref, v].std()     # standardise on the controls
for v in ["CTQ_total", "LEC_direct", "PHQ9", "GAD7", "PCL5"]:
    d[v + "_z"] = zc(v)
d["age_z"] = (d.age - d.age.mean()) / d.age.std()


def cohen_d(a, b):
    a, b = a.dropna(), b.dropna()
    s = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(a) + len(b) - 2))
    dd = (a.mean() - b.mean()) / s
    se = np.sqrt((len(a) + len(b)) / (len(a) * len(b)) + dd ** 2 / (2 * (len(a) + len(b))))
    return dd, dd - 1.96 * se, dd + 1.96 * se, len(a), len(b), 2 * norm.sf(abs(dd / se))


def adjusted(v, grp, refg, cov=("CTQ_total_z", "LEC_direct_z", "age_z", "female")):
    """group effect on v, in standard deviations of the reference group, at equal exposure, age and gender"""
    s = d[(d.g.isin([grp, refg]) if grp != "core" else (d.core | (d.g == refg)))].copy()
    s["y"] = (s[v] - d.loc[d.g == refg, v].mean()) / d.loc[d.g == refg, v].std()
    s["grp"] = (s.g != refg).astype(float)
    m = smf.ols("y ~ grp + " + " + ".join(cov), data=s.dropna(subset=["y", *cov])).fit()
    return m.params["grp"], *m.conf_int().loc["grp"], m.pvalues["grp"], int(m.nobs)


fmt = lambda r: f"{r[0]:+.2f} [{r[1]:+.2f}, {r[2]:+.2f}]"
L = ["# 04 - trauma analyses\n",
     "Cells: Cohen d or adjusted difference in reference-group SD units, [95% CI]. Positive = higher in the aphantasia group.\n",
     d[d.g != "other"].groupby("g").agg(n=("pid", "size"), LEC_answered=("LEC_direct", "count"),
                                        PCL_answered=("PCL5", "count")).to_markdown(),
     f"\nCore congenital (Psi-Q visual items all 0): n = {int(d.core.sum())}.\n"]

# ---------------- A. exposure ----------------
EXP = {"CTQ_total": "Childhood trauma, total", "CTQ_EmoAb": "Emotional abuse", "CTQ_PhyAb": "Physical abuse",
       "CTQ_SexAb": "Sexual abuse", "CTQ_EmoNeg": "Emotional neglect", "CTQ_PhyNeg": "Physical neglect",
       "LEC_direct": "Traumatic events experienced or witnessed"}
rows = []
for v, lab in EXP.items():
    for grp, refg in [("congenital", "control"), ("congenital", "an_imagers")]:
        r = cohen_d(d.loc[d.g == grp, v], d.loc[d.g == refg, v])
        rows.append(dict(part="A_exposure", measure=lab, vs=refg, d=r[0], lo=r[1], hi=r[2], n1=r[3], n2=r[4], p=r[5]))
A = pd.DataFrame(rows)
anyev = d[d.g.isin(["control", "congenital", "an_imagers"])].groupby("g").LEC_direct.apply(lambda s: (s.dropna() > 0).mean())
L += ["## A. Trauma exposure (unadjusted Cohen d)\n",
      A.assign(cell=A.apply(lambda r: fmt((r.d, r.lo, r.hi)), axis=1)).pivot(index="measure", columns="vs", values="cell").to_markdown(),
      f"\nAt least one event experienced or witnessed: {anyev.round(3).to_dict()}\n"]

# ---------------- B. symptom form, adjusted for exposure ----------------
SYM = {"PCL5": "PTSD total", "PCL_no8": "PTSD total without item 8","PCL_intrusion": "Intrusion (B)", "PCL_avoidance": "Avoidance (C)",
       "PCL_negmood": "Negative cognition and mood (D)", "PCL_arousal": "Arousal (E)",
       **{f"PCL_Q{i}": f"PCL item {i}" for i in range(1, 21)},
       "DES_dep": "Depersonalisation/derealisation", "DES_amn": "Amnesia", "DES_abs_unbiased": "Absorption, unbiased items",
       "DES_Q9": "DES 9: no memory of important events", "PHQ9": "Depression", "GAD7": "Anxiety",
       "IAMI_past": "Involuntary past memories (IAMI)", "IAMI_future": "Involuntary future thoughts (IAMI)"}
rows = []
for v, lab in SYM.items():
    for grp, refg in [("congenital", "control"), ("core", "control"), ("congenital", "an_imagers")]:
        r = adjusted(v, grp, refg)
        rows.append(dict(part="B_form", measure=lab, var=v, group=grp, vs=refg, est=r[0], lo=r[1], hi=r[2], p=r[3], n=r[4]))
B = pd.DataFrame(rows)
B.to_csv(f"{ROOT}/output/04_symptom_form.csv", index=False, float_format="%.6f")
L += ["## B. Symptom form at the same exposure (CTQ, LEC, age, gender)\n",
      B.assign(cell=B.apply(lambda r: fmt((r.est, r.lo, r.hi)) + ("*" if r.p < .05 else ""), axis=1),
               col=B.group + " vs " + B.vs).pivot(index="measure", columns="col", values="cell")
      .reindex(list(SYM.values())).to_markdown()]

# ---------------- C. burden: slopes of mood on exposure, by group ----------------
rows = []
for out in ["PHQ9_z", "GAD7_z", "PCL5_z"]:
    for exp in ["CTQ_total_z", "LEC_direct_z"]:
        for refg in ["control", "an_imagers"]:
            s = d[d.g.isin(["congenital", refg])].dropna(subset=[out, exp, "age_z", "female"]).copy()
            s["grp"] = (s.g == "congenital").astype(float)
            m = smf.ols(f"{out} ~ {exp} * grp + age_z + female", data=s).fit()
            sl_ref, inter = m.params[exp], m.params[f"{exp}:grp"]
            ci = m.conf_int().loc[f"{exp}:grp"]
            # equivalence (TOST) of the slope difference within +/- EQ_BOUND control SD; 90% CI <=> two one-sided tests at .05
            se, df_ = m.bse[f"{exp}:grp"], m.df_resid
            ci90 = m.conf_int(alpha=.10).loc[f"{exp}:grp"]
            p_tost = max(t_dist.sf((inter + EQ_BOUND) / se, df_), t_dist.cdf((inter - EQ_BOUND) / se, df_))
            rows.append(dict(outcome=out[:-2], exposure=exp[:-2], vs=refg, slope_ref=sl_ref, slope_aph=sl_ref + inter,
                             diff=inter, lo=ci[0], hi=ci[1], p=m.pvalues[f"{exp}:grp"], n=int(m.nobs),
                             lo90=ci90[0], hi90=ci90[1], p_tost=p_tost, equivalent=p_tost < .05))
C = pd.DataFrame(rows)
C.to_csv(f"{ROOT}/output/04_slopes.csv", index=False, float_format="%.6f")
L += ["\n## C. Slope of symptoms on exposure (control-SD per control-SD); diff = aphantasia minus reference\n",
      f"Equivalence: TOST of the slope difference against +/-{EQ_BOUND} SD (90% CI inside the bounds).\n",
      C.round(3).to_markdown(index=False)]

# ---------------- D. routes: CTQ -> {imagery route, other route} -> depression, per group ----------------
def mediation(s, x, ms, y, B_=2000):
    def est(t):
        X = sm.add_constant(t[[x]]).values
        a = [np.linalg.lstsq(X, t[m].values, rcond=None)[0][1] for m in ms]
        Xb = sm.add_constant(t[[x, *ms]]).values
        b = np.linalg.lstsq(Xb, t[y].values, rcond=None)[0]
        return [a[k] * b[2 + k] for k in range(len(ms))] + [b[1]]
    pt = est(s)
    bs = np.array([est(s.iloc[RNG.integers(0, len(s), len(s))]) for _ in range(B_)])
    return pt, np.percentile(bs, [2.5, 97.5], axis=0), bs


# numbing items (12-14) overlap with depression itself, so the second spec keeps only amnesia (8) and concentration (19)
d["PCL_gaps_conc"] = P([8, 19])
SPECS = {"with numbing": ("PCL_other_route", "memory gaps, numbing, concentration"),
         "without numbing": ("PCL_gaps_conc", "memory gaps, concentration")}
rows = []
for spec, (m2, lab2) in SPECS.items():
    boots, ms, labs = {}, ["PCL_imagery_route", m2], ["indirect via intrusions", f"indirect via {lab2}", "direct"]
    for grp in ["control", "congenital", "an_imagers"]:
        s = d[d.g == grp].dropna(subset=["CTQ_total", *ms, "PHQ9"])
        s = (s[["CTQ_total", *ms, "PHQ9"]] - d.loc[ref, ["CTQ_total", *ms, "PHQ9"]].mean()) / d.loc[ref, ["CTQ_total", *ms, "PHQ9"]].std()
        pt, ci, bs = mediation(s, "CTQ_total", ms, "PHQ9")
        boots[grp] = (pt, bs)
        for k, lab in enumerate(labs):
            rows.append(dict(spec=spec, group=grp, n=len(s), path=lab, est=pt[k], lo=ci[0][k], hi=ci[1][k]))
    for refg in ["control", "an_imagers"]:
        for k, lab in enumerate(labs[:2]):
            diff = boots["congenital"][1][:, k] - boots[refg][1][:, k]
            rows.append(dict(spec=spec, group=f"congenital minus {refg}", n=np.nan, path=lab,
                             est=boots["congenital"][0][k] - boots[refg][0][k],
                             lo=np.percentile(diff, 2.5), hi=np.percentile(diff, 97.5)))
Dm = pd.DataFrame(rows)
Dm.to_csv(f"{ROOT}/output/04_routes.csv", index=False, float_format="%.6f")
L += ["\n## D. Routes from childhood trauma to depression (standardised on controls; bootstrap 95% CI, 2000)\n",
      Dm.round(3).to_markdown(index=False)]

# ---------------- E. memory gaps: PCL-5 item 8 with depression and numbing ----------------
rows = []
for grp in ["control", "congenital", "an_imagers"]:
    s = d[d.g == grp].dropna(subset=["PCL_Q8", "PHQ9", "PCL_numbing", "CTQ_total", "LEC_direct", "age", "female"])
    for y in ["PHQ9", "PCL_numbing", "PCL_intrusion"]:
        m = smf.ols(f"scale({y}) ~ scale(PCL_Q8) + scale(CTQ_total) + scale(LEC_direct) + scale(age) + female", data=s).fit()
        rows.append(dict(group=grp, outcome=y, beta=m.params["scale(PCL_Q8)"], lo=m.conf_int().loc["scale(PCL_Q8)"][0],
                         hi=m.conf_int().loc["scale(PCL_Q8)"][1], p=m.pvalues["scale(PCL_Q8)"], n=int(m.nobs)))
E = pd.DataFrame(rows)
L += ["\n## E. Amnesia for the trauma (PCL-5 item 8) as a predictor, adjusted for exposure, age, gender (standardised beta)\n",
      E.round(3).to_markdown(index=False)]

# ---------------- F. detection by the PCL-5 cut-off ----------------
rows = []
for outv in ["probable_ptsd", "probable_ptsd_nb"]:
    for refg in ["control", "an_imagers"]:
        s = d[d.g.isin(["congenital", refg])].dropna(subset=[outv, "CTQ_total_z", "LEC_direct_z", "PHQ9_z", "GAD7_z", "age_z"]).copy()
        s["grp"] = (s.g == "congenital").astype(float)
        rate = s.groupby("g")[outv].mean()
        for adj in ["grp + age_z + female", "grp + CTQ_total_z + LEC_direct_z + age_z + female",
                    "grp + CTQ_total_z + LEC_direct_z + PHQ9_z + GAD7_z + age_z + female"]:
            m = smf.logit(f"{outv} ~ {adj}", data=s).fit(disp=0)
            rows.append(dict(cutoff_on="PCL-5 total" if outv == "probable_ptsd" else "PCL-5 without items 3 and 8",
                             vs=refg, rate_aph=rate["congenital"], rate_ref=rate[refg], adjusted_for=adj.replace("grp + ", ""),
                             OR=np.exp(m.params["grp"]), lo=np.exp(m.conf_int().loc["grp"][0]),
                             hi=np.exp(m.conf_int().loc["grp"][1]), p=m.pvalues["grp"], n=int(m.nobs)))
F = pd.DataFrame(rows)
F.to_csv(f"{ROOT}/output/04_detection.csv", index=False, float_format="%.6f")
L += ["\n## F. Above the PCL-5 cut-off (>= 33): odds ratio for congenital aphantasia\n", F.round(3).to_markdown(index=False)]

# ---------------- H. DSM-5 symptom rule: which criterion keeps aphantasics out? ----------------
# PCL-5 provisional diagnosis (Weathers et al.): an item counts at >= 2 ("moderately");
# B >= 1 of items 1-5, C >= 1 of 6-7, D >= 2 of 8-14, E >= 2 of 15-20.
cnt = lambda idx: (d[[f"PCL_Q{i}" for i in idx]] >= 2).sum(axis=1).where(d[[f"PCL_Q{i}" for i in idx]].notna().all(axis=1))
d["crit_B"], d["crit_C"] = (cnt(range(1, 6)) >= 1).astype(float), (cnt([6, 7]) >= 1).astype(float)
d["crit_D"], d["crit_E"] = (cnt(range(8, 15)) >= 2).astype(float), (cnt(range(15, 21)) >= 2).astype(float)
for c_ in ["crit_B", "crit_C", "crit_D", "crit_E"]:
    d[c_] = d[c_].where(d.PCL5.notna())
d["dsm_full"] = d[["crit_B", "crit_C", "crit_D", "crit_E"]].min(axis=1)
d["dsm_DE"] = d[["crit_D", "crit_E"]].min(axis=1)                       # the non-imagery clusters only
d["dsm_missing_only_B"] = ((d.crit_B == 0) & (d.crit_C == 1) & (d.dsm_DE == 1)).astype(float).where(d.PCL5.notna())
rows = []
for outv, lab in [("crit_B", "B intrusion"), ("crit_C", "C avoidance"), ("crit_D", "D cognition and mood"),
                  ("crit_E", "E arousal"), ("dsm_full", "all four (provisional PTSD)"), ("dsm_DE", "D and E only"),
                  ("dsm_missing_only_B", "C, D and E met, B missing")]:
    for refg in ["control", "an_imagers"]:
        s = d[d.g.isin(["congenital", refg])].dropna(subset=[outv, "CTQ_total_z", "LEC_direct_z", "age_z"]).copy()
        s["grp"] = (s.g == "congenital").astype(float)
        rate = s.groupby("g")[outv].mean()
        m = smf.logit(f"{outv} ~ grp + CTQ_total_z + LEC_direct_z + age_z + female", data=s).fit(disp=0)
        # same model plus depression and anxiety, to compare with the cut-off result at equal distress (part F)
        s2 = s.dropna(subset=["PHQ9_z", "GAD7_z"])
        m2 = smf.logit(f"{outv} ~ grp + CTQ_total_z + LEC_direct_z + PHQ9_z + GAD7_z + age_z + female", data=s2).fit(disp=0)
        rows.append(dict(criterion=lab, vs=refg, rate_aph=rate["congenital"], rate_ref=rate[refg],
                         OR_at_equal_exposure=np.exp(m.params["grp"]), lo=np.exp(m.conf_int().loc["grp"][0]),
                         hi=np.exp(m.conf_int().loc["grp"][1]), p=m.pvalues["grp"], n=int(m.nobs),
                         OR_equal_exposure_mood=np.exp(m2.params["grp"]), lo_mood=np.exp(m2.conf_int().loc["grp"][0]),
                         hi_mood=np.exp(m2.conf_int().loc["grp"][1]), p_mood=m2.pvalues["grp"], n_mood=int(m2.nobs)))
H = pd.DataFrame(rows)
H.to_csv(f"{ROOT}/output/04_dsm_criteria.csv", index=False, float_format="%.6f")
L += ["\n## H. DSM-5 symptom rule on the PCL-5 (item >= 2): share meeting each criterion, OR at equal exposure\n",
      H.round(3).to_markdown(index=False)]

# ---------------- I. general involuntary memory (IAMI): does it account for the fewer trauma intrusions? ----------------
# If aphantasics have fewer trauma intrusions because they have fewer involuntary memories in general, adding the
# IAMI should shrink the group effect on intrusion (and avoidance) but not on mood and arousal.
# Caveat: the IAMI met only partial scalar invariance against the controls (script 03), so its group mean is not
# fully comparable there; the within-group slopes do not depend on that.
for v in ["IAMI_past", "IAMI_future"]:
    d[v + "_z"] = zc(v)
COV = "CTQ_total_z + LEC_direct_z + age_z + female"
rows = []
for y in ["PCL_intrusion", "PCL_avoidance", "PCL_negmood", "PCL_arousal", "PCL5"]:
    for refg in ["control", "an_imagers"]:
        s = d[d.g.isin(["congenital", refg])].dropna(subset=[y, "IAMI_past_z", "IAMI_future_z", "CTQ_total_z", "LEC_direct_z", "age_z"]).copy()
        s["grp"] = (s.g == "congenital").astype(float)
        s["yz"] = (s[y] - d.loc[d.g == refg, y].mean()) / d.loc[d.g == refg, y].std()
        m0 = smf.ols(f"yz ~ grp + {COV}", data=s).fit()
        m1 = smf.ols(f"yz ~ grp + IAMI_past_z + {COV}", data=s).fit()
        m2 = smf.ols(f"yz ~ grp * IAMI_past_z + {COV}", data=s).fit()
        rows.append(dict(outcome=y, vs=refg, n=int(m0.nobs), grp_without_IAMI=m0.params.grp, grp_with_IAMI=m1.params.grp,
                         grp_with_IAMI_lo=m1.conf_int().loc["grp"][0], grp_with_IAMI_hi=m1.conf_int().loc["grp"][1],
                         grp_with_IAMI_p=m1.pvalues.grp, IAMI_beta=m1.params.IAMI_past_z,
                         slope_ref=m2.params.IAMI_past_z, slope_aph=m2.params.IAMI_past_z + m2.params["grp:IAMI_past_z"],
                         slope_diff=m2.params["grp:IAMI_past_z"], slope_diff_lo=m2.conf_int().loc["grp:IAMI_past_z"][0],
                         slope_diff_hi=m2.conf_int().loc["grp:IAMI_past_z"][1], slope_diff_p=m2.pvalues["grp:IAMI_past_z"]))
I = pd.DataFrame(rows)
I.to_csv(f"{ROOT}/output/04_iami.csv", index=False, float_format="%.6f")
# the same for the odds of passing the PCL-5 cut-off
rows = []
for refg in ["control", "an_imagers"]:
    s = d[d.g.isin(["congenital", refg])].dropna(subset=["probable_ptsd", "IAMI_past_z", "CTQ_total_z", "LEC_direct_z", "PHQ9_z", "GAD7_z", "age_z"]).copy()
    s["grp"] = (s.g == "congenital").astype(float)
    for lab, f_ in [("exposure, mood", f"{COV} + PHQ9_z + GAD7_z"), ("+ IAMI past", f"{COV} + PHQ9_z + GAD7_z + IAMI_past_z")]:
        m = smf.logit(f"probable_ptsd ~ grp + {f_}", data=s).fit(disp=0)
        rows.append(dict(vs=refg, adjusted_for=lab, OR=np.exp(m.params.grp), lo=np.exp(m.conf_int().loc["grp"][0]),
                         hi=np.exp(m.conf_int().loc["grp"][1]), p=m.pvalues.grp, n=int(m.nobs)))
I2 = pd.DataFrame(rows)
L += ["\n## I. Involuntary past memories (IAMI) and the group difference in post-traumatic symptoms\n",
      "grp = congenital minus reference in reference-SD units at equal exposure, without and with the IAMI past-memory "
      "score as a covariate; slopes = IAMI past -> symptom within each group.\n",
      I.round(3).to_markdown(index=False), "\nOdds of a PCL-5 total >= 33:\n", I2.round(3).to_markdown(index=False)]

# ---------------- K. detection among the clinically distressed (proxy for impairment) ----------------
# The survey has no functioning item, so clinical-range depression or anxiety (PHQ-9 >= 10 or GAD-7 >= 10) stands in
# for "needs help". Among those people: how often does the PCL-5 cut-off, or the DSM-5 rule, flag them, by group?
# Clinical-range depression or anxiety is only a proxy for impairment; the survey had no functioning measure.
d["distressed"] = ((d.PHQ9 >= 10) | (d.GAD7 >= 10)).astype(float).where(d.PHQ9.notna() & d.GAD7.notna())
rows = []
for flag, lab in [("probable_ptsd", "PCL-5 total >= 33"), ("dsm_full", "DSM-5 rule")]:
    for refg in ["control", "an_imagers"]:
        s = d[d.g.isin(["congenital", refg]) & (d.distressed == 1)].dropna(subset=[flag, "CTQ_total_z", "LEC_direct_z", "age_z"]).copy()
        s["grp"] = (s.g == "congenital").astype(float)
        rate = s.groupby("g")[flag].mean(); nn = s.groupby("g").size()
        m = smf.logit(f"{flag} ~ grp + CTQ_total_z + LEC_direct_z + PHQ9_z + GAD7_z + age_z + female", data=s).fit(disp=0)
        rows.append(dict(flag=lab, vs=refg, n_aph=nn["congenital"], n_ref=nn[refg], rate_aph=rate["congenital"], rate_ref=rate[refg],
                         OR=np.exp(m.params.grp), lo=np.exp(m.conf_int().loc["grp"][0]), hi=np.exp(m.conf_int().loc["grp"][1]), p=m.pvalues.grp))
K = pd.DataFrame(rows)
K.to_csv(f"{ROOT}/output/04_detection_distressed.csv", index=False, float_format="%.6f")
L += ["\n## K. Detection among people with clinical-range depression or anxiety (PHQ-9 or GAD-7 >= 10)\n",
      "OR adjusted for exposure, PHQ-9, GAD-7, age, gender.\n", K.round(3).to_markdown(index=False)]

# ---------------- L. dissociation against distress: lower level, or weaker response? ----------------
# Does distress turn into dissociation less readily in aphantasia (flatter slope),
# or do aphantasics start lower and respond the same (equal slope, lower intercept)? Outcome: DES depersonalisation
# and derealisation (no biased items), standardised on controls; one model per distress measure.
d["DES_dep_z"] = zc("DES_dep")
rows = []
for x in ["PHQ9_z", "GAD7_z", "CTQ_total_z", "LEC_direct_z", "PCL5_z"]:
    for refg in ["control", "an_imagers"]:
        s = d[d.g.isin(["congenital", refg])].dropna(subset=["DES_dep_z", x, "age_z", "female"]).copy()
        s["grp"] = (s.g == "congenital").astype(float)
        m = smf.ols(f"DES_dep_z ~ grp * {x} + age_z + female", data=s).fit()
        ci = m.conf_int().loc[f"grp:{x}"]
        rows.append(dict(distress=x[:-2], vs=refg, n=int(m.nobs), slope_ref=m.params[x], slope_aph=m.params[x] + m.params[f"grp:{x}"],
                         slope_diff=m.params[f"grp:{x}"], lo=ci[0], hi=ci[1], p=m.pvalues[f"grp:{x}"],
                         level_diff_at_mean=m.params.grp, level_lo=m.conf_int().loc["grp"][0],
                         level_hi=m.conf_int().loc["grp"][1], level_p=m.pvalues.grp))
Ldf = pd.DataFrame(rows)

# ---------------- M. within people who have imagery: is weaker imagery associated with more depersonalisation? ----------------
# Lambert et al. (2001) found this in depersonalisation disorder. Spearman rho between Psi-Q vividness and DES
# depersonalisation/derealisation (no biased items) within each imagery group, raw and partial for depression and anxiety.
# Partial rho = Pearson correlation of the rank residuals after regressing out the ranks of PHQ-9 and GAD-7 (t test on
# n - 4 df). 95% CIs by Fisher's z with the Bonett-Wright standard error for Spearman correlations.
from scipy.stats import spearmanr  # noqa: E402


def rho_ci(r, n, k=0):
    se = np.sqrt((1 + r ** 2 / 2) / (n - 3 - k))
    return np.tanh(np.arctanh(r) - 1.96 * se), np.tanh(np.arctanh(r) + 1.96 * se)


rows = []
for grp in ["control", "an_imagers"]:
    s = d[d.g == grp].dropna(subset=["PSIQ", "DES_dep", "PHQ9", "GAD7"])
    n = len(s)
    rho, p = spearmanr(s.PSIQ, s.DES_dep)
    rk = s[["PSIQ", "DES_dep", "PHQ9", "GAD7"]].rank()
    res = lambda v: sm.OLS(rk[v], sm.add_constant(rk[["PHQ9", "GAD7"]])).fit().resid
    prho = float(np.corrcoef(res("PSIQ"), res("DES_dep"))[0, 1])
    pp = 2 * t_dist.sf(abs(prho) * np.sqrt((n - 4) / (1 - prho ** 2)), n - 4)
    lo, hi = rho_ci(rho, n)
    plo, phi_ = rho_ci(prho, n, k=2)
    rows.append(dict(group=grp, n=n, rho=rho, lo=lo, hi=hi, p=p, partial_rho_mood=prho, partial_lo=plo, partial_hi=phi_,
                     partial_p=pp))
Mdf = pd.DataFrame(rows)
Mdf.to_csv(f"{ROOT}/output/04_imagery_depersonalisation.csv", index=False, float_format="%.6f")

# ---------------- N. trauma amnesia beyond general autobiographical gaps ----------------
# Is the higher PCL-5 item 8 in aphantasia the general memory style (DES item 9: no memory of important life events)
# showing up for the trauma, or something specific to how trauma is remembered? Group effect on item 8 at equal
# exposure, without and with DES 9 as covariate; and the item 8 to depression/numbing association with DES 9 held fixed.
d["DES_Q9_z"] = zc("DES_Q9")
rows = []
for refg in ["control", "an_imagers"]:
    s = d[d.g.isin(["congenital", refg])].dropna(subset=["PCL_Q8", "DES_Q9_z", "CTQ_total_z", "LEC_direct_z", "age_z"]).copy()
    s["grp"] = (s.g == "congenital").astype(float)
    s["y"] = (s.PCL_Q8 - d.loc[d.g == refg, "PCL_Q8"].mean()) / d.loc[d.g == refg, "PCL_Q8"].std()
    for lab, extra in [("exposure", ""), ("exposure + DES 9", " + DES_Q9_z")]:
        m = smf.ols(f"y ~ grp + CTQ_total_z + LEC_direct_z + age_z + female{extra}", data=s).fit()
        rows.append(dict(vs=refg, adjusted_for=lab, n=int(m.nobs), grp=m.params.grp, lo=m.conf_int().loc["grp"][0],
                         hi=m.conf_int().loc["grp"][1], p=m.pvalues.grp,
                         des9_beta=m.params.get("DES_Q9_z", np.nan)))
N1 = pd.DataFrame(rows)
rows = []
for grp in ["control", "congenital"]:
    s = d[d.g == grp].dropna(subset=["PCL_Q8", "DES_Q9", "PHQ9", "PCL_numbing", "CTQ_total", "LEC_direct", "age", "female"])
    for y in ["PHQ9", "PCL_numbing"]:
        m = smf.ols(f"scale({y}) ~ scale(PCL_Q8) + scale(DES_Q9) + scale(CTQ_total) + scale(LEC_direct) + scale(age) + female",
                    data=s).fit()
        rows.append(dict(group=grp, outcome=y, n=int(m.nobs), item8_beta=m.params["scale(PCL_Q8)"],
                         lo=m.conf_int().loc["scale(PCL_Q8)"][0], hi=m.conf_int().loc["scale(PCL_Q8)"][1],
                         p=m.pvalues["scale(PCL_Q8)"], des9_beta=m.params["scale(DES_Q9)"],
                         des9_lo=m.conf_int().loc["scale(DES_Q9)"][0], des9_hi=m.conf_int().loc["scale(DES_Q9)"][1],
                         des9_p=m.pvalues["scale(DES_Q9)"]))
N2 = pd.DataFrame(rows)
# does the association of DES item 9 with depression differ between the groups? interaction in the pooled sample
s = d[d.g.isin(["control", "congenital"])].dropna(subset=["PCL_Q8", "DES_Q9", "PHQ9", "CTQ_total", "LEC_direct", "age", "female"]).copy()
s["grp"] = (s.g == "congenital").astype(float)
m = smf.ols("scale(PHQ9) ~ scale(PCL_Q8) + scale(DES_Q9) * grp + scale(CTQ_total) + scale(LEC_direct) + scale(age) + female",
            data=s).fit()
N3 = pd.DataFrame([dict(term="DES 9 x group (congenital vs control)", estimate=m.params["scale(DES_Q9):grp"],
                        lo=m.conf_int().loc["scale(DES_Q9):grp"][0], hi=m.conf_int().loc["scale(DES_Q9):grp"][1],
                        p=m.pvalues["scale(DES_Q9):grp"], n=int(m.nobs))])
N3.to_csv(f"{ROOT}/output/04_des9_interaction.csv", index=False, float_format="%.6f")
N1.to_csv(f"{ROOT}/output/04_amnesia_beyond_des9.csv", index=False, float_format="%.6f")
N2.to_csv(f"{ROOT}/output/04_amnesia_distress_des9.csv", index=False, float_format="%.6f")
Ldf.to_csv(f"{ROOT}/output/04_dissociation_slopes.csv", index=False, float_format="%.6f")
L += ["\n## L. Depersonalisation/derealisation (DES, unbiased items) on each distress measure, by group\n",
      "slope_diff = congenital minus reference; level_diff_at_mean = group difference at the mean of the distress measure.\n",
      Ldf.round(3).to_markdown(index=False),
      "\n## M. Imagery vividness (Psi-Q) and depersonalisation/derealisation within the imagery groups (Spearman)\n",
      Mdf.round(3).to_markdown(index=False),
      "\n## N. PCL-5 item 8 (trauma amnesia) with and without DES 9 (general autobiographical gaps) as covariate\n",
      "grp = congenital minus reference on item 8, reference-SD units.\n", N1.round(3).to_markdown(index=False),
      "\nItem 8 as predictor of depression and numbing, with DES 9 held fixed (standardised beta):\n",
      N2.round(3).to_markdown(index=False),
      "\nGroup difference in the DES 9 to depression association (interaction, pooled sample):\n",
      N3.round(3).to_markdown(index=False)]

pd.concat([A, B.assign(part="B_form"), C.assign(part="C_slopes"), Dm.assign(part="D_routes"), E.assign(part="E_gaps"),
           F.assign(part="F_detection"), H.assign(part="H_dsm")]).to_csv(f"{ROOT}/output/04_all.csv", index=False, float_format="%.6f")
Path(f"{ROOT}/output/04_trauma_analyses.md").write_text("\n".join(L), encoding="utf-8")
print("\n".join(L))
