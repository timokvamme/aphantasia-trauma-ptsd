"""
Step 2: factor structure of all questionnaire items, and items answered differently in aphantasia.

  A. pooled exploratory factor analysis of all 181 items (Spearman correlations via a rank transform, minimum
     residual extraction, oblimin rotation, number of factors by Horn's parallel analysis)
  D. the structure of the 146 items other than the Psi-Q, fitted separately in controls and in congenital aphantasia,
     compared with Tucker's congruence against 50 random splits of the controls
  F. item bias: each item regressed on the mean of the other items of its factor, group, age and gender, against
     each reference group
  G. pooling check: the pooled analysis repeated after removing every group difference in item means

Input : output/scored.csv (script 01)
Output: output/02_*.csv, output/02_factor_structure.md
Usage : python scripts/02_factor_structure.py [--k N]   (default: parallel analysis picks k)
"""
import argparse
import re
import warnings

import numpy as np
import pandas as pd
from factor_analyzer import FactorAnalyzer
from scipy.linalg import orthogonal_procrustes

warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")   # factor_analyzer uses a renamed sklearn argument
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1].as_posix()
RNG = np.random.default_rng(20260930)

ap = argparse.ArgumentParser()
ap.add_argument("--k", type=int, default=0)
args = ap.parse_args()

d = pd.read_csv(f"{ROOT}/output/scored.csv", low_memory=False)
PREFIX = {"PSIQ": r"^PSIQ_[A-G]_Q\d+$", "DES": r"^DES_Q\d+$", "MAIA": r"^MAIA_\d+$", "IAMI": r"^IAMI_Q\d+$",
          "CTQ": r"^CTQ_[A-E]_Q\d+$", "PHQ": r"^PHQ_Q\d+$", "GAD": r"^GAD_Q\d+$", "PCL": r"^PCL_Q\d+$"}
inst = {c: k for k, p in PREFIX.items() for c in d.columns if re.match(p, c)}
items = list(inst)
IMG_FREE = [c for c in items if inst[c] != "PSIQ"]          # pool without the imagery questionnaire

d["grp"] = d.grp_ext.replace({"cong_floor": "congenital", "cong_above": "congenital"})
cc = d.dropna(subset=items).reset_index(drop=True)          # complete cases on every item
L = [f"# 02 - factor structure of the questionnaire items\n\nItems: {len(items)} "
     f"({', '.join(f'{k} {sum(v == k for v in inst.values())}' for k in PREFIX)}). "
     f"Complete cases: {len(cc)} of {len(d)}.\n", cc.grp_ext.value_counts().to_frame("n").to_markdown()]


def ranks(df):
    """Rank transform then z-score: Pearson on this is Spearman on the raw items."""
    r = df.rank()
    return (r - r.mean()) / r.std(ddof=0)


def parallel(X, reps=30, q=95):
    """Horn's parallel analysis on the correlation matrix eigenvalues."""
    ev = np.linalg.eigvalsh(np.corrcoef(X.T))[::-1]
    n, p = X.shape
    rnd = np.array([np.linalg.eigvalsh(np.corrcoef(RNG.standard_normal((n, p)).T))[::-1] for _ in range(reps)])
    return ev, int(np.argmax(ev <= np.percentile(rnd, q, axis=0)))


def efa(X, k, cols):
    fa = FactorAnalyzer(n_factors=k, method="minres", rotation="oblimin")
    fa.fit(X)
    return fa, pd.DataFrame(fa.loadings_, index=cols, columns=[f"F{i + 1}" for i in range(k)])


def label(load):
    """Name each factor by the instrument that supplies most of its salient (|loading| >= .40) items."""
    out = {}
    for f in load.columns:
        sal = load.index[load[f].abs() >= .40]
        top = pd.Series([inst[i] for i in sal]).value_counts()
        out[f] = f"{f}:" + "+".join(f"{k}{v}" for k, v in top.head(3).items()) if len(top) else f"{f}:none"
    return out


def cohen_d(a, b):
    a, b = a.dropna(), b.dropna()
    s = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(a) + len(b) - 2))
    dd = (a.mean() - b.mean()) / s
    se = np.sqrt((len(a) + len(b)) / (len(a) * len(b)) + dd ** 2 / (2 * (len(a) + len(b))))
    return dd, dd - 1.96 * se, dd + 1.96 * se


# ---------------- A. pooled EFA, all items ----------------
X = ranks(cc[items]).values
ev, k_pa = parallel(X)
k = args.k or k_pa
fa, load = efa(X, k, items)
names = label(load)
load_named = load.rename(columns=names)
load_named.insert(0, "instrument", [inst[i] for i in items])
load_named.to_csv(f"{ROOT}/output/02_loadings_pooled.csv", float_format="%.6f")
pd.Series(ev[:40], name="eigenvalue").to_csv(f"{ROOT}/output/02_eigenvalues.csv")
phi = pd.DataFrame(fa.phi_, index=names.values(), columns=names.values())
phi.to_csv(f"{ROOT}/output/02_factor_correlations.csv", float_format="%.6f")
L += [f"\n## A. Pooled EFA\n\nParallel analysis suggests {k_pa} factors; fitted k = {k}. "
      f"First eigenvalues: {np.round(ev[:k + 3], 2).tolist()}. Variance explained: {fa.get_factor_variance()[2][-1]:.2f}.\n",
      "Salient items (|loading| >= .40) per factor, by instrument:\n",
      pd.crosstab(pd.Series({i: inst[i] for i in items}, name="instrument"),
                  load.abs().idxmax(axis=1).map(names).where(load.abs().max(axis=1) >= .40, "none")
                  .rename("factor")).to_markdown(),
      "\nFactor correlations:\n", phi.round(2).to_markdown()]

# the imagery factor = the one carrying most PSIQ items
prim = load.abs().idxmax(axis=1)
f_img = prim[[i for i in items if inst[i] == "PSIQ"]].value_counts()
img_factors = list(f_img.index)                             # every factor with a PSIQ-primary item

ctrl, cong = cc[cc.grp == "control"], cc[cc.grp == "congenital"]

# ---------------- D. same structure in controls and congenital aphantasia? ----------------
# PSIQ is at floor in aphantasia, so compare on the imagery-free pool. Tucker phi after Procrustes rotation.
Xc, Xa = ranks(ctrl[IMG_FREE]), ranks(cong[IMG_FREE])
ok = [c for c in IMG_FREE if Xc[c].notna().all() and Xa[c].notna().all()]      # drop zero-variance items
kd = min(parallel(Xc[ok].values)[1], parallel(Xa[ok].values)[1])
lc, la = efa(Xc[ok].values, kd, ok)[1], efa(Xa[ok].values, kd, ok)[1]
R, _ = orthogonal_procrustes(la.values, lc.values)
la_rot = pd.DataFrame(la.values @ R, index=ok, columns=lc.columns)
tucker = lambda a, b: float((a * b).sum() / np.sqrt((a ** 2).sum() * (b ** 2).sum()))
cong_f = pd.Series({label(lc)[f]: tucker(lc[f].values, la_rot[f].values) for f in lc.columns}, name="tucker_phi")
# null for phi: split the controls 50 times into a random subsample of the aphantasia group's size (a) and the
# remainder (b). The null compares a with b. Because b is smaller than the full control group, the same split also
# compares the aphantasia group with b, so that both comparisons have the same group sizes (size-matched).
mean_phi = lambda l_ref, l_new: np.mean([tucker(l_ref.values[:, j], (l_new.values @ orthogonal_procrustes(
    l_new.values, l_ref.values)[0])[:, j]) for j in range(kd)])
null, matched = [], []
for _ in range(50):
    idx = RNG.permutation(len(ctrl))
    a_, b_ = ctrl.iloc[idx[:len(cong)]], ctrl.iloc[idx[len(cong):]]
    l1, l2 = efa(ranks(a_[ok]).fillna(0).values, kd, ok)[1], efa(ranks(b_[ok]).fillna(0).values, kd, ok)[1]
    null.append(mean_phi(l1, l2))
    matched.append(mean_phi(la, l2))
cong_f.to_csv(f"{ROOT}/output/02_congruence.csv", float_format="%.6f")
pd.DataFrame({"mean_phi": null, "mean_phi_aphantasia_vs_same_controls": matched}).to_csv(
    f"{ROOT}/output/02_congruence_null.csv", index=False, float_format="%.4f")
L += [f"\n## D. Structure in controls (n = {len(ctrl)}) vs congenital aphantasia (n = {len(cong)})\n\n"
      f"{len(ok)} imagery-free items, k = {kd}. Mean Tucker phi = {cong_f.mean():.3f}. "
      f"Control random-split null (n = {len(cong)} vs {len(ctrl) - len(cong)}): mean {np.mean(null):.3f}, "
      f"5th percentile {np.percentile(null, 5):.3f}. Size-matched (aphantasia vs the same {len(ctrl) - len(cong)} controls): "
      f"mean {np.mean(matched):.3f}; difference from the null {np.mean(np.array(matched) - np.array(null)):+.3f} "
      f"(5th to 95th percentile {np.percentile(np.array(matched) - np.array(null), 5):+.3f} to "
      f"{np.percentile(np.array(matched) - np.array(null), 95):+.3f}).\n", cong_f.round(3).to_frame().to_markdown()]

# ---------------- F. item bias scan: same latent level, different item score? ----------------
# For every imagery-free item: item ~ rest score of its own factor + group + age + gender, in rank-z units, among
# participants with complete items who reported being women or men (gender is a binary covariate here).
# Two contrasts. (1) congenital vs Prolific controls: imagery and recruitment source are confounded.
# (2) congenital vs Aphantasia Network members WITH imagery: same source, so bias here is imagery.
# Linear item-bias regression on rank-z scores, with two purification passes of the anchor (rest) items.
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

fam = prim[IMG_FREE]


def dif_scan(ref):
    s = cc[cc.grp.isin([ref, "congenital"]) & cc.age.notna() & cc.gender.isin(["female", "male"])]
    Zs = ranks(s[IMG_FREE])
    g = (s.grp == "congenital").astype(float).values
    cov = np.column_stack([g, (s.age - s.age.mean()) / s.age.std(), (s.gender == "female").astype(float)])
    anchors = {f: [i for i in IMG_FREE if fam[i] == f] for f in fam.unique()}
    for _ in range(3):                                       # pass 1 all items, then 2 purification passes
        out = {}
        for i in IMG_FREE:
            rest = [a for a in anchors[fam[i]] if a != i]
            if len(rest) < 3 or Zs[i].isna().any():
                continue
            m = sm.OLS(Zs[i].values, sm.add_constant(np.column_stack([Zs[rest].mean(axis=1), cov]))).fit()
            out[i] = (m.params[2], m.bse[2], m.pvalues[2])
        anchors = {f: [i for i in v if i not in out or abs(out[i][0]) < .25] or v for f, v in anchors.items()}
    r_ = pd.DataFrame(out, index=["dif", "se", "p"]).T
    r_["q"] = multipletests(r_.p, method="fdr_bh")[1]
    # latent difference: all items of the factor vs the unbiased anchors only
    lat = {names[f]: (cohen_d(Zs.loc[g == 1, [i for i in IMG_FREE if fam[i] == f]].mean(axis=1),
                              Zs.loc[g == 0, [i for i in IMG_FREE if fam[i] == f]].mean(axis=1))[0],
                      cohen_d(Zs.loc[g == 1, v].mean(axis=1), Zs.loc[g == 0, v].mean(axis=1))[0], len(v),
                      sum(fam == f)) for f, v in anchors.items()}
    return r_, pd.DataFrame(lat, index=["d_all_items", "d_anchors_only", "n_anchors", "n_items"]).T, len(s)


F1, lat1, n1 = dif_scan("control")
F2, lat2, n2 = dif_scan("an_with_imagery")
F = F1.add_suffix("_vs_prolific").join(F2.add_suffix("_vs_an_imagers"))
F.insert(0, "factor", fam.map(names))
F.insert(0, "instrument", [inst[i] for i in F.index])
F["biased_both"] = (F.q_vs_prolific < .05) & (F.q_vs_an_imagers < .05) & \
                   (np.sign(F.dif_vs_prolific) == np.sign(F.dif_vs_an_imagers)) & \
                   (F[["dif_vs_prolific", "dif_vs_an_imagers"]].abs().min(axis=1) >= .25)
F.to_csv(f"{ROOT}/output/02_item_bias.csv", float_format="%.6f")
lat = lat1.add_suffix("_vs_prolific").join(lat2.add_suffix("_vs_an_imagers"))
lat.to_csv(f"{ROOT}/output/02_latent_difference_anchors.csv", float_format="%.6f")
L += [f"\n## F. Item bias (congenital vs Prolific controls, n = {n1}; congenital vs AN imagers, n = {n2})\n\n"
      f"dif = item difference in SD units at the same rest score, age and gender adjusted; negative = aphantasia lower.\n"
      f"r(dif across the two contrasts, {len(F)} items) = {F.dif_vs_prolific.corr(F.dif_vs_an_imagers):.2f}.\n"
      f"Items biased in both contrasts (BH q < .05, same sign, |dif| >= .25): {int(F.biased_both.sum())}.\n",
      F.groupby("instrument").agg(items=("dif_vs_prolific", "size"), biased_both=("biased_both", "sum"),
                                  mean_abs_dif_prolific=("dif_vs_prolific", lambda x: x.abs().mean()),
                                  mean_abs_dif_an=("dif_vs_an_imagers", lambda x: x.abs().mean())).round(2).to_markdown(),
      "\nBiased in both contrasts:\n",
      F[F.biased_both].sort_values("dif_vs_prolific")[["instrument", "dif_vs_prolific", "dif_vs_an_imagers"]]
      .round(2).to_markdown(),
      "\nGroup difference on each factor, all items vs unbiased anchors only:\n", lat.round(2).to_markdown()]

# ---------------- G. pooling check: are the factors made by group mean differences? ----------------
# Rank-standardise every item WITHIN each of the 8 groups of script 01 (complete cases), which removes every
# between-group difference in item means (and in item spread), then refit with the pooled k.
Xw = cc.groupby("grp_ext")[items].transform(lambda s: (s.rank() - s.rank().mean()) / (s.rank().std(ddof=0) or 1)).fillna(0)
kw = parallel(Xw.values)[1]
lw = efa(Xw.values, k, items)[1]
Rw, _ = orthogonal_procrustes(lw.values, load.values)
lw_rot = pd.DataFrame(lw.values @ Rw, index=items, columns=load.columns)
phi_w = pd.Series({names[f_]: tucker(load[f_].values, lw_rot[f_].values) for f_ in load.columns}, name="tucker_phi")
phi_w.to_csv(f"{ROOT}/output/02_pooling_check.csv", float_format="%.6f")
L += [f"\n## G. Pooling check\n\nItems rank-standardised within each of the {cc.grp_ext.nunique()} groups, so no group mean "
      f"difference is left. Parallel analysis: {kw} factors (pooled: {k_pa}). Tucker phi between the pooled and the "
      f"within-group loadings, k = {k}: mean {phi_w.mean():.3f}, minimum {phi_w.min():.3f}.\n"
      f"DES 14 loading on the imagery factor: pooled {load.loc['DES_Q14', img_factors[0]]:.2f}, "
      f"within-group {lw_rot.loc['DES_Q14', img_factors[0]]:.2f}.\n", phi_w.round(3).to_frame().to_markdown()]

# self-check: the rank transform must reproduce Spearman, and the imagery questionnaire must find its own factor
assert np.allclose(np.corrcoef(X[:, :3].T), cc[items[:3]].corr(method="spearman").values, atol=1e-6)
assert f_img.iloc[0] >= 20, "PSIQ items did not form a factor"

Path(f"{ROOT}/output/02_factor_structure.md").write_text("\n".join(L), encoding="utf-8")
print("\n".join(L))
