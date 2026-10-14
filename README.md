# People with aphantasia report fewer PTSD intrusions but similar distress after trauma: analysis code

Chen Jielu Sunny Chan, Timo L. Kvamme and Juha Silvanto.

This repository holds the code that rebuilds every number, table and figure of the paper, and `expected/` holds the
published results, so a run can be compared with ours file by file.

**Data.** Participants were told that their responses would be reported only in the aggregate, so the participant-level
data are not public. The anonymised data file (`data/trauma_survey_anonymised.csv`) is available from the authors on
reasonable request. `data/codebook.csv` describes every column, and `data/README.md` describes how the file was
anonymised. Put the file in `data/` to run the analysis.

## Quick start

1. Install Python 3.12 and the packages in `requirements.txt`: `pip install -r requirements.txt`.
2. Install R 4.4 and run `Rscript install_r_packages.R` (lavaan, semTools, knitr). Only the measurement invariance models use R.
3. Run `python run_all.py`. It takes about 5 minutes plus 10 to 20 minutes for the invariance models.
   `python run_all.py --skip-invariance` skips R and reuses the invariance results already in `output/`, or copies
   them from `expected/` when there are none (the invariance tables then match `expected/` by construction).
4. The last step prints `REPRODUCTION OK` when every result table matches `expected/` to six significant digits.

Results go to `output/` (tables) and `figures/` (the paper's figures under the paper's names).

## Layout

| Folder or file | Contents |
| --- | --- |
| `data/` | the codebook and a description of the anonymised data file, which is available on request |
| `scripts/` | the analysis, run in order by `run_all.py` |
| `expected/` | the published results: every table produced by the scripts, as `.csv` (full precision) and `.md` (readable) |
| `figures/` | Figures 1 to 4 and Supplementary Figures S1 and S2 |
| `tools/` | how the shared data file was made from the raw exports, its identifier check, and the codebook generator |

## The scripts

| Step | Script | What it does | Main outputs |
| --- | --- | --- | --- |
| 1 | `01_build_dataset.py` | scores every questionnaire, counts LEC-5 events, defines the groups | `scored.csv`, `01_descriptives.md` |
| 2 | `02_factor_structure.py` | exploratory factor analysis of all 181 items; the structure in controls and in congenital aphantasia (Tucker congruence against random splits of the controls); item bias; pooling check | `02_factor_structure.md`, `02_loadings_pooled.csv`, `02_congruence*.csv`, `02_item_bias.csv`, `02_pooling_check.csv` |
| 3 | `03_invariance.R` | ordinal multi-group confirmatory factor models per questionnaire (configural, metric, scalar, partial scalar) against each reference group | `03_invariance_fit_*.csv`, `03_invariance_freed_*.csv`, `03_latent_means_*.csv` |
| 4 | `04_trauma_analyses.py` | every trauma analysis in the paper, parts A to N (listed at the top of the script; letters G and J are not used) | `04_trauma_analyses.md` and one `.csv` per part |
| 5 | `05_tables.py` | Tables 1 and 2, Supplementary Tables S1 and S2, and the mediation estimates | `05_tables.md`, `05_table_*.csv` |
| 6 | `06_figure_factor_structure.py` | Supplementary Figure S1 | `figures/Supplementary_Figure_S1.png` |
| 7 | `07_figure_invariance.py` | Supplementary Figure S2 | `figures/Supplementary_Figure_S2.png` |
| 8 | `08_figures_trauma.py` | Figures 1 to 4 | `figures/Figure_1.png` to `Figure_4.png` |
| 9 | `09_check_reproduction.py` | compares every result table with `expected/` | prints `REPRODUCTION OK` |

## Where each result in the paper comes from

| Paper section | Output |
| --- | --- |
| Participant characteristics (table) and internal consistency (supplementary table) | `05_table_participants.csv`, `05_table_alpha.csv` |
| Trauma exposure, depression and anxiety | `04_all.csv` (part `A_exposure`), `04_symptom_form.csv` |
| Measurement equivalence of the questionnaires | `02_factor_structure.md`, `02_congruence*.csv`, `02_pooling_check.csv`, `03_invariance_fit_*.csv`, `02_item_bias.csv`, `05_table_biased_items.csv` |
| Intrusion and avoidance, trauma amnesia, core subgroup | `04_symptom_form.csv` |
| Dissociation and involuntary memory | `04_symptom_form.csv`, `04_dissociation_slopes.csv`, `04_imagery_depersonalisation.csv` |
| Trauma and distress: slopes and equivalence tests | `04_slopes.csv` (columns `lo90`, `hi90`, `p_tost`, `equivalent`) |
| Indirect associations of childhood trauma with depression | `04_routes.csv`, `05_table_mediation.csv` |
| Trauma amnesia, depression and numbing | `04_all.csv` (part `E_gaps`), `04_amnesia_beyond_des9.csv`, `04_amnesia_distress_des9.csv` |
| Screening cut-off and DSM-5 criteria | `04_detection.csv`, `04_dsm_criteria.csv`, `05_table_dsm_criteria.csv`, `04_detection_distressed.csv` |
| Involuntary memories and the group difference in intrusions | `04_iami.csv`, `04_trauma_analyses.md` (part I) |

## Methods in brief

Groups (`scripts/01_build_dataset.py`): Prolific controls have a mean Psi-Q above 4 on its 0 to 10 scale. Congenital
aphantasia is Aphantasia Network members who report never having had imagery and have a mean Psi-Q of 4 or below; the
core subgroup rated every Psi-Q visual item at 0. Network members with a mean Psi-Q above 4 form the second reference
group. Group differences in symptoms are regression coefficients at equal childhood trauma (CTQ) and traumatic events
(LEC-5), adjusted for age and gender, in standard deviations of the reference group of each comparison. Slopes of
distress on exposure and indirect effects are in standard deviations of the controls. Slope differences are tested for
equivalence with two one-sided tests against bounds of plus and minus 0.20 SD.

Random seeds are fixed in every script that resamples (parallel analysis, random splits, bootstrap), so a rerun gives
identical numbers.

## Licence and citation

Code: MIT (`LICENSE`). The data are available from the authors on reasonable request. Please cite the paper and this
code (`CITATION.cff`): https://doi.org/10.5281/zenodo.23375587
