# Data

## Requesting the data

The participant-level data are not public: participants were told that their responses would be reported only in the
aggregate. The anonymised data file is available on reasonable request from the corresponding authors:

- Timo L. Kvamme, timo@cfin.au.dk
- Juha Silvanto, juhasilvanto@um.edu.mo

Please say briefly what you want to use the data for. Put the file you receive in this folder as
`trauma_survey_anonymised.csv`, then run `python run_all.py` from the repository root.

## The data file

`trauma_survey_anonymised.csv` has one row per participant (1,533 rows: 495 recruited through Prolific, 1,038 through
the Aphantasia Network) and 291 columns. `codebook.csv` describes every column.

## What the file contains

| Columns | Content |
| --- | --- |
| `pid` | study code: `pr` (Prolific) or `an` (Aphantasia Network) plus a running number |
| `src` | recruitment source, `prolific` or `an` |
| `age`, `gender_code` | age in years; 1 female, 2 male, 3 another gender |
| `imagery_change` | lifetime change in mental imagery (1 no change, 2 decreased, 3 increased, 4 never had imagery, 5 other) |
| `PSIQ_*`, `DES_*`, `MAIA_*`, `IAMI_*`, `CTQ_*`, `PHQ_*`, `GAD_*`, `PCL_*`, `LEC_*` | item responses in the survey's numeric codes (see `codebook.csv`; `scripts/01_build_dataset.py` rescales and scores them) |

## How it was anonymised

The raw survey exports are not shared. `tools/anonymise.py` builds this file from them with an explicit allowlist of
columns, so everything else is left out: IP addresses, e-mail fields, response IDs, timestamps, the time taken to
complete the survey, country and region, all free-text answers (occupation, the description of the worst event, when
and why imagery changed), and a dream questionnaire that the paper does not use. Respondents who reported an age
under 18 were removed in a pre-cleaning step, so the file holds adults only. Three stray text entries in item columns
became missing values, which is how the analysis treats them anyway. The rows keep the order of the export, so the analysis reproduces the reported
results exactly; run `python tools/check_data.py` to check the file for identifying content.

Item wording is not included, because several questionnaires are copyrighted. Item numbers follow the published
versions: Psi-Q (Andrade et al., 2014), 28-item DES (Bernstein & Putnam, 1986; Carlson & Putnam, 1993), MAIA-2 (Mehling et al., 2018), IAMI
(Berntsen et al., 2015), CTQ short form (Bernstein et al., 2003), PHQ-9 (Kroenke et al., 2001), GAD-7 (Spitzer et
al., 2006), PCL-5 (Weathers, Litz, et al., 2013) and LEC-5 (Weathers, Blake, et al., 2013).

## Groups

Groups are defined in `scripts/01_build_dataset.py` from the source, the mean Psi-Q score (0 to 10) and the answer on
lifetime imagery change:

| Group | Rule | n |
| --- | --- | --- |
| Prolific controls | Prolific, mean Psi-Q above 4 | 471 |
| Congenital aphantasia | Aphantasia Network, never had imagery, mean Psi-Q 4 or below | 340 |
| (core subgroup) | of these, every Psi-Q visual item at 0 | 266 |
| Network members with imagery | Aphantasia Network, mean Psi-Q above 4 | 261 |
| not analysed | Prolific with mean Psi-Q 4 or below (13); no Psi-Q data (181); Network members with mean Psi-Q 4 or below who did not report lifelong absence of imagery or did not answer the question (267) | 461 |
