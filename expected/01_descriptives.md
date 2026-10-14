# 01 - dataset build

Rows: {'an': 1038, 'prolific': 495}, total 1533

## Groups (Psi-Q on the 0-10 metric, cut-off 4)

| grp_ext         |   n |   psiq |   age |   phq |   gad |   pcl |
|:----------------|----:|-------:|------:|------:|------:|------:|
| acquired        |  67 |   1.29 | 50.66 |  9.19 |  6.37 | 22.24 |
| an_low_nochange | 152 |   1.18 | 52.25 |  6.62 |  5.13 | 19.12 |
| an_low_other    |  48 |   1.23 | 51.38 |  6.66 |  4.22 | 19.22 |
| an_with_imagery | 261 |   7.14 | 48.3  |  8.15 |  6.78 | 25.8  |
| cong_above      | 198 |   1.25 | 51.12 |  6.96 |  5.31 | 19.84 |
| cong_floor      | 142 |   0    | 52.85 |  6.62 |  4.73 | 16.46 |
| control         | 471 |   7.79 | 48.81 |  6.49 |  5.61 | 23.42 |
| no_psiq         | 181 | nan    | 51.42 |  5.65 |  4.82 | 18.11 |
| prolific_low    |  13 |   2.75 | 41.62 |  8.15 |  5.77 | 24.31 |


## Complete questionnaire blocks by source

| src      |   ok_MAIA |   ok_DES |   ok_CTQ |   ok_PHQ |   ok_GAD |   ok_PCL |   ok_IAMI |   ok_PSIQ |
|:---------|----------:|---------:|---------:|---------:|---------:|---------:|----------:|----------:|
| an       |       867 |      878 |      886 |     1028 |      851 |      872 |       859 |       868 |
| prolific |       481 |      487 |      489 |      491 |      487 |      488 |       486 |       484 |