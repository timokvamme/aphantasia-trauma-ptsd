# 03 - measurement invariance (ordinal multi-group CFA, WLSMV, Wu & Estabrook 2016 identification)

Criterion per step: dCFI >= -.010 and dRMSEA <= .015.

## Fit



|contrast        |instrument | n_ref| n_cong|model          | chisq.scaled| df.scaled| cfi.scaled| rmsea.scaled|  srmr|  d_cfi| d_rmsea|
|:---------------|:----------|-----:|------:|:--------------|------------:|---------:|----------:|------------:|-----:|------:|-------:|
|an_with_imagery |DES        |   244|    303|configural     |     1419.199|       700|      0.924|        0.061| 0.085|     NA|      NA|
|an_with_imagery |DES        |   244|    303|thresholds     |     1485.389|       753|      0.923|        0.060| 0.085| -0.001|  -0.002|
|an_with_imagery |DES        |   244|    303|metric         |     1506.968|       780|      0.924|        0.058| 0.085|  0.001|  -0.001|
|an_with_imagery |DES        |   244|    303|scalar         |     1929.231|       807|      0.882|        0.071| 0.094| -0.042|   0.013|
|an_with_imagery |DES        |   244|    303|partial_scalar |     1575.235|       798|      0.918|        0.060| 0.093| -0.005|   0.001|
|an_with_imagery |MAIA       |   237|    300|configural     |     3066.231|      1202|      0.932|        0.076| 0.080|     NA|      NA|
|an_with_imagery |MAIA       |   237|    300|thresholds     |     3170.096|      1280|      0.931|        0.074| 0.080| -0.001|  -0.002|
|an_with_imagery |MAIA       |   237|    300|metric         |     3159.563|      1309|      0.932|        0.073| 0.080|  0.001|  -0.002|
|an_with_imagery |MAIA       |   237|    300|scalar         |     3207.015|      1338|      0.931|        0.072| 0.080| -0.001|   0.000|
|an_with_imagery |IAMI       |   238|    295|configural     |     1569.434|       338|      0.958|        0.117| 0.068|     NA|      NA|
|an_with_imagery |IAMI       |   238|    295|thresholds     |     1623.342|       376|      0.957|        0.112| 0.068| -0.001|  -0.005|
|an_with_imagery |IAMI       |   238|    295|metric         |     1595.416|       394|      0.959|        0.107| 0.068|  0.002|  -0.005|
|an_with_imagery |IAMI       |   238|    295|scalar         |     1565.920|       412|      0.960|        0.103| 0.068|  0.002|  -0.004|
|an_with_imagery |CTQ        |   245|    302|configural     |     1018.300|       530|      0.987|        0.058| 0.068|     NA|      NA|
|an_with_imagery |CTQ        |   245|    302|thresholds     |     1068.276|       574|      0.987|        0.056| 0.068|  0.000|  -0.002|
|an_with_imagery |CTQ        |   245|    302|metric         |     1062.938|       594|      0.987|        0.054| 0.069|  0.001|  -0.002|
|an_with_imagery |CTQ        |   245|    302|scalar         |     1051.007|       612|      0.988|        0.051| 0.068|  0.001|  -0.003|
|an_with_imagery |PHQ        |   259|    340|configural     |      244.508|        54|      0.969|        0.109| 0.069|     NA|      NA|
|an_with_imagery |PHQ        |   259|    340|thresholds     |      257.902|        62|      0.968|        0.103| 0.069| -0.001|  -0.006|
|an_with_imagery |PHQ        |   259|    340|metric         |      247.633|        70|      0.971|        0.092| 0.069|  0.003|  -0.011|
|an_with_imagery |PHQ        |   259|    340|scalar         |      225.224|        78|      0.976|        0.080| 0.069|  0.005|  -0.013|
|an_with_imagery |GAD        |   235|    299|configural     |      139.903|        28|      0.992|        0.123| 0.055|     NA|      NA|
|an_with_imagery |GAD        |   235|    299|thresholds     |      151.652|        35|      0.991|        0.112| 0.055|  0.000|  -0.011|
|an_with_imagery |GAD        |   235|    299|metric         |      138.128|        41|      0.993|        0.094| 0.055|  0.001|  -0.018|
|an_with_imagery |GAD        |   235|    299|scalar         |      139.713|        47|      0.993|        0.086| 0.055|  0.000|  -0.008|
|an_with_imagery |PCL        |   238|    307|configural     |      875.556|       328|      0.963|        0.078| 0.063|     NA|      NA|
|an_with_imagery |PCL        |   238|    307|thresholds     |      926.592|       365|      0.962|        0.075| 0.063| -0.001|  -0.003|
|an_with_imagery |PCL        |   238|    307|metric         |      912.309|       381|      0.964|        0.072| 0.064|  0.002|  -0.004|
|an_with_imagery |PCL        |   238|    307|scalar         |      944.019|       397|      0.963|        0.071| 0.065| -0.001|   0.000|


## Items whose thresholds and intercept were freed



|contrast        |instrument | n_items| n_freed|freed                  |still_fails |
|:---------------|:----------|-------:|-------:|:----------------------|:-----------|
|an_with_imagery |DES        |      28|       3|DES_Q14 DES_Q9 DES_Q18 |FALSE       |
|an_with_imagery |MAIA       |      37|       0|                       |FALSE       |
|an_with_imagery |IAMI       |      20|       0|                       |FALSE       |
|an_with_imagery |CTQ        |      25|       0|                       |FALSE       |
|an_with_imagery |PHQ        |       9|       0|                       |FALSE       |
|an_with_imagery |GAD        |       7|       0|                       |FALSE       |
|an_with_imagery |PCL        |      20|       0|                       |FALSE       |


## Latent mean of congenital aphantasia (SD units of the reference group)



|contrast        |instrument |factor      | d_scalar| d_partial|
|:---------------|:----------|:-----------|--------:|---------:|
|an_with_imagery |DES        |dis         |    -0.57|     -0.53|
|an_with_imagery |MAIA       |noticing    |    -0.41|     -0.41|
|an_with_imagery |MAIA       |notdistract |     0.01|      0.01|
|an_with_imagery |MAIA       |notworry    |     0.32|      0.32|
|an_with_imagery |MAIA       |attention   |    -0.42|     -0.42|
|an_with_imagery |MAIA       |emotional   |    -0.50|     -0.50|
|an_with_imagery |MAIA       |selfreg     |    -0.41|     -0.41|
|an_with_imagery |MAIA       |listening   |    -0.72|     -0.72|
|an_with_imagery |MAIA       |trusting    |    -0.08|     -0.08|
|an_with_imagery |IAMI       |future      |    -1.09|     -1.09|
|an_with_imagery |IAMI       |past        |    -1.12|     -1.12|
|an_with_imagery |CTQ        |emoabuse    |    -0.28|     -0.28|
|an_with_imagery |CTQ        |phyabuse    |    -0.18|     -0.18|
|an_with_imagery |CTQ        |sexabuse    |    -0.24|     -0.24|
|an_with_imagery |CTQ        |emoneglect  |    -0.03|     -0.03|
|an_with_imagery |CTQ        |phyneglect  |    -0.18|     -0.18|
|an_with_imagery |PHQ        |dep         |    -0.25|     -0.25|
|an_with_imagery |GAD        |anx         |    -0.41|     -0.41|
|an_with_imagery |PCL        |intrusion   |    -0.74|     -0.74|
|an_with_imagery |PCL        |avoidance   |    -0.56|     -0.56|
|an_with_imagery |PCL        |negmood     |    -0.40|     -0.40|
|an_with_imagery |PCL        |arousal     |    -0.42|     -0.42|
