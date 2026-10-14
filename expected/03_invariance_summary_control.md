# 03 - measurement invariance (ordinal multi-group CFA, WLSMV, Wu & Estabrook 2016 identification)

Criterion per step: dCFI >= -.010 and dRMSEA <= .015.

## Fit



|contrast |instrument | n_ref| n_cong|model          | chisq.scaled| df.scaled| cfi.scaled| rmsea.scaled|  srmr|  d_cfi| d_rmsea|
|:--------|:----------|-----:|------:|:--------------|------------:|---------:|----------:|------------:|-----:|------:|-------:|
|control  |DES        |   466|    303|configural     |     1732.525|       700|      0.971|        0.062| 0.063|     NA|      NA|
|control  |DES        |   466|    303|thresholds     |     1810.930|       753|      0.970|        0.061| 0.063| -0.001|  -0.001|
|control  |DES        |   466|    303|metric         |     1813.464|       780|      0.971|        0.059| 0.063|  0.001|  -0.002|
|control  |DES        |   466|    303|scalar         |     2821.730|       807|      0.944|        0.081| 0.070| -0.027|   0.022|
|control  |DES        |   466|    303|partial_scalar |     2173.314|       801|      0.962|        0.067| 0.067| -0.009|   0.008|
|control  |MAIA       |   464|    300|configural     |     4279.451|      1202|      0.925|        0.082| 0.083|     NA|      NA|
|control  |MAIA       |   464|    300|thresholds     |     4421.533|      1280|      0.924|        0.080| 0.083| -0.002|  -0.002|
|control  |MAIA       |   464|    300|metric         |     4442.826|      1309|      0.924|        0.079| 0.083|  0.000|  -0.001|
|control  |MAIA       |   464|    300|scalar         |     4558.387|      1338|      0.922|        0.079| 0.083| -0.002|   0.000|
|control  |IAMI       |   466|    295|configural     |     1913.994|       338|      0.949|        0.111| 0.064|     NA|      NA|
|control  |IAMI       |   466|    295|thresholds     |     2018.241|       378|      0.946|        0.107| 0.064| -0.002|  -0.004|
|control  |IAMI       |   466|    295|metric         |     1944.798|       396|      0.949|        0.102| 0.064|  0.003|  -0.005|
|control  |IAMI       |   466|    295|scalar         |     1939.408|       414|      0.950|        0.099| 0.064|  0.001|  -0.003|
|control  |CTQ        |   469|    302|configural     |     1264.594|       530|      0.990|        0.060| 0.066|     NA|      NA|
|control  |CTQ        |   469|    302|thresholds     |     1324.288|       574|      0.990|        0.058| 0.066|  0.000|  -0.002|
|control  |CTQ        |   469|    302|metric         |     1312.667|       594|      0.990|        0.056| 0.067|  0.000|  -0.002|
|control  |CTQ        |   469|    302|scalar         |     1336.064|       612|      0.990|        0.055| 0.066|  0.000|  -0.001|
|control  |PHQ        |   469|    340|configural     |      264.259|        54|      0.982|        0.098| 0.057|     NA|      NA|
|control  |PHQ        |   469|    340|thresholds     |      285.526|        62|      0.981|        0.095| 0.057| -0.001|  -0.004|
|control  |PHQ        |   469|    340|metric         |      276.082|        70|      0.983|        0.085| 0.057|  0.001|  -0.009|
|control  |PHQ        |   469|    340|scalar         |      290.086|        78|      0.982|        0.082| 0.057| -0.001|  -0.003|
|control  |GAD        |   467|    299|configural     |      162.685|        28|      0.995|        0.112| 0.041|     NA|      NA|
|control  |GAD        |   467|    299|thresholds     |      178.986|        35|      0.995|        0.104| 0.041|  0.000|  -0.008|
|control  |GAD        |   467|    299|metric         |      168.273|        41|      0.995|        0.090| 0.041|  0.001|  -0.014|
|control  |GAD        |   467|    299|scalar         |      174.954|        47|      0.995|        0.084| 0.040|  0.000|  -0.006|
|control  |PCL        |   467|    307|configural     |     1220.207|       328|      0.967|        0.084| 0.059|     NA|      NA|
|control  |PCL        |   467|    307|thresholds     |     1281.165|       365|      0.966|        0.081| 0.059| -0.001|  -0.003|
|control  |PCL        |   467|    307|metric         |     1261.155|       381|      0.967|        0.077| 0.059|  0.001|  -0.003|
|control  |PCL        |   467|    307|scalar         |     1374.775|       397|      0.963|        0.080| 0.061| -0.004|   0.003|


## Items whose thresholds and intercept were freed



|contrast |instrument | n_items| n_freed|freed          |still_fails |
|:--------|:----------|-------:|-------:|:--------------|:-----------|
|control  |DES        |      28|       2|DES_Q14 DES_Q9 |FALSE       |
|control  |MAIA       |      37|       0|               |FALSE       |
|control  |IAMI       |      20|       0|               |FALSE       |
|control  |CTQ        |      25|       0|               |FALSE       |
|control  |PHQ        |       9|       0|               |FALSE       |
|control  |GAD        |       7|       0|               |FALSE       |
|control  |PCL        |      20|       0|               |FALSE       |


## Latent mean of congenital aphantasia (SD units of the reference group)



|contrast |instrument |factor      | d_scalar| d_partial|
|:--------|:----------|:-----------|--------:|---------:|
|control  |DES        |dis         |    -0.27|     -0.23|
|control  |MAIA       |noticing    |    -0.52|     -0.52|
|control  |MAIA       |notdistract |    -0.22|     -0.22|
|control  |MAIA       |notworry    |     0.55|      0.55|
|control  |MAIA       |attention   |    -0.62|     -0.62|
|control  |MAIA       |emotional   |    -0.70|     -0.70|
|control  |MAIA       |selfreg     |    -0.62|     -0.62|
|control  |MAIA       |listening   |    -1.02|     -1.02|
|control  |MAIA       |trusting    |    -0.49|     -0.49|
|control  |IAMI       |future      |    -1.31|     -1.31|
|control  |IAMI       |past        |    -1.07|     -1.07|
|control  |CTQ        |emoabuse    |     0.13|      0.13|
|control  |CTQ        |phyabuse    |    -0.27|     -0.27|
|control  |CTQ        |sexabuse    |     0.10|      0.10|
|control  |CTQ        |emoneglect  |     0.42|      0.42|
|control  |CTQ        |phyneglect  |     0.04|      0.04|
|control  |PHQ        |dep         |     0.05|      0.05|
|control  |GAD        |anx         |    -0.18|     -0.18|
|control  |PCL        |intrusion   |    -0.75|     -0.75|
|control  |PCL        |avoidance   |    -0.83|     -0.83|
|control  |PCL        |negmood     |    -0.15|     -0.15|
|control  |PCL        |arousal     |    -0.02|     -0.02|
