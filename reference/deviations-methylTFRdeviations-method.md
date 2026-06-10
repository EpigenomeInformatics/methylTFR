# deviations

Extract bias corrected deviations from methylTFRdeviations object.

## Usage

``` r
# S4 method for class 'methylTFRdeviations'
deviations(x)
```

## Arguments

- x:

  methylTFRdeviations object.

## Value

A matrix of bias corrected deviations.

## Examples

``` r
# Load the data
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
# Get deviations
deviations(tc_mem)
#>             Tc-Mem_OP_S5_Long_D1.bedGraph.bed
#> FOXF2                             0.002128691
#> FOXD1                            -0.068104043
#> IRF2                             -0.099855724
#> MZF1(var.2)                      -0.091823810
#> MAX::MYC                         -0.207980740
#> PPARG                            -0.044047525
#> PAX6                              0.002634395
#> PBX1                             -0.008565336
#> RORA                             -0.058758695
#> RORA(var.2)                      -0.078493539
#>             Tc-Mem_OP_S4_Long_D60.bedGraph.bed
#> FOXF2                              0.001343341
#> FOXD1                             -0.042951266
#> IRF2                              -0.122406382
#> MZF1(var.2)                       -0.068451047
#> MAX::MYC                          -0.224920075
#> PPARG                             -0.041315791
#> PAX6                              -0.031466427
#> PBX1                              -0.028415769
#> RORA                              -0.040940878
#> RORA(var.2)                       -0.058993445
#>             Tc-Mem_OP_S4_Long_D28.bedGraph.bed
#> FOXF2                             -0.007352279
#> FOXD1                             -0.025370314
#> IRF2                              -0.114251930
#> MZF1(var.2)                       -0.070261983
#> MAX::MYC                          -0.231861827
#> PPARG                             -0.078707801
#> PAX6                              -0.032462078
#> PBX1                              -0.008190409
#> RORA                              -0.099863662
#> RORA(var.2)                       -0.043749070
#>             Tc-Mem_OP_S4_Long_D1.bedGraph.bed Tc-Mem_OP_S3_High_D1.bedGraph.bed
#> FOXF2                            -0.002867150                     -1.509199e-05
#> FOXD1                            -0.053054377                     -6.417042e-02
#> IRF2                             -0.090927407                     -1.067922e-01
#> MZF1(var.2)                      -0.077510986                     -6.353562e-02
#> MAX::MYC                         -0.201681226                     -2.157996e-01
#> PPARG                            -0.060603258                     -6.165151e-02
#> PAX6                             -0.003925531                      2.580385e-03
#> PBX1                             -0.024302408                     -3.844766e-02
#> RORA                             -0.071378314                     -8.644321e-02
#> RORA(var.2)                      -0.040281824                     -3.292322e-02
```
