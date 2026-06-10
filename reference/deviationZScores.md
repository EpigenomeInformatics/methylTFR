# deviationZScores

Function to get deviation Z-scores from a methylTFRdeviations object.

## Usage

``` r
deviationZScores(x)
```

## Arguments

- x:

  A methylTFRdeviations object.

## Value

A matrix of deviation Z-scores.

## Examples

``` r
# Load example data
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
# Get deviation Z-scores
deviationZScores(tc_mem)
#>             Tc-Mem_OP_S5_Long_D1.bedGraph.bed
#> FOXF2                              0.09952745
#> FOXD1                             -3.37320266
#> IRF2                              -8.81278614
#> MZF1(var.2)                       -7.93557376
#> MAX::MYC                          -6.92619776
#> PPARG                             -8.53704072
#> PAX6                               0.15511514
#> PBX1                              -0.44512022
#> RORA                             -12.11118430
#> RORA(var.2)                       -8.16113084
#>             Tc-Mem_OP_S4_Long_D60.bedGraph.bed
#> FOXF2                               0.05099083
#> FOXD1                              -2.15313860
#> IRF2                              -10.60528416
#> MZF1(var.2)                        -4.63134222
#> MAX::MYC                           -6.83283764
#> PPARG                              -4.52570166
#> PAX6                               -1.98421035
#> PBX1                               -1.55884533
#> RORA                               -7.95318382
#> RORA(var.2)                        -4.79944817
#>             Tc-Mem_OP_S4_Long_D28.bedGraph.bed
#> FOXF2                               -0.2966509
#> FOXD1                               -1.0020442
#> IRF2                                -9.5352471
#> MZF1(var.2)                         -4.2448393
#> MAX::MYC                            -6.9823143
#> PPARG                               -8.1898404
#> PAX6                                -1.9610534
#> PBX1                                -0.2753343
#> RORA                               -16.2778446
#> RORA(var.2)                         -3.1748887
#>             Tc-Mem_OP_S4_Long_D1.bedGraph.bed Tc-Mem_OP_S3_High_D1.bedGraph.bed
#> FOXF2                              -0.1072407                     -6.671444e-04
#> FOXD1                              -2.5315756                     -3.059927e+00
#> IRF2                               -7.1715554                     -1.072176e+01
#> MZF1(var.2)                        -5.4315341                     -4.469167e+00
#> MAX::MYC                           -6.1736131                     -8.222368e+00
#> PPARG                              -6.7031345                     -9.327438e+00
#> PAX6                               -0.2066197                      1.935868e-01
#> PBX1                               -1.3246348                     -2.053792e+00
#> RORA                              -11.4897174                     -2.082799e+01
#> RORA(var.2)                        -3.8143368                     -2.919634e+00
```
