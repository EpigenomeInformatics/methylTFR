# deviationZScores

Extract deviation Z-scores from methylTFRdeviations object.

## Usage

``` r
# S4 method for class 'methylTFRdeviations'
deviationZScores(x)
```

## Arguments

- x:

  methylTFRdeviations object.

## Value

A matrix of deviation Z-scores.

## Examples

``` r
# Load the data
load(system.file("extdata", "tc_mem.rda", package = "methylTFR"))
# Get deviation Z-scores
deviationZScores(tc_mem)
#>              Tc-Mem_01  Tc-Mem_02  Tc-Mem_03   Tc-Mem_04   Tc-Mem_05
#> FOXF2        1.0583028  0.9668552 -1.2630023 -0.20849376 -0.55366194
#> FOXD1       -0.4069572  1.4598517 -1.2034865  0.42422298 -0.27363103
#> IRF2         0.2548902 -1.5752563 -0.3040243  0.95684258  0.66754780
#> MZF1(var.2)  0.6955303 -0.3860781  1.3858082 -0.77776242 -0.91749802
#> MAX::MYC     0.1218820 -0.3340680  1.3863916  0.20727822 -1.38148377
#> PPARG        0.1407703  0.9391524 -0.2081849 -1.59234460  0.72060678
#> PAX6         1.7187724 -0.1670604 -0.1387236 -0.69235173 -0.72063664
#> PBX1         0.4381136 -1.1895353 -0.6219780 -0.04237231  1.41577195
#> RORA         1.3467369 -0.9808149  0.5332867 -0.96710750  0.06789881
#> RORA(var.2) -0.4332507  1.4127288 -1.2038808 -0.30131047  0.52571313
```
