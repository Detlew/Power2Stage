README
================

  - [Power2Stage](#power2stage)
      - [Supported Methods](#supported-methods)
          - [Simulation-based](#simulation-based)
              - [‘Type 1’](#type-1)
              - [‘Type 2’](#type-2)
              - [Blinded Sample Size Re-estimation in the
                Interim](#blinded-sample-size-re-estimation-in-the-interim)
              - [Group Sequential Design](#group-sequential-design)
          - [Inverse-Normal Combination](#inverse-normal-combination)
      - [Functions](#functions)
          - [Main](#main)
          - [Helpers](#helpers)
      - [Examples](#examples)
          - [Method B](#method-b)
          - [Method C](#method-c)
          - [Inverse-Normal Combination](#inverse-normal-combination-1)
      - [Speed Comparisons](#speed-comparisons)
      - [Installation](#installation)

[![cran
checks](https://cranchecks.info/badges/summary/PowerTOST)](https://cran.r-project.org/web/checks/check_results_Power2Stage.html)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/Power2Stage?color=blue)](https://r-pkg.org/pkg/Power2Stage)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/last-month/Power2Stage?color=green)](https://r-pkg.org/pkg/Power2Stage)

# Power2Stage

The package contains functions to obtain the operational characteristics
(power, type I error, percentage of studies proceeding to the second
stage, average and quantiles of total sample sizes) of bioequivalence
studies in adaptive sequential 2-Stage Designs (TSD) via simulations.  
Built 2019-09-22 with R 3.6.1.

## Supported Methods

### Simulation-based

Since the many letters denoting the methods given by various authors
might be confusing, [I
classified](https://doi.org/10.1007/s00228-015-1806-2) the methods as
two ‘types’:

  - **‘Type 1’**  
    An adjusted *α* is used *both* in the interim as well as in the
    final analysis of pooled data.
  - **‘Type 2’**  
    Whether an unadjusted or an adjusted *α* is used depends on interim
    power. An adjusted *α* is used in the final analysis of pooled data.

It should be noted that the adjusted alphas do not necessarily have to
be the same in both stages. Below a summary of conditions used in the
decision schemes of the published methods.

<small>[TOC ↩](#readme)</small>

#### ‘Type 1’

  - [Potvin *et al.*](https://doi.org/10.1002/pst.294) (2008) ‘Method
    B’: *α* 0.0294 (*θ*<sub>0</sub> 0.95, target power 0.80).
  - [Fuglsang](https://doi.org/10.1208/s12248-013-9475-5) (2013) ‘Method
    B’: *α* 0.0284 (*θ*<sub>0</sub> 0.95, target power 0.90).
  - [Karalis](https://doi.org/10.1016/j.ijpharm.2013.08.013) (2013)
    ‘TSD-2’: *α* 0.0294 (*θ*<sub>0</sub> = PE, target power 0.80).
  - [Fuglsang](https://doi.org/10.1208/s12248-014-9571-1) (2014) ‘Method
    B’ (parallel design): *α* 0.0294 (*θ*<sub>0</sub> 0.95, target power
    0.80).
  - [Zheng *et al.*](https://doi.org/10.1002/pst.1672) (2015) ‘MSDBE’:
    *α*<sub>1</sub> 0.01, *α*<sub>2</sub> 0.04.
  - [Xu *et al.*](https://doi.org/10.1002/pst.1721) (2016) ‘Method E’:
    (*θ*<sub>0</sub> 0.95, target power 0.80, *n<sub>max</sub>* 42).  
      - For *CV* 10–30%  
        *α*<sub>1</sub> 0.0294, *α*<sub>2</sub> 0.0357, futility rule on
        CI {0.9374, 1/0.9374}.
      - For *CV* 30–55%  
        *α*<sub>1</sub> 0.0254, *α*<sub>2</sub> 0.0363, futility rule on
        CI {0.9305, 1/0.9305}.
  - [Molins *et al.*](https://doi.org/10.1002/sim.7452) (2017) ‘Type 1
    modified Potvin B’: *α* 0.0301 (*θ*<sub>0</sub> 0.95, target power
    0.80, min. *n<sub>2</sub>* = 1.5*n<sub>1</sub>*, *n<sub>max</sub>*
    150).

#### ‘Type 2’

  - [Potvin *et al.*](https://doi.org/10.1002/pst.294) (2008) ‘Method
    C’: *α* 0.0294 (*θ*<sub>0</sub> 0.95, target power 0.80).
  - [Montague *et al.*](https://doi.org/10.1002/pst.483) (2011) ‘Method
    D’: *α* 0.0280 (*θ*<sub>0</sub> 0.90, target power 0.80).
  - [Fuglsang](https://doi.org/10.1208/s12248-013-9475-5) (2013) ‘Method
    C/D’:  
    *α* 0.0274 (*θ*<sub>0</sub> 0.95, target power 0.90).  
    *α* 0.0269 (*θ*<sub>0</sub> 0.90, target power 0.90).
  - [Karalis and Macheras](https://doi.org/10.1007/s11095-013-1026-3)
    (2013) ‘TSD’: *α* 0.0294 (*θ*<sub>0</sub> = PE, target power 0.80).
  - [Karalis](https://doi.org/10.1016/j.ijpharm.2013.08.013) (2013)
    ‘TSD-1’: *α* 0.0280 (*θ*<sub>0</sub> = PE, target power 0.80).
  - [Xu *et al.*](https://doi.org/10.1002/pst.1721) (2016) ‘Method F’:
    (*θ*<sub>0</sub> 0.95, target power 0.80, *n<sub>max</sub>* 180).  
      - For *CV* 10–30%  
        *α*<sub>1</sub> 0.0248, *α*<sub>2</sub> 0.0364, futility rule on
        CI {0.9492, 1/0.9492}.
      - For *CV* 30–55%  
        *α*<sub>1</sub> 0.0259, *α*<sub>2</sub> 0.0349, futility rule on
        CI {0.9350, 1/0.9350}.
  - [Molins *et al.*](https://doi.org/10.1002/sim.7452) (2017) ‘Type 2
    modified Potvin C’: *α* 0.0280 (*θ*<sub>0</sub> 0.95, target power
    0.80, min. *n<sub>2</sub>* = 1.5*n<sub>1</sub>*, *n<sub>max</sub>*
    150).

#### Blinded Sample Size Re-estimation in the Interim

[Golkowski *et al.*](https://doi.org/10.1002/pst.1617) (2014).

#### Group Sequential Design

[Kieser and Rauch](https://doi.org/10.1002/sim.6487) (2015).

### Inverse-Normal Combination

[König *et al.*](https://doi.org/10.13140/RG.2.1.5190.0967) (2014),
[Kieser and Rauch](https://doi.org/10.1002/sim.6487) (2015), [Wassmer
and Brannath](https://doi.org/10.1007/978-3-319-32562-0) (2016), [Maurer
*et al.*](https://doi.org/10.1002/sim.7614) (2018).

<small>[TOC ↩](#readme)</small>

## Functions

### Main

Defaults employed if not specified in the function call:

| function          | `theta0` | `target power` | `usePE` | `Nmax` | `max.n` | `fCrit` | `fClower` |
| ----------------- | :------: | :------------: | :-----: | :----: | :-----: | :-----: | :-------: |
| `power.tsd()`     |  `0.95`  |     `0.80`     | `FALSE` | `Inf`  |    –    |    –    |     –     |
| `power.tsd.fC()`  |  `0.95`  |     `0.80`     | `FALSE` |   –    |  `Inf`  | `"PE"`  |  `0.80`   |
| `power.tsd.KM()`  |  `0.95`  |     `0.80`     |    –    | `150`  |    –    |    –    |     –     |
| `power.tsd.ssr()` |  `0.95`  |     `0.80`     | `FALSE` |   –    |  `Inf`  |    –    |     –     |
| `power.tsd.GS()`  |  `0.95`  |       –        |    –    |   –    |    –    | `"PE"`  |  `0.80`   |
| `power.tsd.in()`  |  `0.95`  |     `0.80`     | `FALSE` |   –    |  `Inf`  | `"CI"`  |  `0.95`   |
| `power.tsd.p()`   |  `0.95`  |     `0.80`     | `FALSE` | `Inf`  |    –    |    –    |     –     |

All functions are for a 2×2×2 crossover design except `power.tsd.p()`,
which is for a two-group parallel design.  
If `usePE=TRUE` the point estimate in the interim is used in sample size
estimation of the second stage.  
If the estimated total sample size exceeds `max.n` the second stage will
be forced to `max.n-n1` (*i.e.*, it is *not* a futility criterion).  
The method used for interim power and sample size estimation is
specified by the argument `pmethod`. It defaults to `"nct"`
(approximation by the noncentral *t*-distribution) except in
`power.tsd.GS()`, where the total sample size is already fixed.  
The BE limits are specified by the arguments `theta1` and `theta2`
(default to 0.80 and 1.25).  
The number of simulations is specified by the argument `nsims`. It
defaults to 1e5 if simulating power and to 1e6 if simulating the empiric
type I error (*i.e.*, `theta0` set to the value of `theta1` or
`theta2`).

**Futility Criteria in the Interim**

  - `Nmax`: The study will stop if the estimated total sample size
    exceeds `Nmax`.
  - `fCrit` (`"PE"` or `"CI"`): The study will stop if outside `fClower`
    and `1/fClower`.
      - `"PE"`: `fClower` defaults to 0.80.
      - `"CI"`: `fClower` defaults to 0.925 (except in function
        `power.tsd.in()`, where it defaults to 0.95).

<small>[TOC ↩](#readme)</small>

### Helpers

  - `sampleN2.TOST()`  
    Estimates the sample size of stage 2 to achieve at least the target
    power.
  - `interim.tsd.in()`  
    Interim analysis based on the Inverse-Normal Combination method.
  - `final.tsd.in()`  
    Final analysis based on the Inverse-Normal Combination method.

<small>[TOC ↩](#readme)</small>

## Examples

Before running the examples attach the library.

``` r
library(Power2Stage)
```

If not noted otherwise, defaults are employed.

### Method B

Power estimation by the shifted central *t*-distribution.

``` r
power.tsd(CV = 0.20, n1 = 12, pmethod = "shifted")
# TSD with 2x2 crossover 
# Method B: alpha (s1/s2) = 0.0294 0.0294 
# Target power in power monitoring and sample size est. = 0.8
# Power calculation via shifted central t approx. 
# CV1 and GMR = 0.95 in sample size est. used
# No futility criterion
# BE acceptance range = 0.8 ... 1.25
# 
# CV = 0.2; n(stage 1) = 12; GMR = 0.95
# 
# 1e+05 sims at theta0 = 0.95 (p(BE) = 'power').
# p(BE)    = 0.84454
# p(BE) s1 = 0.41333
# Studies in stage 2 = 56.45%
# 
# Distribution of n(total)
# - mean (range) = 20.7 (12 ... 82)
# - percentiles
#  5% 50% 95% 
#  12  18  40
```

Explore the empiric type I error at the upper BE-limit.

``` r
power.tsd(CV = 0.20, n1 = 12, pmethod = "shifted",
          theta0 = 1.25)[["pBE"]]
# [1] 0.046352
```

<small>[TOC ↩](#readme)</small>

### Method C

Power estimation by the shifted central *t*-distribution.

``` r
power.tsd(method = "C", CV = 0.20, n1 = 12, pmethod = "shifted")
# TSD with 2x2 crossover 
# Method C: alpha0 = 0.05, alpha (s1/s2) = 0.0294 0.0294 
# Target power in power monitoring and sample size est. = 0.8
# Power calculation via shifted central t approx. 
# CV1 and GMR = 0.95 in sample size est. used
# No futility criterion
# BE acceptance range = 0.8 ... 1.25
# 
# CV = 0.2; n(stage 1) = 12; GMR = 0.95
# 
# 1e+05 sims at theta0 = 0.95 (p(BE) = 'power').
# p(BE)    = 0.8496
# p(BE) s1 = 0.42656
# Studies in stage 2 = 53.7%
# 
# Distribution of n(total)
# - mean (range) = 20.6 (12 ... 82)
# - percentiles
#  5% 50% 95% 
#  12  18  40
```

Slightly better than ‘Method B’ in terms of power in both stages and
fewer studies are expected to proceed to the second stage.

Explore the empiric type I error at the upper BE-limit (1 milion
simulations).

``` r
power.tsd(method = "C", CV = 0.20, n1 = 12, pmethod = "shifted",
          theta0 = 1.25)[["pBE"]]
# [1] 0.051238
```

Slight inflation of the type I error (although considered negligible by
the authors). However, more adjustment (adjusted *α* 0.0280) controls
the type I error.

``` r
power.tsd(method = "C", alpha = rep(0.0280, 2), CV = 0.20,
          n1 = 12, pmethod = "shifted", theta0 = 1.25)[["pBE"]]
# [1] 0.049903
```

<small>[TOC ↩](#readme)</small>

### Inverse-Normal Combination

Data given by Potvin *et al.* in Example 2: 12 subjects in stage 1, PE
1.0876, CV 0.18213, all defaults of the function used.

``` r
interim.tsd.in(GMR1 = 1.0876, CV1 = 0.18213, n1= 12)
# TSD with 2x2 crossover
# Inverse Normal approach
#  - Maximum combination test with weights for stage 1 = 0.5 0.25 
#  - Significance levels (s1/s2) = 0.02635 0.02635 
#  - Critical values (s1/s2) = 1.93741 1.93741 
#  - BE acceptance range = 0.8 ... 1.25
#  - Observed point estimate from stage 1 is not used for SSR
#  - With conditional error rates and conditional estimated target power
# 
# Interim analysis after first stage
# - Derived key statistics:
#   z1 = 3.10000, z2 = 1.70344,
#   Repeated CI = (0.92491, 1.27891)
#   Median unbiased estimate = NA
# - No futility criterion met
# - Test for BE not positive (not considering any futility rule)
# - Calculated n2 = 6
# - Decision: Continue to stage 2 with 6 subjects
```

The second stage should be initiated with 6 subjects. Note that with
`interim.tsd.in(..., fCrit="No", ssr.conditional="no")` 8 subjects would
be required like in the Methods of Potvin *et al.*

The second stage is performed in 8 subjects, PE 0.9141, CV 0.25618.

``` r
final.tsd.in(GMR1 = 1.0876, CV1 = 0.18213, n1 = 12,
             GMR2 = 0.9141, CV2 = 0.25618, n2 = 8)
# TSD with 2x2 crossover
# Inverse Normal approach
#  - Maximum combination test with weights for stage 1 = 0.5 0.25 
#  - Significance levels (s1/s2) = 0.02635 0.02635 
#  - Critical values (s1/s2) = 1.93741 1.93741 
#  - BE acceptance range = 0.8 ... 1.25
# 
# Final analysis after second stage
# - Derived key statistics:
#   z1 = 2.87952, z2 = 2.60501,
#   Repeated CI = (0.87690, 1.17356)
#   Median unbiased estimate = 1.0135
# - Decision: BE achieved
```

The study passed with a (repeated) CI of 87.69–117.36%. Although
slightly more conservative, same conclusion like based on the 94.12% CI
of 88.45–116.38% reported by Potvin *et al.*

<small>[TOC ↩](#readme)</small>

## Speed Comparisons

Performed on a double Xeon E3-1245v3 3.4 GHz, 8 MB cache, 16 GB RAM, R
3.6.1 64 bit on Windows 7.

‘Method B’ (*CV* 0.20, *n*<sub>1</sub> 12).

    #   method   power seconds
    #  shifted 0.84454    1.09
    #      nct 0.84266    1.61
    #    exact 0.84260   31.98

Despite being the fastest, the shifted central *t*-distribution should
only be used in order to compare with published methods. The noncentral
*t*-distribution is a good compromise between speed and accuracy and
hence, the default in all functions. The exact method based on Owen’s
Q-function is time-consuming and therefore, not recommended in
validating a custom method in a narrow grid of
*n*<sub>1</sub>/*CV*-combinations. However, in designing a new study it
is the method of choice.

Blinded sample size re-estimation (*α* 0.03505, *CV* 0.239,
*n*<sub>1</sub> 10, target power 0.90), 1 million simulations for the
empiric type I error.

    #   method      TIE seconds
    #       ls 0.049054    3.67
    #  shifted 0.046106   12.85
    #      nct 0.046319   18.24
    #    exact 0.046319  429.10

The crude large sample approximation (`pmethod="ls"`) should only be
used to compare with the published method.

<small>[TOC ↩](#readme)</small>

## Installation

You can install the released version of Power2Stage from
[CRAN](https://CRAN.R-project.org) with …

``` r
package <- "Power2Stage"
inst    <- package %in% installed.packages()
if (length(package[!inst]) > 0) install.packages(package[!inst])
```

… and the development version from [GitHub](https://github.com/) with …

    # install.packages("devtools")
    devtools::install_github("Detlew/Power2Stage")

<small>[TOC ↩](#readme)</small>
