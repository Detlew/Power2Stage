# Power2Stage 0.5.2.9000
Published on GitHub 2019-09-14.

## Minor changes

  * NEWS.md instead of NEWS.
  * Stop `interim.2stage.in()` if n1 <3 and final.tsd.in() if n1 | n2 <3 (avoids errors and wranings).
  * References in man-pages of `*.tsd.in()` corr.

# Power2Stage 0.5.2 (easter egg)
On CRAN 2019-04-21.

Maintenance release with mainly bug fixes.

## Bug fixes

  * Calculation of median unbiased estimate in `final.tsd.in()`.
  * `power.tsd.GS()` regarding variable `hw`.
  * `interim.tsd.in()` w.r.t to power > target triggering futulity fixed (THX to mittyri).

## Minor changes

  * Typos and minor wording and formatting corrections (man pages and `*.tds.in()`).

# Power2Stage 0.5.1 (belated easter egg)
On CRAN 2018-04-03.

## Major changes

  * New functions `interim.2stage.in()` and `final.2stage.in()` to perform interim and final analysis of 2 stage designs based on the Standard Combination or Maximum Combination Test (also includes a print class).

## Minor changes

  * All functions `xyz.2stage.ab()` now have an alias `xyz.tsd.ab()`, e.g. insteadof `power.2stage.fC()` one may use `power.tsd.fC()`. The versions with `.2stage.` in their names will be removed in later versions.
  * The futility rule for a CI in `power.2stage.GS()` is now based on the 90% CI to be in line with the other functions using such a futility rule.
  * Default of argument `fCrit` in `power.2stage.GS()` changed to `"CI"`.
  * Updated defaults of `power.2stage.in()` in order to be consistent with paper of Maurer et al (2018).
  * Aliases `power.tsd.xyz` introduced (not public yet).
  * Experimental function `power.tsd.2m()` to deal with 2 PK metrics in the TSD decision scheme Potvin "B" (not public yet).

# Power2Stage 0.4-6
On CRAN 2017-10-27.

## Major changes

  * New function `sampleN2.TOST()` to estimate `n2` in 2x2 and parallel designs based on `df=n-3` in stage 2.
  * df for sample size re-estimation changed to `n-3` for Potvin like TSD schemes (pooled evaluation with additionally stage term).
  * As above for `power.2Stage.ssr()`.
  * Function `power.2stage.in()` to calculate the operational characteristics of TSDs evaluated via p-value combination using the inverse normal approach. Contributed by Benjamin Lang.

## Minor changes

  * GitHub URL and bug reports URL added.

# Power2Stage 0.4-5
On CRAN 2017-01-18.

Maintenance release to reflect changes in power calculation of package PowerTOST.

## Bug fixes

  * Bug in `power.2stage.p()` and `power.2stage.pAF()` removed which caused negative residual variances in case of `test="anova"` and imprecise power values in all other cases for certain scenarios of n1, CV.

## Minor changes

  * DOI of references in man pages added.
  * Large sample approximation for the sample size adapted to deal with asymmetric BE acceptance ranges, one limit = Inf allowed.
  * Small cosmetic changes in code w.r.t min.n2

# Power2Stage 0.4-4
Internal only 2016-01-21 (not released via CRAN).

## Major changes

  * Method `"B0"` (re)introduced in `power.2stage()` and `power.2stage.fC()` which employs the decision scheme of the so-called MSDBE (Modified Sequential Design for BE) of Zheng et al.
  * Method `"B"` in `power.2stage.p()` and `power.2stage.pAF()` also adapted todecision scheme E of Xu et al.

## Minor changes

  * Figures with TSD decision schemes added in `/doc` subdirectory.

# Power2Stage 0.4-3
On CRAN 2015-11-24.

## Major changes

  * `power2stage.fC()` has a new argument `max.n` which constraints the total sample size.
  * Method `"B"` in `power2stage.fC()` and `power.2stage()` redefined to deal with unequal `alpha[1]` and `alpha[2]` in the same manner as method E in Xu et al.This gives small differences compared to previous calculations with unequal alphas, f.i. Haybittle-Peto alphas.
  * The CI futility criterion is now based on the 90% CI according to Xu et al. It was formerly the `1-2*alpha[1]` CI.
  * Default futility criterion changed to `"CI"` with a lower bound 0.925.
  * Defunct `print` argument removed from all `power.2stage.xyz()` functions.

# Power2Stage 0.4-2
On CRAN 2015-07-13.

## Bug fixes

  * Bug removed in `power.2stage.ssr()` which prevented correct sample-size re-estimation if power was calculated via `pmethod="ls"`.

## Major changes

  * New argument `usePE` in `power.2stage.ssr()` to use the point estimate from the interim analysis in sample-size re-estimation (makes only sense if `blind=FALSE`).

# Power2Stage 0.4-1
On CRAN 2015-06-19.

## Major changes

  * Major rewrite of the `power.2stage.xyz()` functions wich now return an S3 object of class `"pwrtsd"`. Output is now done more R-like via the S3 print method for class `"pwrtsd"`. Therefore, the `print` argument is defunc and will be removed in the next version.
  * The `power.2stage.xyz()` functions (except `power.2stage.GS()`) return a component `ntable`, a 'table' object for summarizing the discrete distribution of `ntotal`. The `nhist` component was removed.
  * The default of `fClower` in case of `fCrit="PE"` changed to 0.8 in function `power.2stage.fC()`.

## Minor changes

  * The argument `detail` has now the default `FALSE`. Thus be patient if you simulate for alpha with 1 mio sim's.
  * The `nsims` argument, if missing, is now set to 1E6 if you simulate for alpha (i.e. with `theta0` at border or outside acceptance range `theta1 ... theta2`) and to 1E5 otherwise.

# Power2Stage 0.3-2
Released to alpha testers only 2015-06-10.

## Major changes

  * The `power.2stage.yy()` functions (except `power.2stage.GS()`) now return a component `nhist` with class `"histogram"` which can be used with `plot()` to visualize the distribution of `Ntotal`. Suggested by H. Schuetz.
  * Handling of non-integer degrees of freedom in `power.2stage.p()` in case of Welch's test changed (mo more truncation).
  * Power calculation method `"shifted"` implemented in `power.2stage.fC()`.

## Minor changes

  * Slight improvements and typo fixes in man pages.

# Power2Stage 0.3-1
On CRAN 2015-01-24 (dedicated to my brother Stefan's 50 birthday).

## Major changes

  * Revision to reflect internal changes in upcoming PowerTOST v1.2-06 'raw' power functions.
  * Sample-size estimation routine (vectorized!) reworked which gives a considerable run-time boost.
  * Internal change of start value of sample size estimation to Zhang's formula (may change the extremal value of n and its mean if `usePE=TRUE`).
  * `power.2stage.p()` accepts now unbalanced stage 1 and uses the correct power of Welch's test in the power monitoring steps and in sample size estimation.
  * Former function `power.2stage.p()` now available as `power.2stage.pAF()`  which performs the calculations exactly as described in Fuglsang's paper.

# Power2Stage 0.2-2
On CRAN 2014-12-08.

## Bug fixes

  * Two nasty bugs removed.

# Power2Stage 0.2-1
Released to alpha-testers only 2014-12-06 (dedicated to my daughter Antje's 40 birthday).

## Major changes

  * Deprecated function `power.2stage.Bf()` removed.
  * Functions `power.2stage()` and `power.2stage.fC()` have now an argument `min.n2` to restrict the sample size of stage 2 to a lower limit.
  * Calculation of % studies in stage 2 re-defined to studies having total n>n1 (was in the past studies which had to be evaluated with alpha2; affects only asymmetric alpha settings).
  * Function `power.2stage.ssr()` for (blinded) sample size re-estimation added (2-stage design without BE decision at interim).

## Minor changes

  * Code streamlining, enhancements and unification of output of power functions.

# Power2Stage 0.1-5
On CRAN 2014-10-09.

## Major changes

  * `power.2stage.Bf()` reworked into the new function `power.2stage.fC()` to include Potvin method "C" and to include an additional futility criterion for the CI as well. `power.2stage.Bf()` is deprecated and will not be removed in one of the next versions. Use `power.2stage.fC()` instead.

# Power2Stage 0.1-4
On CRAN 2014-07-24.

## Bug fixes

  * Bug in `power.2stage.GS()` removed which prevented the convergence of power values with increasing number of sim's.

# Power2Stage 0.1-3
On CRAN 2014-07-02.

## Bug fixes

  * Internal change of start value of sample size search to avoid failed searches if variability is high and `theta0` close to 1.

## Major changes

  * Power calculation method `"shifted"` added in `power.2stage()` to the provide comparison with Potvin et al.

# Power2Stage 0.1-2
Not released to the public (special version to alpha testers)

# Power2Stage 0.1-1
On CRAN 2014-05-08.

Released to alpha testers 2014-04-24.

## Major changes

  * New function `power.2stage.p()` for calculation of power in BE studies with sequential (2-stage) designs in 2 parallel groups acc. to Fuglsang.

## Minor changes

  * Internal code streamlining.

# Power2Stage 0.0-8
On CRAN 2014-04-11.

## Bug fixes

  * Bug in `power.2stage.KM()` removed if using `method="B"` (Karalis TSD-2). Thanks to Helmut Schuetz!

## Minor changes

  * Man page of `power.2stage.KM()` changed to the correct definition of the Karalis/Macheras TSD or Karalis TSD-1, TSD-2 respectively.
# Power2Stage 0.0-7
On CRAN 2014-02-13.
Maintenance release to reflect changes in PowerTOST V1.1-10.

## Major changes

  * New function `power.2stage.GS()` for non-adaptive group sequential (2-stage) BE studies.

# Power2Stage 0.0-6
On CRAN 2014-01-02.

## Minor changes

  * Examples adapted to complain with CRAN policy "Examples should run for no more than a few seconds each". Few seconds means below 5 sec as I learned. Attention! Number of sim's are too low to get meaningful results. Minimum number of sim's should be 1E5 for 'power', 1E6 for 'alpha'.

# Power2Stage 0.0-5
First attempt to release via CRAN 2013-12-27.

# Power2Stage 0.0-4
Not released, internal version only 2013-09-12.

## Major changes

  * Default method in `power.2stage.KM()` changed to `"C"`.
  * Arguments in `power.2stage.KM()` removed which are not necessary for the TSDs described in the Karalis & Macheras and in Karalis papers (see references).
  * `power2.2stage()` renamed to `power.2stage.KM()`.
  * New function `power.2stage.Bf()` implemented which evaluates a 2-stage design derived from Potvin method B with a futility criterion for the PE of stage 1.

# Power2Stage 0.0-3
Released to alpha testers 2013-09-03.

## Major changes

  * New function `power2.2stage()` implemented which uses PE and MSE of stage 1 also for the power calculation steps if `usePE=TRUE`.

# Power2Stage 0.0-2
Released to alpha testers 2013-07-10.

## Bug fixes

  * Bug in output removed stating that the sample size is estimated with PE & mse from stage 1 even if argument `usePE=FALSE`. THX to Helmut.

## Major changes

  * Default for `npct` argument changed to `c(0.05, 0.5, 0.95)` in accordance with the Potvin papers.
  * Cases with n(total)>Nmax now counted at stage 1 since I feel this is more logical.
  * Internally stage vector introduced and output based on that.

## Minor changes

  * More checks of input added. THX to Helmut Schuetz.

# Power2Stage 0.0-1
First release to alpha testers 2013-07-09.
