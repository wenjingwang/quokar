---
title: "GSoC 2018 project by Wenjing"
---

This project improves the functionality for diagnosing outliers in quantile regression in the R package `quokar`. First, the implementation of the diagnosing algorithm is optimized in Baysian quantile regression framework, resulting in a significant speed gain. Second, an elemental sets diagnosing method is added, optimizing the leverage diagnostic in previous package. Third, an regression depth estimator is added, providing new diagnosinig method based on robust regression. This project page gives a summary of the work done during Google Summer of Code 2018.

This project is done under guidance of Dianne Cook and Kris Boudt. 

Installing the most current master branch of `quokar` available on Github can be done follows


```r
require(devtools)
install_github("wenjingwang/quokar")
```

## Part 1

The first part of this project can be found in [this](https://github.com/wenjingwang/quokar/pull/3) pull request. Here the function `bayesKL` has been rewrite to improve the speed. The previous version of this function calculate the pairwise comparison of Kullback–Leibler divergence for all observations in quantile regression model, which results in the computational complexity of O(2^(n - 1) * n). In this new version of function `bayesKL`, we improve the pairwise comparison algrithom by first traverse all observations then implement comparison which simplify the computational complexity to O(n).

Furthermore, we revise function `ALDqr_GCD` and `ALDqr_QD` to adapt to the maximum likelihood estimator of quantile regression under asymmetric Laplace distribution (Benites et al, 2016). The calculation of maximum likelihood estimator of quantile regression under asymmetric Laplace distribution is implemented by EM algorithm. These two functions are renamed as `GCDqr` and `QDqr` to diagnose outliers in regression model on each quantile based on the mle estimator or maximum likelihood function. 


In addition, we revise function `bayesKL` and `bayesProb` to adapt to the Bayesian quantile regression estimator under asymmetric Laplace distribution (Santos \& Bolfarine, 2016). The Bayesian quantile regression estimator involves MCMC algorithm which add additional input parameter `ndraw` to function `bayesKL` and `bayesProb`. These two functions are restructured to implement outlier diagnosing based on Baysian quantile regression estimator by comparing the mean probability or Kullback–Leibler divergence of the posterior distribution of the latent variable which is affected by each observation under bayesian framework. 

The data examples of the functions mentioned above can be downloaded here.

For now, this part still has to be merged with the master branch, to load the code in order to use the `GCDqr`, `QDqr`, `bayesKL` and `bayesProb`, do:


```r
require("githubinstall")
gh_install_packages("wenjingwang/quokar", ref = "improve_speed")
```

## Part 2

Next, we implement Elemental Sets method for quantile regression leverage diagnostic based on Rangnai & Vuuren (2014). We propose to use multiple-case high leverage diagnositc for quantile regressions using elemental regressions which is based only on the minimum number of observations p to estimate the parameters of the given model.

We calculate the leverage statistics for each observation by using the relationship among quantile regression, elemental regression and OLS. To split the observations into elemental sets, we use the L1 simplex algorithm for estimating the L1-quantile estimator of  quantile regression. Based on the elemental sets and leverage statsitics from OLS regression, we calculate the leverage statistics for QR model.

This Elmental Sets method performs better than robust distance method for diagnosing leverage points in previous `frame_distance` function in `quokar`. Robust distance method in `frame_distance` function only offer diagnostic for the whole data set while not for regressions on each quantile.

This part of work can be found in [this](https://github.com/wenjingwang/quokar/pull/4) pull request. For now, this part still has to be merged with the master branch, to load the code in order to use the `high_leverage`, do:


```r
require("githubinstall")
gh_install_packages("wenjingwang/quokar", ref = "high_leverage_qr")
```

More data example can be found [here](https://github.com/wenjingwang/gsoc-R/blob/master/R/First_phrase_report.R).

## Part 3

Finally, the third stage of this project implements the regression depth estimator. An estimator based on Rousseeuw & Hubert (1999) and Debruyne et al. (2008). The L1-quantile estimator used in linear quantile regression model is resistant to vertical outliers (observations that are outlying in y given x). However, it can be heavily influenced by leverage points (observations outlying in x-space). Slightest amount of contamination can have a disastrous effect on the resulting estimates. 

L1-quantiles take less information in the residuals of the model into consideration. This estimator only depends on the sign of the residuals and not on the exact value of the response variable which may result in the unrobustness. The quantile regression depth estimator provides robust estimation. Based on the robust estimator of quantile regression model, we can easily detect observations outlying from others by comparing the observed value and fitted value. We provide `plot` method for `qrdepth` class object to visualize outliers in quantile regression based on depth estimator. 

To illustrate the different performance of L1-quantile estimator and quantile regression depth estimator, we provide the following data example:



<img src="final_website_files/figure-html/unnamed-chunk-5-1.png" width="40%" style="display: block; margin: auto;" />

More data example can be found [here](https://github.com/wenjingwang/gsoc-R/blob/master/R/Second_phrase_report.R).

This part of work can be found in [this](https://github.com/wenjingwang/quokar/pull/5) pull request. For now, this part still has to be merged with the master branch, to load the code in order to use the `qrdepth`, do:


```r
require("githubinstall")
gh_install_packages("wenjingwang/quokar", ref = "qrdepth")
```

## Further work

The aim of this project is to improve the functionality of outlier diagnostic methods within R package `quokar`. This project includes computating speed improvements of the previous functions in `quokar`, also integrates the implementation of new diagnosing methods.

I will continue to update the package when new diagnosing methods for quantile regression are proposed in the literature. 

## References

Santos B, Bolfarine H. On Bayesian quantile regression and outliers[J]. arXiv preprint arXiv:1601.07344, 2016.

Benites L E, Lachos V H, Vilca F E. Case-deletion diagnostics for Quantile regression using the asymmetric Laplace distribution[J]. arXiv preprint arXiv:1509.05099, 2015.

Ranganai E, Van Vuuren J O, De Wet T. Multiple case high leverage diagnosis in regression quantiles[J]. Communications in Statistics-Theory and Methods, 2014, 43(16): 3343-3370.

Debruyne M, Hubert M, Portnoy S, et al. Censored depth quantiles[J]. Computational statistics & data analysis, 2008, 52(3): 1604-1614.

Rousseeuw P J, Hubert M. Regression depth[J]. Journal of the American Statistical Association, 1999, 94(446): 388-402.







