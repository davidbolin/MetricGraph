# MetricGraph <a href="https://davidbolin.github.io/MetricGraph/"><img src="/man/figures/logo.png" align="right" height="138" /></a>

[![R-CMD-check](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check.yml)
[![R-CMD-check-windows](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check-windows.yml/badge.svg)](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check-windows.yml)

## Overview 

`MetricGraph` is an R package used for working with data and random fields on metric graphs, such as street or river networks. The main functionality is contained in the `metric_graph` class, which is used for specifying metric graphs, adding data to them, visualization, and other basic functions that are needed for working with data and random fields on metric graphs. The package also implements various Gaussian fields on metric graphs, and in particular the Whittle--Matérn fields introduced in the references below. 

Basic statistical tasks such as likelihood evaluation and prediction is implemented for these Gaussian fields in `MetricGraph`. Further, the package also contains interfaces to [R-INLA][ref5] and [inlabru][ref6] that facilitates using those packages for full Bayesian inference of general Latent Gaussian Models (LGMs) that includes Whittle-Matérn fields on metric graphs. 

We refer to the [package homepage][ref7] for detailed tutorials on all different aspects of the package.

# Installation instructions #
The latest stable version can be installed by using the command
```r
remotes::install_github("davidbolin/metricgraph", ref = "stable")
```
in R. The development version can be installed using the command
```r
remotes::install_github("davidbolin/metricgraph", ref = "devel")
```

# References #
D. Bolin, A. Simas, J. Wallin (2022) [Gaussian Whittle-Matérn fields on metric graphs][ref1]. ArXiv:2205.06163

D. Bolin, M. Kovács, V. Kumar, A. Simas (2023) [Regularity and numerical approximation of fractional elliptic differential equations on compact metric graphs][ref2]. ArXiv:2302.03995

D. Bolin, A. Simas, J. Wallin (2023) [Markov properties of Gaussian random fields on compact metric graphs][ref3]. ArXiv:2304.03190

D. Bolin, A. Simas, J. Wallin (2023) [Statistical inference for Gaussian Whittle-Matérn fields on metric graphs][ref4]. ArXiv:2304.10372

# Repository branch workflows #
The package version format for released versions is `major.minor.bugfix`. All regular development should be performed on the `devel` branch or in a feature branch, managed with `git flow feature`. Ideally, all the changes should be made on the `devel` branch. The `devel` version of the package should contain unit tests and examples for all important functions. Several functions may depend on `INLA`. Examples and tests for such functions might create problems when submitting to CRAN. To solve this problem, we created some Github Actions scripts that get the examples and tests depending on `INLA` on the `devel` branch and adapt to versions that will not fail on CRAN. Therefore, the best way to handle these situations is to avoid as much as possible to do any push to the `stable` branch. The idea is to update the `stable` branch by merges following the workflow that will be described below. 
The examples that depend on `INLA` should have the following structure:

```
#' \donttest{ #devel version
#' library(INLA)
#' 
#' # The contents of the example
#'
#' #devel.tag
#' }
```

The tests that depend on `INLA` should have the following structure:

```
test_that("Description of the test", {
testthat::skip_on_cran()
  if (!requireNamespace("INLA", quietly=TRUE))
    testthat::skip(message = 'INLA package is not installed. (see www.r-inla.org/download-install)')
  
  old_threads <- INLA::inla.getOption("num.threads")
  INLA::inla.setOption(num.threads = "1:1")
  
  # The contents of the test
  
  INLA::inla.setOption(num.threads = old_threads)
})
```

On the `devel` branch, the vestion number is `major.minor.bugfix.9000`, where the first three components reflect the latest released version with changes present in the `default` branch. Bugfixes should be applied via the `git flow bugfix` and `git flow hotfix` methods, as indicated below. For `git flow` configuration, use `master` as the stable master branch, `devel` as the develop branch, and `v` as the version tag prefix. Hotfixes directly `stable` should be avoided whenever possible to minimize conflicts on merges. See [the `git flow` tutorial](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) for more information.

For non `master` and `devel` branches that collaborators need access to (e.g. release branches, feature branches, etc, use the `git flow publish` mechanism).


  * Prepare a new stable release with CRAN submission:
```
git flow release start major.(minor+1).0

# In R, the following updates the version number in DESCRIPTION and NEWS:
usethis::use_version("minor") 
## At this point, see the CRAN submission section below.
git flow release finish 'VERSION'
# In the stable merge, accept all incoming changes.
# Push the changes and do adjustments if needed.
# After handling the stable branch, merge back with devel.
# In R, the following updates the dev version number in DESCRIPTION and NEWS:
usethis::use_dev_version() 
```
  * Do a hotfix (branch from stable branch; use bugfix for release branch bugfixes):
```
git flow hotfix start hotfix_branch_name
## Do the bugfix, update the verison number major.minor.(bugfix+1), and commit
## Optionally, do CRAN submission
git flow hotfix finish hotfix_branch_name
## Resolve merge conflicts, if any
```
  * CRAN submission
```
## Perform CRAN checks on the stable branch
## If unsuccessful then do bugfixes with increasing bugfix version, until ok
## Submit to CRAN
## If not accepted then do more bugfixes and repeat
```

[ref1]: https://arxiv.org/abs/2205.06163 "Gaussian Whittle-Matérn fields on metric graphs"
[ref2]: https://arxiv.org/abs/2302.03995 "Regularity and numerical approximation of fractional elliptic differential equations on compact metric graphs"
[ref3]: https://arxiv.org/abs/2304.03190 "Markov properties of Gaussian random fields on compact metric graphs"
[ref4]: https://arxiv.org/abs/2304.10372 "Statistical inference for Gaussian Whittle-Matérn fields on metric graphs"
[ref5]: https://r-inla.org "INLA homepage"
[ref6]: https://sites.google.com/inlabru.org/inlabru "inlabru homepage"
[ref7]: https://davidbolin.github.io/MetricGraph/ "MetricGraph homepage"

