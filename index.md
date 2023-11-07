# MetricGraph #

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version-last-release/MetricGraph)](https://cran.r-project.org/package=MetricGraph)
[![R-CMD-check](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check.yml)
[![R-CMD-check-windows](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check-windows.yml/badge.svg)](https://github.com/davidbolin/MetricGraph/actions/workflows/R-CMD-check-windows.yml)

`MetricGraph` is an R package used for working with data and random fields on metric graphs, such as street or river networks. The main functionality is contained in the `metric_graph` class, which is used for specifying metric graphs, adding data to them, visualization, and other basic functions that are needed for working with data and random fields on metric graphs. The package also implements various Gaussian fields on metric graphs, and in particular the Whittle--Matérn fields introduced in the references below. 

Basic statistical tasks such as likelihood evaluation and prediction is implemented for these Gaussian fields in `MetricGraph`. Further, the package also contains interfaces to [R-INLA][ref5] and [inlabru][ref6] that facilitates using those packages for full Bayesian inference of general Latent Gaussian Models (LGMs) that includes Whittle-Matérn fields on metric graphs. 

To get started with the package, please go to the [MetricGraph: Random Fields on Metric Graphs](articles/MetricGraph.html) vignette. For more comprehensive examples, please see the vignettes in the `vignettes` tab. 


# References #
D. Bolin, A. Simas, J. Wallin (2023) [Gaussian Whittle-Matérn fields on metric graphs][ref1]. Bernoulli. In press.

D. Bolin, M. Kovács, V. Kumar, A. Simas (2023) [Regularity and numerical approximation of fractional elliptic differential equations on compact metric graphs][ref2]. Mathematics of Computation. In press.

D. Bolin, A. Simas, J. Wallin (2023) [Markov properties of Gaussian random fields on compact metric graphs][ref3]. ArXiv:2304.03190

D. Bolin, A. Simas, J. Wallin (2023) [Statistical inference for Gaussian Whittle-Matérn fields on metric graphs][ref4]. ArXiv:2304.10372


# Installation instructions #
The latest CRAN release of the package can be installed directly from CRAN with `install.packages("MetricGraph")`.

It is also possible to install the CRAN version from github by using the command:
```r
remotes::install_github("davidbolin/metricgraph", ref = "cran")
```

The latest stable version can be installed by using the command
```r
remotes::install_github("davidbolin/metricgraph", ref = "stable")
```
in R. The development version can be installed using the command
```r
remotes::install_github("davidbolin/metricgraph", ref = "devel")
```


[ref1]: https://arxiv.org/abs/2205.06163 "Gaussian Whittle-Matérn fields on metric graphs"
[ref2]: https://arxiv.org/abs/2302.03995 "Regularity and numerical approximation of fractional elliptic differential equations on compact metric graphs"
[ref3]: https://arxiv.org/abs/2304.03190 "Markov properties of Gaussian random fields on compact metric graphs"
[ref4]: https://arxiv.org/abs/2304.10372 "Statistical inference for Gaussian Whittle-Matérn fields on metric graphs"
[ref5]: https://r-inla.org "INLA homepage"
[ref6]: https://sites.google.com/inlabru.org/inlabru "inlabru homepage"
[ref7]: https://davidbolin.github.io/MetricGraph/ "MetricGraph homepage"

