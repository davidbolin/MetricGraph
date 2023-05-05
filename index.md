# MetricGraph #

`MetricGraph` is an R package used for working with data and random fields on metric graphs, such as street or river networks. The main functionality is contained in the `metric_graph` class, which is used for specifying metric graphs, adding data to them, visualization, and other basic functions that are needed for working with data and random fields on metric graphs. The package also implements various Gaussian fields on metric graphs, and in particular the Whittle--Matérn fields introduced in the references below. 

Basic statistical tasks such as likelihood evaluation and prediction is implemented for these Gaussian fields in `MetricGraph`. Further, the package also contains interfaces to [R-INLA][ref5] and [inlabru][ref6] that facilitates using those packages for full Bayesian inference of general Latent Gaussian Models (LGMs) that includes Whittle-Matérn fields on metric graphs. 

We refer to the [package homepage][ref7] for detailed tutorials on all different aspects of the package.


# References #
D. Bolin, A. Simas, J. Wallin (2022) [Gaussian Whittle-Matérn fields on metric graphs][ref1]. ArXiv:2205.06163

D. Bolin, M. Kovács, V. Kumar, A. Simas (2023) [Regularity and numerical approximation of fractional elliptic differential equations on compact metric graphs][ref2]. ArXiv:2302.03995

D. Bolin, A. Simas, J. Wallin (2023) [Markov properties of Gaussian random fields on compact metric graphs][ref3]. ArXiv:2304.03190

D. Bolin, A. Simas, J. Wallin (2023) [Statistical inference for Gaussian Whittle-Matérn fields on metric graphs][ref4]. ArXiv:2304.10372


# Installation instructions #
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

