# Gaussian processes on metric graphs

'MetricGraph' is used for creation and manipulation of metric graphs,
such as street or river networks. It also has several functions
thatfacilitates operations and visualizations of data on metric graphs,
and the creation of a large class of random fields and stochastic
partial differential equations on such spaces. The main models are the
Whittle-Matérn fields, which are specified through the fractional
elliptic SPDE \$\$(\kappa^2 - \Delta)^{\alpha/2} (\tau u(s)) = W,\$\$
\\\kappa,\tau\>0\\ and \\\alpha\>1/2\\ are parameters and \\W\\ is
Gaussian white noise. It contains exact implementations of the above
model for \\\alpha=1\\ and \\\alpha=2\\, and contains approximate
implementations, via the finite element method, for any \\\alpha \>
0.5\\. It also implements models based on graph Laplacians and isotropic
covariance functions. Several utility functions for specifying graphs,
computing likelihoods, performing prediction, simulating processes, and
visualizing results on metric graphs are provided. In particular, linear
mixed effects models including random field components can be fitted to
data based on computationally efficient sparse matrix representations.
Interfaces to the R packages 'INLA' and 'inlabru' are also provided,
which facilitate working with Bayesian statistical models on metric
graphs.

## Details

At the heart of the package is the `R6` class `[metric_graph()]`. This
is used for specifying metric graphs, and contains various utility
functions which are needed for specifying Gaussian processes on such
spaces.

Linear mixed effects models are provided (see `[graph_lme]`) and perform
predictions (see `[predict.graph_lme]`). The package also has interfaces
for 'INLA' (see `[graph_spde]`), and it this interface also works with
'inlabru'.

For a more detailed introduction to the package, see the 'MetricGraph'
Vignettes.

## See also

Useful links:

- <https://davidbolin.github.io/MetricGraph/>

- Report bugs at <https://github.com/davidbolin/MetricGraph/issues>

## Author

**Maintainer**: David Bolin <davidbolin@gmail.com>

Authors:

- Alexandre Simas <alexandre.impa@gmail.com>

- Jonas Wallin <jonas.wallin81@gmail.com>
