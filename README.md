
# Project

Tools for finding heterogeneous treatment effects (and means) based on
partitioning the covariate/feature space via full cross-cuts and solved
via greedy search. A typical usage would be analyzing an experiment to
find the high-level subgroups (a coarse partition that is useful to
humans) that differ in their estimated treatment effects.

This package is inspired by, and uses ideas from, [Causal
Tree](https://github.com/susanathey/causalTree) but aims to have the
partition be more interpretable and have better accuracy. It is slower,
though for high-level partitions this is usually not an issue.

This project is currently in an advanced prototype stage. Issues may
still be found in common usage. Please create issues for these!

## Installation

``` r
# install.package(devtools)
devtools::install_github("microsoft/CausalGrid")
```

## Usage

A simple example

``` r
# Generate some fake data
library(CausalGrid)
N= 500 #Number of obsevations
K=3 #Number of features
X = matrix(runif(N*K), ncol=K) #Features for splitting partition
d = rbinom(N, 1, 0.5) #treatment assignment
tau = as.integer(X[,1]>.5)*2-1 #true treatment effect (just heterogeneous across X1)
y = d*tau + rnorm(N, 0, 1) #outcome

# Estimate a partition
est_part0 = fit_estimate_partition(y, X, d, cv_folds=2)
get_desc_df(est_part0) #summary table
#>            X1 N_est param_ests         pval
#> 1 <=0.5017553   139  -1.138928 3.550247e-09
#> 2  >0.5017553   111   1.039283 7.201483e-08
```

## Documentation

<!--- We follow the documentation structure of https://documentation.divio.com/  --->

Documentation can be found online
[here](https://microsoft.github.io/CausalGrid/index.html) (and in the
package) or jump straight to your use case:

-   Getting started with Tutorials: [Beginner
    tutorial](https://microsoft.github.io/CausalGrid/articles/tutorial_beginner.html).

-   How-to guides: [Parallel processing for
    speedups](https://microsoft.github.io/CausalGrid/articles/howto_parallel.html).

-   Explanations/Overviews: [Separate parameters with a single
    partition](https://microsoft.github.io/CausalGrid/articles/overview_separate_params.html).

-   Reference: [API
    docs](https://microsoft.github.io/CausalGrid/reference/index.html),
    [Changelog](https://microsoft.github.io/CausalGrid/news/index.html)

## Contributing

This project welcomes contributions and suggestions. Most contributions
require you to agree to a Contributor License Agreement (CLA) declaring
that you have the right to, and actually do, grant us the rights to use
your contribution. For details, visit
<https://cla.opensource.microsoft.com>.

When you submit a pull request, a CLA bot will automatically determine
whether you need to provide a CLA and decorate the PR appropriately
(e.g., status check, comment). Simply follow the instructions provided
by the bot. You will only need to do this once across all repos using
our CLA.

This project has adopted the [Microsoft Open Source Code of
Conduct](https://opensource.microsoft.com/codeofconduct/). For more
information see the [Code of Conduct
FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or contact
<opencode@microsoft.com> with any additional questions or comments.

## Trademarks

This project may contain trademarks or logos for projects, products, or
services. Authorized use of Microsoft trademarks or logos is subject to
and must follow [Microsoft’s Trademark & Brand
Guidelines](https://www.microsoft.com/en-us/legal/intellectualproperty/trademarks/usage/general).
Use of Microsoft trademarks or logos in modified versions of this
project must not cause confusion or imply Microsoft sponsorship. Any use
of third-party trademarks or logos are subject to those third-party’s
policies.
