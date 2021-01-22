#' CausalGrid: A package for subgroup effects
#'
#' Tools for finding heterogeneous treatment effects (and means) based on
#' partitioning the covariate/feature space via full cross-cuts and solved via
#' greedy search. A typical usage would be analyzing and experiment to find the
#' high-level subgroups (a coarse partition that is useful to humans) that
#' differ in their estimated treatment effects.
#'
#' This package is inspired by, and uses ideas from, \code{Causal
#' Tree} but aims to have the
#' partition be more interpretable and have better accuracy. It is slower,
#' though for high-level partitions this is usually not an issue.
#'
#' Subgroups are constructed as a grid over the features/covariates \code{X}.
#' For example, with 1 feature going from 0 to 1 it may split at values c1, c2,
#' resulting in segments \code{[0,c1], (c1, c2], (c2,1]}. A split at value
#' \code{c} means it splits <= and >. The segments may be of uneven sizes.
#' Splits along several features result in a grid by constructing the Cartesian
#' product of the feature-specific splits. Not all features will necessarily be
#' split or split the same number of times.
#'
#' The main entry point is \code{\link{fit_estimate_partition}}.
#'
#' Randomization: This package should be able to be run with no randomness. With
#' default/simple parameters the following places randomize but can be
#' overridden. \itemize{ 
#' \item  Generating train/est splits. Can be overridden
#' by providing \code{tr_split} 
#' \item  Generating trtr/trcv splits. Can be
#' overridden by providing \code{cv_folds} 
#' \item  Bumping samples. Can be
#' overridden by providing list of samples for \code{bump_samples} 
#' \item
#' Estimation plans: Provide ones ( \code{lm_est(lasso=TRUE,...)} and
#' \code{grid_rf}) use  \code{cv_folds}. User-made ones should too. 
#' }
#'
#' @useDynLib CausalGrid, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats coef formula lm model.matrix p.adjust pt qt quantile sd
#'   vcov var predict rnorm
#' @importFrom utils combn
#' @import caret
#' @import gsubfn
#' @import Rcpp
#' @import assertthat
#' @docType package
#' @name CausalGrid
NULL
#> NULL

#' Source code layout:
#'
#' - Given the multiple estimates regime, most of the data process needs to
#' account for different formats. Wrappers for these are in `utils.R`
#'
#' - Constant-vector checking: Given estimations are costly, we often check if
#' there's any variance in D before estimating. Therefore we have a
#' check-if-constant vector function. We have a compiled C++ function
#' `const_vect()` that improves speed on big data. See
#' https://stackoverflow.com/questions/4752275/. This can make it bit harder to
#' do iterative development with the package (technically should use "Install
#' and Restart" rather than "Load all" and this is actually necessary with
#' parallel processing). Alternatively one can remove `const_vect` and assign
#' `const_vect = const_vectr` to use the R-only version from `utils.R`
#'
#' Data structure notes:
#'
#' - In order to work with both X as matrix and data.frame I used X[,k], but
#' this is messed up with incoming Tibbles so convert those.
#'
#' R checks currently ignoring:
#'
#' - Undefined global functions or variables due to my use of gsubfn, (anyway to
#' disable?)
#'
#' - undocumented ... in main two entry points (what to say in doc?)
#'
#' - Non-standard files/directories found at top level: 'CODE_OF_CONDUCT.md'
#' 'SECURITY.md' 'SUPPORT.md'

#' TODO (keeping here so centralized and can prioritize easier):
#'
#' Correctness:
#'
#' - When getting factors for s ample from a grid, if the new data is more
#' extreme (even for numbers) we throw NAs. Fix.
#'
#' Cleanup:
#'
#' - Encapsulate valid_partition() + bucket-splits with est_plan (deal with what
#' happens with est error). + Doc.
#'
#' - Styler and lintr; https://style.tidyverse.org/
#'
#' - Merge `factor_from_idxs` and `foldlists_to_foldids`
#'
#' Functionality:
#'
#' - Implement honest version of eval metric
#'
#' - Allow for picking paritition with # cells closest to an 'ideal' number
#'
#' - Allow initial splits to be pre-determined.
#'
#' - Cleanup Predict function and allow y_hat (how do other packages distinguish
#' y_hat from d_hat?). Allow mult. te
#'
#' - Like GRF, When considering a split, require that each child node have
#' min.node.size samples with treatment value less than the average, and at
#' least that many samples with treatment value greater than or equal to the
#' average.
#'
#' - summary method? maybe use get_desc_df
#'
#' - Warn if CV picks end-point (and maybe reduce minsize of cv_tr to 2 in the
#' case we need more complex)
#'
#' - Check that tr+est and trtr's each had all the categories in sufficient
#' quantity
#'
#' - Provide a double-selection version of lm_est (think through multiple D
#' case)
#'
#' - update importance weights in change_complexity
#'
#' - Provide partial dependency functions/graphs
#'
#' - graphs: Show a rug (indicators for data points on the x-axis) or a
#' histogram along axis.
#'
#' Usability:
#'
#' - Documentation descriptions (separate from titles). Define the bump samples
#' explicit case.
#'
#' - Have nicer factor labels (especially if split at bottom point, make [T,T]
#' rather than [T,T+1], and redo top (to not have -1))
#'
#' - ?? switch to have `min_size` apply only to `train_folds*splits` rather than
#' smallest (test) fold * splits. User can work around with math.
#'
#' - Filter out "essentially perfect fit: summary may be unreliable" from lm
#'
#' Performance: Low-priority as doing pretty well so far
#'
#' - Is it faster to disable the check for constant D (and just catch the
#' estimation error?)
#'
#' - For each dim, save stats for each existing section and only update the
#' stats for the section being split.
#'
#' - Additionally, Use incremental update formulas to recompute the metrics
#' after moving split point over a few observations (for affected cells). Can
#' look at (fromo)[https://cran.r-project.org/web/packages/fromo/] and
#' (twextras)[https://github.com/twolodzko/twextras/blob/master/R/cumxxx.R]
#'
#' - Pre-compute the allowable range for new splits in each slice given the cell
#' min size.
#'
#' - see if can swap `findInterval` for `cut()` (do I need the labels)
#'
#' Graphing:
#'
#' - Add marginal plots:
#' https://www.r-graph-gallery.com/277-marginal-histogram-for-ggplot2.html
#'
#' Checks: Check all user input types of exported functions
#'
#' Tests: More!

