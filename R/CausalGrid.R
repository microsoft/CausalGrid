#' CausalGrid: A package for subgroup effects
#'
#' Intervals are (a,b], and [a,b] for the lowest. A split at x means <= and >
#' We randomize in generating train/est and trtr/trcv splits. Possibly cv.glmnet and cv.gglasso as well.
#'
#'
#' @useDynLib CausalGrid, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef formula lm model.matrix p.adjust pt qt quantile sd vcov var predict rnorm
#' @importFrom utils combn
#' @import caret
#' @import gsubfn
#' @import Rcpp
#' @import assertthat
#' @docType package
#' @name CausalGrid
NULL
#> NULL


#TODO:
# Correctness:
# - Ensure case where estimation might not have any observations in a cell
# Cleanup:
# - Encapsulate valid_partition() + bucket-splits with est_plan (deal with what happens with est error). + Doc.
# - Styler and lintr; https://style.tidyverse.org/
# - cleanup the _m functions (with their bare counterparts)
# Functionality:
# - Allow for picking paritition with # cells closest to an 'ideal' number
# - Allow for integer types with range <=g to pick those values rather than the quantiles.
# - Allow initial splits to be pre-determined.
# - Cleanup Predict function and allow y_hat (how do other packages distinguish y_hat from d_hat?). Allow mult. te
# - Like GRF, When considering a split, require that each child node have min.node.size samples with treatment value
#   less than the average, and at least that many samples with treatment value greater than or equal to the average.
# - summary method?
# - Warn if CV picks end-point (and maybe reduce minsize of cv_tr to 2 in the case we need more complex)
# - Check that tr+est and trtr's each had all the categories in sufficient quantity
# - Provide a double-selection version of lm_X_est (think through multiple D case)
# - update importance weights in change_complexity
# - Provide partial dependency functions/graphs
# - graphs: Show a rug (indicators for data points on the x-axis) or a histogram along axis.
# Usability:
# - msg for assertions
# - Have nicer factor labels (especially if split at bottom point, make [T,T] rather than [T,T+1], and redo top 
#   (to not have -1))
# - ?? switch to have min_size apply only to train_folds*splits rather than smallest (test) fold * splits. User can 
#   work around with math.
# Performance: Low-priority as doing pretty well so far
# - For each dim, save stats for each existing section and only update the stats for the section being split.
#   - Additionally, Use incremental update formulas to recompute the metrics after moving split point over a 
#     few observations (for affected cells). Can look at (fromo)[https://cran.r-project.org/web/packages/fromo/] and 
#     (twextras)[https://github.com/twolodzko/twextras/blob/master/R/cumxxx.R]
# - Pre-compute the allowable range for new splits in each slice given the cell min size.
# - see if can swap findInterval for cut() (do I need the labels)
# Checks: Check all user input types of exported functions
# Tests: More!
# R check (currently ignoring): License, top-level dev_notes.md, checking dependencies in R code, 
#          Undefined global functions or variables, tests