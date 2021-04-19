# To run in the command-line with load_all: change do_load_all=T, then run the code in the first if(FALSE), subsequent runs just run that last line of the False block
# To run in the command-line with load_all: then run the code in the first if(FALSE), subsequent runs just run that last line of the False block

library(testthat)
loadNamespace(rprojroot)
root_dir <- rprojroot::find_package_root_file() #R cmd check doesn't copy over git and RStudio proj file

if(FALSE) { #Run manually to debug
  debugSource(paste0(rprojroot::find_testthat_root_file(),"/testrun.R"))
}

library(CausalGrid)

set.seed(1337)

context("Test Run")

if(getwd()==root_dir) {
  print(root_dir)
  print(getwd())
  box::use(tests/dgps)
} else{ #testing
  box::use(../dgps)
}

data <- dgps$mix_data_d(n=1000)
breaks_per_dim = list(c(0.5), c(0))

# Does Bumping work -------------------

ret_bmp1 <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, bump_samples=2, bump_complexity=list(doCV=FALSE, incl_comp_in_pick=FALSE))
ret_bmp2 <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, bump_samples=2, bump_complexity=list(doCV=TRUE, incl_comp_in_pick=FALSE))
ret_bmp3 <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, bump_samples=2, bump_complexity=list(doCV=FALSE, incl_comp_in_pick=TRUE))
ret_bmp4 <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, bump_samples=2, bump_complexity=list(doCV=TRUE, incl_comp_in_pick=TRUE))

# Make sure partition is fine with 0 obs ----------------
X_range = get_X_range(data$X)
ex_part = add_partition_split(grid_partition(data$X), partition_split(1, 0.5))
ex_fact = predict(ex_part, matrix(0.1, ncol=2, nrow=2))
test_that("# of levels will be full even if they don't appear in data.", {expect_equal(length(levels(ex_fact)), 2)})

# Just y ---------------

ret1 <- fit_estimate_partition(data$y, data$X, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim)
print(ret1$partition)
test_that("We get OK results (OOS)", {
  expect_equal(ret1$partition$nsplits_by_dim, c(1,1))
})


# Include d ---------------

ret1d <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim)
print(ret1d$partition)
test_that("We get OK results (OOS)", {
  expect_equal(ret1d$partition$nsplits_by_dim, c(1,1))
})
test_any_sign_effect(ret1d, check_negative=T, method="fdr") #
#test_any_sign_effect(ret1d, check_negative=T, method="sim_mom_ineq") #the sim produces treatment effect with 0 std err, so causes problems

ret2d <- dgps$suppress_warnings(fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, ctrl_method="all"), "essentially perfect fit: summary may be unreliable")
print(ret2d$partition)
#TODO: Should I check this?
#test_that("We get OK results (OOS)", {
#  expect_equal(ret2d$partition$nsplits_by_dim, c(1,1))
#})

ret3d <- dgps$suppress_warnings(fit_estimate_partition(data$y, data$X, data$d, cv_folds=3, verbosity=0, breaks_per_dim=breaks_per_dim, ctrl_method="LassoCV"), "essentially perfect fit: summary may be unreliable")
print(ret3d$partition)
#TODO: Should I check this?
#test_that("We get OK results (OOS)", {
#  expect_equal(ret3d$partition$nsplits_by_dim, c(1,1))
#})

ret4d <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, ctrl_method="RF")
print(ret4d$partition)
#TODO: Should I check this?
#test_that("We get OK results (OOS)", {
#  expect_equal(ret4d$partition$nsplits_by_dim, c(1,1))
#})

ret1db <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, bump_samples=2)


ret1dc <- fit_estimate_partition(data$y, data$X, data$d, cv_folds=2, verbosity=0, breaks_per_dim=breaks_per_dim, importance_type="single")

# Test no cols ------------------------------------
temp <- fit_estimate_partition(rnorm(100), matrix(NA, nrow=100, ncol=0), rnorm(100))


# Test the output/verbosity ------------------------------------
X_3 = data$X
X_3$X3 = data$X$X2
pot_break_points_3 = breaks_per_dim
pot_break_points_3[[3]] = breaks_per_dim[[2]]
print("---------------")
ret1dd <- fit_estimate_partition(data$y, X_3, data$d, cv_folds=2, verbosity=2, breaks_per_dim=pot_break_points_3, importance_type="interaction", bump_samples=3)
print("---------------")
ret1dd <- fit_estimate_partition(data$y, X_3, data$d, cv_folds=2, verbosity=1, breaks_per_dim=pot_break_points_3, importance_type="interaction", bump_samples=3)
print("---------------")
ret1dd <- fit_estimate_partition(data$y, X_3, data$d, cv_folds=2, verbosity=0, breaks_per_dim=pot_break_points_3, importance_type="interaction", bump_samples=3)

