# To run in the command-line with load_all: then run the code in the first if(FALSE), subsequent runs just run that last line of the False block

library(testthat)

if(FALSE) { #Run manually to debug
  loadNamespace(rprojroot)
  debugSource(paste0(rprojroot::find_testthat_root_file(),"/test_multi.R"))
}

library(CausalGrid)

set.seed(1337)

context("Testing multiple estimates")

#Define different data shapes and structures (short vs long, non-matrix vs matrix)
#Short variants
ys = rnorm(60)
ds = rnorm(60)
Xs = data.frame(X1=1:60, X2=61:120)
Xsm = as.matrix(Xs)
ysm = matrix(ys)
dsm = matrix(ds)
tr_splits = seq(1, by=2, length.out=30)
cv_foldss = rep(1, 30)
cv_foldss[seq(2, by=3, length.out=10)] = 2
cv_foldss[seq(3, by=3, length.out=10)] = 3

#Long-variants
yl = rnorm(120)
dl = rnorm(120)
Xl = data.frame(X1=1:120, X2=121:240)
Xlm = as.matrix(Xl)
ylm = matrix(yl)
dlm = matrix(dl)
tr_splitl = seq(1, by=2, length.out=60)
cv_foldsl = rep(1, 60)
cv_foldsl[seq(2, by=3, length.out=20)] = 2
cv_foldsl[seq(3, by=3, length.out=20)] = 3

y0 = ys
d0 = ds
X0 = Xs

y0m = ysm
d0m = dsm
X0m = Xsm

y1 = list(ys, yl, ys)
d1 = list(ds, dl, ds)
X1 = list(Xs, Xl, Xs)

y1m = list(ysm, ylm, ysm)
d1m = list(dsm, dlm, dsm)
X1m = list(Xsm, Xlm, Xsm)

y2 = ys
d2 = cbind(d1=ds, d2=1:10, d3=1:5)
X2 = Xs

y2m = ysm
d2m = cbind(d1=ds, d2=1:10, d3=1:5)
X2m = Xsm

y3 = cbind(y1=ys, y2=ys, y3=ys)
d3 = ds
X3 = Xs

y3m = cbind(y1=ys, y2=ys, y3=ys)
d3m = dsm
X3m = Xsm

#Alternate whether: (a) using random or fixed splits, cv or partition_i, non-matrix vs matrix data structures, controls
#print(0)
fit_estimate_partition(y0, X0, d0, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2, ctrl_method="LassoCV")
fit_estimate_partition(y0, X0, d0, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2, ctrl_method="RF")
fit_estimate_partition(y0, X0, d0, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2)
fit_estimate_partition(y0m, X0m, d0m, bump_samples=2)

#print(1)
fit_estimate_partition(y1m, X1m, d1m, partition_i=2, tr_split = list(tr_splits,tr_splitl,tr_splits), cv_folds=list(cv_foldss,cv_foldsl,cv_foldss), bump_samples=2, ctrl_method="LassoCV")
fit_estimate_partition(y1m, X1m, d1m, partition_i=2, tr_split = list(tr_splits,tr_splitl,tr_splits), cv_folds=list(cv_foldss,cv_foldsl,cv_foldss), bump_samples=2, ctrl_method="RF")
fit_estimate_partition(y1m, X1m, d1m, partition_i=2, tr_split = list(tr_splits,tr_splitl,tr_splits), cv_folds=list(cv_foldss,cv_foldsl,cv_foldss), bump_samples=2)
fit_estimate_partition(y1, X1, d1, bump_samples=2)

#print(2)
fit_estimate_partition(y2, X2, d2, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2, ctrl_method="LassoCV")
fit_estimate_partition(y2, X2, d2, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2, ctrl_method="RF")
fit_estimate_partition(y2, X2, d2, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2)
fit_estimate_partition(y2m, X2m, d2m, bump_samples=2)

#print(3)
fit_estimate_partition(y3m, X3m, d3m, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2, ctrl_method="LassoCV")
fit_estimate_partition(y3m, X3m, d3m, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2, ctrl_method="RF")
fit_estimate_partition(y3m, X3m, d3m, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2)
fit_estimate_partition(y3, X3, d3, bump_samples=2)

#print(4) #outcome-mean
fit_estimate_partition(y1m, X1m, partition_i=2, tr_split = list(tr_splits,tr_splitl,tr_splits), cv_folds=list(cv_foldss,cv_foldsl,cv_foldss), bump_samples=2)
fit_estimate_partition(y3m, X3m, partition_i=2, tr_split = tr_splits, cv_folds=cv_foldss, bump_samples=2)

test_that("We get through the tests", {
  expect_equal(1,1)
})
