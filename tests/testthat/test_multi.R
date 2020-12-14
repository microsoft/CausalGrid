library(testthat)

if(F) {
  devtools::load_all(".", export_all=FALSE, helpers=FALSE)
}
set.seed(1337)

ys = rnorm(100)
ds = rnorm(100)
Xs = matrix(1:200, ncol=2)
ysm = matrix(ys)
dsm = matrix(ds)
tr_splits = seq(1, by=2, length.out=50)
cvsa = seq(1, by=2, length.out=25)
cvsb = seq(2, by=2, length.out=25)
cv_foldss = list(index=list(cvsa, cvsb),
                 indexOut=list(cvsb, cvsa))

yl = rnorm(120)
dl = rnorm(120)
Xl = matrix(1:240, ncol=2)
ylm = matrix(yl)
dlm = matrix(dl)
tr_splitl = seq(1, by=2, length.out=60)
cvla = seq(1, by=2, length.out=30)
cvlb = seq(2, by=2, length.out=30)
cv_foldsl = list(index=list(cvla, cvlb),
                 index=list(cvlb, cvla))

y0 = ys
d0 = ds
X0 = Xs

y0m = ysm
d0m = dsm
X0m = Xs

y1 = list(ys, yl, ys)
d1 = list(ds, dl, ds)
X1 = list(Xs, Xl, Xs)

y1m = list(ysm, ylm, ysm)
d1m = list(dsm, dlm, dsm)
X1m = list(Xs, Xl, Xs)

y2 = ys
d2 = cbind(d1=ds, d2=1:10, d3=1:5)
X2 = Xs

y2m = ysm
d2m = cbind(d1=ds, d2=1:10, d3=1:5)
X2m = Xs

y3 = cbind(y1=ys, y2=ys, y3=ys)
d3 = ds
X3 = Xs

y3m = cbind(y1=ys, y2=ys, y3=ys)
d3m = dsm
X3m = Xs

print(0)
#fit_estimate_partition(y0, X0, d0, partition_i=2, tr_split = tr_splits)
#fit_estimate_partition(y0, X0, d0, bump_B=2, cv_folds=cv_foldss)
#fit_estimate_partition(y0m, X0, d0m, bump_B=2)

print(1)
fit_estimate_partition(y1, X1, d1, partition_i=2, tr_split = list(tr_splits,tr_splitl,tr_splits))
fit_estimate_partition(y1, X1, d1, bump_B=2, cv_folds=list(cv_foldss,cv_foldsl,cv_foldss))
fit_estimate_partition(y1m, X1, d1m, bump_B=2)

print(2)
fit_estimate_partition(y2, X2, d2, partition_i=2, tr_split = tr_splits)
fit_estimate_partition(y2, X2, d2, bump_B=2, cv_folds=cv_foldss)
fit_estimate_partition(y2m, X2, d2m, bump_B=2)

print(3)
fit_estimate_partition(y3, X3, d3, partition_i=2, tr_split = tr_splits)
fit_estimate_partition(y3, X3, d3, bump_B=2, cv_folds=cv_foldss)
fit_estimate_partition(y3m, X3, d3m, bump_B=2)

expect_equal(1,1)