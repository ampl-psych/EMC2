matchfun <- function(d)d$S==d$lR
# design an "average and difference" contrast matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))

# specify design
design_LBABE <- design(data = forstmann,model=LBA,
                formula=list(v~1,sv~1,B~E+lR,A~1,t0~1),
                constants=c(sv=log(1)))

pmean <- c(v=1, B=log(.5),B_Eneutral=log(1.5),
           B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
psd <- c(v=1,
         B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
prior_LBABE <- prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)
emc <- make_emc(forstmann,design_LBABE,type="standard",  prior=prior_LBABE,
                          compress = FALSE)

test_that("make_emc", {
  expect_snapshot(str(emc, give.attr = F))
})
