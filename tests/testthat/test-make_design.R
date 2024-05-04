test_that("make_design", {
  expect_snapshot(
    design <- make_design(data = data.frame(forstmann, CO = 1:nrow(forstmann)),
                model=LBA,matchfun=function(d)d$S==d$lR,
                formula=list(v~lM,sv~lM,B~E+lR,t0~E2 + CO),
                contrasts=list(v = list(lM=matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")))),
                constants=c(sv=log(1)),
                functions = list(
                  E2 = function(d) factor(d$E!="speed",labels=c("speed","nonspeed"))
                ))
  )
})
