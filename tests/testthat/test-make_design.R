design_data <- design(data = data.frame(forstmann, CO = 1:nrow(forstmann)),
            model=LBA,matchfun=function(d)d$S==d$lR,
            formula=list(v~lM,sv~lM,B~E+lR,t0~E2 + CO),
            contrasts=list(v = list(lM=matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")))),
            constants=c(sv=log(1)),
            functions = list(
              E2 = function(d) factor(d$E!="speed",labels=c("speed","nonspeed"))
            ), report_p_vector=FALSE)


design_custom <- design(factors = list(S = c("left", "right"),
                                         subjects = 1:3),
                          Rlevels = c("left", "right"), model = LNR,
                          formula =list(m~0+S,s~1, t0~1),
                          constants=c(s=log(1)), report_p_vector=FALSE)

test_that("design builds coefficients from contrasts, covariates and factors", {
  expect_named(sampled_pars(design_data),
               c("v", "v_lMd", "sv_lMTRUE", "B", "B_Eneutral", "B_Eaccuracy",
                 "B_lRright", "t0", "t0_E2nonspeed", "t0_CO"))
  expect_equal(design_data$constants, c(sv = 0, A = -Inf))
  expect_equal(design_data$Clist$v$lM[, 1], c(-.5, .5))
  expect_equal(as.character(design_data$Ffunctions$E2(forstmann)),
               ifelse(forstmann$E == "speed", "speed", "nonspeed"))
  expect_named(sampled_pars(design_custom), c("m_Sleft", "m_Sright", "t0"))
  expect_equal(design_custom$Ffactors$subjects, 1:3)
  expect_equal(design_custom$Rlevels, c("left", "right"))
})
