test_full <- plot_defective_density(forstmann, factors = c("E"),
                               correct_fun = function(d) d$R == d$S)

test_subject <- plot_defective_density(forstmann, factors = c("E", "S"), subject = 1,
                               mfcol = TRUE, layout = c(2,3))
expect_doppelganger("density_full", test_full)
expect_doppelganger("density_subject", test_subject)
