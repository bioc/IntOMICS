test_that("output from bn_module is OK", {
  if(interactive())
  {
    data("OMICS_mod_res")
    expected <- bn_module(burn_in = 100, 
                          thin = 20,
                          OMICS_mod_res = OMICS_mod_res,
                          minseglen = 2)
  
    expect_is(expected, "MCMC_sapling_res")
  
    expect_is(estimated_beta(expected), "numeric")
    expect_is(estimated_len(expected), "numeric")
    expect_is(B_prior_mat_weighted(expected), "matrix")
    expect_is(CPDAGs_sim1(expected), "list")
    expect_is(CPDAGs_sim2(expected), "list")
    expect_is(beta_tuning(expected), "matrix")
    expect_is(rms(expected), "numeric")
  
    expect_length(estimated_beta(expected), 1)
    expect_length(estimated_len(expected), 1)
    expect_equal(dim(B_prior_mat_weighted(expected)), 
                 dim(OMICS_mod_res$B_prior_mat))
    expect_equal(unique(mapply(is, CPDAGs_sim1(expected))), "bn")
    expect_equal(unique(mapply(is, CPDAGs_sim2(expected))), "bn")
    expect_equal(rownames(beta_tuning(expected)), c("value", "len"))
  } 
})
