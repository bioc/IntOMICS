test_that("output from bn_module is OK", {
  data("OMICS_mod_res")
  expected <- bn_module(burn_in = 100, 
                        thin = 20,
                        OMICS_mod_res = OMICS_mod_res,
                        minseglen = 2)
  
  expect_is(expected, "MCMC_sapling_res")
  
  expect_is(expected@estimated_beta, "numeric")
  expect_is(expected@estimated_len, "numeric")
  expect_is(expected@B_prior_mat_weighted, "matrix")
  expect_is(expected@CPDAGs_sim1, "list")
  expect_is(expected@CPDAGs_sim2, "list")
  expect_is(expected@beta_tuning, "matrix")
  expect_is(expected@rms, "numeric")
  
  expect_length(expected@estimated_beta, 1)
  expect_length(expected@estimated_len, 1)
  expect_equal(dim(expected@B_prior_mat_weighted), 
               dim(OMICS_mod_res$B_prior_mat))
  expect_equal(unique(mapply(is, expected@CPDAGs_sim1)), "bn")
  expect_equal(unique(mapply(is, expected@CPDAGs_sim2)), "bn")
  expect_equal(rownames(expected@beta_tuning), c("value", "len"))
  
})
