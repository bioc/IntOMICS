test_that("output from edge_weights is OK", {
  data(list=c("BN_mod_res"))
  expected <- edge_weights(mcmc_res = BN_mod_res, burn_in = 10000, 
                           thin = 500, edge_freq_thres = 0.3) 
  expect_is(expected, "data.frame")
  expect_equal(colnames(expected), c("edge", "from.x",  "to.x", "strength.x",
      "direction.x", "from.y", "to.y", "strength.y", "direction.y", "strength"))
})