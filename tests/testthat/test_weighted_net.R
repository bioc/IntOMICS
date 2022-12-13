test_that("output from weighted_net is OK", {
  data(list=c("OMICS_mod_res", "BN_mod_res", "gene_annot", "TFtarg_mat", "PK"))
  res_weighted <- edge_weights(mcmc_res = BN_mod_res, burn_in = 10000, 
    thin = 500, edge_freq_thres = 0.3) 
  expected <- weighted_net(cpdag_weights = res_weighted, 
    gene_annot = gene_annot, PK = PK, OMICS_mod_res = OMICS_mod_res,
    gene_ID = "gene_symbol", TFtargs = TFtarg_mat,
    B_prior_mat_weighted = B_prior_mat_weighted(BN_mod_res))
 
  expect_is(expected, "list")
  expect_named(expected, c("edge_list", "node_palette", "node_list", 
    "net_weighted", "borders_GE", "borders_CNV", "borders_METH"))
  
  expect_is(expected$edge_list, "matrix")
  expect_equal(colnames(expected$edge_list), c("from", "to", "weight", "edge_type",
    "edge"))
  expect_lt(length(expected$node_palette), 30)
  expect_is(expected$node_list, "matrix")
  expect_equal(colnames(expected$node_list), c("label", "color"))
  expect_is(expected$net_weighted, "igraph")
  expect_is(expected$borders_GE, "numeric")
})
