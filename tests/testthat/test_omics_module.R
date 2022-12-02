test_that("output from omics_module is OK", {
  data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", 
              "gene_annot"))
  expected <- omics_module(omics = omics, 
                           PK = PK, 
                           layers_def = layers_def, 
                           TFtargs = TFtarg_mat,
                           annot = annot, 
                           gene_annot = gene_annot)
  
  expect_is(expected, "list")
  expect_named(expected, c("pf_UB_BGe_pre", "B_prior_mat", "annot", "omics", 
                           "layers_def", "omics_meth_original"))
  expect_is(expected$pf_UB_BGe_pre, "list")
  expect_is(expected$B_prior_mat, "matrix")
  expect_is(expected$annot, "list")
  expect_is(expected$omics, "list")
  expect_is(expected$layers_def, "data.frame")
  
  expect_named(expected$pf_UB_BGe_pre, c("partition_func_UB",
                                         "parents_set_combinations", 
                                         "energy_all_configs_node", 
                                         "BGe_score_all_configs_node"))
  expect_named(expected$omics, expected$layers_def$omics, 
               ignore.order = TRUE)
  
  if(is(omics,'MultiAssayExperiment'))
  {
    expect_equal(dim(expected$omics_meth_original), 
                 dim(t(assay(omics[["meth"]]))))
    expect_equal(dim(expected$B_prior_mat), 
                 rep(sum(mapply(ncol, expected$omics)),2))
  } else if(is(omics,'list')) {
    expect_equal(dim(expected$omics_meth_original), 
                 dim(omics$meth))
    features <- 0
    for(i in seq(1,length(omics)))
    {
      features <- features + nrow(assay(omics[[i]]))
    }
    expect_equal(dim(expected$B_prior_mat), 
                 rep(features,2))
  }
})