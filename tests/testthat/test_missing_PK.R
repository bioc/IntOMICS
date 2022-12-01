test_that("Output from omics_module", {
  data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", 
              "gene_annot"))
  expect_message(
    omics_module(omics = omics, 
                 layers_def = layers_def, 
                 TFtargs = TFtarg_mat,
                 annot = annot, 
                 gene_annot = gene_annot),
    'Please, consider adding PK to increase the IntOMICS prediction accuracy.')
})