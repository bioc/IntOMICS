#' omics_module
#' @description
#' `omics_module` data preprocessing + B_prior_mat definition + 
#'  partition function upper bound estimation + 
#'  all possible parent sets per node definition + 
#'  BGe score computation for all possible parent sets
#'
#' @param omics MultiAssayExperiment or named list containing the gene 
#' expression (possibly copy number variation and methylation data). If using 
#' named list, be aware rownames (samples) match across all objects.
#' @param PK data.frame with known interactions.
#' @param layers_def data.frame containing the modality ID, corresponding 
#' layer in BN and maximal number of parents from given layer to GE nodes.
#' @param TFtargs matrix containing the direct interactions between TFs
#' (columns) and their targets (rows).
#' @param annot named list containing the associated methylation probes 
#' of given gene.
#' @param gene_annot data.frame containing the entrez ID and corresponding 
#' gene symbol for conversion.
#' @param lm_METH logical asking whether to use linear regression to filter
#' methylation data (default=TRUE).
#' @param r_squared_thres numeric vector to define the R^2 used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.3).
#' @param p_val_thres numeric vector to define the p-value used as a threshold
#' of significance in linear regression if lm_METH=TRUE (default=0.05).
#' @param TFBS_belief numeric vector to define the belief concerning the TF 
#' and its target interaction (default=0.75).
#' @param nonGE_belief numeric vector to define the belief concerning
#' interactions of features except GE-GE (default=0.5).
#' @param woPKGE_belief numeric vector to define the belief concerning GE-GE
#' interactions without prior knowledge (default=0.5).
#' @importFrom SummarizedExperiment assay
#' @importFrom methods is
#' 
#' @examples
#' data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics", 
#' "gene_annot"), package="IntOMICS")
#' OMICS_mod_res <- omics_module(omics = omics, PK = PK, 
#'     layers_def = layers_def, TFtargs = TFtarg_mat, annot = annot, 
#'     gene_annot = gene_annot, r_squared_thres = 0.3, lm_METH = TRUE)
#'
#' @return List of 6 elements needed to init MCMC simulation 
#' @keywords internal         
#' @export
omics_module <- function(omics, PK=NULL, layers_def, TFtargs=NULL, annot=NULL, 
lm_METH = TRUE, r_squared_thres = 0.3, p_val_thres = 0.05, TFBS_belief = 0.75, 
nonGE_belief = 0.5, woPKGE_belief = 0.5, gene_annot)
{
  if (length(omics)!=nrow(layers_def)) {
    stop('Number of modalities differs in "omics" and "layers_def".')
  }
  
  if(!is(PK,"NULL"))
  {
    if(!is(PK,'data.frame') | !all(colnames(PK) %in% c("src_entrez",
      "dest_entrez", "edge_type")) | 
      !all(regexpr("EID:",c(PK$src_entrez,PK$dest_entrez), fixed = TRUE)==1))
    {
        message('Invalid input "PK". Must be data.frame with colnames 
            c("src_entrez","dest_entrez","edge_type") and gene names must 
            be in EID:XXXX format indicating Entrez IDs.')
    }
  }
  
  if(!is(layers_def,'data.frame') | !all(colnames(layers_def) %in% 
                                       c("omics", "layer", "fan_in_ge")) | 
     !all(layers_def$omics %in% names(omics)))
  {
    message('Invalid input "layers_def". Must be data.frame with colnames 
            c("omics", "layer", "fan_in_ge") and layers_def$omics 
            must fit names(omics).')
  }
  
  if(!is(TFtargs,"NULL"))
  {
    if(!is(TFtargs,'matrix') | 
     !all(grepl("EID:",unlist(dimnames(TFtargs)), fixed = TRUE)))
    {
      message('Invalid input "TFtargs". Must be matrix and dimnames must 
            be in EID:XXXX format indicating Entrez IDs.')
    }
  }
   
  if(!is(annot,"NULL"))
  {
    if(!is(annot,'list') | is(names(annot),'NULL') )
    {
      message('Invalid input "annot". Must be named list, 
            names corresponding to gene symbols.')
    }
    if(!all(names(annot) %in% gene_annot$gene_symbol))
    {
      message('Invalid input "annot". All names(annot) must be present in 
            gene_annot$gene_symbol.')
    }
  } else {
    if("meth" %in% names(omics))
    message('Input "annot" is missing even if "meth" 
        modality is available.')
  }
  
  if(!is(gene_annot, 'data.frame') | !all(colnames(gene_annot) %in% 
                                       c("entrezID","gene_symbol")) | 
     !all(grepl("EID:",gene_annot$entrezID, fixed = TRUE)))
  {
    message('Invalid input "gene_annot". Must be data.frame with colnames 
            c("entrezID","gene_symbol") and entrezID must 
            be in EID:XXXX format indicating Entrez IDs.')
  }
  
  if(TFBS_belief==woPKGE_belief)
  {
    message('GE-GE interactions without PK have the same belief as TFs-targets 
            interactions. TFs-targets interactions will be considered as GE-GE 
            interactions without prior knowledge.')
  }
  
  if(!is(r_squared_thres,'numeric') | !is(p_val_thres,'numeric') | 
     !is(TFBS_belief,'numeric') | !is(nonGE_belief,'numeric') |
     !is(woPKGE_belief,'numeric') | length(r_squared_thres)>1 |
     length(p_val_thres)>1 | length(TFBS_belief)>1 |
     length(nonGE_belief)>1 | length(woPKGE_belief)>1)
  {
    message('Invalid input. "r_squared_thres", "p_val_thres", "TFBS_belief", 
          "nonGE_belief", and "woPKGE_belief" must be numeric of length 1.')  
  }
  
  if(!is(lm_METH,'logical'))
  {
    message('Invalid input "lm_METH", must be logical.')  
  }
  
  layers_def <- layers_def[order(layers_def$layer, decreasing = TRUE),]
  
  if(is(omics,'MultiAssayExperiment'))
  {
    omics <- omics_to_list(omics = omics, gene_annot = gene_annot, 
                           layers_def = layers_def)
  } else if(is(omics,'list') & !is(names(omics),'NULL')) {
    
    omics <- omics[layers_def$omics[order(layers_def$layer, 
                                          decreasing = TRUE)]]
  } else {
    message('Invalid input "omics". Must be MultiAssayExperiment or 
            named list.')
  }
  
  if(length(unique(mapply(nrow,omics)))!=1) {
    message('Invalid input "omics". 
            Number of samples differ across modalities.')
  }
  
  if(!is(annot,"NULL"))
  {
    names(annot) <- gene_annot$entrezID[
    match(names(annot),gene_annot$gene_symbol)]
  }
  
  B <- b_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
      annot = annot, lm_METH = lm_METH, r_squared_thres = r_squared_thres, 
      p_val_thres = p_val_thres, TFtargs = TFtargs, 
      TFBS_belief = TFBS_belief, nonGE_belief = nonGE_belief, 
      woPKGE_belief = woPKGE_belief)
  pf_UB_res <- pf_ub_est(omics = B$omics, layers_def = layers_def, 
      B_prior_mat = B$B_prior_mat, annot = B$annot)
  return(list(pf_UB_BGe_pre = pf_UB_res, B_prior_mat = B$B_prior_mat, 
      annot = B$annot, omics = B$omics, layers_def = layers_def,
      omics_meth_original = B$omics_meth_original))
}
