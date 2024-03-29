#' Wnt signalling pathway
#'
#' A dataset containing known direct interactions between 7 genes.
#'
#' @format A data.frame with 6 rows and 3 variables:
#' \describe{
#'   \item{src_entrez}{the parent node}
#'   \item{dest_entrez}{the child node}
#'   \item{edge_type}{the edge from parent node to child node is present or missing}
#' }
#' @source \url{https://www.kegg.jp/entry/map04310}
"PK"

#' transcription factors and their known targets
#'
#' A dataset containing the direct interactions between TFs and their targets.
#'
#' @format A matrix with 22452 rows and 181 variables:
#' \describe{
#'   columns refer to TFs and rows to their targets
#' }
#' @source \url{https://maayanlab.cloud/Harmonizome/dataset/ENCODE+Transcription+Factor+Targets}
"TFtarg_mat"

#' Gene ID conversion table
#'
#' A data.frame containing the entrez ID and corresponding gene symbol.
#'
#' @format A data.frame with 8 rows and 2 variables:
#' \describe{
#'   \item{entrezID}{Entrez ID}
#'   \item{gene_symbol}{gene symbol}
#' }
"gene_annot"

#' Genes and associated methylation probes
#'
#' A named list containing the associated methylation probes of given gene.
#'
#' @format A named list with 5 components - each component corresponds to one gene:
#' \describe{
#' each component of the list is a character vector with probe names associated with given gene
#' }
#' @source \url{https://www.cancer.gov/tcga}
"annot"

#' Omics data
#'
#' A named list containing the gene expression, copy number variation and methylation data.
#'
#' @format A named list with 3 components - each component corresponds to one omics data:
#' \describe{
#' each component of the list is a matrix with samples in rows and features in columns
#' }
#' @source \url{https://www.cancer.gov/tcga}
"omics"

#' Gene ID conversion table
#'
#' A data.frame containing the modality ID, corresponding layer in BN and maximal number of parents from given layer to GE nodes.
#'
#' @format A data.frame with 3 rows and 3 variables:
#' \describe{
#'   \item{omics}{modality}
#'   \item{layer}{layer ID}
#'   \item{fan_in_ge}{maximal number of parents from given layer to single GE node}
#' }
"layers_def"

#' IntOMICS MCMC simulation result
#'
#' The output from IntOMICS::BN_module function. A named list containing results from the MCMC sampling (resulting sample is thinned and converted into corresponding CPDAGs)
#'
#' @format A named list with 3 components - each component corresponds to one omics data:
#' \describe{
#'   \item{B_prior_mat_weighted}{IntOMICS estimated empirical biological knowledge}
#'   \item{sampling.phase_res}{results from the conventional MCMC sampling - two independent simulations}
#'   \item{beta_tuning}{result from the automatically tuned MCMC algorithm}
#' }
"BN_mod_res"
