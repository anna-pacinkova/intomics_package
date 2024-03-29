% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BN_module_func.R
\name{first_adapt_phase}
\alias{first_adapt_phase}
\title{1st adaption phase}
\usage{
first_adapt_phase(
  omics,
  B_prior_mat,
  energy_all_configs_node,
  len,
  layers_def,
  prob_mbr,
  BGe_score_all_configs_node,
  parent_set_combinations,
  annot
)
}
\arguments{
\item{omics}{named list containing the gene expression
(possibly copy number variation and methylation data).
Each component of the list is a matrix with samples in rows and
features in columns.}

\item{B_prior_mat}{a biological prior matrix.}

\item{energy_all_configs_node}{list of nodes energy for all possible
parent set configurations.}

\item{len}{numeric vector initial width of the sampling interval
for hyperparameter beta.}

\item{layers_def}{data.frame containing the modality ID, corresponding
layer in BN and maximal number of parents from given layer to GE nodes.}

\item{prob_mbr}{numeric vector probability of the MBR step.}

\item{BGe_score_all_configs_node}{list of nodes BGe score
for all possible parent set configurations.}

\item{parent_set_combinations}{list of all possible parent set
configuration for all nodes available.}

\item{annot}{named list containing the associated methylation
probes of given gene.}
}
\value{
List of 1 element: first adaption phase result
}
\description{
\code{first_adapt_phase} 1st adaption phase of the adaptive MCMC:
the variance of the proposal distribution is changed to achieve
the MC acceptance rate of 0.44.
}
\examples{
data(list=c("PK", "TFtarg_mat", "annot", "layers_def", 
   "omics"), package="IntOMICS")
OMICS_mod_res <- OMICS_module(omics = omics, PK = PK, 
   layers_def = layers_def, TFtargs = TFtarg_mat,
   annot = annot, r_squared_thres = 0.3, 
   lm_METH = TRUE)
first_adapt_phase(omics = OMICS_mod_res$omics, 
   B_prior_mat = OMICS_mod_res$B_prior_mat, len = 5, 
   energy_all_configs_node = 
   OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
   layers_def = OMICS_mod_res$layers_def, prob_mbr = 0.07,
   BGe_score_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
   parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
   annot = OMICS_mod_res$annot)

}
