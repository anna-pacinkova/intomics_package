% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BN_module_func.R
\name{transient_phase}
\alias{transient_phase}
\title{transient phase}
\usage{
transient_phase(
  first.adapt.phase_net,
  omics,
  B_prior_mat,
  layers_def,
  energy_all_configs_node,
  prob_mbr,
  BGe_score_all_configs_node,
  parent_set_combinations,
  annot
)
}
\arguments{
\item{first.adapt.phase_net}{list output of the first.adapt.phase
function.}

\item{omics}{named list containing the gene expression
(possibly copy number variation and methylation data).
Each component of the list is a matrix with samples in rows and
features in columns.}

\item{B_prior_mat}{a biological prior matrix.}

\item{layers_def}{data.frame containing the modality ID, corresponding
layer in BN and maximal number of parents from given layer to GE nodes.}

\item{energy_all_configs_node}{list of nodes energy for all possible
parent set configurations.}

\item{prob_mbr}{numeric vector probability of the MBR step.}

\item{BGe_score_all_configs_node}{list of nodes BGe score for all
possible parent set configurations.}

\item{parent_set_combinations}{list of all possible parent set
configuration for all nodes available.}

\item{annot}{named list containing the associated methylation
probes of given gene.}
}
\value{
List of 1 element: first adaption phase and transient phase result
}
\description{
\code{transient_phase} This phase verify if the chain is moving towards
the mode of target distribution.
}
\examples{
data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics"),
    package="IntOMICS")
OMICS_mod_res <- OMICS_module(omics = omics, PK = PK, 
    layers_def = layers_def, TFtargs = TFtarg_mat, annot = annot, 
    r_squared_thres = 0.3, lm_METH = TRUE)
first.adapt.phase_net <- first_adapt_phase(len = 5, 
    omics = OMICS_mod_res$omics, B_prior_mat = OMICS_mod_res$B_prior_mat, 
    energy_all_configs_node = 
    OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
    layers_def = OMICS_mod_res$layers_def, prob_mbr = 0.07, 
    BGe_score_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
    parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
    annot = OMICS_mod_res$annot)
transient_phase(first.adapt.phase_net = first.adapt.phase_net, 
    omics = OMICS_mod_res$omics, B_prior_mat = OMICS_mod_res$B_prior_mat, 
    layers_def = OMICS_mod_res$layers_def, 
    energy_all_configs_node = 
    OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
    prob_mbr = 0.07, annot = OMICS_mod_res$annot,
    BGe_score_all_configs_node = 
    OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
    parent_set_combinations = 
    OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations) 

}
