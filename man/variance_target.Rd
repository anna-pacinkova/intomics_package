% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BN_module_func.R
\name{variance_target}
\alias{variance_target}
\title{Second adaption phase variance tuning}
\usage{
variance_target(
  transient.phase_net,
  constant,
  fin,
  B_prior_mat,
  omics,
  parent_set_combinations,
  BGe_score_all_configs_node,
  layers_def,
  prob_mbr,
  annot
)
}
\arguments{
\item{transient.phase_net}{list output of the variance_target or
transient.phase function.}

\item{constant}{numeric vector used to multiply the beta_sd to determine
the variance of the distribution of the hyperparameter beta.}

\item{fin}{numeric vector iteration to stop.}

\item{B_prior_mat}{a biological prior matrix.}

\item{omics}{named list containing the gene expression (possibly copy number
variation and methylation data).
Each component of the list is a matrix with samples in rows and
features in columns.}

\item{parent_set_combinations}{list of all possible parent set configuration
for all nodes available.}

\item{BGe_score_all_configs_node}{list of nodes BGe score for all possible
parent set configurations.}

\item{layers_def}{data.frame containing the modality ID, corresponding layer
in BN and maximal number of parents from given layer to GE nodes.}

\item{prob_mbr}{numeric vector probability of the MBR step.}

\item{annot}{named list containing the associated methylation probes
of given gene.}
}
\value{
Large List of 3 elements: second adaptive phase result
with possible MCMC mixing; acceptance rate of hyperparameter beta;
SD of hyperparameter beta
}
\description{
\code{variance_target} This phase identifies the proposal distribution that
has a similar covariance structure with the target distribution.
This is part of second_adapt_phase.
}
\examples{
data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics"),
    package="IntOMICS")
OMICS_mod_res <- OMICS_module(omics = omics, PK = PK, annot = annot, 
    layers_def = layers_def, TFtargs = TFtarg_mat, r_squared_thres = 0.3,
    lm_METH = TRUE)
first.adapt.phase_net <- first_adapt_phase(omics = OMICS_mod_res$omics, 
    B_prior_mat = OMICS_mod_res$B_prior_mat, 
    energy_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
    len = 5, layers_def = OMICS_mod_res$layers_def, prob_mbr = 0.07,
    BGe_score_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
    parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
    annot = OMICS_mod_res$annot)
transient.phase_net <- transient_phase(annot = OMICS_mod_res$annot, 
    first.adapt.phase_net = first.adapt.phase_net, 
    omics = OMICS_mod_res$omics, prob_mbr = 0.07, 
    B_prior_mat = OMICS_mod_res$B_prior_mat, 
    layers_def = OMICS_mod_res$layers_def, 
    energy_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node,
    BGe_score_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
    parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations) 
variance_target(transient.phase_net = transient.phase_net, 
    constant = 1.586667, fin = 200, B_prior_mat = OMICS_mod_res$B_prior_mat,
    parent_set_combinations = OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations, 
    BGe_score_all_configs_node = OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node, 
    layers_def = OMICS_mod_res$layers_def, omics = OMICS_mod_res$omics, 
    prob_mbr = 0.07, annot = OMICS_mod_res$annot)

}
