% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BN_module_func.R
\name{BN_module}
\alias{BN_module}
\title{#' BN module}
\usage{
BN_module(
  burn_in,
  thin,
  seed1,
  seed2,
  OMICS_mod_res,
  minseglen,
  len = 5,
  prob_mbr = 0.07
)
}
\description{
\code{BN_module} Performs automatically tuned MCMC sampling from posterior distribution together with conventional MCMC sampling using empirical biological prior matrix to sample network structures from posterior distribution.
}
