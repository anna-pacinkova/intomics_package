% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BN_module_func.R
\name{epsilon}
\alias{epsilon}
\title{Epsilon}
\usage{
epsilon(net, B_prior_mat)
}
\arguments{
\item{net}{adjacency matrix of given network.}

\item{B_prior_mat}{a biological prior matrix.}
}
\value{
Numeric vector of length 1: epsilon of given adjacency matrix
(needed to compute energy of given adjacency matrix)
}
\description{
\code{epsilon} This function returns the epsilon value for each variable/node
of the network.
The sum of the epsilons of all variables/nodes in the network gives us
the energy of given network.
}
\examples{

data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics"),
     package="IntOMICS")
OMICS_mod_res <- OMICS_module(omics = omics, PK = PK, annot = annot, 
     layers_def = layers_def, lm_METH = TRUE,
     TFtargs = TFtarg_mat, r_squared_thres = 0.3)
adjacency_matrix <- OMICS_mod_res$B_prior_mat
adjacency_matrix[,] <- 0
epsilon(net = adjacency_matrix, B_prior_mat = OMICS_mod_res$B_prior_mat)

}
