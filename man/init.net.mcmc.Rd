% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NET_manipulation_func.R
\name{init.net.mcmc}
\alias{init.net.mcmc}
\title{Random initial network}
\usage{
init.net.mcmc(omics, layers_def, B_prior_mat)
}
\arguments{
\item{omics}{named list containing the gene expression (possibly copy number
variation and methylation data). Each component of the list is a matrix
with samples in rows and features in columns.}

\item{layers_def}{data.frame containing the modality ID, corresponding layer
in BN and maximal number of parents from given layer to GE nodes.}

\item{B_prior_mat}{a biological prior matrix.}
}
\value{
List of 2 elements: random adjacency network and empty network
}
\description{
\code{init.net.mcmc} This function is used to sample random initial network.
The edges are sampled only between GE nodes.
}
\examples{
data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics"),
package="IntOMICS")
B <- B_prior_mat(omics = omics, PK = PK, annot = annot, lm_METH = TRUE,
     layers_def = layers_def, r_squared_thres = 0.3, p_val_thres = 0.05,
     TFtargs = TFtarg_mat, TFBS_belief = 0.75, nonGE_belief = 0.5, 
     woPKGE_belief = 0.5)
init.net.mcmc(omics = B$omics, layers_def = layers_def, 
     B_prior_mat = B$B_prior_mat)

}
