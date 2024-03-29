% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NET_manipulation_func.R
\name{sample.chain}
\alias{sample.chain}
\title{Random initial network edge generation}
\usage{
sample.chain(empty_net, omics_ge)
}
\arguments{
\item{empty_net}{adjacency matrix of an empty network/graph
(all values are 0).}

\item{omics_ge}{matrix with gene expression data (samples in rows and
features in columns).}
}
\value{
BN object with conditional probabilities
}
\description{
\code{sample.chain} This function is used to sample random initial network.
The edges are sampled only between GE nodes.
}
\examples{
data(list=c("PK", "TFtarg_mat", "annot", "layers_def", "omics"),
package="IntOMICS")
B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, 
     annot = annot, lm_METH = TRUE, r_squared_thres = 0.3,
     p_val_thres = 0.05, TFtargs = TFtarg_mat, TFBS_belief = 0.75, 
     nonGE_belief = 0.5, woPKGE_belief = 0.5)
empty.net <- matrix(0, nrow = sum(mapply(ncol,B$omics)), ncol =
sum(mapply(ncol,B$omics)), dimnames = list(unlist(mapply(colnames,B$omics)),
unlist(mapply(colnames,B$omics))))
sample.chain(empty_net = empty.net, 
     omics_ge = B$omics[[layers_def$omics[1]]])

}
