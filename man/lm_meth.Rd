% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OMICS_module_func.R
\name{lm_meth}
\alias{lm_meth}
\title{Linear regression GE~METH}
\usage{
lm_meth(ge_mat, meth_mat, gene, meth_probes, r_squared_thres, p_val_thres)
}
\arguments{
\item{ge_mat}{matrix of gene expression with samples in rows and
features in columns.}

\item{meth_mat}{matrix of DNA methylaton with samples in rows and
features in columns.}

\item{gene}{character vector with given node name.}

\item{meth_probes}{character vector methylation probes associated
with a gene.}

\item{r_squared_thres}{numeric vector to define the R^2 used as a threshold
of significance in linear regression if lm_METH=TRUE (default=0.3).}

\item{p_val_thres}{numeric vector to define the p-value used as a threshold
of significance in linear regression if lm_METH=TRUE (default=0.05).}
}
\value{
Character vector with methylation probes
}
\description{
\code{lm_meth} The linear regression model for a dependent variable GE and
explanatory variable METH. Returns METH with significant coefficient,
R^2 > threshold and R~Gaussian residuals.
}
\examples{
data(list=c("annot", "omics"), package="IntOMICS")
lm_meth(ge_mat = omics$ge, meth_mat = omics$meth, 
    gene = "ENTREZID:7482", meth_probes = annot[["ENTREZID:7482"]], 
    r_squared_thres = 0.3, p_val_thres = 0.05)

}
