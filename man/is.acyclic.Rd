% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NET_manipulation_func.R
\name{is.acyclic}
\alias{is.acyclic}
\title{Acyclic network identification.}
\usage{
is.acyclic(g)
}
\arguments{
\item{g}{adajcency matrix of given network/graph.}
}
\value{
boolean of length 1
}
\description{
\code{is.acyclic} This function is from bnstruct R package. Check if the directed
graph is acyclic.
}
\examples{
is.acyclic(matrix(c(1,rep(0,20)), nrow=3))

}
