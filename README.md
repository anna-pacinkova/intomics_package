# IntOMICS: an R package for integrative analysis of multi-omics data to infer regulatory networks.

IntOMICS is an efficient integrative framework based on Bayesian networks.
IntOMICS systematically analyses gene expression, DNA methylation, copy number variation and biological prior knowledge to infer regulatory networks. 
IntOMICS complements the missing biological prior knowledge by so-called empirical biological knowledge, estimated from the available experimental data. 
An automatically tuned MCMC algorithm ([Yang and Rosenthal, 2017](http://probability.ca/jeff/ftpdir/jinyoung1.pdf)) estimates model parameters and the empirical biological knowledge.
Conventional MCMC algorithm with additional Markov blanket resampling step ([Su and Borsuk, 2016](https://jmlr.org/papers/volume17/su16a/su16a.pdf)) infers resulting regulatory network structure consisting of three types of nodes: GE nodes refer to gene expression levels, CNV nodes refer to associated copy number variations, and METH nodes refer to associated DNA methylation probe(s).
Regulatory networks derived from IntOMICS provide deeper insights into the complex flow of genetic information. 
IntOMICS is a powerful resource for exploratory systems biology and can provide valuable insights into the complex mechanisms of biological processes 
that has a vital role in personalised medicine.

IntOMICS takes as input (i) gene expression matrix GE with m samples and n1 genes, (ii) associated copy number variation matrix CNV (m x n2), 
(iii) associated DNA methylation matrix of beta-values METH (m x n3) sampled from the same individuals, and (iv) the biological prior knowledge 
with information on known interactions among molecular features. 
An automatically tuned MCMC algorithm (Yang and Rosenthal, 2009) estimates parameters and empirical biological knowledge. 
Conventional MCMC algorithm with additional Markov blanket resampling step (Su and Borsuk, 2016) is used to infer resulting regulatory network structure 
consisting of three types of nodes: GE nodes refer to gene expression levels, CNV nodes refer to copy number variations, and METH nodes refer to DNA methylation. 
The resulting regulatory network structure contains the edge weights $w_i$ representing the empirical frequency of given edge over samples of network structures from two independent MCMC simulations.

<p align="center">
  <img src="vignettes/figures/IntOMICS_framework_METH_empB.png" width="460" height="300" alt="IntOMICS framework">
</p>
     
For further details, see manuscript [Pacinkova \& Popovici, 2022](https://assets.researchsquare.com/files/rs-1291540/v1_covered.pdf?c=1643735189). 


## Installation
<pre><code>install_github("anna-pacinkova/intomics_package")  
library(devtools)</code></pre>


Comprehensive tutorial: vignettes/IntOMICS_vignette.Rmd
