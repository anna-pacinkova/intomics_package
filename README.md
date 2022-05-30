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
  <img src="vignettes/figures/IntOMICS_framework_METH_empB.png" width="300" height="450" alt="IntOMICS framework">
</p>
     
For further details, see manuscript [Pacinkova \& Popovici, 2022](https://assets.researchsquare.com/files/rs-1291540/v1_covered.pdf?c=1643735189). 


# Installation

<pre><code>
install_github("anna-pacinkova/intomics_package")  
library(devtools)
</code></pre>

# Usage

This tutorial will show you how to use the IntOMICS package with a toy example.
The example dataset is from [the TCGA data portal](https://portal.gdc.cancer.gov/): 30 colon cancer samples (COAD) with microsatellite instability (MSI).
We choose the set of 7 genes from the [KEGG Colorectal cancer pathway](https://www.genome.jp/pathway/hsa05210).


## Part 1: Input data loading and preprocessing

<pre><code>
library(knitr)
library(IntOMICS)
library(bestNormalize)
library(foreach)
library(bnlearn)
library(RColorBrewer)
library(png) 
library(matrixStats)
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(bnstruct))
</code></pre>

IntOMICS framework takes as input:  
  
* gene expression matrix $GE$ ($m$ x $n_1$) with microarray intensities or RNA-seq count data transformed into a continuous domain ($m$ samples and $n_1$ features)
  
* associated copy number variation matrix $CNV$ ($m$ x $n_2$) with continuous segment mean values derived for each gene ($n_2 \leq n_1$),
  
* associated DNA methylation matrix of beta-values $METH$ ($m$ x $n_3$),
  
* data.frame including all known interactions between molecular features (information from public available databases such as KEGG (Ogata et al., 1999) or REACTOME (Wu \& Haw, 2017)). However, any other source of prior knowledge can be used.  
  
All data matrices are sampled from the same individuals.  

Available omics data in the example TCGA COAD MSI dataset are gene expression (GE) of 7 genes + copy number variation (CNV) of 7 genes + beta value of 115 DNA methylation (METH) probes:

Comprehensive tutorial: vignettes/IntOMICS_vignette.Rmd
