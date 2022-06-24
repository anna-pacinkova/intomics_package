# IntOMICS: an R package for integrative analysis of multi-omics data to infer regulatory networks.

IntOMICS is an efficient integrative framework based on Bayesian networks.
IntOMICS systematically analyses gene expression (GE), DNA methylation (METH), copy number variation (CNV) and biological prior knowledge (B) to infer regulatory networks. 
IntOMICS complements the missing biological prior knowledge by so-called empirical biological knowledge (empB), estimated from the available experimental data. 
An automatically tuned MCMC algorithm ([Yang and Rosenthal, 2017](http://probability.ca/jeff/ftpdir/jinyoung1.pdf)) estimates model parameters and the empirical biological knowledge.
Conventional MCMC algorithm with additional Markov blanket resampling (MBR) step ([Su and Borsuk, 2016](https://jmlr.org/papers/volume17/su16a/su16a.pdf)) infers resulting regulatory network structure consisting of three types of nodes: GE nodes refer to gene expression levels, CNV nodes refer to associated copy number variations, and METH nodes refer to associated DNA methylation probe(s).
Regulatory networks derived from IntOMICS provide deeper insights into the complex flow of genetic information. 
IntOMICS is a powerful resource for exploratory systems biology and can provide valuable insights into the complex mechanisms of biological processes 
that has a vital role in personalised medicine.

IntOMICS takes as input: 
* gene expression matrix, 
* associated copy number variation matrix sampled from the same individuals, 
* associated DNA methylation matrix of beta-values sampled from the same individuals, and 
* the biological prior knowledge with information on known interactions among molecular features. 

The resulting regulatory network structure contains the edge weights $w_i$ representing the empirical frequency of given edge over samples of network structures from two independent MCMC simulations.

<p align="center">
  <img src="vignettes/figures/IntOMICS_framework_METH_empB_modules.png" width="400" height="680" alt="IntOMICS framework">
</p>
     
For further details about the IntOMICS algorithm, its performance and benchmark analysis, see manuscript [Pacinkova \& Popovici, 2022](https://assets.researchsquare.com/files/rs-1291540/v1_covered.pdf?c=1643735189). 


# Installation
```ruby
library(remotes)
install_github("anna-pacinkova/intomics_package")  
library(IntOMICS)
```

# Usage

This tutorial will show you how to use the IntOMICS package with a toy example.
The example dataset consisting of processed gene expression and DNA methylation (Illumina Infinium HumanMethylation450 BeadChip) is from [the TCGA data portal](https://portal.gdc.cancer.gov/): 30 colon cancer samples (COAD) with microsatellite instability (MSI). The copy number variation of the associated genes from TCGA-COAD samples were downloaded from [the Broad Institute GDAC Firehose](https://gdac.broadinstitute.org/).
We choose the set of 8 genes from the [KEGG Colorectal cancer pathway](https://www.genome.jp/pathway/hsa05210).


## Part 1: Input data
```ruby
library(IntOMICS)
# load required libraries
library(bnlearn)
library(bnstruct)
library(rlist)
library(matrixStats)
library(parallel)
library(foreach)
library(ggraph)
library(dplyr)
library(sm)
library(tibble)
library(RColorBrewer)
library(rstatix)
library(gridExtra)
library(bestNormalize)
library(igraph)
library(colorspace)
library(gplots)
```

IntOMICS framework takes as input:  
  
* gene expression matrix $GE$ ($m$ x $n_1$) with microarray intensities or RNA-seq count data transformed into a continuous domain ($m$ samples and $n_1$ features)
  
* associated copy number variation matrix $CNV$ ($m$ x $n_2$) with continuous segment mean values derived for each gene ($m$ samples and $n_2$ features, $n_2 \leq n_1$),
  
* associated DNA methylation matrix of beta-values $METH$ ($m$ samples and $n_3$ features, $m$ x $n_3$),
  
* data.frame including all known interactions between molecular features (information from public available databases such as KEGG ([Ogata et al., 1999](https://academic.oup.com/nar/article/27/1/29/1238108?login=true)) or REACTOME ([Wu \& Haw, 2017](https://link.springer.com/protocol/10.1007/978-1-4939-6783-4_11))). Any other source of prior knowledge can be used. (IntOMICS is designed to run even if prior knowledge is not available. However, we highly recommend to use it.)

* logical matrix with known transcription factors (TFs) and its targets (information from public available databases such as ENCODE ([ENCODE Project Consortium, 2004](https://pubmed.ncbi.nlm.nih.gov/15499007/)). Any other source can be used. (IntOMICS is designed to run even if any TF-target interaction is not available.)

All data matrices are sampled from the same individuals.  

Available omics data saved in a names list ```omics``` in the example TCGA COAD MSI dataset are gene expression (GE) of 7 genes + copy number variation (CNV) of 7 genes + beta value of 115 DNA methylation (METH) probes:
```ruby
omics$ge[1:5,1:5]
```
```diff
#>              ENTREZID:7482 ENTREZID:2535 ENTREZID:1857 ENTREZID:2932 ENTREZID:8312
#> TCGA-A6-5661    5.04929458    -0.8516422      7.140874      6.444206      2.913619
#> TCGA-AD-5900    1.35701918     1.4948136      5.588249      4.698387      5.800814
#> TCGA-CM-4743   -0.10079261    -0.4195763      5.701181      5.145362      5.086909
#> TCGA-G4-6586    0.07274478     0.5322971      5.537289      4.935165      5.591844
#> TCGA-F4-6570   -0.21499097     1.2462155      5.794873      5.420462      5.361082
```
These values correspond to normalised RNA-seq data. 
However, the user is not limited to this platform. Another assay, such as microarray data, can be used. The column names of ```omics$ge``` matrix must be entrez ID in the format ENTREZID:XXXX.

```ruby
omics$cnv[1:5,1:5]
```
```diff
#>              entrezid:7482 entrezid:2535 entrezid:1857 entrezid:2932 entrezid:8312
#> TCGA-A6-5661         0.000         0.005         0.000         0.000        -0.010
#> TCGA-AD-5900         0.000         0.001         0.009         0.002        -0.011
#> TCGA-CM-4743         0.001        -0.014        -0.003        -0.003        -0.016
#> TCGA-G4-6586         0.003         0.002        -0.001        -0.001         0.005
#> TCGA-F4-6570         0.002         0.018         0.007         0.007         0.002
```
These copy number values represent segment mean values equal to $log_2(\frac{copy-number}{2})$.
The column names of ```omics$cnv``` matrix must be entrez ID in the format entrezid:XXXX.
In the ```omics$cnv``` matrix, define only columns with available CNV data.

```ruby
omics$meth[1:5,1:5]
```
```diff
#>              cg08739433 cg11806528 cg14013390 cg08741842 cg10111629
#> TCGA-A6-5661  0.5277487  0.4595131  0.9353342  0.6480896  0.9147693
#> TCGA-AD-5900  0.5883163  0.5817354  0.8950541  0.5856380  0.7314535
#> TCGA-CM-4743  0.4499933  0.4935660  0.9525694  0.7287386  0.8809010
#> TCGA-G4-6586  0.5763992  0.5834095  0.9489423  0.7939009  0.8992803
#> TCGA-F4-6570  0.5128817  0.5656870  0.9377459  0.6671709  0.8519091
```
These values represent DNA methylation beta values. The column names of the ```omics$meth``` matrix are probe IDs.  

IntOMICS is designed to infer regulatory networks even if the copy number variation or DNA methylation data (or both) are not available. In such a case, omics must be a named list with a single element (matrix with gene expression). 


If methylation data are available, we have to provide an annotation:

```ruby
str(annot)
```
```diff
#> List of 5
#> $ ENTREZID:2932: chr [1:2] "cg14479617" "cg12588208"
#> $ ENTREZID:8312: chr [1:2] "cg27308245" "cg26525629"
#> $ ENTREZID:7482: chr "cg04571584"
#> $ ENTREZID:1499: chr [1:3] "cg04180460" "cg06626556" "cg02247160"
#> $ ENTREZID:2535: chr "cg09339219"
```
```annot``` is a named list. Each component of the list is a character vector and corresponds to probe IDs associated with a given gene. Names of the annot must be again in the format ENTREZID:XXXX.  

To generate comprehensive figures with gene IDs, we need to provide a gene annotation table:
```ruby
gene_annot
```
```diff
#>         entrezID gene_symbol   alias_KEGG
#> 14 ENTREZID:7482       WNT2B          Wnt
#> 21 ENTREZID:2535        FZD2     Frizzled
#> 34 ENTREZID:1857        DVL3          DVL
#> 35 ENTREZID:2932       GSK3B     GSK3beta
#> 38 ENTREZID:8312       AXIN1         AXIN
#> 40 ENTREZID:1499      CTNNB1 beta-catenin
#> 36  ENTREZID:324         APC          APC
#> 43 ENTREZID:6934      TCF7L2      TCF/LEF
```
```gene_annot``` is Gene ID conversion table with "entrezID" and "gene_symbol" column names. Gene symbols are used for the final regulatory network visualisation.  

And finally, the prior knowledge from any source chosen by the user:
```ruby
PK
```
```diff
#>      src_entrez   dest_entrez edge_type
#> 1 ENTREZID:7482 ENTREZID:2535   present
#> 2 ENTREZID:2535 ENTREZID:1857   present
#> 3 ENTREZID:1857 ENTREZID:2932   present
#> 4  ENTREZID:324 ENTREZID:1499   present
#> 5  ENTREZID:324 ENTREZID:1857   present
#> 6 ENTREZID:8312  ENTREZID:324   present
```
```PK``` is the data.frame with biological prior knowledge about known interactions between features. Column names are "src_entrez" (the parent node), "dest_entrez" (the child node) and "edge_type" (the prior knowledge about the direct interaction between parent and child node; the allowed values are "present" or "missing").

```ruby
TFtarg_mat[1:5,1:5]
```
```diff
#>                 ENTREZID:1820 ENTREZID:466 ENTREZID:1386 ENTREZID:467 ENTREZID:571
#> ENTREZID:1                  1            0             0            0            0
#> ENTREZID:503538             1            0             0            0            0
#> ENTREZID:29974              1            0             0            0            0
#> ENTREZID:2                  1            0             0            0            0
#> ENTREZID:144568             0            0             1            0            0
```
```TFtarg_mat``` is the logical matrix with known transcription factors (TFs) and its targets. TFs are listed in columns, corresponding targets are listed in rows.


## Part 2: Data preprocessing

The first step is to define the biological prior matrix and estimate the upper bound of the partition function needed to define the prior distribution of network structures.
We also need to define all possible parent set configurations for each node.  For each parent set configuration, we compute the energy (needed to define the prior distribution of network structures) and the BGe score (needed to determine the posterior probability of network structures).
These functionalities are available through the ```OMICS_module``` function.  
We can use linear regression to filter irrelevant DNA methylation probes. We set the parameter "lm_METH = TRUE".
We can also specify the threshold for the R^2 to choose DNA methylation probes with significant coefficient using argument "r_squared_thres" (default = 0.3), or p-value using "p_val_thres" (default = 0.05).  
There are several other arguments: "woPKGE_belief" (default = 0.5) refers to the belief concerning GE-GE interactions without prior knowledge, "nonGE_belief" (default = 0.5) refers to the belief concerning the belief concerning interactions of features except GE (e.g. CNV-GE, METH-GE), "TFBS_belief" refers to the belief concerning the TF and its target interaction (default = 0.75).
Note that all interactions with belief equal to "woPKGE_belief" in biological prior knowledge will be updated in empirical biological knowledge.

```ruby
OMICS_mod_res <- OMICS_module(omics = omics, 
                              PK = PK, 
                              layers_def = layers_def, 
                              TFtargs = TFtarg_mat,
                              annot = annot, 
                              r_squared_thres = 0.1, 
                              lm_METH = TRUE)
```

This function returns several outputs:
```ruby
names(OMICS_mod_res)
```

```diff
#> "pf_UB_BGe_pre"       "B_prior_mat"         "annot"               "omics"               "layers_def"          "omics_meth_original"
```

1. ```OMICS_mod_res$pf_UB_BGe_pre``` is a list that contains:
- ```OMICS_mod_res$pf_UB_BGe_pre$partition_func_UB``` the upper bound of the partition function for hyperparameter 
$\beta = 0$,
- ```OMICS_mod_res$pf_UB_BGe_pre$parents_set_combinations``` all possible parent set configuration for given node,
- ```OMICS_mod_res$pf_UB_BGe_pre$energy_all_configs_node``` energy for given parent set configurations,
- ```OMICS_mod_res$pf_UB_BGe_pre$BGe_score_all_configs_node``` BGe score for given parent set configurations.  
2. ```OMICS_mod_res$B_prior_mat``` is a biological prior matrix.
3. ```OMICS_mod_res$annot``` contains DNA methylation probes that passed the filter.  
4. ```OMICS_mod_res$omics``` is a list with gene expression, copy number variation and normalised methylation data (possibly filtered if we use "lm_METH = TRUE").
5. ```OMICS_mod_res$omics_meth_original``` the original methylation data.


## Part 3: MCMC simulation

Now, we can use the automatically tuned MCMC algorithm (Yang and Rosenthal, 2017) to estimate model parameters and empirical biological knowledge and the conventional MCMC algorithm with additional Markov blanket resampling step (Su and Borsuk, 2016) to infer regulatory network structure consisting of three types of nodes: GE, CNV and METH nodes. 
This step can be time-consuming (you can skip it and use the pre-computed result -> R object ```BN_mod_res```).
```ruby
BN_mod_res <- BN_module(burn_in = 100000, 
                        thin = 500, 
                        seed1 = 101,
                        seed2 = 505,
                        OMICS_mod_res = OMICS_mod_res,
                        minseglen = 50000)
```
There are two optional arguments: "len" specifies the initial width of the sampling interval for hyperparameter $\beta$. However, this parameter will be tuned during the adaptive phases of the MCMC algorithm. "prob_mbr" specifies the probability of the MBR step (default = TRUE). We strongly recommend to use the default setting (for further details on how this argument affects MCMC scheme results, see [Su and Borsuk, 2016](https://jmlr.org/papers/volume17/su16a/su16a.pdf)).

Let's check the outputs of ```BN_module``` function:
```ruby
names(BN_mod_res)
```
```diff
#> "B_prior_mat_weighted" "sampling.phase_res"   "beta_tuning"
```
1. ```BN_mod_res$B_prior_mat_weighted``` is empirical biological knowledge matrix. Interactions from the biological prior knowledge and TFs-target interactions are constant (if "TFBS_belief" is not equal to "woPKGE_belief").
2. ```BN_mod_res$sampling.phase_res``` is a named list with results from the last sampling phase of the algorithm:
- ```BN_mod_res$sampling.phase_res$rms``` trace of root mean square used for c_rms measure to evaluate the convergence of MCMC simulation,
- ```BN_mod_res$sampling.phase_res$mcmc_sim_part_res``` contains CPDAGs of two independent MCMC simulations (for details see [Pacinkova \& Popovici, 2022](https://assets.researchsquare.com/files/rs-1291540/v1_covered.pdf?c=1643735189)).
4. ```BN_mod_res$beta_tuning``` is a named list for each iteration of adaptive phases that contains:
- ```BN_mod_res$beta_tuning$value``` trace of hyperparameter beta tuning,
- ```BN_mod_res$beta_tuning$len``` trace of width of the sampling interval for hyperparameter $\beta$.

## Part 4: MCMC diagnostics

Trace plots provide an important tool for assessing mixing of a Markov chain and should be inspected carefully.
Now we have to create a directory to store trace plots. 
Once it is created, we can run the trace_plots function, which generates:  
  
1. beta_values.svg: trace plot of beta values (we want to explore the sample space many times and avoid flat bits - the chain stays in the same state for too long; in this case, beta value fluctuates around single value, so the Markov chain is mixing well)  
<p align="center">
<img src="vignettes/figures/beta_values.svg" width="400" height="400">
</p>

2. post_prob_edges.svg: consistency of edges posterior probabilities in two independent MCMC simulations (scatter plot of the edge weights confidence using two independent MCMC runs; the convergence is determined by the spread of the points around the y=x line; in this case, the edge weights seems to be consistent in two independent simulations)   
<p align="center">
<img src="vignettes/figures/post_prob_edges.svg" width="400" height="400">
</p>

3. convergence_RMS.svg: the c<sub>rms</sub> strength for the convergence evaluation (summarizes the spread of the points around the line y=x in post_prob_edges.svg, for details see ([Agostinho et al., 2015](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0734-6) and [Pacinkova \& Popovici, 2022](https://assets.researchsquare.com/files/rs-1291540/v1_covered.pdf?c=1643735189))).  
<p align="center">
<img src="vignettes/figures/convergence_RMS.svg" width="400" height="400">
</p>

The parameter "edge_freq_thres" determines the quantile of all edge weights used to filter only reliable edges (default = NULL, all edges will be considered as present). For illustration, we use quite low edge weights filter to capture both CNV-GE and METH-GE interactions. We recommend to use some edge weights filtering, such as 0.75 quantile of all edge weights in the resulting networks using "edge_freq_thres = 0.75".
The parameter "gene_ID" determines the IDs used in the final network. There are two options: "gene_symbol" (default) or "entrezID".
If you have changed the default value of the "TFBS_belief" argument in the ```OMICS_module``` function, you have to use the same argument in ```trace_plots``` function.
```ruby
res_weighted <- trace_plots(mcmc_res = BN_mod_res, 
                            figures_dir = "figures/MSI/", 
                            burn_in = 100000, 
                            thin = 500, 
                            gene_annot = gene_annot, 
                            PK = PK, 
                            OMICS_mod_res = OMICS_mod_res,
                            gene_ID = "gene_symbol", 
                            edge_freq_thres = 0.4,
                            TFtargs = TFtarg_mat)
```

## Part 5: IntOMICS resulting network structure

We can plot the resulting regulatory network inferred by IntOMICS using ```ggraph``` function from the [ggraph](https://github.com/thomasp85/ggraph) package:

```ruby
ggraph(res_weighted$net_weighted, layout = 'dh') + 
  geom_edge_link(aes(end_cap = circle(node2.degree + 7, "pt"), edge_color = edge, label = weight),
                 arrow = arrow(angle = 20, length = unit(0.1, "inches"),
                              ends = "last", type = "closed"))+
  geom_node_point(aes(color = factor(color)), size = 10) +
  scale_colour_manual(values = res_weighted$node_palette, guide = "none")+
  geom_node_text(aes(label = label),family="serif")
```
<p align="center">
<img src="vignettes/figures/fin_net.svg" width="600" height="600">
</p>

Edges highlighted in blue are known from the biological prior knowledge. 
The edge labels reflects its empirical frequency over the final set of CPDAGs.
GE node names are in upper case, CNV node names are in lower case, METH node names are the same as DNA methylation probe names in ```omics$meth``` matrix.  

Node colours legend:
```ruby
legend_custom(net = res_weighted)
```
<p align="center">
<img src="vignettes/figures/fin_net_legend.svg" width="400" height="400">
</p>
Node colour scales are given by GE/CNV/METH values of all features from the corresponding input data matrix.  

We can also change the edge labels to inspect the empirical prior knowledge inferred by IntOMICS using the argument "edge_weights = empB" (default = "MCMC_freq"):
```ruby
res_weighted <- trace_plots(mcmc_res = BN_mod_res, 
                            figures_dir = "figures/MSI/", 
                            burn_in = 100000, 
                            thin = 500, 
                            gene_annot = gene_annot, 
                            PK = PK, 
                            OMICS_mod_res = OMICS_mod_res, 
                            gene_ID = "gene_symbol", 
                            edge_freq_thres = 0.4, 
                            edge_weights = "empB",
                            TFtargs = TFtarg_mat)

ggraph(res_weighted$net_weighted, layout = 'dh') + 
  geom_edge_link(aes(end_cap = circle(node2.degree + 7, "pt"), edge_color = edge, label = weight),
                 arrow = arrow(angle = 20, length = unit(0.1, "inches"),
                              ends = "last", type = "closed"))+
  geom_node_point(aes(color = factor(color)), size = 10) +
  scale_colour_manual(values = res_weighted$node_palette, guide = "none")+
  geom_node_text(aes(label = label),family="serif")
```
<p align="center">
<img src="vignettes/figures/fin_net_empB.svg" width="600" height="600">
</p>

Function ```empB_heatmap``` can be used to check the difference between empirical biological knowledge and biological prior knowledge of GE-GE interactions:
```ruby
empB_heatmap(mcmc_res = BN_mod_res, 
             OMICS_mod_res = OMICS_mod_res, 
             gene_annot = gene_annot, 
             TFtargs = TFtarg_mat)
```
<p align="center">
<img src="vignettes/figures/empB_matrix.svg" width="450" height="450">
</p>
Interactions with constant biological knowledge are highlighted in gray.

Interesting could be also density of the edge weights inferred by IntOMICS, that can be plot using ```ggplot``` function. First of all, we have to use the ```trace_plots``` function without the edge weights filtering:
```ruby
res_weighted <- trace_plots(mcmc_res = BN_mod_res, 
                            figures_dir = "figures/MSI/", 
                            burn_in = 100000, 
                            thin = 500, 
                            gene_annot = gene_annot, 
                            PK = PK, 
                            OMICS_mod_res = OMICS_mod_res, 
                            gene_ID = "gene_symbol",  
                            TFtargs = TFtarg_mat)
                            
df <- data.frame(edge_weight = as.numeric(res_weighted$edge_list[,"weight"]))
q3 <- quantile(df$edge_weight,0.75)
p <- ggplot(df, aes(x=edge_weight, y = ..scaled..)) +
  geom_density(col="dodgerblue")
p+geom_vline(xintercept=q3, size=0.5, color="black", linetype = "dashed")
```
<p align="center">
<img src="vignettes/figures/edge_weights.svg" width="300" height="300">
</p>

If you find a bug or have a comment let us know, please, via an e-mail: ana.pacinkova@gmail.com


