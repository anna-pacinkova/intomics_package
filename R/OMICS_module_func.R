#' OMICS_module
#' @description
#' `OMICS_module` data preprocessing + B_prior_mat definition + partition function upper bound estimation + all possible parent sets per node definition + BGe score computation for all possible parent sets
#' @export
OMICS_module <- function(omics, PK, layers_def, annot=NULL, r_squared_thres = 0.3, lm_METH = TRUE)
{
  
  # sort according to the layer number
  layers_def <- layers_def[order(layers_def$layer, decreasing = TRUE),]
  omics <- omics[layers_def$omics[order(layers_def$layer, decreasing = TRUE)]]
  
  # The biological prior matrix definition:
  B <- B_prior_mat(omics = omics, PK = PK, layers_def = layers_def, annot = annot, r_squared_thres = r_squared_thres, lm_METH = lm_METH)
  
  # The upper bound of the partition function with beta_min = 0 (the resulting partition_func_UB is log transformed!)
  pf_UB_res <- pf_UB_est(omics = B$omics, layers_def = layers_def, B_prior_mat = B$B_prior_mat, annot = B$annot)
  
  return(list(pf_UB_BGe_pre = pf_UB_res,
              B_prior_mat = B$B_prior_mat, 
              annot = B$annot, 
              omics = B$omics,
              layers_def = layers_def,
              omics_meth_original = B$omics_meth_original))
}

#' Biological prior matrix definition
#' @export
B_prior_mat <- function(omics, PK, layers_def, annot, r_squared_thres, lm_METH)
{
  
  if(any(regexpr("ENTREZID:",colnames(omics[[layers_def$omics[1]]]))<0))
  {
  	stop("Gene names in GE matrix are not in the correct form. Please, use ENTREZID:XXXX.")
  } # end if(regexpr("ENTREZID:",colnames(omics[[layers_def$omics[1]]]))<0)
  
  if(!all(sort(names(omics))==sort(layers_def$omics)))
  {
    stop("Names of the list 'omics' does not match omics names in the 'layers_def'.")
  } # end if(!all(order(names(omics))==order(layers_def$omics)))
  
  # the first value indicates the belief of present nodes, the second value indicates the belief of absent nodes
  pk_belief <- c(1,0)

  ## interactions from the last layer
  features <- colnames(omics[[layers_def$omics[1]]])

  # 0.5 denotes that we have not any information about the presence or absence of an edge
  B_layer_max <- matrix(0.5, ncol=length(features), nrow=length(features), dimnames=list(features,features))
  
  # B_ij = 1 we have prior evidence that there is a directed edge pointing from node i to node j
  PK_present <- PK[PK$edge_type=="present",]
  PK_absent <- PK[PK$edge_type=="absent",]
    
  for(i in c(1:nrow(B_layer_max)))
  {
    
    if (sum(PK_present$src_entrez==rownames(B_layer_max)[i]) != 0)
    {
      B_layer_max[rownames(B_layer_max)[i],intersect(PK_present[PK_present$src_entrez==rownames(B_layer_max)[i],"dest_entrez"],rownames(B_layer_max))] <- pk_belief[1]
    } # end if sum...
    
    if (sum(PK_absent$src_entrez==rownames(B_layer_max)[i]) != 0)
    {
      B_layer_max[rownames(B_layer_max)[i],intersect(PK_absent[PK_absent$src_entrez==rownames(B_layer_max)[i],"dest_entrez"],rownames(B_layer_max))] <- pk_belief[2]
    } # end if sum...
    
  } # end for i...
  B <- B_layer_max
  diag(B) <- 0
  new_annot <- list()
  
  omics_meth_original <- matrix(nrow = 0, ncol = 0)
  ## interactions from other layers
  if(length(layers_def$omics)>1)
  {
    for(j in c(2:length(layers_def$omics)))
    {
      features_lower <- colnames(omics[[layers_def$omics[j]]])
      if(any(regexpr("entrezid:",colnames(omics[[layers_def$omics[j]]]))>0)==TRUE)
      {
        B_layer_lower <- matrix(0, ncol=ncol(B),
                                nrow=length(features_lower), 
                                dimnames=list(features_lower[match(colnames(B), toupper(features_lower), nomatch=0)],colnames(B)))
        # CNV data: all possible interactions are set to 0.5
        diag(B_layer_lower) <- 0.5
        B_layer1 <- matrix(0, ncol=length(features_lower), nrow=nrow(B)+length(features_lower), 
                           dimnames=list(c(rownames(B), features_lower), rownames(B_layer_lower)))
        B <- cbind(rbind(B,B_layer_lower),B_layer1)
      } else {
        annot <- annot[intersect(names(annot),colnames(omics[[layers_def$omics[1]]]))]
        annot <- lapply(annot, FUN=function(s) intersect(s,colnames(omics[[layers_def$omics[j]]])))
        
        omics_meth_original <- omics[[layers_def$omics[j]]]
        # Methylation data transformation:
        omics[[layers_def$omics[j]]] <- apply(omics[[layers_def$omics[j]]], 2, FUN =function(column) orderNorm(column)$x.t)
        
        if(lm_METH)
        {
          new_annot <- lapply(seq_along(annot), function(list) lm_meth(omics[[layers_def$omics[1]]], omics[[layers_def$omics[j]]], names(annot)[[list]], annot[[list]], r_squared_thres))
          names(new_annot) <- names(annot)
        } else {
          new_annot <- annot
        } # end if(lm_METH)
        
        if(sum(!mapply(is.null,new_annot))!=0)
        {
         new_annot <- new_annot[!mapply(is.null,new_annot)]
         if(any(mapply(length,new_annot)>10)){print("WARNING: There are more than 10 methylation probes per gene. Consider some methylation probes filtering due to computational complexity.")}
         
         features_lower <- unlist(new_annot)
         B_layer_lower <- matrix(0, ncol=ncol(B), 
                                 nrow=length(features_lower), 
                                 dimnames=list(features_lower,colnames(B)))
         # METH data: all possible interactions from annot are set to 0.5
         for(a in c(1:length(new_annot)))
         {
           B_layer_lower[intersect(features_lower,new_annot[[a]]),names(new_annot)[a]] <- 0.5
         } # end for a
         B_layer1 <- matrix(0, ncol=length(features_lower), nrow=nrow(B)+length(features_lower), 
                            dimnames=list(c(rownames(B), features_lower), rownames(B_layer_lower)))
         B <- cbind(rbind(B,B_layer_lower),B_layer1)
         omics[[layers_def$omics[j]]] <- omics[[layers_def$omics[j]]][,unlist(new_annot)]
        } else {
          new_annot <- list()
          omics[[layers_def$omics[j]]] <- matrix(nrow = 0, ncol = 0)
        } # end if else (sum(!mapply(is.null,new_annot))==0)
        # there are some parents of METH_candidates[i] from the IntOMICS network and the METH value has more than 1 factor level
      } # end if else (any(regexpr("entrezid:",colnames(omics[[layers_def$omics[j]]]))>0)==TRUE)
    } # end for j
  } # end if(length(layers_def$omics)>1)
  return(list(B_prior_mat = B, annot = new_annot, omics = omics, omics_meth_original = omics_meth_original))
}

#' BGe score
#' @description
#' `DAGcorescore` The log of the BGe score simplified as much as possible. This function is from BiDAG package.
#' @export
DAGcorescore <- function(j,parentnodes,n,param) {
  
    TN <- param$TN
    awpN <- param$awpN
    scoreconstvec <- param$scoreconstvec
    
    lp <- length(parentnodes) #number of parents
    awpNd2 <- (awpN-n+lp+1)/2
    A <- TN[j,j]

    switch(as.character(lp),
           "0"={# just a single term if no parents
             corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
           },
           
           "1" = {# no need for matrices
             D <- TN[parentnodes,parentnodes]
             logdetD <- log(D)
             B <- TN[j,parentnodes]
             logdetpart2 <- log(A-B^2/D)
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           },
           
           "2" = {
             D <- TN[parentnodes,parentnodes]
             
             dettwobytwo <- function(D) {
               D[1,1]*D[2,2]-D[1,2]*D[2,1]
             }
             
             detD <- dettwobytwo(D)
             logdetD <- log(detD)
             B <- TN[j,parentnodes]
             logdetpart2 <- log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           },
           
           {# otherwise we use cholesky decomposition to perform both
             D<-as.matrix(TN[parentnodes,parentnodes])
             choltemp<-chol(D)
             logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
             B<-TN[j,parentnodes]
             logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           })
  
  return(corescore)
}

#' Node energy function
#' @description
#' `energy_function_node_specific`  For each node returns its energy over all parent set configurations, the empty parent set is included.
#' @export
energy_function_node_specific <- function(all_parents_config, B_prior_mat, int_node)
{
  epsilon <- apply(all_parents_config,2,FUN=function(z)
    if(is.na(z[1]))
    {
      sum(B_prior_mat[,int_node])
    } else {
      sum(1-B_prior_mat[z,int_node]) + sum(B_prior_mat[-match(z,rownames(B_prior_mat)),int_node])
    })
  return(epsilon)
}

#' Flatten a list of lists
#' @description
#' `flattenlist` https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists
#' @export
flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  } # end if else
}

#' Linear regression GE~METH
#' @description
#' `lm_meth` The linear regression model for a dependent variable GE and explanatory variable METH. Returns METH with significant coefficient, R^2 > threshold and R~Gaussian residuals.
#' @export
lm_meth <- function(ge_mat, meth_mat, gene, meth_probes, r_squared_thres)
{
  meth_probes_sig <- c()
  if(length(meth_probes)>0)
  {
    for(f in c(1:length(meth_probes)))
    {
      res <- lm(ge_mat[,gene] ~ meth_mat[,meth_probes[f]])
      cond1 <- summary(res)$coefficients[2,"Pr(>|t|)"] < 0.05
      cond2 <- summary(res)$r.squared > r_squared_thres
      cond3 <- shapiro.test(summary(res)$resid)$p.value > 0.1
      if(cond1 & cond2 & cond3)
      {
        meth_probes_sig <- c(meth_probes_sig, meth_probes[f])
      } # end if(cond1 & cond2 & cond3)
    } # end for f
  } # end if(length(meth_probes)>0)
  
  return(meth_probes_sig)
}

#' Partition function upper bound
#' @description
#' `pf_UB_est` Partition function upper bound estimation with beta = 0. For each node returns energy over all possible parent set configurations and BGe score.
#' @export
pf_UB_est <- function(omics, B_prior_mat, layers_def, annot)
{
  ## PART I: the upper bound of the partition function
  # define all parents configurations
  comb_all <- foreach (i=1:ncol(omics[[layers_def$omics[1]]])) %do% {
    int_node <- colnames(omics[[layers_def$omics[1]]])[i]
    comb_some <- list()
    
    ## GE all parents set combinations
    for(rep in c(1:layers_def$fan_in_ge[1]))
    {
      comb_some[[rep]] <- combn(colnames(omics[[layers_def$omics[1]]])[-i],rep)
    } # end for rep
    
    # add empty parent sets
    comb_some[[layers_def$fan_in_ge[1]+1]] <- matrix(NA,1,1)
    
    if(length(layers_def$omics)>1)
    {
      ## add CNV for each parents set combination including empty parent sets
      modalities <- layers_def$omics[-1]
      if(any(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0)==TRUE, omics[modalities])) & tolower(int_node) %in% unlist(lapply(omics[modalities],colnames)))
      {
        comb_some[seq(length(comb_some)+1,length.out = length(comb_some))] <- lapply(comb_some,FUN=function(list) rbind(list,tolower(int_node)))
      } # end if(any(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0)==TRUE, omics[modalities])))
      
      ## is there METH?
      if(any(mapply(omics,FUN=function(list) any(regexpr("entrezid:",colnames(list), ignore.case = TRUE)<0))) & length(annot[[int_node]])>0)
      {
        modality <- names(which(mapply(omics,FUN=function(list) any(regexpr("entrezid:",colnames(list), ignore.case = TRUE)<0))==TRUE))
        max_fan_in <- max(layers_def$fan_in_ge[layers_def$omics==modality],length(annot[[int_node]]), na.rm = TRUE)
        
        ## only METH all parents set combinations
        comb_some_meth <- list()
        for(rep in c(1:max_fan_in))
          {
            comb_some_meth[[rep]] <- combn(annot[[int_node]],rep)
        } # end for rep

        ## add METH for each parents set combination including empty parent sets
        comb_some_meth_add <- list()
        for(meth_pr in c(1:length(comb_some_meth)))
        {
          comb_some_meth_add <- c(comb_some_meth_add, foreach(a=1:ncol(comb_some_meth[[meth_pr]])) %do% {
            lapply(comb_some, FUN=function(par_def) apply(par_def,2,FUN=function(par_def_col) c(par_def_col, comb_some_meth[[meth_pr]][,a])))
          }) # end foreach(a=1:ncol(comb_some_meth[[meth_pr]])) %do%
        } # end for meth_pr
        
        comb_some_meth_add <- flattenlist(comb_some_meth_add)
        comb_some <- c(comb_some, comb_some_meth_add)
      } # if if(any(mapply(omics,FUN=function(list) any(regexpr("entrezid:",colnames(list), ignore.case = TRUE)<0))))
    } # end if(length(layers_def$omics)>1)
    
    # now remove NAs (delete an empty parent set too)
    comb_some <- lapply(comb_some,na.omit)
    # add again an empty parent set
    comb_some[[1]] <- cbind(comb_some[[1]], NA)
    comb_some <- comb_some[mapply(comb_some,FUN=function(x) nrow(x))!=0]
    
    parents_config <- list()
    for(l in c(1:max(mapply(comb_some,FUN=function(x) nrow(x)))))
    {
      parents_config[[l]] <- do.call(cbind, comb_some[mapply(comb_some,FUN=function(x) nrow(x))==l])
    } # end for l
    parents_config
  } # end foreach (i=1:ncol(omics[[layers_def$omics[1]]]))
  names(comb_all) <- colnames(omics[[layers_def$omics[1]]])
  
  if(length(layers_def$omics)>1)
  {
    ## add also CNV/METH nodes and their empty parent set configurations
    comb_all_others <- vector(mode = "list", length = sum(mapply(ncol,omics[setdiff(layers_def$omics,layers_def$omics[1])])))
    comb_all_others <- lapply(comb_all_others, FUN=function(list) list <- matrix(NA))
    names(comb_all_others) <- unlist(mapply(colnames,omics[setdiff(layers_def$omics,layers_def$omics[1])]))
    comb_all <- c(comb_all, comb_all_others)
  } # end if(length(layers_def$omics)>1)
  
  
  energy_all_configs_node <- list()
  ## GE
  for(i in c(1:ncol(omics[[layers_def$omics[1]]])))
  {
    energy_all_configs_node[[i]] <- unlist(lapply(comb_all[[i]],FUN=function(list) energy_function_node_specific(list, B_prior_mat, names(comb_all)[i])))
  } # end for i
  ## CNV/METH
  # because we force the CNV/METH nodes to be orphaned, we can only sum the values in the B_prior_mat (now it is 0 because of 0 pk_belief)
  if(length(layers_def$omics)>1)
  {
    for(i in c((ncol(omics[[layers_def$omics[1]]])+1):(ncol(omics[[layers_def$omics[1]]])+sum(mapply(ncol,omics[setdiff(layers_def$omics,layers_def$omics[1])])))))
    {
      energy_all_configs_node[[i]] <- sum(B_prior_mat[,names(comb_all)[i]])
    } # end for i
  } # end if(length(layers_def$omics)>1)

  partition_func_UB <- sum(log(mapply(energy_all_configs_node,FUN=function(x) sum(exp(-0*x)))))
  
  ## PART II: BGe score for each parent set configuration
  data <- do.call(cbind, omics[mapply(nrow,omics)>0])
  myScore <- scoreparameters_BiDAG_BGe(n = ncol(data), data = data)
  
  # for each parent set of given node
  n <- ncol(myScore$data)
  BGe_score_list <- list()
  ## GE
  for(i in c(1:ncol(omics[[layers_def$omics[1]]])))
  {
    BGe_score_list[[i]] <- lapply(comb_all[[i]], FUN=function(list) apply(list, 2, FUN=function(column) 
      if(is.na(column[1]))
      {
        DAGcorescore(names(comb_all)[i], integer(length = 0), n = myScore$n, param = myScore)
      } else {
        DAGcorescore(names(comb_all)[i], column, n = myScore$n, param = myScore)
      }))
  } # end for(i in c(1:ncol(omics[[layers_def$omics[1]]])))
  ## CNV/METH
  if(length(layers_def$omics)>1)
  {
    for(i in c((ncol(omics[[layers_def$omics[1]]])+1):(ncol(omics[[layers_def$omics[1]]])+sum(mapply(ncol,omics[setdiff(layers_def$omics,layers_def$omics[1])])))))
    {
      BGe_score_list[[i]] <- matrix(DAGcorescore(names(comb_all)[i], integer(length = 0), n = myScore$n, param = myScore))
    } # end for i
  } # end if(length(layers_def$omics)>1)
  
  names(BGe_score_list) <- names(comb_all)
  
  return(list(partition_func_UB = partition_func_UB, energy_all_configs_node = energy_all_configs_node, parents_set_combinations = comb_all, BGe_score_all_configs_node = BGe_score_list))
}

#' Range between 0 and 1
#' @description
#' `range_01` This function re-scales a numeric vector so that it ranges between 0 and 1.
#' @export
range_01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' BGe score parameters
#' @description
#' `scoreparameters_BiDAG_BGe` Returns parameters needed for calculation of the BGe score. This function is from BiDAG package.
#' @export
scoreparameters_BiDAG_BGe <- function (n, 
                                       data, 
                                       bgepar = list(am = 1, aw = NULL))
{

  if (anyNA(data)) {
    print("Warning: Dataset contains missing data (covariance matrix computation: complete.obs parameter - missing values are handled by casewise deletion)")
  }

  if (all(is.character(colnames(data)))) {
    nodeslabels <- colnames(data)
  } else {
    nodeslabels <- sapply(c(1:n), function(x) paste("v", x, sep = ""))
  }

  colnames(data) <- nodeslabels
  initparam <- list()
  initparam$labels <- nodeslabels
  initparam$type <- "bge"
  initparam$DBN <- FALSE
  initparam$weightvector <- NULL
  initparam$data <- data

  initparam$bgnodes <- NULL
  initparam$static <- NULL
  initparam$mainnodes <- c(1:n)
  
  initparam$bgn <- 0
  initparam$n <- n
  initparam$nsmall <- n
  initparam$labels.short <- initparam$labels
  initparam$logedgepmat <- NULL

  N <- nrow(data)
  covmat <- cov(data, use = "complete.obs") * (N - 1)
  means <- colMeans(data, na.rm = TRUE)
  bgepar$aw <- n + bgepar$am + 1
  
  initparam$am <- bgepar$am
  initparam$aw <- bgepar$aw
  initparam$N <- N
  initparam$means <- means
  mu0 <- numeric(n)
  T0scale <- bgepar$am * (bgepar$aw - n - 1)/(bgepar$am + 1)
  T0 <- diag(T0scale, n, n)
  initparam$TN <- T0 + covmat + ((bgepar$am * N)/(bgepar$am + N)) * (mu0 - means) %*% t(mu0 - means)
  initparam$awpN <- bgepar$aw + N
  constscorefact <- -(N/2) * log(pi) + (1/2) * log(bgepar$am/(bgepar$am +  N))
  initparam$muN <- (N * means + bgepar$am * mu0)/(N + bgepar$am)
  initparam$SigmaN <- initparam$TN/(initparam$awpN - n - 1)
  initparam$scoreconstvec <- numeric(n)
  for (j in (1:n)) {
    awp <- bgepar$aw - n + j
    initparam$scoreconstvec[j] <- constscorefact - lgamma(awp/2) + 
      lgamma((awp + N)/2) + ((awp + j - 1)/2) * log(T0scale)
  }
  attr(initparam, "class") <- "scoreparameters"
  initparam
}
