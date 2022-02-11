#' Trace plots of MCMC simulation & resulting network definition
#' @description
#' `trace_plots` Create trace plots of MCMC simulation and filter low reliable edges based on the edge_freq_thres parameter. Defines the resulting network structure and determines the color scale for each modality.
#' @export
trace_plots <- function(mcmc_res, burn_in, thin, figures_dir, gene_annot, PK, OMICS_module_res, edge_freq_thres = NULL, gene_ID)
{
  if(!dir.exists(figures_dir)){dir.create(figures_dir)}

  # Trace plot of beta values
  df1 <- data.frame(beta = mapply(mcmc_res$beta_tuning,FUN=function(x) x$value), k=1:length(mapply(mcmc_res$beta_tuning,FUN=function(x) x$value)), accept = 1)

  # RMS strength (threshold)
  rms_strength <- abs(diff(mcmc_res$sampling.phase_res$rms))
  strength_threshold <- quantile(rms_strength, 0.75, na.rm = TRUE)

  # custom.strength estimates the strength of each arc as its empirical frequency over a set of networks:
  cpdags1 <- unique(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed1$cpdags[(burn_in/thin+1):length(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed1$cpdags)])
  cpdags2 <- unique(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed2$cpdags[(burn_in/thin+1):length(mcmc_res$sampling.phase_res$mcmc_sim_part_res$seed2$cpdags)])

  cpdag_weights1 <- custom.strength(cpdags1, nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
  cpdag_weights2 <- custom.strength(cpdags2, nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
  cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
  cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]

  cpdag_weights1$edge <- paste(cpdag_weights1$from, cpdag_weights1$to, sep="_")
  cpdag_weights2$edge <- paste(cpdag_weights2$from, cpdag_weights2$to, sep="_")
  cpdag_weights <- merge(cpdag_weights1, cpdag_weights2, by = "edge")
  cpdag_weights$strength <- round(rowMeans(cbind(cpdag_weights$strength.x, cpdag_weights$strength.y)),2)
  
  if(!is.null(edge_freq_thres))
  {
    strength_quant <- quantile(x = cpdag_weights$strength, probs = edge_freq_thres)
    cpdag_weights <- cpdag_weights[cpdag_weights$strength >= strength_quant,]
  }

  total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))

  svg(paste(figures_dir,"beta_values.svg",sep="/"))
  plot(df1$beta ~ df1$k, type = "l", col= "darkblue",
       main = "Beta values of adaptive MCMC",
       xlab = "iteration",
       ylab = "beta")
  dev.off()

  svg(paste(figures_dir,"post_prob_edges.svg",sep="/"))
  plot(total$strength.x ~ total$strength.y,
       main = "Consistency of edges posterior probabilities",
       xlab="MCMC run 2",
       ylab = "MCMC run 1")
  abline(0,1, col="orange")
  dev.off()

  svg(paste(figures_dir,"convergence_RMS.svg",sep="/"))
  plot(rms_strength, main="Convergence RMS strength (C.RMS.str)", pch = 18, col="gray30")
  abline(h=strength_threshold, col="#E69F00", lwd = 1.5)
  text(label = paste("3rd quartile of C.RMS.str = ", round(strength_threshold,3),sep=""), x = 100, y = strength_threshold+0.015, col="#E69F00")
  dev.off()

  if(gene_ID=="entrezID")
  {
    edge_list <- matrix(data = c(cpdag_weights$from.x,
                                        cpdag_weights$to.x,
                                        cpdag_weights$strength,
                                        rep(NA,length(cpdag_weights$strength)),
                                        rep(NA,length(cpdag_weights$strength))),
                               nrow = length(cpdag_weights$strength),
                               dimnames = list(c(), c("from", "to", "weight", "edge_type", "edge")))
    # needs to be sorted because of colors in the final figure
    node_list <- unique(c(edge_list[,"from"], edge_list[,"to"]))
    edge_list[,"edge"] <- paste(edge_list[,"from"], edge_list[,"to"], sep="_")
    return_list <- edge_types(PK = PK, gene_annot = gene_annot, edge_list = edge_list, node_list = node_list, OMICS_module_res = OMICS_module_res)
  } else if (gene_ID=="gene_symbol") {
    from <- as.character(gene_annot$gene_symbol[match(cpdag_weights$from.x, gene_annot$entrezID)])
    from[is.na(from)] <- cpdag_weights$from.x[is.na(from)]
    from[regexpr("entrezid",from)>0] <- tolower(as.character(gene_annot$gene_symbol[match(toupper(from[regexpr("entrezid",from)>0]), gene_annot$entrezID)]))
    to <- as.character(gene_annot$gene_symbol[match(cpdag_weights$to.x, gene_annot$entrezID)])
    edge_list <- matrix(data = c(from,
                                 to,
                                 cpdag_weights$strength,
                                 rep(NA,length(cpdag_weights$strength)),
                                 rep(NA,length(cpdag_weights$strength))),
                        nrow = length(cpdag_weights$strength),
                        dimnames = list(c(), c("from", "to", "weight", "edge_type", "edge")))
    # needs to be sorted because of colors in the final figure
    node_list <- unique(c(edge_list[,"from"], edge_list[,"to"]))
    edge_list[,"edge"] <- paste(edge_list[,"from"], edge_list[,"to"], sep="_")
    return_list <- edge_types(PK = PK, gene_annot = gene_annot, edge_list = edge_list, node_list = node_list, OMICS_module_res = OMICS_module_res)
  } else {
    print('gene_ID argument must be either "entrezID" or "gene_symbol"')
  } # end if else if else (gene_ID=="entrezID")

  net_weighted <- graph_from_edgelist(return_list$edge_list[,c("from","to")])
  V(net_weighted)$color <- return_list$node_list[match(as_ids(V(net_weighted)),return_list$node_list[,"label"]),"color"]
  palette <- return_list$node_palette
  names(palette) <- 1:length(palette)
  palette <- palette[unique(V(net_weighted)$color)]
  V(net_weighted)$label <- return_list$node_list[match(as_ids(V(net_weighted)),return_list$node_list[,"label"]),"label"]
  E(net_weighted)$edge <- return_list$edge_list[match(sub("|","_",as_ids(E(net_weighted)), fixed = TRUE),return_list$edge_list[,"edge"]),"edge_type"]
  E(net_weighted)$weight <- return_list$edge_list[match(sub("|","_",as_ids(E(net_weighted)), fixed = TRUE),return_list$edge_list[,"edge"]),"weight"]
  
  # arrow in the network plot
  V(net_weighted)$degree <- degree(net_weighted, mode = "in")
  V(net_weighted)$degree <- normalise(V(net_weighted)$degree, to = c(3, 11))
  
  return(list(edge_list = return_list$edge_list,
              node_list = return_list$node_list,
              borders_GE = return_list$borders_GE,
              borders_CNV = return_list$borders_CNV,
              borders_METH = return_list$borders_METH,
              node_palette = palette,
              net_weighted = net_weighted))
}

#' Arrow of directed edges tuning
#' @description
#' `normalise` This function is from the ambient package. It is used to determine the position of directed edge arrows.
#' @export
normalise <- function (x, from = range(x), to = c(0, 1)) {
  x <- (x - from[1])/(from[2] - from[1])
  if (!identical(to, c(0, 1))) {
    x <- x * (to[2] - to[1]) + to[1]
  }
  x
}

#' Resulting network definition
#' @description
#' `edge_types` Defines the resulting network structure and determines the color scale for each modality. This is part of trace_plots.
#' @export
edge_types <- function(PK, gene_annot, edge_list, node_list, OMICS_module_res)
{
  omics <- OMICS_module_res$omics
  layers_def <- OMICS_module_res$layers_def
  omics_meth_original <- OMICS_module_res$omics_meth_original
  
  if(any(regexpr("ENTREZID:",node_list)>0))
  {
    PK <- paste(PK$src_entrez, PK$dest_entrez, sep="_")
    
    # unique edges for IntOMICS are gray ("#999999"):
    edge_list[match(setdiff(edge_list[,"edge"],PK),edge_list[,"edge"]),"edge_type"] <- "new"
    # edges present in both IntOMICS and PK are green ("#009E73"):
    edge_list[match(intersect(edge_list[,"edge"],PK),edge_list[,"edge"]), "edge_type"] <- "PK"
    
    # GE node colors
    ge_cols <- brewer.pal(9, "Blues")
    ge_common <- intersect(unique(node_list),colnames(omics[[layers_def$omics[1]]]))
    omics_ge_gs <- as.matrix(omics[[layers_def$omics[1]]][,ge_common])
    colnames(omics_ge_gs) <- ge_common
    
    borders_ge_b1 <- unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]<=0], seq(from=min(omics[[layers_def$omics[1]]]), to=0, length.out=5), include.lowest = T)),","),FUN=function(l) l[1]))
    borders_ge_b1[1] <- sub("[","(",borders_ge_b1[1], fixed = TRUE)
    borders_ge_b1 <- as.numeric(sub("(","",borders_ge_b1, fixed = TRUE))
    borders_ge_b2 <- unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]>0], seq(from=0, to=max(omics[[layers_def$omics[1]]]), length.out=6), include.lowest = T)),","),FUN=function(l) l[1]))
    borders_ge_b2[1] <- sub("[","(",borders_ge_b2[1], fixed = TRUE)
    borders_ge_b2 <- as.numeric(sub("(","",borders_ge_b2, fixed = TRUE))
    borders_ge_b <- c(borders_ge_b1,borders_ge_b2)
    
    borders_ge_t1 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]<=0], seq(from=min(omics[[layers_def$omics[1]]]), to=0, length.out=5), include.lowest = T)),","),FUN=function(l) l[2]))))
    borders_ge_t2 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]>0], seq(from=0, to=max(omics[[layers_def$omics[1]]]), length.out=6), include.lowest = T)),","),FUN=function(l) l[2]))))
    borders_ge_t <- c(borders_ge_t1,borders_ge_t2)
    borders <- sort(unique(c(borders_ge_b,borders_ge_t)))
    expr_group <- cut(colMeans(omics_ge_gs), breaks = borders, include.lowest = T, labels = FALSE)
    names(expr_group) <- colnames(omics_ge_gs)
    node_list <- matrix(data = c(node_list,
                                 as.numeric(expr_group[match(node_list, names(expr_group))])),
                        nrow = length(node_list),
                        dimnames = list(c(), c("label", "color")))
    ind_cols <- paste(paste("(", paste(borders[-length(borders)], borders[-1]),sep=""),"]",sep="")
    borders_cnv <- NULL
    borders_meth <- NULL
    
    if(any(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))
    {
      # CNV node colors
      cnv_cols <- brewer.pal(11, "PiYG")
      cnv_common <- intersect(node_list[,"label"][regexpr("entrez",node_list[,"label"])>0],colnames(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]))
      omics_cnv_gs <- as.matrix(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][,cnv_common])
      
      borders_cnv_b1 <- unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]<=0], seq(from=min(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), to=0, length.out=6), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_cnv_b1[1] <- sub("[","(",borders_cnv_b1[1], fixed = TRUE)
      borders_cnv_b1 <- as.numeric(sub("(","",borders_cnv_b1, fixed = TRUE))
      borders_cnv_b2 <- unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]>0], seq(from=0, to=max(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), length.out=7), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_cnv_b2[1] <- sub("[","(",borders_cnv_b2[1], fixed = TRUE)
      borders_cnv_b2 <- as.numeric(sub("(","",borders_cnv_b2, fixed = TRUE))
      borders_cnv_b <- c(borders_cnv_b1,borders_cnv_b2)
      
      borders_cnv_t1 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]<=0], seq(from=min(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), to=0, length.out=6), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_cnv_t2 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]>0], seq(from=0, to=max(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), length.out=7), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_cnv_t <- c(borders_cnv_t1,borders_cnv_t2)
      borders_cnv <- sort(unique(c(borders_cnv_b,borders_cnv_t)))
      
      cnv_group <- cut(colMeans(omics_cnv_gs, na.rm = TRUE), breaks = borders_cnv, include.lowest = T, labels = FALSE) + length(ge_cols)
      names(cnv_group) <- colnames(omics_cnv_gs)
      node_list[regexpr("entrez",node_list[,"label"])>0,"color"] <- as.numeric(cnv_group[match(node_list[regexpr("entrez",node_list[,"label"])>0,"label"], names(cnv_group))])
      
      ge_cols <- c(ge_cols, cnv_cols)
    } # end if(any(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))
    
    if(any(mapply(omics,FUN=function(list) any(regexpr("entrezid:",colnames(list), ignore.case = TRUE)<0))))
    {
      # METH node colors
      meth_cols <- brewer.pal(9, "YlOrRd")
      meth_common <- intersect(node_list[,"label"],colnames(omics_meth_original))
      omics_meth_gs <- as.matrix(omics_meth_original[,meth_common])
      
      borders_meth_b1 <- unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original<=0.5], seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5, length.out=5), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_meth_b1[1] <- sub("[","(",borders_meth_b1[1], fixed = TRUE)
      borders_meth_b1 <- as.numeric(sub("(","",borders_meth_b1, fixed = TRUE))
      borders_meth_b2 <- unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original>0.5], seq(from=0.5, to=max(omics_meth_original, na.rm = TRUE), length.out=6), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_meth_b2[1] <- sub("[","(",borders_meth_b2[1], fixed = TRUE)
      borders_meth_b2 <- as.numeric(sub("(","",borders_meth_b2, fixed = TRUE))
      borders_meth_b <- c(borders_meth_b1,borders_meth_b2)
      
      borders_meth_t1 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original<=0.5], seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5, length.out=5), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_meth_t2 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original>0.5], seq(from=0.5, to=max(omics_meth_original, na.rm = TRUE), length.out=6), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_meth_t <- c(borders_meth_t1,borders_meth_t2)
      borders_meth <- sort(unique(c(borders_meth_b,borders_meth_t)))
      
      meth_group <- cut(colMeans(omics_meth_gs, na.rm = TRUE), breaks = borders_meth, include.lowest = T, labels = FALSE) + length(ge_cols)
      names(meth_group) <- colnames(omics_meth_gs)
      
      node_list[!is.na(match(node_list[,"label"],colnames(omics_meth_original))),"color"] <- as.numeric(meth_group[match(node_list[,"label"][!is.na(match(node_list[,"label"],colnames(omics_meth_original)))],names(meth_group))])

      ge_cols <- c(ge_cols, meth_cols)
    } # end if(any(mapply(omics,FUN=function(list) any(regexpr("entrezid:",colnames(list), ignore.case = TRUE)<0))))
    
  } else {
    PK_src_dest <- as.character(gene_annot$gene_symbol[match(PK$src_entrez,gene_annot$entrezID)])
    PK_src_dest[regexpr("entrezid",PK_src_dest)>0] <- tolower(as.character(gene_annot$gene_symbol[match(toupper(PK_src_dest[regexpr("entrezid",PK_src_dest)>0]),gene_annot$entrezID)]))
    PK_src_dest[is.na(PK_src_dest)] <- PK$src_entrez[is.na(PK_src_dest)]
    
    PK <- paste(PK_src_dest,
                as.character(gene_annot$gene_symbol[match(PK$dest_entrez,gene_annot$entrezID)]), sep="_")
    
    # unique edges for IntOMICS are gray ("#999999"):
    edge_list[match(setdiff(edge_list[,"edge"],PK),edge_list[,"edge"]),"edge_type"] <- "new"
    # edges present in both IntOMICS and PK are green ("#009E73"):
    edge_list[match(intersect(edge_list[,"edge"],PK),edge_list[,"edge"]), "edge_type"] <- "PK"

    # GE node colors
    ge_cols <- brewer.pal(9, "Blues")
    ge_common <- intersect(gene_annot$entrezID[match(unique(node_list),gene_annot$gene_symbol)],colnames(omics[[layers_def$omics[1]]]))
    omics_ge_gs <- as.matrix(omics[[layers_def$omics[1]]][,ge_common])
    colnames(omics_ge_gs) <- gene_annot$gene_symbol[match(ge_common,gene_annot$entrezID)]
    
    borders_ge_b1 <- unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]<=0], seq(from=min(omics[[layers_def$omics[1]]]), to=0, length.out=5), include.lowest = T)),","),FUN=function(l) l[1]))
    borders_ge_b1[1] <- sub("[","(",borders_ge_b1[1], fixed = TRUE)
    borders_ge_b1 <- as.numeric(sub("(","",borders_ge_b1, fixed = TRUE))
    borders_ge_b2 <- unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]>0], seq(from=0, to=max(omics[[layers_def$omics[1]]]), length.out=6), include.lowest = T)),","),FUN=function(l) l[1]))
    borders_ge_b2[1] <- sub("[","(",borders_ge_b2[1], fixed = TRUE)
    borders_ge_b2 <- as.numeric(sub("(","",borders_ge_b2, fixed = TRUE))
    borders_ge_b <- c(borders_ge_b1,borders_ge_b2)
    
    borders_ge_t1 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]<=0], seq(from=min(omics[[layers_def$omics[1]]]), to=0, length.out=5), include.lowest = T)),","),FUN=function(l) l[2]))))
    borders_ge_t2 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[layers_def$omics[1]]][omics[[layers_def$omics[1]]]>0], seq(from=0, to=max(omics[[layers_def$omics[1]]]), length.out=6), include.lowest = T)),","),FUN=function(l) l[2]))))
    borders_ge_t <- c(borders_ge_t1,borders_ge_t2)
    borders <- sort(unique(c(borders_ge_b,borders_ge_t)))
    expr_group <- cut(colMeans(omics_ge_gs), breaks = borders, include.lowest = T, labels = FALSE)
    names(expr_group) <- colnames(omics_ge_gs)
    node_list <- matrix(data = c(node_list,
                                 as.numeric(expr_group[match(node_list, names(expr_group))])),
                        nrow = length(node_list),
                        dimnames = list(c(), c("label", "color")))
    ind_cols <- paste(paste("(", paste(borders[-length(borders)], borders[-1]),sep=""),"]",sep="")
    borders_cnv <- NULL
    borders_meth <- NULL
    
    if(any(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))
    {
      # CNV node colors
      cnv_cols <- brewer.pal(11, "PiYG")
      cnv_common <- intersect(tolower(gene_annot$entrezID[match(toupper(node_list[is.na(node_list[,"color"]),"label"]),gene_annot$gene_symbol)]),colnames(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]))
      omics_cnv_gs <- as.matrix(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][,cnv_common])
      colnames(omics_cnv_gs) <- node_list[node_list[,"label"]==tolower(node_list[,"label"]),"label"][gene_annot$entrezID[match(toupper(node_list[node_list[,"label"]==tolower(node_list[,"label"]),"label"]),gene_annot$gene_symbol)] %in% toupper(colnames(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]))]
      
      borders_cnv_b1 <- unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]<=0], seq(from=min(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), to=0, length.out=6), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_cnv_b1[1] <- sub("[","(",borders_cnv_b1[1], fixed = TRUE)
      borders_cnv_b1 <- as.numeric(sub("(","",borders_cnv_b1, fixed = TRUE))
      borders_cnv_b2 <- unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]>0], seq(from=0, to=max(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), length.out=7), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_cnv_b2[1] <- sub("[","(",borders_cnv_b2[1], fixed = TRUE)
      borders_cnv_b2 <- as.numeric(sub("(","",borders_cnv_b2, fixed = TRUE))
      borders_cnv_b <- c(borders_cnv_b1,borders_cnv_b2)
      
      borders_cnv_t1 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]<=0], seq(from=min(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), to=0, length.out=6), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_cnv_t2 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]][omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]>0], seq(from=0, to=max(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]], na.rm = TRUE), length.out=7), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_cnv_t <- c(borders_cnv_t1,borders_cnv_t2)
      borders_cnv <- sort(unique(c(borders_cnv_b,borders_cnv_t)))
      
      cnv_group <- cut(colMeans(omics_cnv_gs, na.rm = TRUE), breaks = borders_cnv, include.lowest = T, labels = FALSE) + length(ge_cols)
      names(cnv_group) <- colnames(omics_cnv_gs)
      node_list[node_list[,"label"]==tolower(node_list[,"label"]),"color"][gene_annot$entrezID[match(toupper(node_list[node_list[,"label"]==tolower(node_list[,"label"]),"label"]),gene_annot$gene_symbol)] %in% toupper(colnames(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]))] <- as.numeric(cnv_group[match(node_list[node_list[,"label"]==tolower(node_list[,"label"]),"label"][gene_annot$entrezID[match(toupper(node_list[node_list[,"label"]==tolower(node_list[,"label"]),"label"]),gene_annot$gene_symbol)] %in% toupper(colnames(omics[[names(which(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))]]))], names(cnv_group))])
      
      ge_cols <- c(ge_cols, cnv_cols)
    } # end if(any(mapply(FUN=function(mod) any(regexpr("entrezid:",colnames(mod))>0), omics)==TRUE))
    
    if(any(mapply(omics,FUN=function(list) any(regexpr("entrezid:",colnames(list), ignore.case = TRUE)<0))))
    {
      # METH node colors
      meth_cols <- brewer.pal(9, "YlOrRd")
      meth_common <- intersect(node_list[,"label"],colnames(omics_meth_original))
      omics_meth_gs <- omics_meth_original[,meth_common]
      
      borders_meth_b1 <- unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original<=0.5], seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5, length.out=5), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_meth_b1[1] <- sub("[","(",borders_meth_b1[1], fixed = TRUE)
      borders_meth_b1 <- as.numeric(sub("(","",borders_meth_b1, fixed = TRUE))
      borders_meth_b2 <- unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original>0.5], seq(from=0.5, to=max(omics_meth_original, na.rm = TRUE), length.out=6), include.lowest = T)),","),FUN=function(l) l[1]))
      borders_meth_b2[1] <- sub("[","(",borders_meth_b2[1], fixed = TRUE)
      borders_meth_b2 <- as.numeric(sub("(","",borders_meth_b2, fixed = TRUE))
      borders_meth_b <- c(borders_meth_b1,borders_meth_b2)
      
      borders_meth_t1 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original<=0.5], seq(from=min(omics_meth_original, na.rm = TRUE), to=0.5, length.out=5), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_meth_t2 <- as.numeric(sub("]","",unlist(lapply(strsplit(levels(cut(omics_meth_original[omics_meth_original>0.5], seq(from=0.5, to=max(omics_meth_original, na.rm = TRUE), length.out=6), include.lowest = T)),","),FUN=function(l) l[2]))))
      borders_meth_t <- c(borders_meth_t1,borders_meth_t2)
      borders_meth <- sort(unique(c(borders_meth_b,borders_meth_t)))
      
      meth_group <- cut(colMeans(omics_meth_gs, na.rm = TRUE), breaks = borders_meth, include.lowest = T, labels = FALSE) + length(ge_cols)
      names(meth_group) <- colnames(omics_meth_gs)
      
      node_list[!is.na(match(node_list[,"label"],colnames(omics_meth_original))),"color"] <- as.numeric(meth_group[match(node_list[,"label"][!is.na(match(node_list[,"label"],colnames(omics_meth_original)))],names(meth_group))])
      ge_cols <- c(ge_cols, meth_cols)
    } # end if(any(mapply(omics,FUN=function(list) any(regexpr("entrezid:",colnames(list), ignore.case = TRUE)<0))))
  } # end if else(any(regexpr("ENTREZID:",node_list)>0))
  return(list(edge_list = edge_list, node_list = node_list, borders_GE = borders, borders_CNV = borders_cnv, borders_METH = borders_meth, node_palette = ge_cols))
}

#' Node color
#' @description
#' `legend_custom` Determines the color scale for each modality.
#' @export
legend_custom <- function(net)
{
  xl <- 1.2;yb <- 1;xr <- 1.3;yt <- 2
  par(oma=c(0,0,0,0))
  plot(NA,type="n",ann=FALSE,xlim=c(0.94,2),ylim=c(1.2,1.71),xaxt="n",yaxt="n",bty="n")
  text(x = 0.95,y = 1.25,labels = "GE")
  rect(head(seq(yb,yt,(yt-yb)/9),-1), xl, tail(seq(yb,yt,(yt-yb)/9),-1), xr, col=brewer.pal(9, "Blues"))
  text(x = unique(c(head(seq(yb,yt,(yt-yb)/9),-1), tail(seq(yb,yt,(yt-yb)/9),-1))),y = 1.32,labels = round(net$borders_GE,2))
  if(!is.null(net$borders_CNV))
  {
    text(x = 0.94,y = 1.45,labels = "CNV")
    rect(head(seq(yb, yt, (yt-yb)/11), -1), xl+0.2, tail(seq(yb ,yt , (yt-yb)/11), -1), xr+0.2, col=brewer.pal(11, "PiYG"))
    text(x = unique(c(head(seq(yb,yt,(yt-yb)/11),-1), tail(seq(yb,yt,(yt-yb)/11),-1))),y = 1.52,labels = round(net$borders_CNV,2))
  } # end if(!is.null(net$borders_CNV))
  if(!is.null(net$borders_METH))
  {
    rect(head(seq(yb, yt, (yt-yb)/9), -1), xl+0.4, tail(seq(yb ,yt , (yt-yb)/9), -1), xr+0.4, col=brewer.pal(9, "YlOrRd"))
    text(x = unique(c(head(seq(yb,yt,(yt-yb)/9),-1), tail(seq(yb,yt,(yt-yb)/9),-1))),y = 1.72,labels = round(net$borders_METH,2))
    text(x = 0.94,y = 1.65,labels = "METH")
  } # end if(!is.null(net$borders_METH))
}
