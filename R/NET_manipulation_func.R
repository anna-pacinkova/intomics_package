#' Number of reverse edge candidates
#' @description
#' `fan_in_reverse` Determine the number of edges that can be reversed using the fan-in restriction in the largest layer.
#' @export
fan_in_reverse <- function(positions, net_layer_max, layers_def)
{
  net_layer_max[positions["col"],positions["row"]] <- 1
  possible_rev_edges <- sum(net_layer_max[,positions["col"]]) <= layers_def$fan_in_ge[1]
  return(possible_rev_edges)
}

#' Random initial network
#' @description
#' `init.net.mcmc` This function is used to sample random initial network. The edges are sampled only between GE nodes.
#' @export
init.net.mcmc <- function(omics, seed, layers_def)
{
  empty.net <- matrix(0, nrow = sum(mapply(ncol,omics)), ncol = sum(mapply(ncol,omics)),dimnames = list(unlist(mapply(colnames,omics)),
                                                                                                        unlist(mapply(colnames,omics))))
  init.net <- suppressWarnings(sample.chain(empty_net = empty.net, omics_ge = omics[[layers_def$omics[1]]], seed=seed))
  # Let's create an initial network with edges ONLY between GE nodes. 

  rownames(init.net@dag) <- rownames(empty.net)
  colnames(init.net@dag) <- rownames(empty.net)
  
  # if the initial network is not acyclic, we have to define the new one
  while(!is.acyclic(init.net@dag) | any(colSums(init.net@dag[colnames(omics[[layers_def$omics[1]]]),colnames(omics[[layers_def$omics[1]]])]) > layers_def$fan_in_ge[1]))
  {
    init.net <- sample.chain(empty_net = empty.net, omics_ge = omics[[layers_def$omics[1]]])
    rownames(init.net@dag) <- rownames(empty.net)
    colnames(init.net@dag) <- rownames(empty.net)
  }
  
  source.net <- list(adjacency = init.net@dag, nbhd.size = c(), proposal.distr = c(), energy = c(), prior = c(), BGe = c(), likelihood_part = c(), likelihood = c(), acceptance = c(), edge_move = c())
  return(list(source.net = source.net, empty.net = empty.net))
}

#' Acyclic network identification.
#' @description
#' `is.acyclic` This function is from bnstruct R package. Check if the directed graph is acyclic.
#' @export
is.acyclic <- function(g)
{
  rem <- rep(FALSE,nrow(g))
  while( !all(rem) ) # still some edges to remove
  {
    leaves <- (rowSums(g) == 0)
    if( !any(leaves & !rem) )
      return(FALSE)
    g[,leaves] <- 0L
    rem <- rem | leaves
  }
  return(TRUE)
}

#' Neighborhood size
#' @description
#' `neighborhood_size` This function is determines number of network structures that can be reached from the current network structure.
#' @export
neighborhood_size <- function(net, layers_def, B_prior_mat, omics)
{
  ### remove an edge
  remove.edge.size <- sum(net)
  
  ### reverse an edge
  ## choose GE-GE edges (can be reversed)
  layer_max <- colnames(omics[[layers_def$omics[1]]])
  net_layer_max <- net[layer_max,layer_max]
  reverse_edge_pos <- which(net_layer_max==1, arr.ind = TRUE)
  reverse.edge.size <- sum(apply(reverse_edge_pos,1,FUN=function(x) fan_in_reverse(positions = x, net_layer_max = net_layer_max, layers_def = layers_def)))

  ### add an edge
  ## choose all possible CNV-GE/METH-GE edges that are missing
  layer_lower <- unlist(lapply(omics[layers_def$omics[-1]],colnames))
  net_layer_lower <- net[layer_lower,layer_max]
  B_prior_mat_layer_lower <- B_prior_mat[layer_lower,layer_max]
  
  # if there is some prior knowledge about CNV-GE interaction (not only 0.5), the condition could be: B_prior_mat_layer_lower>0
  add.edge.size <- sum(net_layer_lower[B_prior_mat_layer_lower > 0]==0)
  
  ## choose all possible GE-GE edges that are missing (fan in restriction)
  add.edge.size <- add.edge.size + sum(net_layer_max[,colSums(net_layer_max) <= (layers_def$fan_in_ge[1]-1)]==0)
  
  nbhd_size <- remove.edge.size + reverse.edge.size + add.edge.size
  return(nbhd_size)
}

#' Random initial network edge generation
#' @description
#' `sample.chain` This function is used to sample random initial network. The edges are sampled only between GE nodes.
#' @export
sample.chain <- function(empty_net, omics_ge, seed)
{
  dataset_BND <- BNDataset(data = empty_net,
                           discreteness = rep('d',ncol(empty_net)),
                           variables = c(colnames(empty_net)),
                           node.sizes = rep(2,ncol(empty_net)), starts.from=0)
  # The number of values a variable can take is called its cardinality
  # node.sizes: vector of variable cardinalities (for discrete variables)
  # discretenes: now I use 'd', even if the values are continuous
  
  net <- BN(dataset_BND)
  net.dag <- dag(net)
  if(missing(seed))
  {
    n <- ncol(omics_ge)
    chain <- sample(n,n)
    for( i in 2:n )
    {
      net.dag[chain[i-1],chain[i]] <- 1
    }
    dag(net) <- net.dag
    return( suppressMessages(learn.params(net,dataset_BND)) )
  } else {
    set.seed(seed)
    n <- ncol(omics_ge)
    chain <- sample(n,n)
    for( i in 2:n )
    {
      net.dag[chain[i-1],chain[i]] <- 1
    }
    dag(net) <- net.dag
    return( suppressMessages(learn.params(net,dataset_BND)) )
  }
}

#' Source network for MCMC simulation
#' @description
#' `source_net_def` This function is used to create the initial network with its features necessary for MCMC simulation.
#' @export
source_net_def <- function(init.net.mcmc.output, omics, parent_set_combinations, BGe_score_all_configs_node, B_prior_mat, layers_def, energy_all_configs_node, len)
{
  beta_init <- runif(1, min = 0, max = 10)
  beta.source <- list(value = beta_init, prior = c())
  source.net <- init.net.mcmc.output$source.net
  ### source.net parameters necessary for MC3 method: 
  ## Marginal likelihood P(D|G) using ln(BGe score):
  source.net$BGe <- BGe_score(adjacency_matrix = source.net$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
  source.net$nbhd.size <- neighborhood_size(net = source.net$adjacency, layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
  ## Prior probability
  source.net$energy <- sum(epsilon(net = source.net$adjacency, B_prior_mat = B_prior_mat))
  partition_func_UB_beta_source <- sum(mapply(energy_all_configs_node,FUN=function(x) logSumExp(-beta.source$value*x)))
  # partition_func_UB is already log transformed
  # the output of this function is also log transformed (natural logarithm)!
  source.net$prior <- (-beta.source$value*source.net$energy) - partition_func_UB_beta_source
  source.net$likelihood_part <- source.net$BGe + source.net$prior
  
  beta.source$prior <- source.net$prior
  beta.source$len <- len
  acceptance_saved <- vector("numeric")
  acceptance_beta_saved <- vector("numeric")
  method_choice_saved <- vector("numeric")
  nets <- list()
  nets[[1]] <- source.net
  betas <- list()
  betas[[1]] <- beta.source
  return(list(source.net=source.net, 
              partition_func_UB_beta_source=partition_func_UB_beta_source, 
              beta.source = beta.source, 
              acceptance_saved = acceptance_saved, 
              acceptance_beta_saved = acceptance_beta_saved, 
              method_choice_saved = method_choice_saved, 
              energy_all_configs_node = energy_all_configs_node,
              nets = nets,
              B_prior_mat = B_prior_mat,
              betas = betas))
}

#' Markov Blanket Resampling
#' @description
#' `MBR` This function performs the markov blanket resampling method according to Su and Borsuk, 2016.
#' @export
MBR <- function(source_net_adjacency, layers_def, omics, BGe_score_all_configs_node, parent_set_combinations) 
{
  # Randomly select a node Xi in current graph G (MBR method is only for GE nodes)
  selected_node <- sample(colnames(omics[[layers_def$omics[1]]]),1)
  ## PART I: propose the candidate graph using MBR
  dag_tmp <- source_net_adjacency
  
  # store the current parent set
  current_parent_set <- names(which(dag_tmp[,selected_node]==1))
  if(length(current_parent_set)==0)
  {
    current_parent_set <- NA
  } # end if(length(current_parent_set)==0)
  
  # store the children of selected_node
  children_selected_node <- names(which(dag_tmp[selected_node,]==1))
  
  # delete all edges pointing into selected_node and into its children (except of selected_node itself)
  dag_tmp[,selected_node] <- 0
  dag_tmp[-which(rownames(dag_tmp)==selected_node),children_selected_node] <- 0
  
  # store the set of all selected_node's descendant nodes in the current graph
  dag_tmp_bn <- empty.graph(rownames(dag_tmp),1)
  amat(dag_tmp_bn) <- dag_tmp
  descendants_selected_node <- bnlearn::descendants(x = dag_tmp_bn, node = selected_node)
  
  # possible parent sets of selected_node that contain a node from descendants_selected_node and its parents
  selected_node_parents_scores <- parent_sets_sum_scores_X(parent_set_combinations = parent_set_combinations, selected_node = selected_node, descendants = descendants_selected_node, parent_set = current_parent_set, BGe_score_all_configs_node = BGe_score_all_configs_node)
  # add the new parent set to dag_tmp
  if(!is.na(selected_node_parents_scores$new_parent_set[1]))
  {
    amat(dag_tmp_bn)[selected_node_parents_scores$new_parent_set, selected_node] <- 1
  }
  # for j from 1 to J (J is number of selected_node's children) in a specific order for each MBR move:
  # find and store the set of all child's descendant nodes in the current graph
  if(length(children_selected_node)>0)
  {
    child_order <- sample(1:length(children_selected_node),length(children_selected_node))
    selected_node_children_parents_scores <- parent_sets_sum_scores_childrenX(parent_set_combinations = parent_set_combinations, selected_node = selected_node, children_selected_node = children_selected_node, child_order = child_order, dag_tmp_bn = dag_tmp_bn, new_parent_set = TRUE, source_net_adjacency = source_net_adjacency, BGe_score_all_configs_node = BGe_score_all_configs_node)
    dag_tmp_bn <- selected_node_children_parents_scores$dag_tmp_bn
  } # if(length(children_selected_node)>0)
  
  candidate_net_adjacency <- amat(dag_tmp_bn)
  
  ## PART II: complementary inverse move
  # delete all edges pointing into selected_node and into its children (except of selected_node itself) in candidate network
  amat(dag_tmp_bn)[,selected_node] <- 0
  amat(dag_tmp_bn)[-which(rownames(dag_tmp)==selected_node),children_selected_node] <- 0
    
  # mark pre-computed scores of possible parent sets of selected_node which contain a node from descendants of selected_node and which contain its parents in the candidate dag
  selected_node_parents_scores_candidate <- parent_sets_sum_scores_X(parent_set_combinations = parent_set_combinations, selected_node = selected_node, descendants = descendants_selected_node, parent_set = selected_node_parents_scores$new_parent_set, BGe_score_all_configs_node = BGe_score_all_configs_node)
    
  # add the new parent set selected node to dag_G_0_candidate resulting in dag_G_1_candidate
  amat(dag_tmp_bn)[current_parent_set, selected_node] <- 1
    
  if(length(children_selected_node)>0)
  {
    selected_node_children_parents_scores_candidate <- parent_sets_sum_scores_childrenX(parent_set_combinations = parent_set_combinations, selected_node = selected_node, children_selected_node = children_selected_node, child_order = child_order, dag_tmp_bn = dag_tmp_bn, new_parent_set = FALSE, source_net_adjacency = source_net_adjacency, BGe_score_all_configs_node = BGe_score_all_configs_node)
    r_source_candidate <- (selected_node_parents_scores$sum_score_unmarked + sum(selected_node_children_parents_scores$sum_score_unmarked)) - (selected_node_parents_scores_candidate$sum_score_unmarked + sum(selected_node_children_parents_scores_candidate$sum_score_unmarked))
  } else {
    r_source_candidate <- (selected_node_parents_scores$sum_score_unmarked) - (selected_node_parents_scores_candidate$sum_score_unmarked)
  }

  return(list(adjacency = candidate_net_adjacency, nbhd.size = c(), proposal.distr = c(), energy = c(), prior = c(), BGe = c(), likelihood_part = c(), likelihood = c(), acceptance = r_source_candidate, edge_move = c()))
}

#' MBR sum of children scores
#' @description
#' `parent_sets_sum_scores_childrenX` This function determines the sum of BGe scores of given node's children.
#' @export
parent_sets_sum_scores_childrenX <- function(parent_set_combinations, selected_node, children_selected_node, child_order, dag_tmp_bn, new_parent_set, source_net_adjacency, BGe_score_all_configs_node)
{
  sum_score_unmarked <- c()
  
  if(new_parent_set)
  {
    for(j in child_order)
    {
      descendants <- bnlearn::descendants(x = dag_tmp_bn, node = children_selected_node[j])
      # possible parent sets of selected_node child that contain a node from descendants and that does not contain selected_node
      BGe_marked <- lapply(parent_set_combinations[[children_selected_node[j]]], FUN=function(list) apply(list, 2, FUN=function(column) length(intersect(column, descendants))>0 | !any(column==selected_node)))
      # empty parent set is also possible (it is always possible - it does not contain selected_node)
      BGe_marked[[1]][is.na(parent_set_combinations[[children_selected_node[j]]][[1]])] <- TRUE
      names(BGe_marked) <- paste(as.character(1:length(BGe_marked)),"_",sep="")
      # sample new parent set of selected_node's children_selected_node
      BGe_marked_compressed <- lapply(BGe_marked,FUN=function(list) which(list==FALSE))
      possible_parent_sets_ind <- unlist(BGe_marked_compressed, use.names = TRUE)
      if(length(possible_parent_sets_ind)==0)
      {
        new_parent_set <- NA
        sum_score_unmarked[j] <- NA
      } else if(length(possible_parent_sets_ind)==1)
      {
        sum_score_unmarked[j] <- unlist(Map(function(pos, scores) scores[!pos], BGe_marked, BGe_score_all_configs_node[[children_selected_node[j]]]))
        ind <- as.numeric(unlist(lapply(strsplit(names(possible_parent_sets_ind),"_"),FUN=function(list) list[1])))
        new_parent_set <- parent_set_combinations[[children_selected_node[j]]][[ind]][,possible_parent_sets_ind]
      } else {
        score_unmarked <- unlist(Map(function(pos, scores) scores[!pos], BGe_marked, BGe_score_all_configs_node[[children_selected_node[j]]]))
        new_parent_set_ind <- sample(x = 1:length(possible_parent_sets_ind), size = 1, prob = range_01(score_unmarked - sum(score_unmarked)))
        ind <- as.numeric(unlist(lapply(strsplit(names(possible_parent_sets_ind[new_parent_set_ind]),"_"),FUN=function(list) list[1])))
        new_parent_set <- parent_set_combinations[[children_selected_node[j]]][[ind]][,possible_parent_sets_ind[new_parent_set_ind]]
        sum_score_unmarked[j] <- logSumExp(score_unmarked)
      } # end if(length(possible_parent_sets_ind)==0)
      # add the new parent set to dag_tmp
      amat(dag_tmp_bn)[new_parent_set, children_selected_node[j]] <- 1
    } # end for j
    return(list(new_parent_set = new_parent_set, sum_score_unmarked = sum_score_unmarked, dag_tmp_bn = dag_tmp_bn, BGe_marked = BGe_marked))
  } else {
    for(j in child_order)
    {
      descendants <- bnlearn::descendants(x = dag_tmp_bn, node = children_selected_node[j])
      # possible parent sets of selected_node child that contain a node from descendants_child and that does not contain selected_node
      BGe_marked <- lapply(parent_set_combinations[[children_selected_node[j]]], FUN=function(list) apply(list, 2, FUN=function(column) length(intersect(column, descendants))>0 | !any(column==selected_node)))
      # empty parent set is also possible (it is always possible - it does not contain selected_node)
      BGe_marked[[1]][is.na(parent_set_combinations[[children_selected_node[j]]][[1]])] <- TRUE
      names(BGe_marked) <- paste(as.character(1:length(BGe_marked)),"_",sep="")
      # store the sum of unmarked scores
      BGe_marked_compressed <- lapply(BGe_marked,FUN=function(list) which(list==FALSE))
      possible_parent_sets_ind <- unlist(BGe_marked_compressed, use.names = TRUE)
      if(length(possible_parent_sets_ind)==0)
      {
        sum_score_unmarked[j] <- NA
      } else if(length(possible_parent_sets_ind)==1)
      {
        sum_score_unmarked[j] <- unlist(Map(function(pos, scores) scores[!pos], BGe_marked, BGe_score_all_configs_node[[children_selected_node[j]]]))
      } else {
        sum_score_unmarked[j] <- logSumExp(unlist(Map(function(pos, scores) scores[!pos], BGe_marked, BGe_score_all_configs_node[[children_selected_node[j]]])))
      } # end if(length(possible_parent_sets_ind)==0)
      # add the original parent set to the candidate dag
      amat(dag_tmp_bn)[names(which(source_net_adjacency[,children_selected_node[j]]==1)), children_selected_node[j]] <- 1
    } # end for j
    return(list(sum_score_unmarked = sum_score_unmarked, dag_tmp_bn = dag_tmp_bn, BGe_marked = BGe_marked))
  } # end if(new_parent_set)
  
}

#' MBR sum of scores
#' @description
#' `parent_sets_sum_scores_X` This function determines the sum of BGe scores of given node's parents.
#' @export
parent_sets_sum_scores_X <- function(parent_set_combinations, selected_node, descendants, parent_set, BGe_score_all_configs_node)
{
  BGe_marked <- lapply(parent_set_combinations[[selected_node]], FUN=function(list) apply(list, 2, FUN=function(column) length(intersect(column, descendants))>0 | length(intersect(column, parent_set))==length(parent_set)))
  names(BGe_marked) <- paste(as.character(1:length(BGe_marked)),"_",sep="")
  # sample new parent set of selected_node
  BGe_marked_compressed <- lapply(BGe_marked,FUN=function(list) which(list==FALSE))
  possible_parent_sets_ind <- unlist(BGe_marked_compressed, use.names = TRUE)
  if(length(possible_parent_sets_ind)==0)
  {
    new_parent_set <- NA
    sum_score_unmarked <- 0
  } else if(length(possible_parent_sets_ind)==1)
  {
    score_unmarked <- unlist(Map(function(pos, scores) scores[!pos], BGe_marked, BGe_score_all_configs_node[[selected_node]]))
    ind <- as.numeric(unlist(lapply(strsplit(names(possible_parent_sets_ind),"_"),FUN=function(list) list[1])))
    new_parent_set <- parent_set_combinations[[selected_node]][[ind]][,possible_parent_sets_ind]
    sum_score_unmarked <- score_unmarked
  } else {
    score_unmarked <- unlist(Map(function(pos, scores) scores[!pos], BGe_marked, BGe_score_all_configs_node[[selected_node]]))
    new_parent_set_ind <- sample(x = 1:length(possible_parent_sets_ind), size = 1, prob = range_01(score_unmarked - sum(score_unmarked)))
    ind <- as.numeric(unlist(lapply(strsplit(names(possible_parent_sets_ind[new_parent_set_ind]),"_"),FUN=function(list) list[1])))
    new_parent_set <- parent_set_combinations[[selected_node]][[ind]][,possible_parent_sets_ind[new_parent_set_ind]]
    sum_score_unmarked <- logSumExp(score_unmarked)
  } # end if(length(possible_parent_sets_ind)==0)
  return(list(new_parent_set = new_parent_set, sum_score_unmarked = sum_score_unmarked, BGe_marked = BGe_marked))
}

#' Markov Chain conventional single edge proposal move
#' @description
#' `MC3` This function samples a conventional single edge proposal move.
#' @export
MC3 <- function(source_net_adjacency, omics, layers_def, B_prior_mat, beta.source, partition_func_UB_beta_source, parent_set_combinations, BGe_score_all_configs_node, annot)
{
  dag.G.0 <- source_net_adjacency
  mc3_edge <- NA
  # Randomly select a node in current graph G
  from <- sample(rownames(source_net_adjacency),1)
  
  if(regexpr("entrezid", from)>0)
  {
    to <- toupper(from)
    if(dag.G.0[from,to]==0)
    {
      dag.G.0[from,to] <- 1
      mc3_edge <- "add"
    } else {
      dag.G.0[from,to] <- 0
      mc3_edge <- "delete"
    } # end if else (dag.G.0[from,to]==0)
    
  } else if(regexpr("ENTREZID", from)>0){
    ge_nodes <- colnames(omics[[layers_def$omics[1]]])
    if(sum(source_net_adjacency[from,])==0)
    {
      # choose only those GE nodes with number of parents < layers_def$fan_in_ge[1]
      GE_opts <- names(which(colSums(dag.G.0[ge_nodes,setdiff(ge_nodes,from)]) < layers_def$fan_in_ge[1]))
      
      if(identical(GE_opts, character(0)))
      {
        mc3_edge <- "delete"
        # randomly choose an edge which will be deleted
        edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
        dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
      } else {
        to <- sample(setdiff(GE_opts, from), 1)
        dag.G.0[from,to] <- 1
        mc3_edge <- "add"
      }# end if else (identical(GE_opts, character(0)))
      
      while(!is.acyclic(dag.G.0))
      {
        dag.G.0[from,to] <- 0
        GE_opts <-setdiff(GE_opts, to)
        
        if(identical(GE_opts, character(0)))
        {
          mc3_edge <- "delete"
          # randomly choose an edge which will be deleted
          edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
          dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
        } else {
          to <- sample(setdiff(GE_opts, from), 1)
          dag.G.0[from,to] <- 1
          mc3_edge <- "add"
        }# end if else (identical(GE_opts, character(0)))
      } # end while(!is.acyclic(dag.G.0))
      
    } else { # end if(sum(source_net_adjacency[from,])==0)
      
      if(sum(dag.G.0[ge_nodes,from])==layers_def$fan_in_ge[1])
      {
        mc3_edge <- sample(c("add","delete"),1)
        if(mc3_edge=="delete")
        {
          dag.G.0[from, sample(names(which(dag.G.0[from,ge_nodes]==1)), 1)] <- 0
        }
        if(mc3_edge=="add")
        {
          GE_opts <- names(which(dag.G.0[from,setdiff(ge_nodes,from)]==0))
          GE_opts <- intersect(GE_opts, names(which(colSums(dag.G.0[ge_nodes,GE_opts])<layers_def$fan_in_ge[1])))
          
          if(identical(GE_opts, character(0)) | is.null(GE_opts))
          {
            mc3_edge <- "delete"
            # randomly choose an edge which will be deleted
            edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
            dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
          } else{
            to <- sample(GE_opts, 1)
            dag.G.0[from, to] <- 1
          }# end if else (identical(GE_opts, character(0)))
          
          while(!is.acyclic(dag.G.0))
          {
            dag.G.0[from, to] <- 0
            GE_opts <- setdiff(GE_opts, to)
            
            if(identical(GE_opts, character(0)))
            {
              mc3_edge <- "delete"
              # randomly choose an edge which will be deleted
              edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
              dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
            } else {
              to <- sample(GE_opts, 1)
              dag.G.0[from, to] <- 1
            }# end if else (identical(GE_opts, character(0)))
          } # end while(!is.acyclic(dag.G.0))
        } # end if(mc3_edge=="add")
      } else {
        mc3_edge <- sample(c("add","reverse","delete"),1)
        if(mc3_edge=="delete")
        {
          to <- sample(names(which(dag.G.0[from,ge_nodes]==1)), 1)
          dag.G.0[from, to] <- 0
        }
        if(mc3_edge=="add")
        {
          GE_opts <- names(which(dag.G.0[from,setdiff(ge_nodes,from)]==0))
          GE_opts <- intersect(GE_opts, names(which(colSums(dag.G.0[ge_nodes, GE_opts, drop=FALSE]) < layers_def$fan_in_ge[1])))
          
          if(identical(GE_opts, character(0)) | is.null(GE_opts))
          {
            mc3_edge <- "delete" 
            # randomly choose an edge which will be deleted
            edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
            dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
          } else{
            to <- sample(GE_opts, 1)
            dag.G.0[from, to] <- 1
          }# end if else (identical(GE_opts, character(0)) | is.null(GE_opts))
          
          while(!is.acyclic(dag.G.0))
          {
            dag.G.0[from, to] <- 0
            GE_opts <- setdiff(GE_opts, to)
            
            if(identical(GE_opts, character(0)))
            {
              mc3_edge <- "delete"
              # randomly choose an edge which will be deleted
              edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
              dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
            } else {
              to <- sample(GE_opts, 1)
              dag.G.0[from, to] <- 1
            }# end if else (identical(GE_opts, character(0)))
          } # end while(!is.acyclic(dag.G.0))
        } # end if(mc3_edge=="add")
        if(mc3_edge=="reverse")
        {
          GE_opts <- names(which(dag.G.0[from,ge_nodes]==1))
          to <- sample(GE_opts, 1)
          dag.G.0[from, to] <- 0
          dag.G.0[to, from] <- 1
          while(!is.acyclic(dag.G.0))
          {
            dag.G.0[from, to] <- 1
            dag.G.0[to, from] <- 0
            GE_opts <- setdiff(GE_opts, to)
            
            if(identical(GE_opts, character(0)))
            {
              mc3_edge <- sample(c("add","delete"),1)
              if(mc3_edge=="delete")
              {
                dag.G.0[from, sample(names(which(dag.G.0[from,ge_nodes]==1)), 1)] <- 0
              } else {
                GE_opts <- names(which(dag.G.0[from,setdiff(ge_nodes,from)]==0))
                GE_opts <- intersect(GE_opts, names(which(colSums(dag.G.0[ge_nodes,GE_opts])<layers_def$fan_in_ge[1])))
                
                if(identical(GE_opts, character(0)))
                {
                  mc3_edge <- "delete"
                  # randomly choose an edge which will be deleted
                  edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
                  dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
                } else {
                  to <- sample(GE_opts, 1)
                  dag.G.0[from, to] <- 1 
                }# end if else (identical(GE_opts, character(0)))
                
                while(!is.acyclic(dag.G.0))
                {
                  dag.G.0[from, to] <- 0
                  GE_opts <- setdiff(GE_opts, to)
                  
                  if(identical(GE_opts, character(0)))
                  {
                    mc3_edge <- "delete"
                    # randomly choose an edge which will be deleted
                    edge_opts <- which(dag.G.0[ge_nodes,ge_nodes]==1)
                    dag.G.0[ge_nodes, ge_nodes][sample(edge_opts,1)] <- 0
                  } else {
                    to <- sample(GE_opts, 1)
                    dag.G.0[from, to] <- 1
                  }# end if else (identical(GE_opts, character(0)))
                } # end while(!is.acyclic(dag.G.0))
              } # end if else (mc3_edge=="delete")
            } else {
              to <- sample(GE_opts, 1)
              dag.G.0[from, to] <- 0
              dag.G.0[to, from] <- 1
            } # end if else (identical(GE_opts, character(0)))
          } # end while(!is.acyclic(dag.G.0))
        } # end if(mc3_edge=="reverse")
      } # end if else sum(dag.G.0[ge_nodes,from])==layers_def$fan_in_ge[1]
    } # end if else (sum(source_net_adjacency[from,])==0)
  } else {
    GE_opts <- names(which(mapply(FUN=function(char) any(char==from),annot)==TRUE))
    fan_in_meth <- layers_def$fan_in_ge[names(which(mapply(FUN=function(l) any(colnames(l)==from),omics)==TRUE))]
    if(is.na(fan_in_meth))
    {
      to <- sample(GE_opts, 1)
      if(dag.G.0[from, to]==1)
      {
        dag.G.0[from, to] <- 0
        mc3_edge <- "delete"
      } else{
        dag.G.0[from, to] <- 1
        mc3_edge <- "add"
      } # end if else (dag.G.0[from, to]==1)
    } else {
      ge_nodes <- colnames(omics[[layers_def$omics[1]]])
      GE_opts <- intersect(GE_opts, names(which(colSums(dag.G.0[ge_nodes,GE_opts])<fan_in_meth)))
      if(identical(GE_opts, character(0)) | is.null(GE_opts))
      {
        mc3_edge <- "delete"
        # randomly choose an edge which will be deleted
        edge_opts <- which(dag.G.0==1)
        dag.G.0[sample(edge_opts,1)] <- 0
       } else {
        to <- sample(GE_opts, 1)
        if(dag.G.0[from, to]==1)
        {
          dag.G.0[from, to] <- 0
          mc3_edge <- "delete"
        } else{
          dag.G.0[from, to] <- 1
          mc3_edge <- "add"
        } # end if else (dag.G.0[from, to]==1)
      }# end if else (identical(GE_opts, character(0)) | is.null(GE_opts))
    } # end if else (is.na(layers_def$fan_in_ge[names(which(mapply(FUN=function(l) any(colnames(l)==from),omics)==TRUE))]))
  }# end if else (regexpr("entrezid", from)>0)

  nbhd.size <- neighborhood_size(net = dag.G.0, layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
  
  ## Prior probability P(G)
  energy <- sum(epsilon(net = dag.G.0, B_prior_mat = B_prior_mat))
  prior <- (-beta.source$value*energy) - partition_func_UB_beta_source
  
  ### Marginal likelihood P(D|G) using BGe score:
  BGe <- BGe_score(adjacency_matrix = dag.G.0, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
  
  likelihood_part <- BGe + prior
  
  return(list(adjacency = dag.G.0, nbhd.size = nbhd.size, proposal.distr = c(), energy = energy, prior = prior, BGe = BGe, likelihood_part = likelihood_part, likelihood = c(), acceptance = c(), edge_move = mc3_edge))
}
