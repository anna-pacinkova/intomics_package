#' #' BN module
#' @description
#' `BN_module` Performs automatically tuned MCMC sampling from posterior distribution together with conventional MCMC sampling using empirical biological prior matrix to sample network structures from posterior distribution.
#' @export
BN_module <- function(burn_in, thin, seed1, seed2, OMICS_module_res, minseglen, len = 5, prob_mbr = 0.07)
{
  energy_all_configs_node <- OMICS_module_res$pf_UB_BGe_pre$energy_all_configs_node
  BGe_score_all_configs_node <- OMICS_module_res$pf_UB_BGe_pre$BGe_score_all_configs_node
  parent_set_combinations <- OMICS_module_res$pf_UB_BGe_pre$parents_set_combinations
  
  ##########################
  ### 1st adaption phase ###
  ### we change the variance of the proposal distribution to achieve the MC acceptance rate of 0.44 
  first.adapt.phase_net <- first_adapt_phase(seed1 = seed1, 
                                             omics = OMICS_module_res$omics,
                                             B_prior_mat = OMICS_module_res$B_prior_mat,
                                             energy_all_configs_node = energy_all_configs_node,
                                             len = len,
                                             layers_def = OMICS_module_res$layers_def,
                                             annot = OMICS_module_res$annot,
                                             prob_mbr = prob_mbr, 
                                             BGe_score_all_configs_node = BGe_score_all_configs_node,
                                             parent_set_combinations = parent_set_combinations)

  #######################
  ### transient phase ###
  ### to check if the chain is moving towards the mode of target distribution 
  transient.phase_net <- transient_phase(first.adapt.phase_net = first.adapt.phase_net, 
                  omics = OMICS_module_res$omics,
                  B_prior_mat = OMICS_module_res$B_prior_mat, 
                  layers_def = OMICS_module_res$layers_def,
                  annot = OMICS_module_res$annot,
                  energy_all_configs_node = energy_all_configs_node, 
                  prob_mbr = prob_mbr, 
                  BGe_score_all_configs_node = BGe_score_all_configs_node,
                  parent_set_combinations = parent_set_combinations)

  ##########################
  ### 2nd adaption phase ###
  ### modified Adaptive Metropolis algorithm to find the proposal distribution that has a similar covariance structure with the target distribution
  second.adapt.phase_net <- second_adapt_phase(transient.phase_net = transient.phase_net, 
                                         omics = OMICS_module_res$omics,
                                         B_prior_mat = OMICS_module_res$B_prior_mat, 
                                         energy_all_configs_node = energy_all_configs_node, 
                                         layers_def = OMICS_module_res$layers_def,
                                         prob_mbr = prob_mbr, 
                                         BGe_score_all_configs_node = BGe_score_all_configs_node,
                                         parent_set_combinations = parent_set_combinations,
                                         annot = OMICS_module_res$annot)

  ######################
  ### sampling phase ###
  ### Now we apply 2 MCMC simulations and check the RMS value (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0734-6) to stop the simulation
  ### After the burn-in period, we discard the values from the first half of this phase
  sampling.phase_net <- sampling_phase(second.adapt.phase_net = second.adapt.phase_net, seed1 = seed1, seed2 = seed2, thin = thin,
                                       omics = OMICS_module_res$omics, layers_def = OMICS_module_res$layers_def, prob_mbr = prob_mbr, minseglen = minseglen, burn_in = burn_in)

  # save only beta values and resulting cpdags
  sampling.phase_net$mcmc_sim_part_res$seed1 <- sampling.phase_net$mcmc_sim_part_res$seed1[c("betas","cpdags")]
  sampling.phase_net$mcmc_sim_part_res$seed2 <- sampling.phase_net$mcmc_sim_part_res$seed2[c("betas","cpdags")]
  beta_tuning <- second.adapt.phase_net$betas
  
  return(list(B_prior_mat_weighted = second.adapt.phase_net$B_prior_mat_weighted, 
              sampling.phase_res = sampling.phase_net, 
              beta_tuning = beta_tuning))
}

#' 1st adaption phase
#' @description
#' `first_adapt_phase` 1st adaption phase of the adaptive MCMC: the variance of the proposal distribution is changed to achieve the MC acceptance rate of 0.44.
#' @export
first_adapt_phase <- function(seed1, omics, B_prior_mat, energy_all_configs_node, len, layers_def, prob_mbr, BGe_score_all_configs_node, parent_set_combinations, annot) 
{
  init.net1 <- init.net.mcmc(omics = omics, seed = seed1, layers_def = layers_def)
  first.adapt.phase_net <- source_net_def(init.net.mcmc.output = init.net1, 
                                omics = omics, 
                                parent_set_combinations = parent_set_combinations,
                                BGe_score_all_configs_node = BGe_score_all_configs_node,
                                B_prior_mat = B_prior_mat,
                                layers_def = layers_def,
                                energy_all_configs_node = energy_all_configs_node,
                                len = len)
  first.adapt.phase_net$seed <- seed1
  set.seed(first.adapt.phase_net$seed)
  
  # For every 100 iteration, we calculate the beta acceptance rate of the past 100 iterations 
  # and we add = 0.05 to log(len) if the acceptance rate is higher than 0.44, and subtract from log(len) if the acceptance rate is lower or equal than 0.44.
  # we do this until the acceptance rate for beta falls between 0.28 and 0.60.
  first.adapt.phase_net <- acceptance_check(first.adapt.phase_net = first.adapt.phase_net, round_check = 100, last_iter_check = 100, prob_mbr = prob_mbr, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node, omics = omics, annot = annot)
  
  # We run 100 more iterations with same len, which have made the acceptance rate to fall between 0.28 and 0.60, 
  # and monitor the acceptance rate for the past 200 iterations.
  # If the beta acceptance rate falls outside of 0.28 and 0.60, then we adjust log(len) for every 100 iterations and monitor the acceptance rate for the past 200 iterations until the acceptance rate comes between 0.28 and 0.60.
  first.adapt.phase_net <- acceptance_check(first.adapt.phase_net = first.adapt.phase_net, round_check = 100, last_iter_check = 200, prob_mbr = prob_mbr, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node, omics = omics, annot = annot)

  # We run 200 more iterations with same len, which have made the acceptance rate to fall between 0.28 and 0.60, 
  # and monitor the acceptance rate for the past 400 iterations. 
  # If the beta acceptance rate falls outside of 0.28 and 0.60, then we adjust log(len) for every 200 iterations and monitor the acceptance rate for the past 400 iterations until the acceptance rate comes between 0.28 and 0.60.
  first.adapt.phase_net <- acceptance_check(first.adapt.phase_net = first.adapt.phase_net, round_check = 200, last_iter_check = 400, prob_mbr = prob_mbr, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node, omics = omics, annot = annot)
  # We stop the chain and save len.
  return(first.adapt.phase_net)
}

#' transient phase 
#' @description
#' `transient_phase` This phase verify if the chain is moving towards the mode of target distribution.
#' @export
transient_phase <- function(first.adapt.phase_net, omics, B_prior_mat, layers_def, energy_all_configs_node, prob_mbr, BGe_score_all_configs_node, parent_set_combinations, annot) 
{
  beta_1st_adapt <- first.adapt.phase_net$betas[[length(first.adapt.phase_net$betas)]]$len
  source.net <- first.adapt.phase_net$nets[[length(first.adapt.phase_net$nets)]]
  beta.source <- first.adapt.phase_net$betas[[length(first.adapt.phase_net$betas)]]
  start <- length(first.adapt.phase_net$nets)
 
  # We employ a standard MCMC algorithm with proposals determined by the 1st adaption phase.
  # To check if the chain is moving towards the mode of target distribution, for every 200 iteration, the values generated for beta are averaged
  # With 5 different averages, a linear model is fitted to see if there is any trend in the beta coordinate of the chain.
  for(i in (start+1):(start+1000))
  {
    ## method choice:
    first.adapt.phase_net$method_choice_saved[i] <- sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
    if(first.adapt.phase_net$method_choice_saved[i]=="MC3")
    {
      ### candidate.net definition: 
      candidate.net <- MC3(source_net_adjacency = source.net$adjacency,
                           layers_def =  layers_def, 
                           B_prior_mat = B_prior_mat, 
                           beta.source = beta.source, 
                           partition_func_UB_beta_source = first.adapt.phase_net$partition_func_UB_beta_source, 
                           omics = omics, 
                           parent_set_combinations = parent_set_combinations, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           annot = annot)
      ## Proposal distribution Q(..|..)
      # Q(source.net|candidate.net)
      # Q(candidate.net|source.net)
      candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
      source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
      candidate.net$likelihood <- candidate.net$likelihood_part + source.net$proposal.distr
      source.net$likelihood <- source.net$likelihood_part + candidate.net$proposal.distr
      first.adapt.phase_net$acceptance_saved[i] <- candidate.net$likelihood - source.net$likelihood
      
      candidate_edge <- which(candidate.net$adjacency!=source.net$adjacency, arr.ind = TRUE)

      u <- log(runif(1))
      if (u < first.adapt.phase_net$acceptance_saved[i])
      {
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      } # end if (u < first.adapt.phase_net$acceptance_saved[i])
      first.adapt.phase_net$nets[[i]] <- source.net
      
    } else {
      
      candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                           layers_def = layers_def, 
                           omics = omics, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           parent_set_combinations = parent_set_combinations)
      first.adapt.phase_net$acceptance_saved[i] <- candidate.net$acceptance
      ### source.net parameters necessary for MC3 method: 
      ## Marginal likelihood P(D|G) using ln(BGe score):
      candidate.net$BGe <-BGe_score(adjacency_matrix = candidate.net$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
      candidate.net$nbhd.size <- neighborhood_size(net = candidate.net$adjacency, layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
      ## Prior probability
      candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency, B_prior_mat = B_prior_mat))
      # partition_func_UB is already log transformed
      # the output of this function is also log transformed (natural logarithm)!
      candidate.net$prior <- (-beta.source$value*candidate.net$energy) - first.adapt.phase_net$partition_func_UB_beta_source
      candidate.net$likelihood_part <- candidate.net$BGe + candidate.net$prior
      
      u <- log(runif(1))
      if (u < first.adapt.phase_net$acceptance_saved[i])
      {
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      } # end if (u < acceptance_saved[i])
      first.adapt.phase_net$nets[[i]] <- source.net
      
      partition_func_UB_beta_source <- sum(mapply(first.adapt.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.source$value*x)))
      
    } # end if(method.choice=="MC3")
    
    ### MH for beta with fixed network structure ###
    beta.candidate <- list(value = rnorm(1, mean = beta.source$value, sd = beta_1st_adapt), prior = c(), len = beta_1st_adapt)
    if(beta.candidate$value < 0.5)
    {
      beta.candidate$value <- 0.5
    } # end if(beta.candidate$value < 0.5)
    
    partition_func_UB_beta_candidate <- sum(mapply(first.adapt.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.candidate$value*x)))
    beta.candidate$prior <- (-beta.candidate$value*source.net$energy) - partition_func_UB_beta_candidate
    
    first.adapt.phase_net$acceptance_beta_saved[i] <- beta.candidate$prior - beta.source$prior
    u_beta <- log(runif(1))
    
    if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
    {
      beta.source <- beta.candidate
      first.adapt.phase_net$partition_func_UB_beta_source <- partition_func_UB_beta_candidate
    } # end if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
    first.adapt.phase_net$betas[[i]] <- beta.source
  } # end for(i in (start+1):(start+1000))
  
  beta_means <- colMeans(matrix(mapply(tail(first.adapt.phase_net$betas,1000), FUN=function(list) list$value),nrow=200))
  reg_dat <- data.frame(beta_means = tail(beta_means,5), iter = 1:5)
  model <- lm(beta_means ~ iter, data = reg_dat)
  p.val <- summary(model)$coefficients[1,4]
  
  # We use a regression method to make sure the chain values are moving to only one direction, neither increasing nor deceasing. 
  # If a regression method confirms that the chain values show a linear trend, we presume that the chain is still moving to a local mode. 
  # The p-value for the slope coefficient is used to determine whether there is any linear trend. 
  # If p-value for every coordinate is greater than 0.1, the chain gets stopped and this phase ends.
  while(p.val < 0.1)
  {
    for(i in (length(first.adapt.phase_net$nets)+1):(length(first.adapt.phase_net$nets)+200))
    {
      ## method choice:
      first.adapt.phase_net$method_choice_saved[i] <- sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
      if(first.adapt.phase_net$method_choice_saved[i]=="MC3")
      {
        ### candidate.net definition: 
        candidate.net <- MC3(source_net_adjacency = source.net$adjacency,
                             layers_def =  layers_def, 
                             B_prior_mat = B_prior_mat, 
                             beta.source = beta.source, 
                             partition_func_UB_beta_source = first.adapt.phase_net$partition_func_UB_beta_source, 
                             omics = omics, 
                             parent_set_combinations = parent_set_combinations, 
                             BGe_score_all_configs_node = BGe_score_all_configs_node, 
                             annot = annot)
        ## Proposal distribution Q(..|..)
        # Q(source.net|candidate.net)
        # Q(candidate.net|source.net)
        candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
        source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
        candidate.net$likelihood <- candidate.net$likelihood_part + source.net$proposal.distr
        source.net$likelihood <- source.net$likelihood_part + candidate.net$proposal.distr
        first.adapt.phase_net$acceptance_saved[i] <- candidate.net$likelihood - source.net$likelihood
        
        candidate_edge <- which(candidate.net$adjacency!=source.net$adjacency, arr.ind = TRUE)

        u <- log(runif(1))
        if (u < first.adapt.phase_net$acceptance_saved[i])
        {
          source.net <- candidate.net
          beta.source$prior <- source.net$prior
        } # end if (u < first.adapt.phase_net$acceptance_saved[i])
        first.adapt.phase_net$nets[[i]] <- source.net
        
      } else {
        
        candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                             layers_def = layers_def, 
                             omics = omics, 
                             BGe_score_all_configs_node = BGe_score_all_configs_node, 
                             parent_set_combinations = parent_set_combinations)
        first.adapt.phase_net$acceptance_saved[i] <- candidate.net$acceptance
        ### source.net parameters necessary for MC3 method: 
        ## Marginal likelihood P(D|G) using ln(BGe score):
        candidate.net$BGe <-BGe_score(adjacency_matrix = candidate.net$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
        candidate.net$nbhd.size <- neighborhood_size(net = candidate.net$adjacency, layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
        ## Prior probability
        candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency, B_prior_mat = B_prior_mat))
        # partition_func_UB is already log transformed
        # the output of this function is also log transformed (natural logarithm)!
        candidate.net$prior <- (-beta.source$value*candidate.net$energy) - first.adapt.phase_net$partition_func_UB_beta_source
        candidate.net$likelihood_part <- candidate.net$BGe + candidate.net$prior
        
        u <- log(runif(1))
        if (u < first.adapt.phase_net$acceptance_saved[i])
        {
          source.net <- candidate.net
          beta.source$prior <- source.net$prior
        } # end if (u < acceptance_saved[i])
        first.adapt.phase_net$nets[[i]] <- source.net
        
        partition_func_UB_beta_source <- sum(mapply(first.adapt.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.source$value*x)))
        
      } # end if(method.choice=="MC3")
      
      ### MH for beta with fixed network structure ###
      beta.candidate <- list(value = rnorm(1, mean = beta.source$value, sd = beta_1st_adapt), prior = c(), len = beta_1st_adapt)
      if(beta.candidate$value < 0.5)
      {
        beta.candidate$value <- 0.5
      } # end if(beta.candidate$value < 0.5)
      
      partition_func_UB_beta_candidate <- sum(mapply(first.adapt.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.candidate$value*x)))
      beta.candidate$prior <- (-beta.candidate$value*source.net$energy) - partition_func_UB_beta_candidate
      
      first.adapt.phase_net$acceptance_beta_saved[i] <- beta.candidate$prior - beta.source$prior
      u_beta <- log(runif(1))
      
      if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
      {
        beta.source <- beta.candidate
        first.adapt.phase_net$partition_func_UB_beta_source <- partition_func_UB_beta_candidate
      } # end if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
      first.adapt.phase_net$betas[[i]] <- beta.source
    } # end for(i in (length(first.adapt.phase_net$nets)+1):(length(first.adapt.phase_net$nets)+200))
    
    beta_means <- colMeans(matrix(mapply(tail(first.adapt.phase_net$betas,1000), FUN=function(list) list$value),nrow=200))
    reg_dat <- data.frame(beta_means = beta_means, iter = 1:5)
    model <- lm(beta_means ~ iter, data = reg_dat)
    p.val <- summary(model)$coefficients[1,4]
  } # end while(p.val < 0.1)
  
  return(first.adapt.phase_net)
}

#' Second adaption phase   
#' @description
#' `second_adapt_phase` This phase identifies the proposal distribution that has a similar covariance structure with the target distribution.
#' @export
second_adapt_phase <- function(transient.phase_net, omics, layers_def, B_prior_mat, energy_all_configs_node, prob_mbr, BGe_score_all_configs_node, parent_set_combinations, annot) 
{
  # constant <- 2.38/1.5
  # for details see the original paper
  second.adapt.phase_net <- variance_target(transient.phase_net = transient.phase_net, 
                                            constant = 1.586667, 
                                            fin = 200, 
                                            B_prior_mat = B_prior_mat, 
                                            omics = omics, 
                                            parent_set_combinations = parent_set_combinations, 
                                            BGe_score_all_configs_node = BGe_score_all_configs_node, 
                                            layers_def = layers_def,
                                            prob_mbr = prob_mbr,
                                            annot = annot)
  i <- 1
  while(second.adapt.phase_net$acceptance.rate_betas < 0.02)
  {
    i <- i+1
    constant <- 2.38/(1.5^i)
    second.adapt.phase_net <- variance_target(transient.phase_net = transient.phase_net, 
                                              constant = constant, 
                                              fin = 200, 
                                              B_prior_mat = B_prior_mat, 
                                              omics = omics, 
                                              parent_set_combinations = parent_set_combinations, 
                                              BGe_score_all_configs_node = BGe_score_all_configs_node, 
                                              layers_def = layers_def,
                                              prob_mbr = prob_mbr,
                                              annot = annot)
    
  } # end while(second.adapt.phase_net$acceptance.rate_betas < 0.02)

  # Now calculate the average squared jumping distance for beta (beta_n - beta_n-1)^2.
  # With 5 different averages we fit a linear model to see if there is any trend in squared jumping distance for each coordinate.
  # If the average squared jumping distance stops to increase, the mixing stops to improve and we stop this phase.
  constant <- 2.38/(1.5^i)
  squared.jump_second.adapt.phase_net <- squared_jumping(second.adapt.phase_net = second.adapt.phase_net$variance.target_net, 
                                                         constant = constant,
                                                         beta_sd = second.adapt.phase_net$beta_sd,
                                                         fin = (nrow(B_prior_mat)^2)*5, 
                                                         B_prior_mat = B_prior_mat, 
                                                         omics = omics, 
                                                         parent_set_combinations = parent_set_combinations, 
                                                         BGe_score_all_configs_node = BGe_score_all_configs_node, 
                                                         layers_def = layers_def,
                                                         prob_mbr = prob_mbr,
                                                         annot = annot)

  betas_check <- mapply(tail(squared.jump_second.adapt.phase_net$betas,1001), FUN=function(list) list$value)
  betas_check <- colMeans(matrix((betas_check[-1] - betas_check[-1001])^2,nrow=200))
  reg_dat <- data.frame(beta_means = betas_check, iter = 1:5)
  model <- lm(beta_means ~ iter, data = reg_dat)
  squared.jump_second.adapt.phase_net$p.val <- summary(model)$coefficients[2,4]
 
  while(squared.jump_second.adapt.phase_net$p.val < 0.1)
  {
    squared.jump_second.adapt.phase_net <- squared_jumping(second.adapt.phase_net = squared.jump_second.adapt.phase_net, 
                                                           constant = constant,
                                                           beta_sd = second.adapt.phase_net$beta_sd,
                                                           fin = 200, 
                                                           B_prior_mat = B_prior_mat, 
                                                           omics = omics, 
                                                           parent_set_combinations = parent_set_combinations, 
                                                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                                                           layers_def = layers_def,
                                                           prob_mbr = prob_mbr,
                                                           annot = annot)
    
    betas_check <- mapply(tail(squared.jump_second.adapt.phase_net$betas,1001), FUN=function(list) list$value)
    betas_check <- colMeans(matrix((betas_check[-1] - betas_check[-1001])^2,nrow=200))
    reg_dat <- data.frame(beta_means = betas_check, iter = 1:5)
    model <- lm(beta_means ~ iter, data = reg_dat)
    squared.jump_second.adapt.phase_net$p.val <- summary(model)$coefficients[2,4]
  } # end while(squared.jump_second.adapt.phase_net$p.val < 0.1)

  squared.jump_second.adapt.phase_net$constant <- constant
  squared.jump_second.adapt.phase_net$beta_sd <- second.adapt.phase_net$beta_sd
  B_prior_mat_weighted <- c(B_prior_mat)
  conditions <- c(B_prior_mat)==0.5 & c(squared.jump_second.adapt.phase_net$iter_edges[,,1])>0
  B_prior_mat_weighted[conditions] <- c(squared.jump_second.adapt.phase_net$iter_edges[,,2])[conditions] / c(squared.jump_second.adapt.phase_net$iter_edges[,,1])[conditions]
  squared.jump_second.adapt.phase_net$B_prior_mat_weighted <- matrix(B_prior_mat_weighted, nrow=nrow(B_prior_mat), dimnames = list(rownames(B_prior_mat), colnames(B_prior_mat)))
  squared.jump_second.adapt.phase_net$partition_func_UB <- pf_UB_est(omics = omics, layers_def = layers_def, B_prior_mat = squared.jump_second.adapt.phase_net$B_prior_mat_weighted, annot = annot)
  return(squared.jump_second.adapt.phase_net)
}

#' Acceptance rate   
#' @description
#' `acceptance_check` This phase verify if the acceptance is in range of 0.28 and 0.6.
#' @export
acceptance_check <- function(first.adapt.phase_net, round_check, last_iter_check, prob_mbr, layers_def, parent_set_combinations, BGe_score_all_configs_node, omics, annot)
{
  source.net <- first.adapt.phase_net$nets[[length(first.adapt.phase_net$nets)]]
  beta.source <- first.adapt.phase_net$betas[[length(first.adapt.phase_net$betas)]]
  i <- length(first.adapt.phase_net$nets)
  acceptance.rate_betas <- 0

  while(!(acceptance.rate_betas > 0.28 & acceptance.rate_betas < 0.60))
  {
    i <- i+1
    ## method choice:
    first.adapt.phase_net$method_choice_saved[i] <- sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
    if(first.adapt.phase_net$method_choice_saved[i]=="MC3")
    {
      ### candidate.net definition: 
      candidate.net <- MC3(source_net_adjacency = source.net$adjacency,
                           layers_def = layers_def, 
                           B_prior_mat = first.adapt.phase_net$B_prior_mat, 
                           beta.source = beta.source, 
                           partition_func_UB_beta_source = first.adapt.phase_net$partition_func_UB_beta_source, 
                           omics = omics, 
                           parent_set_combinations = parent_set_combinations, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node,
                           annot = annot)
      ## Proposal distribution Q(..|..)
      # Q(source.net|candidate.net)
      # Q(candidate.net|source.net)
      candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
      source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
      candidate.net$likelihood <- candidate.net$likelihood_part + source.net$proposal.distr
      source.net$likelihood <- source.net$likelihood_part + candidate.net$proposal.distr
      first.adapt.phase_net$acceptance_saved[i] <- candidate.net$likelihood - source.net$likelihood
      
      u <- log(runif(1))
      if (u < first.adapt.phase_net$acceptance_saved[i])
      {
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      }
      first.adapt.phase_net$nets[[i]] <- source.net
    } else {
      
      candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                           layers_def = layers_def, 
                           omics = omics, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           parent_set_combinations = parent_set_combinations)
      first.adapt.phase_net$acceptance_saved[i] <- candidate.net$acceptance
      ### source.net parameters necessary for MC3 method: 
      ## Marginal likelihood P(D|G) using ln(BGe score):
      candidate.net$BGe <-BGe_score(adjacency_matrix = candidate.net$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
      candidate.net$nbhd.size <- neighborhood_size(net = candidate.net$adjacency, layers_def = layers_def, B_prior_mat = first.adapt.phase_net$B_prior_mat, omics = omics)
      ## Prior probability
      candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency, B_prior_mat = first.adapt.phase_net$B_prior_mat))
      # partition_func_UB is already log transformed
      # the output of this function is also log transformed (natural logarithm)!
      candidate.net$prior <- (-beta.source$value*candidate.net$energy) - first.adapt.phase_net$partition_func_UB_beta_source
      candidate.net$likelihood_part <- candidate.net$BGe + candidate.net$prior
      
      u <- log(runif(1))
      if (u < first.adapt.phase_net$acceptance_saved[i])
      {
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      } # end if (u < acceptance_saved[i])
      first.adapt.phase_net$nets[[i]] <- source.net
  
    } # end if(method.choice=="MC3")
    
    ### MH for beta with fixed network structure ###
    beta.candidate <- list(value = rnorm(n = 1, mean = beta.source$value, sd = beta.source$len), prior = c(), len = beta.source$len)
    if(beta.candidate$value < 0.5)
    {
      beta.candidate$value <- 0.5
    } # end if(beta.candidate$value < 0.5)
    
    partition_func_UB_beta_candidate <- sum(mapply(first.adapt.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.candidate$value*x)))
    beta.candidate$prior <- (-beta.candidate$value*source.net$energy) - partition_func_UB_beta_candidate
    
    first.adapt.phase_net$acceptance_beta_saved[i] <- beta.candidate$prior - beta.source$prior
    u_beta <- log(runif(1))
    
    if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
    {
      beta.source <- beta.candidate
      first.adapt.phase_net$partition_func_UB_beta_source <- partition_func_UB_beta_candidate
    } # end if (u_beta < first.adapt.phase_net$acceptance_beta_saved[i])
    
    if(numbers::mod(length(first.adapt.phase_net$nets), round_check)==0)
    {
      acceptance.trace_betas <- unlist(lapply(tail(first.adapt.phase_net$betas, last_iter_check),FUN=function(list) list$prior))
      acceptance.trace_betas <- c(1,acceptance.trace_betas[1:(length(acceptance.trace_betas)-1)] - acceptance.trace_betas[2:length(acceptance.trace_betas)])
      acceptance.trace_betas[acceptance.trace_betas!=0] <- 1
      acceptance.rate_betas <- sum(acceptance.trace_betas==1)/length(acceptance.trace_betas)
      # modify len if necessary
      if(acceptance.rate_betas > 0.44)
      {
        beta.source$len <- exp(log(beta.source$len) + 0.05)
      } else {
        beta.source$len <- exp(log(beta.source$len) - 0.05)
      } # end if else (acceptance.rate_betas > 0.44)
    } # end if(numbers::mod(i,round_check)==0)
    first.adapt.phase_net$betas[[i]] <- beta.source
  } # end while(!(acceptance.rate_betas > 0.28 & acceptance.rate_betas < 0.60))
  return(first.adapt.phase_net)
}

#' BGe score  
#' @description
#' `BGe_score` Computes the BGe score of given network using precomputed sets of possible parents.
#' @export
BGe_score <- function(adjacency_matrix, omics, layers_def, parent_set_combinations, BGe_score_all_configs_node)
{
  score_nodes <- 0
  ## GE
  #for(node in rownames(adjacency_matrix)[1:ncol(omics[[layers_def$omics[1]]])])
  for(node in rownames(adjacency_matrix)[1:5])
  {
    parents <- names(which(adjacency_matrix[,node]==1))
    if(length(parents)>0)
    {
      parents_ind <- which(apply(parent_set_combinations[[node]][[length(parents)]],2,FUN=function(column) all(column==parents))==TRUE)
      score_nodes <- score_nodes + BGe_score_all_configs_node[[node]][[length(parents)]][parents_ind]
    } else {
    } # end if(length(parents)>0)
  } # end for(node in rownames(adjacency_matrix))
  ## CNV
  if(length(layers_def$omics)>1)
  {
    for(node in rownames(adjacency_matrix)[(ncol(omics[[layers_def$omics[1]]])+1):nrow(adjacency_matrix)])
    {
      score_nodes <- score_nodes + BGe_score_all_configs_node[[node]]
    } # end for(node in rownames(adjacency_matrix)[(ncol(omics[[layers_def$omics[1]]])+1):nrow(adjacency_matrix)])
  } # end if(ncol(omics$cnv)>0)
  return(BGe_score_net = score_nodes)
}

#' Epsilon  
#' @description
#' `epsilon` This function returns the epsilon value for each variable/node of the network. The sum of the epsilons of all variables/nodes in the network gives us the energy of given network.
#' @export
epsilon <- function(net, B_prior_mat)
{
  epsilon <- rep(NA,nrow(net))
  
  for(i in c(1:nrow(net)))
  {
    iter_feature <- rownames(net)[i]
    parent <- names(net[,iter_feature])[net[,iter_feature]==1]
    noparent <- setdiff(rownames(net),parent)
    epsilon[i] <- sum(1-B_prior_mat[parent,iter_feature]) + sum(B_prior_mat[noparent,iter_feature])
  } # end for i
  
  return(epsilon)
}

#' Sampling phase 
#' @description
#' `mcmc.simulation_sampling.phase` This function performs the final sampling of network structures with estimated hyperparameters. It if part of sampling_phase function.
#' @export
mcmc.simulation_sampling.phase <- function(first, last, sim_init, prob_mbr, B_prior_mat, omics, parent_set_combinations, BGe_score_all_configs_node, layers_def, len, thin, counting, energy_all_configs_node)
{
  source.net <- sim_init$nets[[length(sim_init$nets)]]
  beta.source <- sim_init$betas[[length(sim_init$betas)]]
  set.seed(sim_init$seed+counting)
  for (i in first:last)
  {
    ## method choice:
    sim_init$method_choice_saved[i] <- sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
    if(sim_init$method_choice_saved[i]=="MC3")
    {
      ### candidate.net definition: 
      candidate.net <- MC3(source_net_adjacency = source.net$adjacency,
                           layers_def =  layers_def, 
                           B_prior_mat = B_prior_mat, 
                           beta.source = beta.source, 
                           partition_func_UB_beta_source = sim_init$partition_func_UB_beta_source, 
                           omics = omics, 
                           parent_set_combinations = parent_set_combinations, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           annot = annot)
      ## Proposal distribution Q(..|..)
      # Q(source.net|candidate.net)
      # Q(candidate.net|source.net)
      candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
      source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
      candidate.net$likelihood <- candidate.net$likelihood_part + source.net$proposal.distr
      source.net$likelihood <- source.net$likelihood_part + candidate.net$proposal.distr
      sim_init$acceptance_saved[i] <- candidate.net$likelihood - source.net$likelihood

    } else {
      candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                           layers_def = layers_def, 
                           omics = omics, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           parent_set_combinations = parent_set_combinations)
      sim_init$acceptance_saved[i] <- candidate.net$acceptance
      ### source.net parameters necessary for MC3 method: 
      ## Marginal likelihood P(D|G) using ln(BGe score):
      candidate.net$BGe <-BGe_score(adjacency_matrix = candidate.net$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
      candidate.net$nbhd.size <- neighborhood_size(net = candidate.net$adjacency, layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
      ## Prior probability
      candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency, B_prior_mat = B_prior_mat))
      # partition_func_UB is already log transformed
      # the output of this function is also log transformed (natural logarithm)!
      candidate.net$prior <- (-beta.source$value*candidate.net$energy) - sim_init$partition_func_UB_beta_source
      candidate.net$likelihood_part <- candidate.net$BGe + candidate.net$prior

    } # end if(method.choice=="MC3")
    
    u <- log(runif(1))
    
    if (u < sim_init$acceptance_saved[i])
    {
      source.net <- candidate.net
    } # end if (u < acceptance_saved[i])
    sim_init$nets[[i]] <- source.net
  
    if(i==last)
    {
      sim_init$cpdags[[length(sim_init$cpdags)+1]] <- empty.graph(rownames(sim_init$nets[[i]]$adjacency))
      amat(sim_init$cpdags[[length(sim_init$cpdags)]]) <- sim_init$nets[[i]]$adjacency
      sim_init$cpdags[[length(sim_init$cpdags)]] <- cpdag(sim_init$cpdags[[length(sim_init$cpdags)]])
    } # end if(i==last)
  } # end for (i in first:last)
  
  return(sim_init)
}

#' Sampling phase 
#' @description
#' `sampling_phase` Now we apply 2 MCMC simulations and check the RMS value. After the burn-in period, we discard the values from the first half of this phase.
#' @export
sampling_phase <- function(second.adapt.phase_net, seed1, seed2, omics, layers_def, prob_mbr, thin, minseglen, burn_in) 
{
  rms <- c()
  seeds_res <- list(seed1=list(),seed2=list())
  # seed1 network is the last one from the 2nd adaption phase
  seeds_res$seed1$nets <- tail(second.adapt.phase_net$nets,1)
  seeds_res$seed1$nets[[1]]$nbhd.size <- neighborhood_size(net = seeds_res$seed1$nets[[1]]$adjacency, layers_def = layers_def, B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted, omics = omics)
  seeds_res$seed1$nets[[1]]$proposal.distr <- c()
  seeds_res$seed1$betas <- tail(second.adapt.phase_net$betas,1)
  seeds_res$seed1$betas[[1]]$prior <- seeds_res$seed1$nets[[1]]$prior
  seeds_res$seed1$partition_func_UB_beta_source <- second.adapt.phase_net$partition_func_UB_beta_source
  seeds_res$seed1$nets[[1]]$energy <- sum(epsilon(net = seeds_res$seed1$nets[[1]]$adjacency, B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted))
  seeds_res$seed1$nets[[1]]$prior <- ((-seeds_res$seed1$betas[[1]]$value)*seeds_res$seed1$nets[[1]]$energy) - seeds_res$seed1$partition_func_UB_beta_source
  seeds_res$seed1$nets[[1]]$likelihood_part <- seeds_res$seed1$nets[[1]]$BGe + seeds_res$seed1$nets[[1]]$prior
  seeds_res$seed2$nets[[1]]$likelihood <- c()
  seeds_res$seed2$nets[[1]]$acceptance <- c()
  seeds_res$seed2$nets[[1]]$edge_move <- c()
  seeds_res$seed1$acceptance_saved <- vector("numeric")
  seeds_res$seed1$method_choice_saved <- vector("numeric")
  seeds_res$seed1$layers <- second.adapt.phase_net$layers
  seeds_res$seed1$seed <- seed1
  seeds_res$seed1$cpdags <- list()
  
  # seed2 network is the empty
  seeds_res$seed2 <- seeds_res$seed1
  seeds_res$seed2$nets[[1]]$adjacency[seeds_res$seed2$nets[[1]]$adjacency==1] <- 0
  seeds_res$seed2$nets[[1]]$nbhd.size <- neighborhood_size(net = seeds_res$seed2$nets[[1]]$adjacency, layers_def = layers_def, B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted, omics = omics)
  seeds_res$seed2$nets[[1]]$energy <- sum(epsilon(net = seeds_res$seed2$nets[[1]]$adjacency, B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted))
  seeds_res$seed2$nets[[1]]$prior <- (-seeds_res$seed2$betas[[1]]$value*seeds_res$seed2$nets[[1]]$energy) - seeds_res$seed2$partition_func_UB_beta_source
  seeds_res$seed2$nets[[1]]$BGe <- BGe_score(adjacency_matrix = seeds_res$seed2$nets[[1]]$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = second.adapt.phase_net$partition_func_UB$parents_set_combinations, BGe_score_all_configs_node = second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node)
  seeds_res$seed2$nets[[1]]$likelihood_part <- seeds_res$seed2$nets[[1]]$BGe + seeds_res$seed2$nets[[1]]$prior
  seeds_res$seed2$betas[[1]]$prior <- seeds_res$seed2$nets[[1]]$prior
  seeds_res$seed2$seed <- seed2

  counting <- 0
  mcmc_sim_part_res <- lapply(seeds_res, 
                              FUN=function(list_l) mcmc.simulation_sampling.phase(first = 1,
                                                                                  last = thin,
                                                                                  sim_init = list_l,
                                                                                  prob_mbr = prob_mbr,
                                                                                  B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted,
                                                                                  omics = omics,
                                                                                  parent_set_combinations = second.adapt.phase_net$partition_func_UB$parents_set_combinations,
                                                                                  BGe_score_all_configs_node = second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node,
                                                                                  layers_def = layers_def,
                                                                                  len = seeds_res$seed1$betas[[1]]$len, 
                                                                                  counting = counting,
                                                                                  thin = thin,
                                                                                  energy_all_configs_node = second.adapt.phase_net$partition_func_UB$energy_all_configs_node))
  
  # remove duplicated cpdags
  cpdags1 <- mcmc_sim_part_res$seed1$cpdags
  cpdags2 <- mcmc_sim_part_res$seed2$cpdags
  cpdag_weights1 <- custom.strength(cpdags1, nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
  cpdag_weights2 <- custom.strength(cpdags2, nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
  cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
  cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]
  total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
  N <- nrow(total)
  dist_i <- abs(total$strength.x - total$strength.y)^2 / 2
  rms <- c(rms,sqrt(1/N*sum(dist_i)))
 
  while(length(mcmc_sim_part_res$seed1$nets)<(2*burn_in))
  {
    counting <- counting+1
    mcmc_sim_part_res <- lapply(mcmc_sim_part_res, 
                                FUN=function(list_l) mcmc.simulation_sampling.phase(first = length(list_l$nets)+1,
                                                                          last = length(list_l$nets)+thin,
                                                                          sim_init = list_l, 
                                                                          prob_mbr = prob_mbr,
                                                                          B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted,
                                                                          omics = omics,
                                                                          parent_set_combinations = second.adapt.phase_net$partition_func_UB$parents_set_combinations,
                                                                          BGe_score_all_configs_node = second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node,
                                                                          layers_def = layers_def,
                                                                          len = seeds_res$seed1$betas[[1]]$len, 
                                                                          counting = counting,
                                                                          thin = thin,
                                                                          energy_all_configs_node = second.adapt.phase_net$partition_func_UB$energy_all_configs_node))
    # remove duplicated cpdags
    cpdags1 <- unique(mcmc_sim_part_res$seed1$cpdags)
    cpdags2 <- unique(mcmc_sim_part_res$seed2$cpdags)
    cpdag_weights1 <- custom.strength(cpdags1, nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
    cpdag_weights2 <- custom.strength(cpdags2, nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
    cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
    cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]
    total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
    N <- nrow(total)
    dist_i <- abs(total$strength.x - total$strength.y)^2 / 2
    rms <- c(rms,sqrt(1/N*sum(dist_i)))
  } # end while(length(mcmc_sim_part_res$seed1$nets)<(2*burn_in))
  rms_strength <- abs(diff(rms))
  strength_threshold <- quantile(rms_strength, 0.75, na.rm = TRUE) # 3rd quartile of RMS strength is the threshold
  
  while(any(tail(rms_strength,minseglen/thin)>strength_threshold))
  {
    counting <- counting + 1
    mcmc_sim_part_res <- lapply(mcmc_sim_part_res, 
                                FUN=function(list_l) mcmc.simulation_sampling.phase(first = length(list_l$nets)+1,
                                                                          last = length(list_l$nets)+thin,
                                                                          sim_init = list_l, 
                                                                          prob_mbr = prob_mbr,
                                                                          B_prior_mat = second.adapt.phase_net$B_prior_mat_weighted,
                                                                          omics = omics,
                                                                          parent_set_combinations = second.adapt.phase_net$partition_func_UB$parents_set_combinations,
                                                                          BGe_score_all_configs_node = second.adapt.phase_net$partition_func_UB$BGe_score_all_configs_node,
                                                                          layers_def = layers_def,
                                                                          len = seeds_res$seed1$betas[[1]]$len,
                                                                          counting = counting,
                                                                          thin = thin,
                                                                          energy_all_configs_node = second.adapt.phase_net$partition_func_UB$energy_all_configs_node))
    # remove duplicated cpdags
    cpdags1 <- unique(mcmc_sim_part_res$seed1$cpdags)
    cpdags2 <- unique(mcmc_sim_part_res$seed2$cpdags)
    cpdag_weights1 <- custom.strength(cpdags1, nodes = bnlearn::nodes(cpdags1[[1]]), weights = NULL)
    cpdag_weights2 <- custom.strength(cpdags2, nodes = bnlearn::nodes(cpdags2[[1]]), weights = NULL)
    cpdag_weights1 <- cpdag_weights1[cpdag_weights1$direction>=0.5,]
    cpdag_weights2 <- cpdag_weights2[cpdag_weights2$direction>=0.5,]
    total <- merge(cpdag_weights1, cpdag_weights2, by = c("from","to"))
    N <- nrow(total)
    dist_i <- abs(total$strength.x - total$strength.y)^2 / 2
    rms <- c(rms,sqrt(1/N*sum(dist_i)))
    rms_strength <- abs(diff(rms))
  } # end while(any(tail(rms_strength,minseglen/thin)>strength_threshold))
  
  return(list(mcmc_sim_part_res = mcmc_sim_part_res, rms = rms))
}

#' Squared jumping of adaptive MCMC algorithm
#' @description
#' `squared_jumping` Squared jumping of adaptive MCMC algorithm is used to tune the variance of the beta parameter.
#' @export
squared_jumping <- function(second.adapt.phase_net, constant, beta_sd, fin, B_prior_mat, omics, parent_set_combinations, BGe_score_all_configs_node, layers_def, prob_mbr, annot)
{
  # beta variance after n iterations is found by calculating the variance of all past values generated by the chain from the point the trending stops in the transient phase to beta(n-1)
  # next beta is generated from N~(beta_n,beta_sd*constant)
  source.net <- second.adapt.phase_net$nets[[length(second.adapt.phase_net$nets)]]
  beta.source <- second.adapt.phase_net$betas[[length(second.adapt.phase_net$betas)]]
  start <- length(second.adapt.phase_net$nets)
  second.adapt.phase_net$iter_edges <- array(0, dim=c(dim(B_prior_mat),2), dimnames = list(rownames(B_prior_mat), colnames(B_prior_mat), c("frequency", "acceptance")))

  for(i in (length(second.adapt.phase_net$nets)+1):(length(second.adapt.phase_net$nets)+fin))
  {
    ## method choice:
    second.adapt.phase_net$method_choice_saved[i] <- sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
    if(second.adapt.phase_net$method_choice_saved[i]=="MC3")
    {
      ### candidate.net definition: 
      candidate.net <- MC3(source_net_adjacency = source.net$adjacency,
                           layers_def =  layers_def, 
                           B_prior_mat = B_prior_mat, 
                           beta.source = beta.source, 
                           partition_func_UB_beta_source = second.adapt.phase_net$partition_func_UB_beta_source, 
                           omics = omics, 
                           parent_set_combinations = parent_set_combinations, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           annot = annot)
      ## Proposal distribution Q(..|..)
      # Q(source.net|candidate.net)
      # Q(candidate.net|source.net)
      candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
      source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
      candidate.net$likelihood <- candidate.net$likelihood_part + source.net$proposal.distr
      source.net$likelihood <- source.net$likelihood_part + candidate.net$proposal.distr
      second.adapt.phase_net$acceptance_saved[i] <- candidate.net$likelihood - source.net$likelihood
      
      candidate_edge <- which(candidate.net$adjacency!=source.net$adjacency, arr.ind = TRUE)

      # now save which edge was chosen in this iteration and which operation was done (add/delete/reverse)
      if(candidate.net$edge_move=="reverse")
      {
        second.adapt.phase_net$iter_edges[candidate_edge[1,"row"], candidate_edge[1,"col"], 1] <- second.adapt.phase_net$iter_edges[candidate_edge[1,"row"], candidate_edge[1,"col"], 1] + 1
        second.adapt.phase_net$iter_edges[candidate_edge[2,"row"], candidate_edge[2,"col"], 1] <- second.adapt.phase_net$iter_edges[candidate_edge[2,"row"], candidate_edge[2,"col"], 1] + 1
      } else {
        second.adapt.phase_net$iter_edges[candidate_edge[,"row"], candidate_edge[,"col"], 1] <- second.adapt.phase_net$iter_edges[candidate_edge[,"row"], candidate_edge[,"col"], 1] + 1
      } # end if else (candidate.net$edge_move=="reverse")
      
      u <- log(runif(1))
      if (u < second.adapt.phase_net$acceptance_saved[i])
      {
        # the edge operation was accepted
        if(candidate.net$edge_move=="add")
        {
          second.adapt.phase_net$iter_edges[candidate_edge[,"row"], candidate_edge[,"col"], 2] <- second.adapt.phase_net$iter_edges[candidate_edge[,"row"], candidate_edge[,"col"], 2] + 1
        } else if (candidate.net$edge_move=="reverse")
        {
          second.adapt.phase_net$iter_edges[candidate_edge[source.net$adjacency[candidate_edge]==0,"row"],candidate_edge[source.net$adjacency[candidate_edge]==0,"col"],2] <- second.adapt.phase_net$iter_edges[candidate_edge[source.net$adjacency[candidate_edge]==0,"row"],candidate_edge[source.net$adjacency[candidate_edge]==0,"col"],2] + 1
        } # end if else if (candidate.net$edge_move=="add")
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      } else {
        # the edge operation was rejected
        if(candidate.net$edge_move=="reverse")
        {
          second.adapt.phase_net$iter_edges[candidate_edge[source.net$adjacency[candidate_edge]==1,"row"],candidate_edge[source.net$adjacency[candidate_edge]==1,"col"],2] <- second.adapt.phase_net$iter_edges[candidate_edge[source.net$adjacency[candidate_edge]==1,"row"],candidate_edge[source.net$adjacency[candidate_edge]==1,"col"],2] + 1
        } else if(candidate.net$edge_move=="delete")
        {
          second.adapt.phase_net$iter_edges[candidate_edge[,"row"], candidate_edge[,"col"], 2] <- second.adapt.phase_net$iter_edges[candidate_edge[,"row"], candidate_edge[,"col"], 2] + 1
        }# end if else if(candidate.net$edge_move=="reverse")
      } # end if else(u < second.adapt.phase_net$acceptance_saved[i])
      second.adapt.phase_net$nets[[i]] <- source.net
      
    } else {
      
      candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                           layers_def = layers_def, 
                           omics = omics, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           parent_set_combinations = parent_set_combinations)
      second.adapt.phase_net$acceptance_saved[i] <- candidate.net$acceptance
      ### source.net parameters necessary for MC3 method: 
      ## Marginal likelihood P(D|G) using ln(BGe score):
      candidate.net$BGe <-BGe_score(adjacency_matrix = candidate.net$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
      candidate.net$nbhd.size <- neighborhood_size(net = candidate.net$adjacency, layers_def = layers_def, B_prior_mat = B_prior_mat, omics = omics)
      ## Prior probability
      candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency, B_prior_mat = B_prior_mat))
      # partition_func_UB is already log transformed
      # the output of this function is also log transformed (natural logarithm)!
      candidate.net$prior <- (-beta.source$value*candidate.net$energy) - second.adapt.phase_net$partition_func_UB_beta_source
      candidate.net$likelihood_part <- candidate.net$BGe + candidate.net$prior
      
      u <- log(runif(1))
      if (u < second.adapt.phase_net$acceptance_saved[i])
      {
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      } # end if (u < acceptance_saved[i])
      second.adapt.phase_net$nets[[i]] <- source.net
      
    } # end if(method.choice=="MC3")
    
    ### MH for beta with fixed network structure ###
    beta.candidate <- list(value = rnorm(1, mean = beta.source$value, sd = beta_sd*constant), prior = c(), len = beta_sd*constant)
    if(beta.candidate$value < 0.5)
    {
      beta.candidate$value <- 0.5
    } # end if(beta.candidate$value < 0.5)
    
    partition_func_UB_beta_candidate <- sum(mapply(second.adapt.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.candidate$value*x)))
    beta.candidate$prior <- (-beta.candidate$value*source.net$energy) - partition_func_UB_beta_candidate
    
    second.adapt.phase_net$acceptance_beta_saved[i] <- beta.candidate$prior - beta.source$prior
    u_beta <- log(runif(1))
    
    if (u_beta < second.adapt.phase_net$acceptance_beta_saved[i])
    {
      beta.source <- beta.candidate
      second.adapt.phase_net$partition_func_UB_beta_source <- partition_func_UB_beta_candidate
    } # end if (u_beta < second.adapt.phase_net$acceptance_beta_saved[i])
    second.adapt.phase_net$betas[[i]] <- beta.source
  } # end for(i in (length(second.adapt.phase_net$nets)+1):(length(second.adapt.phase_net$nets)+fin))
  return(second.adapt.phase_net)
}

#' Second adaption phase
#' @description
#' `variance_target` This phase identifies the proposal distribution that has a similar covariance structure with the target distribution. This is part of second_adapt_phase.
#' @export
variance_target <- function(transient.phase_net, constant, fin, B_prior_mat, omics, parent_set_combinations, BGe_score_all_configs_node, layers_def, prob_mbr, annot)
{
  # beta variance after n iterations is found by calculating the variance of all past values generated by the chain from the point the trending stops in the transient phase to beta(n-1)
  beta_sd <- sd(mapply(tail(transient.phase_net$betas,1000), FUN=function(list) list$value)[-1000])
  # next beta is generated from N~(beta_n,beta_sd*constant)
  source.net <- transient.phase_net$nets[[length(transient.phase_net$nets)]]
  beta.source <- transient.phase_net$betas[[length(transient.phase_net$betas)]]
  start <- length(transient.phase_net$nets)

  
  # We employ a standard MCMC algorithm with proposals determined by the 1st adaption/transient/2nd adaption phase.
  # Check if the acceptance rate is too low (less than 0.02)
  for(i in (start+1):(start+fin))
  {
    ## method choice:
    transient.phase_net$method_choice_saved[i] <- sample(x = c("MC3", "MBR"), size = 1, prob = c(1-prob_mbr, prob_mbr))
    if(transient.phase_net$method_choice_saved[i]=="MC3")
    {
      ### candidate.net definition: 
      candidate.net <- MC3(source_net_adjacency = source.net$adjacency,
                           B_prior_mat = B_prior_mat, 
                           beta.source = beta.source, 
                           partition_func_UB_beta_source = transient.phase_net$partition_func_UB_beta_source, 
                           omics = omics, 
                           parent_set_combinations = parent_set_combinations, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           layers_def = layers_def, 
                           annot = annot)
      ## Proposal distribution Q(..|..)
      # Q(source.net|candidate.net)
      # Q(candidate.net|source.net)
      candidate.net$proposal.distr <- log(1/source.net$nbhd.size)
      source.net$proposal.distr <- log(1/candidate.net$nbhd.size)
      candidate.net$likelihood <- candidate.net$likelihood_part + source.net$proposal.distr
      source.net$likelihood <- source.net$likelihood_part + candidate.net$proposal.distr
      transient.phase_net$acceptance_saved[i] <- candidate.net$likelihood - source.net$likelihood
      
      candidate_edge <- which(candidate.net$adjacency!=source.net$adjacency, arr.ind = TRUE)
     
      u <- log(runif(1))
      if (u < transient.phase_net$acceptance_saved[i])
      {
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      } # end if (u < transient.phase_net$acceptance_saved[i])
      transient.phase_net$nets[[i]] <- source.net
      
    } else {
      
      candidate.net <- MBR(source_net_adjacency = source.net$adjacency, 
                           layers_def = layers_def, 
                           omics = omics, 
                           BGe_score_all_configs_node = BGe_score_all_configs_node, 
                           parent_set_combinations = parent_set_combinations)
      transient.phase_net$acceptance_saved[i] <- candidate.net$acceptance
      ### source.net parameters necessary for MC3 method: 
      ## Marginal likelihood P(D|G) using ln(BGe score):
      candidate.net$BGe <-BGe_score(adjacency_matrix = candidate.net$adjacency, omics = omics, layers_def = layers_def, parent_set_combinations = parent_set_combinations, BGe_score_all_configs_node = BGe_score_all_configs_node)
      candidate.net$nbhd.size <- neighborhood_size(net = candidate.net$adjacency, B_prior_mat = B_prior_mat, layers_def = layers_def, omics = omics)
      ## Prior probability
      candidate.net$energy <- sum(epsilon(net = candidate.net$adjacency, B_prior_mat = B_prior_mat))
      # partition_func_UB is already log transformed
      # the output of this function is also log transformed (natural logarithm)!
      candidate.net$prior <- (-beta.source$value*candidate.net$energy) - transient.phase_net$partition_func_UB_beta_source
      candidate.net$likelihood_part <- candidate.net$BGe + candidate.net$prior
      
      u <- log(runif(1))
      if (u < transient.phase_net$acceptance_saved[i])
      {
        source.net <- candidate.net
        beta.source$prior <- source.net$prior
      } # end if (u < acceptance_saved[i])
      transient.phase_net$nets[[i]] <- source.net
      
      partition_func_UB_beta_source <- sum(mapply(transient.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.source$value*x)))

    } # end if(method.choice=="MC3")
    
    ### MH for beta with fixed network structure ###
    beta.candidate <- list(value = rnorm(1, mean = beta.source$value, sd = beta_sd*constant), prior = c(), len = beta_sd*constant)
    if(beta.candidate$value < 0.5)
    {
      beta.candidate$value <- 0.5
    } # end if(beta.candidate$value < 0.5)
    
    partition_func_UB_beta_candidate <- sum(mapply(transient.phase_net$energy_all_configs_node,FUN=function(x) logSumExp(-beta.candidate$value*x)))
    beta.candidate$prior <- (-beta.candidate$value*source.net$energy) - partition_func_UB_beta_candidate
    
    transient.phase_net$acceptance_beta_saved[i] <- beta.candidate$prior - beta.source$prior
    u_beta <- log(runif(1))
    
    if (u_beta < transient.phase_net$acceptance_beta_saved[i])
    {
      beta.source <- beta.candidate
      transient.phase_net$partition_func_UB_beta_source <- partition_func_UB_beta_candidate
    } # end if (u_beta < transient.phase_net$acceptance_beta_saved[i])
    transient.phase_net$betas[[i]] <- beta.source
  } # end for(i in (start+1):(start+200))
  
  acceptance.trace_betas <- unlist(lapply(tail(transient.phase_net$betas, 200),FUN=function(list) list$prior))
  acceptance.trace_betas <- c(1,acceptance.trace_betas[1:(length(acceptance.trace_betas)-1)] - acceptance.trace_betas[2:length(acceptance.trace_betas)])
  acceptance.trace_betas[acceptance.trace_betas!=0] <- 1
  acceptance.rate_betas <- sum(acceptance.trace_betas==1)/length(acceptance.trace_betas)
  
  return(list(variance.target_net = transient.phase_net, acceptance.rate_betas = acceptance.rate_betas, beta_sd = beta_sd))
}
