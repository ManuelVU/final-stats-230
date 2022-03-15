# Metropolis Hastings for a single chain, theta is the vector of parameter
# values organized as alpha, beta, mu-1, nu, theta_r and theta_f

mh_chain <- function(current_chains, tests_rams, tests_fec, cattle, pen, theta, 
                     log_q_current = NULL, log_pi_current = NULL){
  source(file = "src/forward-filter.R")
  source(file = "src/backward-sampling.R")
  source(file = "src/evaluate-target.R")
  
  proposed_chain <- current_chains
  
  prob_current <- forward_filter(id_number = cattle, id_pen = pen, 
                                 pen_states = current_chains, 
                                 external_inf_rate = theta[1], 
                                 within_inf_rate = theta[2],
                                 prob_infect_1 = theta[4],
                                 infectious_period = theta[3] + 1, 
                                 sens_rams = theta[5], sens_fec = theta[6],
                                 meas_rams = tests_rams, meas_fec = tests_fec)
  
  if(missing(log_q_current)){
    source(file = "src/eval-initial-proposal.R")
    log_q_current <- evaluate_q(proposal_pinfected = prob_current[, 2],
                                pen_states = current_chains, 
                                id_number = cattle, id_pen = pen, 
                                external_inf_rate = theta[1], 
                                within_inf_rate = theta[2], 
                                infectious_period = theta[3]+1)
  }
  
  new_chain <- backward_sampling(proposal_pinfected = prob_current[,2],
                                 pen_states = current_chains, 
                                 id_number = cattle, id_pen = pen, 
                                 external_inf_rate = theta[1], 
                                 within_inf_rate = theta[2],
                                 infectious_period = theta[3] + 1)
  
  proposed_chain[, cattle, pen] <- new_chain$proposal_state
  
  log_q_proposal <- new_chain$loglike
  
  log_pi_proposal <- pi_target(chains = proposed_chain, id_pen = pen, 
                               measures_rams = tests_rams, 
                               measures_fec = tests_fec,
                               alpha = theta[1], beta = theta[2], m = theta[3], 
                               nu = theta[4], theta_r = theta[5], 
                               theta_f = theta[6])
  
  if(missing(log_pi_current)){
    log_pi_current <- pi_target(chains = current_chains, id_pen = pen, 
                                measures_rams = tests_rams, 
                                measures_fec = tests_fec,
                                alpha = theta[1], beta = theta[2], m = theta[3], 
                                nu = theta[4], theta_r = theta[5], 
                                theta_f = theta[6])
  }
  
  log_u <- log(runif(n = 1))
  
  a <- min(0, log_q_current - log_q_proposal + log_pi_proposal - log_pi_current)
  
  if(log_u <= a){
    post_sample <- proposed_chain,
    log_q <- log_q_proposal
    log_pi <- log_pi_proposal
    accepted <- 1
  }
  
  else{
    post_sample <- current_chains,
    log_q <- log_q_current
    log_pi <- log_pi_current
    accepted <- 0
  }
  
  results <- list()
  
  reults$mhchain <- post_sample
  results$log_q <- log_q
  results$log_pi <- log_pi
  results$accepted <- accepted
  
  return(results)
  
}