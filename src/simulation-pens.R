# Simulation of multiple pens using the susceptible-infected-susceptible model.
# Pens are considered to be independent. 

# trials: number of time points to be simulated.
# n.cattle: number of chains in each pen.
# n.pens: number of independent pens.
# alpha: external infection rate. 
# beta: within-pen infection rate.
# mean_infection: mean infectious period.
# initial_pinfect: probability of infection on first trial.
# theta_r: sensitivity of the RAMS test (p(r^c_t = 1 | x^c_t = 1)).
# theta_f: sensitivity of the fecal test (p(f^c_t = 1 | x^c_t = 1)).
# prob_test: probability that at any given trial a measure of both tests will be 
#         taken.

# Variables alpha, beta, mean_infection, initial_pinfect, theta_r, theta_f and 
# prob_test have default values set equal to the ones used on the original 
# simulations of the paper.

sis_simulations_mp_pens <- function(trials, n.cattle, n.pens, alpha = 0.009,
                                    beta = 0.01, mean_infection = 9,
                                    initial_pinfect = 0.1, theta_r = 0.8,
                                    theta_f = 0.5, prob_test = 2/7){
  results <- list()
  results$hstates <- array(data = NA, dim = c(trials, n.cattle, n.pens))
  results$resp$rams <- array(data = NA, dim = c(trials, n.cattle, n.pens))
  results$resp$fec <- array(data = NA, dim = c(trials, n.cattle, n.pens))
  for(p in 1:n.pens){
    results$hstates[,,p] <- sis_sim_hidenstates(samples = trials, al = alpha,
                                                bt = beta,
                                                m_inf = mean_infection,
                                                p_initial =  initial_pinfect,
                                                nc = n.cattle)
    tmp <- sis_sim_observations(h_states = results$hstates[,,p],
                                p_test = prob_test,
                                th_r = theta_r, th_f = theta_f)
    
    results$resp$rams[,,p] <- tmp[,,1]
    results$resp$fec[,,p] <- tmp[,,2]
    
  }
  return(results)
}
