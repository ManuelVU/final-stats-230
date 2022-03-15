# Forward filter algorithm, returns the proposal mass function used by the 
# backward sample algorithm to generate a proposal for the chain.

forward_filter <- function(id_number, id_pen, pen_states, external_inf_rate, 
                           within_inf_rate, prob_infect_1, infectious_period, 
                           sens_rams, sens_fec, meas_rams, meas_fec){
  
  trials <- dim(pen_states)[1]
  
  # take the measures of the subject under analysis
  meas_rams <- meas_rams[, id_number, id_pen]
  meas_fec <- meas_fec[, id_number, id_pen]
  
  # matrix to save filter probabilities
  filter_prob <- matrix(data = NA, nrow = trials, ncol = 2)
  
  # matrix to save probosal mass function
  proposal_mass <- matrix(data = NA, nrow = trials, ncol = 2)
  
  total_infected <- apply(X = pen_states[, -id_number, id_pen], MARGIN = 1, 
                          FUN = sum)
  
  # probability of joint event rams/fec conditional on susceptible
  cond_jointprob_susceptible <- c(1, 0, 0, 0)
  
  # probability of joint event rams/fec conditional on infected
  cond_jointprob_infected <- c((1 - sens_rams) * (1 - sens_fec),
                               sens_rams * (1 - sens_fec),
                               (1 - sens_rams) * sens_fec,
                               sens_rams * sens_fec)
  
  # indicator vectors for joint outcome
  vec_ind <- cbind(ifelse(test = meas_rams == 0 & meas_fec == 0, yes = 1, no = 0),
                   ifelse(test = meas_rams == 1 & meas_fec == 0, yes = 1, no = 0),
                   ifelse(test = meas_rams == 0 & meas_fec == 1, yes = 1, no = 0),
                   ifelse(test = meas_rams == 1 & meas_fec == 1, yes = 1, no = 0))
  
  # filter probability at time one for susceptible and infected condition
  filter_prob[1,] <- c((1 - prob_infect_1) * vec_ind[1,] %*% cond_jointprob_susceptible, 
                       prob_infect_1 * vec_ind[1,] %*% cond_jointprob_infected)
  
  # proposal probability mass function at trial 1
  proposal_mass[1,] <- filter_prob[1,]/sum(filter_prob[1,])
  
  for(t in 2:trials){
    
  # model probability of 0|0, 1|0; 
  #                      0|1, 1|1
    transition <- rbind(
      c(exp(-external_inf_rate - within_inf_rate * total_infected[t-1]),
        1 - exp(-external_inf_rate - within_inf_rate * total_infected[t-1])),
      c(1 / infectious_period,
        (infectious_period - 1) / infectious_period))
    
    # part update equation 6
    filter_prob[t, ] <- c(transition[,2] %*% proposal_mass[t-1, ],
                          transition[,1] %*% proposal_mass[t-1, ])
    
    # finish update equation 7
    filter_prob[t, ] <- c(filter_prob[t, 1] *  vec_ind[t,] %*% cond_jointprob_susceptible,
                          filter_prob[t, 2] *  vec_ind[t,] %*% cond_jointprob_infected)
    
    # proposal mass function at trial t
    proposal_mass[t, ] <- filter_prob[t, ] / sum(filter_prob[t,])
    
  }
  return(proposal_mass)
}
