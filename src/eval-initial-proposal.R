# This function is design to evaluate q at the initial proposal only, after 
# future evaluations can be taken from the backward sample algorithm

evaluate_q <- function(proposal_pinfected, pen_states, id_number, id_pen, 
                       external_inf_rate, within_inf_rate, 
                       infectious_period){
  
  trials <- length(proposal_pinfected)
  
  total_infected <- apply(X = pen_states[, -id_number, id_pen], MARGIN = 1, 
                          FUN = sum)
  
  proposal_state <- pen_states[, id_number, id_pen]
  
  proposal_like <- rep(0, trials)
  
  proposal_like[trials] <- ifelse(test = proposal_state[trials] == 1,
                                  yes = proposal_pinfected[trials],
                                  no = 1 - proposal_pinfected[trials])
  
  for(t in (trials - 1):1){
    # transition to proposal_state[trials] from x_t = 1 * proposal_pinfected
    # transition to proposal_state[trials] from x_t = 0 and x_t = 1 * 1-prop_pinf and prop_pinf respectively
    
    # model probability of 0|0, 1|0; 
    #                      0|1, 1|1
    transition <- rbind(
      c(exp(-external_inf_rate - within_inf_rate * total_infected[t]),
        1 - exp(-external_inf_rate - within_inf_rate * total_infected[t])),
      c(1 / infectious_period,
        (infectious_period - 1) / infectious_period))
    
    prob_state_t <-  transition[2,(proposal_state[t+1] + 1)] * proposal_pinfected[t] / 
      sum(c(transition[1,(proposal_state[t+1] + 1)] * (1 - proposal_pinfected[t]),
            transition[2,(proposal_state[t+1] + 1)] * proposal_pinfected[t]))
    
    proposal_like[t] <- ifelse(test = proposal_state[t] == 1,
                               yes = proposal_pinfected[t],
                               no = 1 - proposal_pinfected[t])
    
  }
  
  return(sum(log(proposal_like)))
  
}
