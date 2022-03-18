# This function is used to calculate the ln of the target distribution for a 
# given pen with chains equal to the number of cattle.

pi_target <- function(chains, id_pen, measures_rams, measures_fec, alpha, beta, 
                      m, nu, theta_r, theta_f){
  
  trials <- dim(chains)[1]
  
  cattle <- dim(chains)[2]
  
  # m is equal to the mean infectious period minus 1
  
  priors <- sum(dgamma(x = alpha, shape = 1, rate = 1, log = TRUE), 
                dgamma(x = beta, shape = 1, rate = 1, log = TRUE),
                dgamma(x = m, shape = 0.01, rate = 0.01, log = TRUE),
                dbeta(x = nu, shape1 = 1, shape2 = 1, log = TRUE),
                dbeta(x = theta_r, shape1 = 1, shape2 = 1, log = TRUE),
                dbeta(x = theta_f, shape1 = 1, shape2 = 1, log = TRUE))
  
  # probability of joint event rams/fec conditional on susceptible
  cond_jointprob_susceptible <- c(1, 0, 0, 0)
  
  # probability of joint event rams/fec conditional on infected
  cond_jointprob_infected <- c((1 - theta_r) * (1 - theta_f),
                               theta_r * (1 - theta_f),
                               (1 - theta_r) * theta_f,
                               theta_r * theta_f)
  
  cond_jointprob <- cbind(cond_jointprob_susceptible, cond_jointprob_infected)
  
  transition_log <- matrix(data = NA, nrow = trials, ncol = cattle)
  
  emissions_log <- matrix(data = NA, nrow = trials, ncol = cattle)
  
  for(cc in 1:cattle){
    
    # indicator vectors for joint outcome
    vec_ind <- cbind(ifelse(test = measures_rams[, cc, id_pen] == 0 & measures_fec[, cc, id_pen] == 0, yes = 1, no = 0),
                     ifelse(test = measures_rams[, cc, id_pen] == 1 & measures_fec[, cc, id_pen] == 0, yes = 1, no = 0),
                     ifelse(test = measures_rams[, cc, id_pen] == 0 & measures_fec[, cc, id_pen] == 1, yes = 1, no = 0),
                     ifelse(test = measures_rams[, cc, id_pen] == 1 & measures_fec[, cc, id_pen] == 1, yes = 1, no = 0))
    
    total_infected <- apply(X = chains[, -cc, id_pen], MARGIN = 1, 
                            FUN = sum)
    
    for(tt in 1:trials){
      
      if(tt == 1){
        transition_log[tt, cc] <- log((1 - chains[tt, cc, id_pen]) * (1 - nu) + 
                                       chains[tt, cc, id_pen] * nu)
        
        emissions_log[tt, cc] <-  log(vec_ind[tt,] %*% cond_jointprob[,(chains[tt, cc, id_pen] + 1)])
      }
      
      else{
        # model probability of 0|0, 1|0; 
        #                      0|1, 1|1
        model_transition <- rbind(
          c(exp(-alpha - beta * total_infected[tt-1]), 1 - exp(-alpha - beta * total_infected[tt-1])),
          c(1 / (m+1), m / (m+1)))
        
        transition_log[tt, cc] <- log(model_transition[(chains[tt,cc,id_pen] + 1), (chains[(tt-1), cc, id_pen] + 1)])
        
        emissions_log[tt, cc] <-  log(vec_ind[tt,] %*% cond_jointprob[,(chains[tt, cc, id_pen] + 1)])
      }
    }
    


  }
  
  em <- sum(emissions_log)
  
  ts <- sum(transition_log)
  
  #print(ts)
  
  #print(em)
  
  #print(priors)
  
  return(priors + ts + em)
  
}
