# This function is used to calculate the ln of the target distribution for a 
# given pen with chains equal to the number of cattle.

pi_target <- function(chains, id_pen, measures_rams, measures_fec, alpha, beta, 
                      m, nu, theta_r, theta_f){
  
  trials <- dim(chains)[1]
  
  cattle <- dim(chains)[2]
  
  # m is equal to the mean infectious period minus 1
  
  priors <- sum(dgamma(x = alpha, shape = 1, rate = 1, log = TRUE), 
                dgamma(x = beta, shape = 1, rate = 1, log = TRUE),
                dgamma(x = m + 1, shape = 0.01, rate = 0.01, log = TRUE),
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
        model_transition <- rbind(c(exp(-alpha - beta * total_infected[tt-1]), 
                                    1 - exp(-alpha - beta * total_infected[tt-1])),
                                  c(1 / m, 
                                    m - 1 / m))
        
        transition_log[tt, cc] <- log(model_transition[(chains[(tt-1), cc, id_pen] + 1), (chains[tt,cc,id_pen] + 1)])
        
        emissions_log[tt, cc] <-  log(vec_ind[tt,] %*% cond_jointprob[,(chains[tt, cc, id_pen] + 1)])
      }
    }
  }
  
  em <- sum(emissions_log)
  
  ts <- sum(transition_log)
  
  if(is.nan(priors + ts + em)){
    return(-Inf)
  }
  else{
    return(priors + ts + em)
  }
}


# loglik joint conditional ------------------------------------------------
pi_target2 <- function(chains, id_cow, id_pen, measures_rams, measures_fec, theta_r, 
                      theta_f, nu, alpha, beta, m){
  
  trials <- dim(chains)[1]
  
  cows <- dim(chains)[2]
  
  chain <- chains[, id_cow, id_pen]
  
  y_r <- measures_rams[, id_cow, id_pen]
  
  y_f <- measures_fec[, id_cow, id_pen]
  
  # evaluate log-likelihood of emissions
  
  n_00 <- sum(y_r == 0 & y_f == 0 & chain == 1)
  n_01 <- sum(y_r == 0 & y_f == 1 & chain == 1)
  n_10 <- sum(y_r == 1 & y_f == 0 & chain == 1)
  n_11 <- sum(y_r == 1 & y_f == 1 & chain == 1)
  
  llike_em <- n_00 * log((1 - theta_r) * (1 - theta_f)) +
              n_01 * log((1 - theta_r) * theta_f) +
              n_10 * log(theta_r * (1 - theta_f)) +
              n_11 * log(theta_r * theta_f)
  
  # evaluate log likelihood of transitions
  infected_at1 <- sum(chains[1,,id_pen] == 1)
  
  llike_trans <- c()
  
  llike_trans[1] <- (cows - infected_at1) * (1 - nu) + infected_at1 * nu
  
  for(t in 2:trials){
    n_00 <- sum(chains[t-1,, id_pen] == 0 & chains[t,, id_pen] == 0)
    n_01 <- sum(chains[t-1,, id_pen] == 0 & chains[t,, id_pen] == 1)
    n_10 <- sum(chains[t-1,, id_pen] == 1 & chains[t,, id_pen] == 0)
    n_11 <- sum(chains[t-1,, id_pen] == 1 & chains[t,, id_pen] == 1)
    
    total_infatprev <- sum(chains[t-1, , id_pen] == 1)
    
    llike_trans[t] <- n_00 * (-alpha -beta * total_infatprev) + 
                      n_01 * log(1 - exp(-alpha -beta * total_infatprev)) - 
                      (n_11 + n_10) * log(m[3]+1) + 
                      n_11 * log(m[3])
  }
  
  return(sum(llike_trans, na.rm = TRUE) + llike_em)
}
