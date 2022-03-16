# Hamiltonian Monte Carlo, function is working but is unstable. The parameters
# that seem to work most of the time are: leaps = 30 and stepsize = 0.03


hamiltonian_mc <- function(u, grad_u, current_theta, prior_theta, 
                           hidden_current, leaps, stepsize){
  
  trials <- dim(hidden_current)[1]
  
  pens <- dim(hidden_current)[3]
  
  total_infected <- apply(X = hidden_current, MARGIN = c(1,3), FUN = sum)
  
  transitions_added <- array(data = NA, dim = c(trials - 1, pens, 4))
  
  for(t in 2:trials){
    transitions_added[(t-1),,1] <- apply(X = ifelse(test = hidden_current[t-1, , ] == 0 & 
                                                           hidden_current[t, , ] == 0, 
                                                    yes = 1, no = 0), 
                                         MARGIN = 2, FUN = sum)
    
    transitions_added[(t-1),,2] <- apply(X = ifelse(test = hidden_current[t-1, , ] == 0 & 
                                                           hidden_current[t, , ] == 1, 
                                                    yes = 1, no = 0), 
                                         MARGIN = 2, FUN = sum)
    
    transitions_added[(t-1),,3] <- apply(X = ifelse(test = hidden_current[t-1, , ] == 1 & 
                                                           hidden_current[t, , ] == 0, 
                                                    yes = 1, no = 0), 
                                         MARGIN = 2, FUN = sum)
    
    transitions_added[(t-1),,4] <- apply(X = ifelse(test = hidden_current[t-1, , ] == 1 & 
                                                           hidden_current[t, , ] == 1, 
                                                    yes = 1, no = 0), 
                                         MARGIN = 2, FUN = sum)
  }
  
  q <- current_theta
  
  p <- rnorm(n = length(q), mean = 0, sd = 1)
  
  current_p <- p
  
  p <- p + stepsize * grad_u(cr_theta = q, transitions = transitions_added, infected_totals = total_infected, initial_theta = prior_theta) / 2
  
  for(l in 1:leaps){
    q <- q - stepsize * p
    
    if(l < leaps){
      p <- p - stepsize * grad_u(cr_theta = q, transitions = transitions_added, infected_totals = total_infected, initial_theta = prior_theta)
    }
  }
  
  p <- p - stepsize * grad_u(cr_theta = q, transitions = transitions_added, infected_totals = total_infected, initial_theta = prior_theta) / 2
  
  p <- - p
  
  current_u <- u(cr_theta = current_theta, transitions = transitions_added, infected_totals = total_infected, initial_theta = prior_theta)
  
  current_k <- sum(current_p^2) / 2
  
  proposed_u <- u(cr_theta = q, transitions = transitions_added, infected_totals = total_infected, initial_theta = prior_theta)
  
  proposed_k <- sum(p^2) / 2
  
  mh_ratio <- min(1, exp(current_u - proposed_u + current_k - proposed_k))
  
  if(runif(n = 1) < mh_ratio){
    return(q)
  }
  else{
    return(current_theta)
  }
}
