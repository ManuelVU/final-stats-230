# Hamiltonian Monte Carlo, function is working but is unstable. The parameters
# that seem to work most of the time are: leaps = 30 and stepsize = 0.03


hamiltonian_mc <- function(u, grad_u, current_theta, prior_theta, 
                           hidden_current, leaps, stepsize){
  
  trials <- dim(hidden_current)[1]
  
  pens <- dim(hidden_current)[3]
  
  q <- current_theta
  
  p <- rnorm(n = length(q), mean = 0, sd = 1)
  
  current_p <- p
  
  p <- p - stepsize * grad_u(cr_theta = q, cr_chains = hidden_current, initial_theta = prior_theta) / 2
  
  for(l in 1:leaps){
    q <- q + stepsize * p
    
    if(l != leaps){
      p <- p - stepsize * grad_u(cr_theta = q, cr_chains = hidden_current, initial_theta = prior_theta)
    }
  }
  
  p <- p - stepsize * grad_u(cr_theta = q, cr_chains = hidden_current, initial_theta = prior_theta) / 2
  
  p <- - p
  
  current_u <- u(cr_theta = current_theta, cr_chains = hidden_current, initial_theta = prior_theta)
  
  current_k <- sum(current_p^2) / 2
  
  proposed_u <- u(cr_theta = q, cr_chains = hidden_current, initial_theta = prior_theta)
  
  proposed_k <- sum(p^2) / 2
  
  if(runif(n = 1) < exp(current_u - proposed_u + current_k - proposed_k)){
    return(q)
  }
  else{
    return(current_theta)
  }
}
