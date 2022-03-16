# This function calculates the gradient of the logarithm of the joint posterior 
# distribution, function is used as grad_u in the hamiltonian_mc function.

# transitions shoudl be organize in an array with trials, pens and then 
# organized as follows: x_t-1 | x_t 
# 0|0, 0|1, 1|0 and 1|1.
# transitions will have dimentions trials-1, pens, 4

# infected totals should be organized as a matrix with rows equal to the number 
# of trials and columns equal to the number of pens.

# initial_theta is the matrix of initial values foe the prior distributions, 
# each parameter on a single row organized as alpha, beta, m. Dimentions should
# be 3, 2.

gradient_log_joint_post <- function(cr_theta, transitions, infected_totals,
                                    initial_theta){
  
  trials <- dim(transitions)[1]
  
  n_pens <- dim(infected_totals)[2]
  
  at <- exp(cr_theta[1])
  
  bt <- exp(cr_theta[2])
  
  tmp_alpha <- matrix(data = NA, nrow = trials-1, ncol = n_pens)
  
  tmp_beta <- matrix(data = NA, nrow = trials-1, ncol = n_pens)
  
  tmp_m <- matrix(data = NA, nrow = trials-1, ncol = n_pens)
  
  for(p in 1:n_pens){
    for(t in 1:(trials-1)){
      eratio <- exp(-at - bt * infected_totals[t, p]) / (1 - exp(-at - bt * infected_totals[t, p]))
      
      tmp_alpha[t, p] <- transitions[t, p, 2] * at * eratio - transitions[t, p, 1] * at
      
      tmp_beta[t, p] <- transitions[t, p, 2] * bt * infected_totals[t, p] * eratio - transitions[t, p, 1] * bt * infected_totals[t, p]
      
      tmp_m[t, p] <- transitions[t, p, 4] / (cr_theta[3] + 1)^2 - transitions[t, p, 3] / (cr_theta[3] + 1)^2 
  
    }
  }
  
  d_ln_a <- sum(tmp_alpha) + initial_theta[1, 1] - initial_theta[1, 2] * exp(cr_theta[1])
  d_ln_b <- sum(tmp_beta) + initial_theta[2, 1] - initial_theta[2, 2] * exp(cr_theta[2])
  d_ln_m <- sum(tmp_m) + 1 + (initial_theta[3, 1] - 1) / cr_theta[3] - initial_theta[3, 2]
  
  return(c(d_ln_a, d_ln_b, d_ln_m))
  
}
