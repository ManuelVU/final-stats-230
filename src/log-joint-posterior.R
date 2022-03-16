# This function calculates the logarithm of the joint posterior distribution, 
# function is used as u in the hamiltonian_mc function.

# transitions shoudl be organize in an array with trials, pens and then 
# organized as follows: x_t-1 | x_t 
# 0|0, 0|1, 1|0 and 1|1.
# transitions will have dimentions trials-1, pens, 4

# infected totals should be organized as a matrix with rows equal to the number 
# of trials and columns equal to the number of pens.

# initial_theta is the matrix of initial values foe the prior distributions, 
# each parameter on a single row organized as alpha, beta, m. Dimentions should
# be 3, 2.

log_joint_post <- function(cr_theta, transitions, infected_totals,
                           initial_theta){
  
  trials <- dim(transitions)[1]
  
  n_pens <- dim(infected_totals)[2]
  
  # matrix to temporarily store values of the log joint posterior
  tmp <- matrix(data = NA, nrow = trials - 1, ncol = n_pens)
  
  for(p in 1:n_pens){
    for(t in 1:(trials-1)){
      tmp[t, p] <- transitions[t, p, 1] * (- exp(cr_theta[1]) - exp(cr_theta[2]) * infected_totals[t, p]) +
                   transitions[t, p, 2] * (log(1 - exp(- exp(cr_theta[1]) - exp(cr_theta[2]) * infected_totals[t, p]))) +
                   transitions[t, p, 3] * (1 / (cr_theta[3] + 1)) + 
                   transitions[t, p, 4] * (cr_theta[3] / (cr_theta[3] + 1))
    }
  }
  
  ln_priors <- dgamma(x = exp(cr_theta[1]), shape = initial_theta[1, 1], rate = initial_theta[1, 2]) +
               dgamma(x = exp(cr_theta[2]), shape = initial_theta[2, 1], rate = initial_theta[2, 2]) +
               dgamma(x = exp(cr_theta[3]), shape = initial_theta[3, 1], rate = initial_theta[3, 2])
  
  return(sum(tmp) + ln_priors + sum(cr_theta))
  
}
