# This function calculates the logarithm of the joint posterior distribution, 
# function is used as u in the hamiltonian_mc function.

# transitions should be organize in an array with trials, pens and then 
# organized as follows: x_t-1 | x_t 
# 0|0, 0|1, 1|0 and 1|1.
# transitions will have dimensions trials-1, pens, 4

# infected totals should be organized as a matrix with rows equal to the number 
# of trials and columns equal to the number of pens.

# initial_theta is the matrix of initial values foe the prior distributions, 
# each parameter on a single row organized as alpha, beta, m. Dimensions should
# be 3, 2.

log_joint_post <- function(cr_theta, cr_chains, initial_theta){
  
  at <- exp(cr_theta[1])
  bt <- exp(cr_theta[2])
  
  trials <- dim(cr_chains)[1]
  
  n_pens <- dim(cr_chains)[2]
  
  tmp <- 0
  
  if(cr_theta[3]>0){
    for(p in 1:n_pens){
      
      ind <- cr_chains[,,p]
      
      for(t in 1:(trials-1)){
        
        # transitions 0|0
        n_00 <- length(ind[t,] == 0 & ind[t+1,] == 0)
        
        # transitions 1|0
        n_01 <- length(ind[t,] == 0 & ind[t+1,] == 1)
        
        # transitions 0|1
        n_10 <- length(ind[t,] == 1 & ind[t+1,] == 0)
        
        # transitions 1|1
        n_11 <- length(ind[t,] == 1 & ind[t+1,] == 1)
        
        # number of infected at time t
        infected <- length(ind[t,] == 1)
        
        tmp <- tmp + n_00 * (- at - bt * infected) + 
          n_01 * log(1 - exp(- at - bt * infected)) -
          (n_11 + n_10) * log(cr_theta[3] + 1) + 
          n_11 * log(cr_theta[3])
      }
    }
    
    ln_priors <- dgamma(x = at, shape = initial_theta[1, 1], rate = initial_theta[1, 2], log = TRUE) + log(at) +
      dgamma(x = bt, shape = initial_theta[2, 1], rate = initial_theta[2, 2], log = TRUE) + log(bt) + 
      dgamma(x = cr_theta[3], shape = initial_theta[3, 1], rate = initial_theta[3, 2], log = TRUE)
    
    return(tmp + ln_priors)
  }
  else{
    return(log(tmp))
  }
}
