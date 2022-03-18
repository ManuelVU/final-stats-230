# This function calculates the gradient of the logarithm of the joint posterior 
# distribution, function is used as grad_u in the hamiltonian_mc function.

# transitions should be organize in an array with trials, pens and then 
# organized as follows: x_t-1 | x_t 
# 0|0, 0|1, 1|0 and 1|1.
# transitions will have dimensions trials-1, pens, 4

# infected totals should be organized as a matrix with rows equal to the number 
# of trials and columns equal to the number of pens.

# initial_theta is the matrix of initial values foe the prior distributions, 
# each parameter on a single row organized as alpha, beta, m. Dimensions should
# be 3, 2.

gradient_log_joint_post <- function(cr_theta, cr_chains, initial_theta){
  
  trials <- dim(cr_chains)[1]
  
  n_pens <- dim(cr_chains)[2]
  
  at <- exp(cr_theta[1])
  
  bt <- exp(cr_theta[2])
  
  grad_alpha <- 0
  
  grad_beta <- 0
  
  grad_m <- 0
  
  if(cr_theta[3] > 0){
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
        
        eratio <- exp(-at - bt * infected) / (1 - exp(-at - bt * infected))
        
        grad_alpha <- grad_alpha + n_01 * at * eratio - n_00 * at
        
        grad_beta <- grad_beta + n_01 * bt * infected * eratio - n_00 * bt * infected
        
        grad_m <- n_11 / cr_theta[3] - (n_11 + n_10) / (cr_theta[3] + 1) 

      }
    }
    
    d_ln_a <- grad_alpha + initial_theta[1, 1] - initial_theta[1, 2] * at
    d_ln_b <- grad_beta + initial_theta[2, 1] - initial_theta[2, 2] * bt
    d_ln_m <- grad_m  + (initial_theta[3, 1] - 1) / cr_theta[3] - initial_theta[3, 2]
    
    return(c(d_ln_a, d_ln_b, d_ln_m))
  }
  else{
    return(c(log(grad_alpha), log(grad_beta), log(grad_m)))
  }
}
