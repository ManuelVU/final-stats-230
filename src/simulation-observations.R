# Simulation of outcomes given states for the susceptible-infected-susceptible 
# model, the number of chains is given by the number of columns in h_states.
# Measures of two variables are taken from all chains at the same time. The 
# probability that a measure will be taken is set as p_test.

# h_states: matrix containing simulated values of the hidden states, with trials
#           indicated by rows and chains by the columns.
# th_r: sensitivity of the RAMS test (p(r^c_t = 1 | x^c_t = 1)).
# th_f: sensitivity of the fecal test (p(f^c_t = 1 | x^c_t = 1)).
# p_test: probability that at any given trial a measure of both tests will be 
#         taken.

sis_sim_observations <- function(h_states, th_r, th_f, p_test){
  observations <- array(NA, dim = c(dim(h_states)[1], dim(h_states)[2], 2))
  sample_times <- which(rbinom(n = dim(observations)[1], size = 1, prob = p_test) == 1)
  for(tt in sample_times){
    for(cattle in 1:dim(observations)[2]){
      observations[tt, cattle, ] <- c(rbinom(n = 1, size = 1,
                                             prob = th_r * h_states[tt,cattle]),
                                      rbinom(n = 1, size = 1,
                                             prob = th_f * h_states[tt,cattle]))
    }
  }
  return(observations)
}
