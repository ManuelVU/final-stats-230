# Simulation of hidden states from the susceptible-infected-susceptible model 
# using multiple chains (n.cattle). The possible states are susceptible 
# (x_t = 0) and infected (x_t = 1).

# samples: number of trials.
# al: external infection rate.
# bt: within-pen infection rate.
# m_inf: mean infectious period.
# p_initial: probability of infection on first trial
# nc: number of cattle on a given pen.

sis_sim_hidenstates <- function(samples, nc, al, bt, m_inf,
                                p_initial){
  infect_status <- matrix(NA,nrow = samples, ncol = nc)
  infect_status[1,] <- rbinom(n = nc, size = 1, prob = p_initial)
  for(t in 2:samples){
    for(k in 1:nc){
      infect_status[t,k] <- rbinom(n = 1, size = 1,
                                   prob = ifelse(test = infect_status[t-1,k] == 0,
                                                 yes = 1 - exp(-al - bt * sum(infect_status[(t-1),-k])),
                                                 no = (m_inf - 1) / m_inf))
    }
  }
  return(infect_status)
}

