# Simulating data from the susceptible-infected-susceptible model. The 
# parameters for this simulation have been set by default to the ones used in 
# the original paper.

# load functions
source(file = "src/simulation-hidden-states.R")
source(file = "src/simulation-observations.R")
source(file = "src/simulation-pens.R")  

sis_data <- sis_simulations_mp_pens(trials = 99, n.cattle = 8, n.pens = 20, prob_test = 1)

save(sis_data, file = "data/simulations.Rdata")
