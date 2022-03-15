# Simulating data from the susceptible-infected-susceptible model. The 
# parameters for this simulation have been set by default to the ones used in 
# the original paper.

# set working directory
setwd("~/gitprojects/stats-230/final")

# load functions
source(file = "simulation-hidden-states.R")
source(file = "simulation-observations.R")
source(file = "simulation-pens.R")  

sis_data <- sis_simulations_mp_pens(trials = 99, n.cattle = 8, n.pens = 20)

save(sis_data, file = "data/simulations.Rdata")
