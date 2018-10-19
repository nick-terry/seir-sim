# seir-sim
A discrete-event simulation of the SEIR epidemic model on a contact network

The SEIR model is an important extension of the SIR model for the spread of infectious disease. Where SIR models represent the population of study by three states (Susceptible,Infected,Removed), SEIR models add a fourth state, Exposed, which represents a latent period during which the disease incubates. This extension is an important generalization to capture the dynamics of diseases such as Chicken Pox and Ebola.

seir-sim uses a constant recovery rate algorithm to update the network each time step. The constant recovery rate algorithm is preferable to the constant disease duration variation, as its results are more readily comparable to analytical solutions of the ODE SEIR model (Holme, 2014). 

I used this simulation to evaluate the effectiveness some different algorithms for targeted immunization.

![none](https://github.com/nick-terry/seir-sim/blob/master/example_plot.png "A plot of population levels of each state over time for an SEIR simulation")
