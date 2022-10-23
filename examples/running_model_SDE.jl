using Revise
using RetinalChaos
using Plots

import RetinalChaos.SDEModel #import the ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters
reload_parameters()

#Step 3: determine the timespan
tmin = 0.0
tmax = 60e3

#Step 4: set up the problem
probSDE = SDEProblem(SDEModel, u0, (tmin, tmax), parameters)

#Step 5: Solve the problem
@time sol = solve(probSDE, SOSRI()); #So far the best method is SOSRI
plot(sol, idxs=[v, W], layout=(2, 1))