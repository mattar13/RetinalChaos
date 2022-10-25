using Revise
using RetinalChaos
using Plots

import RetinalChaos.u0 #import the 
import RetinalChaos.parameters
reload_parameters()

#Step 5: Solve the problem
probSDE = RetinalChaos.loadSDE(u0, parameters)
@time sol = solve(probSDE, SOSRI()); #So far the best method is SOSRI
plot(sol, idxs=[v, W], layout=(2, 1))