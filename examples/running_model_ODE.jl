using Revise
using RetinalChaos
using Plots

import RetinalChaos.ODEModel #import the ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters

#Step 3: determine the timespan
tmin = 0.0
tmax = 300e3

#Step 4: set up the problem
prob = ODEProblem(ODEModel, u0, (tmin, tmax), parameters)

#Step 5: Solve the problem
@time sol = solve(prob);

#Once you have the solution you can do anything like plotting or other math
plot(sol, idxs=[v, I_Na, I_Ca, I_K], layout=(4, 1), lw=2.0, c=:red)