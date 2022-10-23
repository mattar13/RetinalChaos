using Revise
using RetinalChaos
using Plots

#Step 1: Import the:
import RetinalChaos.ODEModel # the ODEModel
import RetinalChaos.u0 #the initial condutions 
import RetinalChaos.parameters #the parameters

#Step 2: determine the timespan
tmin = 0.0
tmax = 300e3

#Step 3: set up the problem
prob = ODEProblem(ODEModel, u0, (tmin, tmax), parameters)

#Step 4: Solve the problem
@time sol = solve(prob);

#Once you have the solution you can do anything like plotting or other math
plot(sol, idxs=[v, I_Na, I_Ca, I_K], layout=(4, 1), lw=2.0, c=:red)