# RetinalChaos
A system of dynamical equations simulating the phenomena of retinal waves

Julia code for a system of differential equations used for simulating cholinergic retinal waves

## To install:
1) Install the most up to date Julia: https://julialang.org/ 

2) Run Julia

3) Press "]" and create and activate a new environment: 

     generate Modelling

     activate Modelling

4) Install the package

     add https://github.com/mattar13/RetinalChaos.git

5) Hit backspace to return to Julia REPL mode

## To use
1) Run a ODE: (This is also located in examples/running_model_ODE.jl)
     
     using RetinalChaos
     using Plots

     #Step 1: Import the model parameters:
     import RetinalChaos.ODEModel # the ODEModel
     import RetinalChaos.u0 #the initial condutions 
     import RetinalChaos.parameters #the parameters

     parameters[I_app] = 10.0 #Set the applied current to 100pA
     #Step 2: determine the timespan
     tmin = 0.0
     tmax = 300e3

     #Step 3: set up the problem
     prob = ODEProblem(ODEModel, u0, (tmin, tmax), parameters)

     #Step 4: Solve the problem
     @time sol = solve(prob);

     #Once you have the solution you can do anything like plotting or other math
     plot(sol, idxs=[v, c, b], layout=(4, 1), lw=2.0, c=:red)

Download source files using 