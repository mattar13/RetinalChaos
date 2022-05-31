using Revise
using RetinalChaos

#Step 1: Import the initial conditions
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict

#Step 2: Import the parameters
pars_dict = read_JSON("params\\params.json")
p = pars_dict |> extract_dict

#Step 3: determine the timespan
tspan = (0.0, 300e3)

#Step 4: set up the problem
prob = SDEProblem(T_sde, noise, u0, tspan, p)

#Step 5: Solve the problem
@time sol = solve(prob, SOSRI());