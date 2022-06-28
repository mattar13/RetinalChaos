using Revise
using RetinalChaos
import RetinalChaos: calculate_threshold, get
import RetinalChaos: extract_equilibria, find_equilibria
include("../figures/figure_setup.jl");

#%% Try to fix the equilibria stuff
conds_dict = read_JSON("params\\conds.json")
u0 = conds_dict |> extract_dict

pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = 10.0 #Set initial applied current to 0
pars_dict[:ρi] = 0.0 #remove GABA influence
pars_dict[:ρe] = 0.0 #remove ACh influence
pars_dict[:g_TREK] = 0.0 #Remove the sAHP
p = pars_dict |> extract_dict
tspan = (0.0, 100.0)
prob_eq = ODEProblem(T_ODE, u0, tspan, p)
# Conduct the codim analysis
codim1 = (:I_app)
c1_lims = (45.0, 50.0)
@time c1_map = codim_map(prob_eq, codim1, c1_lims, equilibrium_resolution=10)

#%% Plot Codim Solutions
res = extract_equilibria(c1_map) #Pass back all of the equilibria
points = res[1]
saddle_p = res[2]
stable_p = res[3]
unstable_p = res[4]
unstable_focus_p = res[5]
stable_focus_p = res[6]
#bif_val, bif_eq = find_bifurcation(c1_map)
#saddle_vs = map(x -> x.saddle[1][1], bif_eq)
# Plot 
plot(points, saddle_p, c=:blue)
plot!(points, stable_p, c=:green)
plot!(points, unstable_p, c=:red)
plot!(points, stable_focus_p, c=:red, linestyle=:dash)
plot!(points, unstable_focus_p, c=:green, linestyle=:dash)


#%% Lets look deeper into some issues
eq_val = findall((isnan.(stable_p)) .== 0)[2]
u0_eq = c1_map.equilibria[eq_val].stable[1] #We will use this as out Initial condition

pars_dict = read_JSON("params\\params.json")
pars_dict[:I_app] = c1_map.points[eq_val][1] #Set initial applied current to 0
pars_dict[:ρi] = 0.0 #remove GABA influence
pars_dict[:ρe] = 0.0 #remove ACh influence
pars_dict[:g_TREK] = 0.0 #Remove the sAHP
p = pars_dict |> extract_dict
tspan = (0.0, 1000.0)
prob_eq = ODEProblem(T_ODE, u0_eq, tspan, p)
sol_eq = solve(prob_eq)
eq_obj = find_equilibria(prob_eq)
print(eq_obj)
plot(sol_eq, vars = 1)