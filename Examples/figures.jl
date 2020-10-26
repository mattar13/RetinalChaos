using RetinalChaos
using Plots, Colors
using LaTeXStrings
using Plots.Measures

#Define plotting attributes
font_title = Plots.font("Arial", 24)
font_axis = Plots.font("Arial", 12)
font_legend = Plots.font("Arial", 8)
pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
homedir = pwd()
pars_path = joinpath(homedir, "Settings")
figs_path = joinpath(homedir, "Examples")
#%% Figure 1


#Read JSON files for initial conditions and parameters
u0 = read_JSON(joinpath(pars_path,"conds.json")) |> extract_dict;
#Read JSON files for parameters and set the applied current to 10pA
p_dict = read_JSON(joinpath(pars_path, "params.json")) 
p_dict[:I_app] = 10.0
p = p_dict |> extract_dict
#Set the time span from 0s -> 60s
tspan = (0.0, 60e3)
#Create the problem
prob = ODEProblem(T_ode, u0, tspan, p)
println("Time it took to simulate 60s:")
@time sol = solve(prob); 