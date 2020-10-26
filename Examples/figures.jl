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
println(homedir)
#%% Figure 1

pars_path = joinpath(homedir, "Settings")
println(conds_path)
#Read JSON files for initial conditions and parameters
u0 = read_JSON(join_path(pars_path,"conds.json")) |> extract_dict;
p_dict = read_JSON(joinpath(pars_path, "params.json")) 