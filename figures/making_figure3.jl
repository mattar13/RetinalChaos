#%% Running and analyzing the model using RetinalChaos.jl
using RetinalChaos
import RetinalChaos: Φ, diffuse, ħ
using Dates
using StatsBase, Statistics
using LaTeXStrings
#Setup the fonts and stuff

font_title = font("Arial", 24)
font_axis = font("Arial", 12)
font_legend = font("Arial", 8)
#RetinalChaos.pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
gr(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
#Set up the file root and default parameters
param_root = "params\\"
params_file = joinpath(param_root, "params.json")
conds_file = joinpath(param_root, "conds.json")

#save everything in the figures folder
save_figs = "figures\\"
if isdir(save_figs) == false
    #The directory does not exist, we have to make it 
    mkdir(save_figs)
end