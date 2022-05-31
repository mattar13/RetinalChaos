module RetinalChaos

const version = :master
param_path = joinpath(splitpath(pathof(RetinalChaos))[1:end-2]..., "params")
version_info() = println(version)

#======================================Import all the pre-loaded packages======================================#
using DifferentialEquations #Differential Equations packages
#import DiffEqBase.AbstractODEProblem
#export SDEProblem, ODEProblem, solve, SOSRI
#export EnsembleProblem, EnsembleThreads


#===========================================Loading the Parameters=============================================#
using JSON2, JLD2, BSON #Imports for reading and writing parameters and solutions
include("utilities.jl")
export read_JSON

#=
using Plots: text_box_width
const verbose = false #Adjust the to print out statments relevant to the module import
#Import small functions
import Base.length
import Base.print
import Dates.now
export now
#This is for showing the progress of the wave finding function. Which also should be looked at

using Telegram, Telegram.API, ConfigEnv
using ProgressMeter
if verbose 
     println("[$(now())]: Small functions imported")
end

#Imports if using GPU
using CUDA
export cu, allowscalar

#Include all the plotting utilities

if verbose 
     print("[$(now())]: Plotting utilities imported in: ")
     @time using Plots, Measures     
else
     using Plots, Measures 
end

using ResettableStacks, RandomNumbers, LinearAlgebra #These packages are needed to load the problems

#Imported for dynamical analysis
using ForwardDiff, NLsolve

if verbose
     println("[$(now())]: Extra utilities imported")
end
#For Binomial Nullification of parameters
using Distributions

#These imports are for distributions and statistics. Not necessary for the package, can load based on your needs
using Statistics, StatsBase
using Images, ImageSegmentation
check_version() = println("Version 1.0")
######################UTILITIES######################

if verbose 
     println("[$(now())]: Importing Models")
end
include("models.jl")
export T_ode, T_sde, SOSRI, SOSRA #Load all the DiffEq Interface
#Export functions related to creating the 2D network
export Network, noise
export tar_conds, tar_pars, p_find, u_find

include("utilities.jl")
export data_bytesize #we can use this to predict the approximate bytesize saved by BSON
export extract_dict, read_JSON, write_JSON #Load the parameter loading functions 
export run_model, save_solution, load_solution, convert_to_cpu, animate_solution

# Export functions for dynamical analysis
include("dynamical_analysis.jl")
export ensemble_func
export find_equilibria, find_bifurcation
export codim_map

#include("fitting.jl") 

include("wave_extraction.jl")
export calculate_threshold
export get_timestamps, max_interval_algorithim, extract_interval, timeseries_analysis
export extract_waves
#println("Fininshed Importing")




include("plotting.jl")
export Plots
#Exporting backends
export gr, font
#exporting main plot types
export plot, plot!, heatmap!, heatmap, scatter!, contour, contour!
#export minor plotting types
export hline!, vline!
#Export annotations and additions
export title!, annotate!, grid, stroke, mm
#Export some animations
export @animate, gif, mov
#Export save figure
export savefig
#Import some other plotting utilities
using Colors, ColorSchemes, LaTeXStrings, StatsPlots
export colormatch, colormap, colorschemes

include("logging.jl")
export dotenv, env_location, BotNotify, BotFigure
=#
end