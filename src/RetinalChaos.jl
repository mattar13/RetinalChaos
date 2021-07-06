module RetinalChaos

#Import small functions
import Base.length
import Base.print
import Dates.now
export now
#This is for showing the progress of the wave finding function. Which also should be looked at

using Telegram, Telegram.API, ConfigEnv
using ProgressMeter
println("Small functions imported")

#Imports if using GPU
using CUDA
export cu, allowscalar
CUDA.allowscalar(false)

#Include all the plotting utilities
using Plots, Measures

using DifferentialEquations
#using ModelingToolkit, DiffEqBase
import DiffEqBase.AbstractODEProblem
export SDEProblem, ODEProblem, solve, SOSRI
export EnsembleProblem, EnsembleThreads
#These macros are only useful with ModelingToolkit which currently we are not using
#export @parameters, @variables, @derivatives, @register
#export ODESystem, SDESystem

println("Modelling utilities imported")
#Imports for reading and writing parameters and solutions
using JSON2, JLD2

println("Extra utilities imported")
#Imported for dynamical analysis
using ForwardDiff, LinearAlgebra, NLsolve

#For Binomial Nullification of parameters
using Distributions

#These imports are for distributions and statistics. Not necessary for the package, can load based on your needs
#using Statistics, StatsBase
    
check_version() = println("Version 1.0")
######################UTILITIES######################

println("Importing Models")
include("models.jl")
export T_ode, T_sde, SOSRI, SOSRA #Load all the DiffEq Interface
#Export functions related to creating the 2D network
export Network, noise
export tar_conds, tar_pars, p_find, u_find

include("utilities.jl")
export extract_dict, read_JSON, write_JSON #Load the parameter loading functions 

# Export functions for dynamical analysis
include("dynamical_analysis.jl")
export ensemble_func
export find_equilibria, find_bifurcation
export codim_map

#include("fitting.jl") 

include("wave_extraction.jl")
export calculate_threshold
export get_timestamps, max_interval_algorithim, timescale_analysis
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

end