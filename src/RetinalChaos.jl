module RetinalChaos

#===========================================Import logging materials===========================================#
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using ProgressMeter
export @showprogress

const version = :master
param_path = joinpath(splitpath(pathof(RetinalChaos))[1:end-2]..., "params")
version_info() = println(version)
import Base.length
import Base.print
import Dates.now
export now
#======================================Import all the pre-loaded packages======================================#

using DifferentialEquations #Differential Equations packages
export ODEProblem, SDEProblem, solve
export SOSRI, SOSRA, SROCK1, SRIW1, SKenCarp #Export any algorithims
export SRIW1
export SROCK1 #Needs a dt specification
#import DiffEqBase.AbstractODEProblem
#export SDEProblem, ODEProblem, solve, SOSRI
export EnsembleProblem, EnsembleThreads

using CUDA
export cu, allowscalar

include("models.jl")
export T_ODE, T_SDE, T_PDE
#export GABA_ODE, GABA_SDE, GABA_PDE
export GABA_PDE_gNULL
export t_pars, t_conds
export GABA_pars, GABA_conds
export noise

#===========================================Loading the Parameters=============================================#
using JSON2, JLD2, BSON #Imports for reading and writing parameters and solutions
include("utilities.jl")
export read_JSON, extract_dict, p_find, u_find

using Distributions, StatsBase
export Binomial, Histogram, fit


#==========================================Extracting wave and events==========================================#
include("wave_extraction.jl") #Export functions for wave extraction
export calculate_threshold
export get_timestamps, max_interval_algorithim, extract_interval, timeseries_analysis
export extract_waves


using ForwardDiff, NLsolve, LinearAlgebra #Imported for dynamical analysis
include("dynamical_analysis.jl") # Export functions for dynamical analysis
export ensemble_func
export find_equilibria, find_bifurcation
export codim_map

using RecipesBase #I want to avoid using thisPlots
include("plotting.jl")
#export animate_solution

#=
using Telegram, Telegram.API, ConfigEnv

if verbose 
     print("[$(now())]: Plotting utilities imported in: ")
     @time using Plots, Measures     
else
     using Plots, Measures 
end

using ResettableStacks, RandomNumbers, LinearAlgebra #These packages are needed to load the problems

#These imports are for distributions and statistics. Not necessary for the package, can load based on your needs
using Statistics, StatsBase
using Images, ImageSegmentation
check_version() = println("Version 1.0")
######################UTILITIES######################

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