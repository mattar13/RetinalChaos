module RetinalChaos

#Import small functions
import Base.length
import Base.print
import Dates.now
export now
#This is for showing the progress of the wave finding function. Which also should be looked at
using ProgressMeter
println("Small functions imported")

#Imports if using GPU
using CuArrays

#using OrdinaryDiffEq, StochasticDiffEq, ModelingToolkit
using DifferentialEquations, ModelingToolkit
export SDEProblem, ODEProblem, solve, SOSRI
export EnsembleProblem, EnsembleThreads
#These macros will be useful for extending the model
export @parameters, @variables, @derivatives, @register
export ODESystem, SDESystem

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

#println("Importing Models")
include("models.jl")
export T_ode, T_sde, SOSRI #Load all the DiffEq Interface
#Export functions related to creating the 2D network
export Network, noise
export tar_conds, tar_pars, p_find, u_find

include("utilities.jl")
export extract_dict, read_JSON #Load the parameter loading functions 

# Export functions for dynamical analysis
include("dynamical_analysis.jl")
export ensemble_func

#include("fitting.jl") 

include("wave_extraction.jl")
export calculate_threshold
export get_timestamps, max_interval_algorithim, timescale_analysis
#println("Fininshed Importing")

#Include all the plotting utilities
using Plots
export Plots
include("plotting.jl")
export pyplot, font, Measures
export plot, plot!, grid, @animate #Out of the box, I want to be able to plot
#Import some other plotting utilities
using Colors, ColorSchemes, LaTeXStrings, StatsPlots, Dates
export colormatch, colormap, colorschemes

end