module RetinalChaos



#Import small functions
import Base.length
import Base.print
import Dates.now
export now

#This is for showing the progress of the wave finding function. Which also should be looked at
@time using ProgressMeter

#Imports if using GPU
@time using CuArrays

#These imports deal with modelling and running the models
@time using OrdinaryDiffEq, StochasticDiffEq, ModelingToolkit
export SDEProblem, ODEProblem, solve, SOSRI
#These macros will be useful for extending the model
export @parameters, @variables, @derivatives, @register
export ODESystem, SDESystem

#Imports for reading and writing parameters and solutions
@time using JSON2, JLD2

#Imported for dynamical analysis
@time using ForwardDiff, LinearAlgebra, NLsolve

#For Binomial Nullification of parameters
@time using Distributions

#Imported for saving statistics to excel. Not necessary at this time. Might remove
#using DataFrames, XLSX

#These imports may not be used. They have to do with analysis of waves, but this has become irrelevant. 
#using Images, ImageSegmentation

#These imports are for distributions and statistics. Not necessary for the package, can load based on your needs
#using Statistics, StatsBase
    
check_version() = println("Version 1.0")
######################UTILITIES######################

include("models.jl")
include("utilities.jl")
include("dynamical_analysis.jl")
include("fitting.jl")
include("wave_extraction.jl")
export_plotting() = include("plotting.jl")

# Export functions for dynamical analysis
export EnsembleProblem, EnsembleThreads, ensemble_func
#We are exporting the minimum functions needed to run a 1D simulation
export T_ode, T_sde, SOSRI #Load all the DiffEq Interface
export extract_dict, read_JSON #Load the parameter loading functions 
#Export functions related to creating the 2D network
export Network
export tar_conds, tar_pars, p_find, u_find
println("Fininshed Importing")

end
