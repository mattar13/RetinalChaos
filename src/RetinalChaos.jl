module RetinalChaos

#==Base imports==#
import Base.length
import Base.print
import Dates.now
export now

#===========================================Import logging materials===========================================#
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using ProgressMeter
export @showprogress

#======================================Exports======================================#
using DifferentialEquations, ModelingToolkit

export ODEProblem, SDEProblem, solve
export SOSRI, SOSRA, SROCK1, SRIW1, SKenCarp #Export any algorithims
export SRIW1
export SROCK1 #Needs a dt specification
export EnsembleProblem, EnsembleThreads

#create a function that completes setup

include("open_parameters.jl") #Load all of the parameters
include("auxillary_functions.jl") #Load all of the necessary functions
include("ODE_equations.jl") #Load all model equations
#using CUDA
#export cu, allowscalar

#Import all the auxillary functions
#===========================================Loading the Parameters=============================================#
#using JSON2, JLD2, BSON #Imports for reading and writing parameters and solutions
#include("utilities.jl")
#export read_JSON, extract_dict, p_find, u_find
#export save, load
#export @save

using Distributions, StatsBase
export Binomial, Histogram, fit


#==========================================Extracting wave and events==========================================#
#using ePhys #Export the wave extraction utilities

include("wave_extraction.jl") #Export functions for wave extraction
export calculate_threshold, timeseries_analysis
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


using ForwardDiff
import ForwardDiff as FD
export FD
include("fitting.jl")
export extract_spike_trace, extract_burst_trace, extract_IBI_trace
#export SpikeLoss, BurstLoss, IBILoss
export IntervalLoss
export MeanSquaredErrorSOL, MeanSquaredError
export TimescaleLoss

#=
using Telegram, Telegram.API, ConfigEnv

#These imports are for distributions and statistics. Not necessary for the package, can load based on your needs
using Statistics, StatsBase
using Images, ImageSegmentation
######################UTILITIES######################
=#
end