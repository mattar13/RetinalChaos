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
import ModelingToolkit as MTK
export MTK
using Symbolics
export reload_parameters
export ODEProblem, SDEProblem, solve
export SOSRI, SOSRA, SROCK1, SRIW1, SKenCarp #Export any algorithims
export SRIW1
export SROCK1 #Needs a dt specification
export EnsembleProblem, EnsembleThreads

#Don't explicitly export anything
include("open_parameters.jl") #Load all of the parameters
include("auxillary_functions.jl") #Load all of the necessary functions
include("equations.jl") #Load all model equations

#export all of the parameters so we can edit something
export v, I_Ca, I_Na, I_K, n, m, h, c, a, b, e, i, W #Initial conditions
export g_leak, E_leak, g_Ca, V1, V2, E_Ca, g_K, E_K, g_TREK, g_ACh, k_ACh, E_ACh, g_GABA, k_GABA, E_Cl, I_app, C_m
export V3, V4, τn
export C_0, λ, δ, τc
export α, τa, β, τb, ρe, ρi, τACh, τGABA, VSe, VSi, V0e, V0i
export De, Di
export τw, σ
export g_Na, E_Na
export V7, V8, V9
export V10, V11, V12
export V13, V14, V15
export V16, V17, V18
#using CUDA
#export cu, allowscalar

#Import all the auxillary functions
#===========================================Loading the Parameters=============================================#
using JSON2, JLD2, BSON #Imports for reading and writing parameters and solutions
include("utilities.jl")
export read_JSON, extract_dict, indexof
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

using BifurcationKit
import BifurcationKit as BK
export BK #Use this explicitly for Bifurcation kit utilities
using ForwardDiff, NLsolve, Setfield, LinearAlgebra #Imported for dynamical analysis
include("dynamical_analysis.jl") # Export functions for dynamical analysis
export ensemble_func
export find_equilibria, find_bifurcation
export codim_map
export @lens, norminf

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