module RetinalChaos

using Requires
#==Base imports==#
import Base.length
import Base.print
import Dates.now
export now


#===========================================Import logging materials===========================================#
#using Logging: global_logger
#using TerminalLoggers: TerminalLogger
#global_logger(TerminalLogger())
#using ProgressMeter
#export @showprogress

#======================================Exports======================================#
using DifferentialEquations, ModelingToolkit
import ModelingToolkit as MTK
export MTK
using MethodOfLines, DomainSets, DiffEqOperators #These are for loading the PDEs
export Interval
using Symbolics
export reload_parameters
export ODEProblem, SDEProblem, solve
export SOSRI, SOSRA, SROCK1, SRIW1, SKenCarp #Export any algorithims
export SRIW1, EM
export SROCK1 #Needs a dt specification
export EnsembleProblem, EnsembleThreads
export PresetTimeCallback
#Don't explicitly export anything
include("open_parameters.jl") #Load all of the parameters
include("auxillary_functions.jl") #Load all of the necessary functions
include("equations.jl") #Load all model equations
include("load_models.jl")
#export all of the parameters so we can edit something
export loadODE, loadSDE, loadPDE, loadSPDE
export remake
export x,y,t
export v, n, m, h, c, a, b, e, i, W #Initial conditions for ODEs and SDEs
export v̂, n̂, m̂, ĥ, ĉ, â, b̂, ê, î, Ŵ #Initial conditions for PDEs

export I_Ca, I_Na, I_K
export I_ext
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
#using Plots
#using PyPlot
#using PyCall
#@pyimport matplotlib.animation as anim

using JSON2, JLD2, BSON #Imports for reading and writing parameters and solutions
include("utilities.jl")
export read_JSON, extract_dict, indexof
export save, load
export @save

using Distributions, StatsBase
export Binomial, Histogram, fit


#==========================================Extracting wave and events==========================================#
#using ePhys #Export the wave extraction utilities



#using RecipesBase #I want to avoid using thisPlots


#export animate_solution


#Load all of the old modelling aspects. We can use that one for PDE
#import Plots.@animate #Only import 
include("models.jl")
export T_PDE_w_NA, noise
using Statistics, StatsBase
export std
export run_model

function __init__()
     @require ePhys = "69cbc4a0-077e-48a7-9b45-fa8b7014b5ca" begin
          import ePhys: calculate_threshold, get_timestamps, timeseries_analysis, max_interval_algorithim
          using Images, ImageSegmentation 
          using DataFrames, Query, XLSX
          include("wave_extraction.jl") #Export functions for wave extraction
          export calculate_threshold, timeseries_analysis
          export get_timestamps, max_interval_algorithim, extract_interval, timeseries_analysis
          export WaveSegmentation
     end

     @require BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665" begin
          import BifurcationKit as BK
          export BK #Use this explicitly for Bifurcation kit utilities
          using ForwardDiff, NLsolve, Setfield, LinearAlgebra #Imported for dynamical analysis
          include("dynamical_analysis.jl") # Export functions for dynamical analysis
          export ensemble_func
          export find_equilibria, find_bifurcation
          export codim_map
          export @lens, norminf
          export t_pars, t_conds
     end

     @requires ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210" begin
          import ForwardDiff as FD
          export FD
          include("fitting.jl")
          export extract_spike_trace, extract_burst_trace, extract_IBI_trace
          export SpikeLoss, BurstLoss, IBILoss
          export IntervalLoss, MSELoss
          export MeanSquaredErrorSOL, MeanSquaredError
          export TimescaleLoss
     end
end

end