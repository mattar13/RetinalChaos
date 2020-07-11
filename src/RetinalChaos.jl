module RetinalChaos

using Dates
import Dates.now
using JSON2, JLD
using DifferentialEquations, ParameterizedFunctions
using ModelingToolkit
using LinearAlgebra, ForwardDiff, NLsolve
using Distributions
using Images, ImageSegmentation
using ProgressMeter
#using Logging, TerminalLoggers #For some reason these are causing problems randomly
#global_logger(TerminalLogger());
using Plots
using DataFrames, XLSX
using Loess, StatsBase, Statistics

# We make it so, CuArrays is attempted to be loaded onto a computer. If it cannot, then GPU arrays are disabled
using CuArrays
    
check_version() = println("Version 1.0")
######################UTILITIES######################

include("models.jl")
include("utilities.jl")
include("dynamical_analysis.jl")
include("fitting.jl")
include("wave_extraction.jl")
include("plotting.jl")

export read_JSON
export run_model
export append_modeldata

export model_pars, model_conds
###### The main simulation loop is here#########################################
"""
This function contains everything you need to run a single instance of the model,
    and then save the stats and params.
"""

end
