module RetinalChaos

using Dates
using JSON2, JLD
using DifferentialEquations, ParameterizedFunctions
using ModelingToolkit
using LinearAlgebra, ForwardDiff, NLsolve
using Distributions
using Images, ImageSegmentation
using ProgressMeter
using Logging, TerminalLoggers
global_logger(TerminalLogger());
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
function run_simulation(prob; run_time = 300e3, warmup_time = 60e3, dt = 10.0)
    println("Warming up solution")
    prob = SDEProblem(prob.f, prob.g, prob.u0, (0.0, warmup_time), prob.p);
    sol = solve(
            prob,
            SOSRI(),
            abstol = 0.2,
            reltol = 2e-2,
            maxiters = 1e7,
            progress = true, 
            save_everystep = false
        )
    #Take the last solution to the ODE and use it as the initial conditions
    warmed_up_prob = SDEProblem(prob.f, prob.g, sol[end], (0.0, run_time), prob.p);
    sol = solve(
            warmed_up_prob,
            SOSRI(),
            abstol = 0.2,
            reltol = 2e-2,
            maxiters = 1e7,
            progress = true, 
            saveat = dt,
        )
    return sol
end
end
