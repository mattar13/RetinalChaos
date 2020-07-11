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
function simulation_loop(net::Network, u0::Array{Float64,3}, p::Array{Float64, 1}; 
        root::String = "D:\\ModellingData\\", sim_name::String = "default", 
        wu_time::Float64 = 60e3, 
        tmax = 300e3, dt::Float64 = 1.0, checkpoint::Float64 = 100e3,    #We will run the solution for 300s with a dt of 1ms, checkpointing every 100s
        )
    
    path = joinpath(root, sim_name)

    ic_path = joinpath(path, "previous.jld")
    #This loop checks whether or not the solution has already been created
    success, ic = parse_ic(path, ic_path)
    if success == false #Unsuccessful loading of the initial condition leads to needing to warmup
        println("[$(now())]: Beginning simulation '$sim_name'")
        prob = SDEProblem(net, noise_2D, u0, (0.0, wu_time), p);
        println("[$(now())]: Warming Up solutiuon for $wu_time ms")
        @time sol = solve(
            prob,
            SOSRI(),
            abstol = 0.2,
            reltol = 2e-2,
            maxiters = 1e7,
            progress = true, 
            save_everystep = false
        );
        println("[$(now())]: Saving warmed-up solution")
        jldopen(ic_path, "w") do file
            file["ic"] = sol[end]
        end
        ic = sol[end]
    end
    
    sim_path = joinpath(path, "dataset.jld")
    #The initial two datasets for the sim_path are the time series and the size of the array, which will come in handy for data analysis
    tsteps = collect(0.0:dt:tmax)
    jldopen(sim_path, "w") do file
        file["time"] = tsteps
        file["size"] = (size(ic,1), size(ic,2), length(tsteps))
    end
    #We will simulate each step in chunks based on checkpoint
    for t = 1.0:checkpoint:tmax
        tfin = (t-1)+checkpoint
        prob = SDEProblem(net, noise_2D, ic, (t, tfin), p);
        println("[$(now())]: Simulating timesteps $t -> $tfin -> [dt:$dt]")
        sol = solve(
            prob,
            SOSRI(),
            abstol = 0.2,
            reltol = 2e-2,
            maxiters = 1e7,
            progress = true, 
            saveat = dt
            )
        println("[$(now())]: Saving data points")
        jldopen(sim_path, "r+") do file
            for t in sol.t
                #println("Saving timepoint $t")
                file["$t"] = sol(t)[:,:,1]
            end
        end
        println("[$(now())]: Saving final data point to be reused")
        jldopen(ic_path, "w") do file
            file["ic"] = sol[end]
        end
        ic = sol[end]
    end
end

end
