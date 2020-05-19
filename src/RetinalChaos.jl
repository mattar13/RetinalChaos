module RetinalChaos

using Dates
using JSON2
using DifferentialEquations, ParameterizedFunctions
using ModelingToolkit
using LinearAlgebra, ForwardDiff, NLsolve
using Distributions
using Images, ImageSegmentation
using ProgressMeter
using Plots
using DataFrames, XLSX
using Loess, StatsBase

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
function run_model(p_dict, u_dict, tspan; dt = 10.0, nx = 96, ny = 96, μ = 0.25, gpu = true)
    SACnet = BurstPDE(nx, ny; μ = μ, gpu = gpu)
    u0_mat = cat(map(x -> fill(u_dict[x], (ny, nx)), model_conds)..., dims = 3)
    if gpu
        u0_mat = u0_mat |> cu
        CuArrays.allowscalar(false)
    end
    u0 = extract_dict(u_dict, model_conds)
    p0 = extract_dict(p_dict, model_pars)
    #warm up the model
    println("[$(now())]: Warming up the model for 60s")
    SDE_mat_prob = SDEProblem(SACnet, noise_2D, u0_mat, (0.0, 60e3), p0);
    @time SDE_mat_sol = solve(
        SDE_mat_prob,
        SOSRI(),
        abstol = 0.2,
        reltol = 2e-2,
        maxiters = 1e7,
        progress = true,
        save_everystep = false,
    );
    #get the last solution from the warmup
    println("[$(now())]: Running the model for $(tspan[end]/1000)s")
    u0_new = SDE_mat_sol[end]
    SDE_mat_prob = SDEProblem(SACnet, noise_2D, u0_new, tspan, p0);
    @time SDE_mat_sol = solve(
        SDE_mat_prob,
        SOSRI(),
        abstol = 0.2,
        reltol = 2e-2,
        maxiters = 1e7,
        progress = true,
        saveat = dt,
    );
    println("[$(now())]: Model completed")
    timestamp = now()
    df_params = params_to_datasheet(timestamp, p0, u0)
    #Converting Solution to array
    SDE_sol_arr = zeros(size(SDE_mat_sol)...)
    for i in size(SDE_sol_arr, 4)
        SDE_sol_arr[:,:,:,i] .= Array(SDE_mat_sol(i))
    end
    println("[$(now())]: Running statistics")
    df_stats = @time run_wavestats(timestamp, SDE_sol_arr[:,:,1,:])
    return SDE_sol_arr, df_params, df_stats
end

end
