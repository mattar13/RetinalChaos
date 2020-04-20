module RetinalChaos

using JSON2
using DifferentialEquations, ParameterizedFunctions
using LinearAlgebra, ForwardDiff, NLsolve
using Distributions
using Images, ImageSegmentation
using ProgressMeter

check_version() = println("Version 1.0")
######################UTILITIES######################

include("models.jl")
include("utilities.jl")
include("dynamical_analysis.jl")
include("plotting.jl")
include("fitting.jl")
include("wave_extraction.jl")

###### The main simulation loop is here#########################################
"""
This function contains everything you need to run a single instance of the model,
    and then save the stats and params.
"""
function run_model(p_dict, u_dict, tspan; nx = 96, ny = 96, μ = 0.25)
    SACnet = BurstPDE(nx, ny; μ = μ)
    u0_mat = cat(map(x -> fill(u_dict[x], (ny, nx)), BurstModel.syms)..., dims = 3)
    u0 = extract_dict(u_dict, BurstModel.syms)
    p0 = extract_dict(p_dict, BurstModel.params)
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
        #saveat = 100
    );
    #get the last solution from the warmup
    println("[$(now())]: Warming up the model for $(tspan[end]/1000)s")
    u0_new = SDE_mat_sol[end]
    SDE_mat_prob = SDEProblem(SACnet, noise_2D, u0_new, tspan, p0);
    @time SDE_mat_sol = solve(
        SDE_mat_prob,
        SOSRI(),
        abstol = 0.2,
        reltol = 2e-2,
        maxiters = 1e7,
        progress = true,
        saveat = 10.0,
    );
    println("[$(now())]: Model completed")
    timestamp = now()
    df_params = params_to_datasheet(timestamp, p0, u0)
    SDE_sol_arr = Array(SDE_mat_sol);
    println("[$(now())]: Running statistics")
    df_stats = @time run_wavestats(timestamp, SDE_sol_arr[:,:,1,:])
    return SDE_sol_arr, df_params, df_stats
end

end
