module RetinalChaos

#There is no way around it, if I include functions for plotting, I have to include plots
using Plots
export plot, plot!, grid, @animate #Out of the box, I want to be able to plot
#Import some other plotting utilities
using Colors, ColorSchemes, LaTeXStrings, StatsPlots, Dates
export colormatch, colormap, colorschemes

using Plots.Measures

#Import small functions
import Base.length
import Base.print
import Dates.now
export now

#This is for showing the progress of the wave finding function. Which also should be looked at
using ProgressMeter

#Imports if using GPU
using CuArrays

#These imports deal with modelling and running the models
using OrdinaryDiffEq, StochasticDiffEq, ModelingToolkit
export SDEProblem, ODEProblem, solve, SOSRI
#These macros will be useful for extending the model
export @parameters, @variables, @derivatives, @register
export ODESystem, SDESystem

#Imports for reading and writing parameters and solutions
using JSON2, JLD2

#Imported for dynamical analysis
using ForwardDiff, LinearAlgebra, NLsolve

#Imported for saving statistics to excel. Not necessary at this time. Might remove
#using DataFrames, XLSX

#These imports may not be used. They have to do with analysis of waves, but this has become irrelevant. 
#using Images, ImageSegmentation

#These imports are for distributions and statistics. Not necessary for the package, can load based on your needs
#using Distributions, Statistics, StatsBase

#Since almost every notebook makes use of Plots, and we really don't get anything extra by importing it. I will leave this commented
#using Plots
#import Plots.Measures


    
check_version() = println("Version 1.0")
######################UTILITIES######################

include("models.jl")
include("utilities.jl")
include("dynamical_analysis.jl")
include("fitting.jl")
include("wave_extraction.jl")
include("plotting.jl")

# Export functions for dynamical analysis
export EnsembleProblem, EnsembleThreads, ensemble_func
#We are exporting the minimum functions needed to run a 1D simulation
export T_ode, T_sde, SOSRI #Load all the DiffEq Interface
export extract_dict, read_JSON #Load the parameter loading functions 
#Export functions related to creating the 2D network
export Network
export tar_conds, tar_pars, p_find, u_find


###### The main simulation loop is here#########################################
"""
This function contains everything you need to run a single instance of the model,
    and then save the stats and params.
"""
function simulation_loop(net::Network, u0::Array{Float64,3}, p::Array{Float64, 1}; 
        save_sim::Bool = true, root::String = "D:\\ModellingData\\", sim_name::String = "default", 
        wu_time::Float64 = 60e3, 
        tmax = 120e3, dt::Float64 = 1.0, checkpoint::Float64 = 60e3,    #We will run the solution for 300s with a dt of 1ms, checkpointing every 100s
        )
    #Create the new path as the sim_name (or title)
    path = joinpath(root, sim_name)
    #Create both the initial condition file and the data file
    ic_path = joinpath(path, "previous.jld2")
    sim_path = joinpath(path, "dataset.jld2")
    #This loop checks whether or not the solution has already been created
    success, ic = parse_ic(path, ic_path)
    if success == false #Unsuccessful loading of the initial condition leads to needing to warmup
        println("[$(now())]: Beginning simulation '$sim_name'")
        prob = SDEProblem(net, noise, u0, (0.0, wu_time), p);
        println("[$(now())]: Warming Up solutiuon for $wu_time ms")
        @time sol = solve(prob,SOSRI(),abstol = 0.2, reltol = 2e-2,  maxiters = 1e7, progress = true, save_start = false, save_everystep = false);
        println("[$(now())]: Saving warmed-up solution")
        last_time = 0.0
        ic = sol[end]
        if save_sim
            JLD2.@save ic_path ic last_time
        end
    end
    
    is_dataset, data_size = parse_ic(path, sim_path; data_name = "data_size")
    println(is_dataset)
    #We only want to create a new file if one does not already exist. Otherwise we want to just append
    if is_dataset == false
        println("Dataset not yet created")
        #The initial two datasets for the sim_path are the time series and the size of the array, which will come in handy for data analysis
        tsteps = collect(1.0:dt:tmax)
        data_size = (size(ic,1), size(ic,2), length(tsteps))
        #We will simulate each step in chunks based on checkpoint
        println("Saving timsetps and data size")
        if save_sim
            JLD2.@save sim_path tsteps data_size
        end
    end
        
    for t = 1.0:checkpoint:tmax
        #Load the last time and initial condition from the previous solution  
        JLD2.@load ic_path ic last_time
        if t < last_time
            println("previously saved solution is not complete")
            continue
        else
            tfin = (t-1)+checkpoint
            prob = SDEProblem(net, noise_2D, ic, (t, tfin), p);
            println("[$(now())]: Simulating timesteps $t -> $tfin -> [dt:$dt]")
            sol = solve(prob, SOSRI(), abstol = 0.2, reltol = 2e-2, maxiters = 1e7, progress = true, saveat = dt)
            println("[$(now())]: Saving data points")
            if save_sim
                jldopen(sim_path, "a+") do file
                    for t in sol.t
                        #println("Saving timepoint $t")
                        file["$t"] = sol(t)[:,:,1]
                    end
                end
                println("[$(now())]: Saving final data point to be reused")
                jldopen(ic_path, "w") do file
                    file["ic"] = sol[end]
                    file["last_time"] = tfin
                end
            end
        end     
    end
    ic
end

export simulation_loop

end
