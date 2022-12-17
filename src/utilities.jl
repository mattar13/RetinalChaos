function data_bytesize(arr::Vector{T}; bson_header_bytes = 133) where T <: Real
    data_bytes = sizeof(arr)
    size_bytes = 11 #There is a section in the header that 
    file_size = bson_header_bytes + data_bytes + size_bytes
    return file_size
end

function data_bytesize(arr::Matrix{T}; bson_header_bytes = 133) where T <: Real
    data_bytes = sizeof(arr)
    size_bytes = 11 * ndims(arr) #There is a section in the header that 
    file_size = bson_header_bytes + data_bytes + size_bytes
    return file_size
end

function data_bytesize(arr::Vector{Vector{T}}; bson_header_bytes = 133) where T <: Real
    data_bytes = sizeof(arr[1])
    size_bytes = 111 * length(arr)
    return bson_header_bytes + data_bytes + 111*length(arr)
end


#Reading and writing JSON files
"""
This file writes either a named tuple or a dictionary into a JSON file
"""
function write_JSON(data::T, filename::String) where T
    string_data = JSON2.write(data)
    open(filename, "w") do f
        write(f, string_data)
    end
end

"""
This function reads JSON files into whatever type is_type is.
"""
function read_JSON(::Type{T}, name_file::String) where T <: Dict
    nt = nothing
    open(name_file, "r") do f
        nt = JSON2.read(f, T)
    end
    nt
end

read_JSON(name_file::String) = read_JSON(Dict{Symbol, Float64}, name_file)

"""
Extract a parameter, condition dictionary
USAGE:
#p_dict -> [The parameter dictionary read directly from JSON]
> p0 = extract_dict(p_dict, BurstModel.params)

If your initial conditions are for a network, you can add the extra argument for dims

> In[1]: dims = (64,64);
> In[2]: p0 = extract_dict(p_dict, BurstModel.params, dims)
> Out[2]: 125x125xN AbstractArray:
"""
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}) = map(x -> Float64(dict_item[x]), pars)
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}, dims::Tuple) = cat(map(x -> fill(Float64(dict_item[x]), dims), pars)..., dims = length(dims)+1)

#get the index of the conditions in the list
indexof(num::Symbolics.Num; syms = ODEModel.ps) = findfirst(isequal(num), syms)
indexof(sym::Symbol; syms = ODEModel.ps) = findfirst(isequal(ODEModel.var_to_name[sym]), syms)
"""
This function can do a few things

    1) It can take a dictionary and extract it in the order of the dict_names
    2) It can read through a dictionary and search for certain keys
"""
function extract_dict(dict_item::Dict{Symbol,T}, key_names; dims = (1)) where {T<:Real}
    #This is a parameter and we should extract it
    items = T[]
    for key in key_names
        push!(items, dict_item[key])
    end
    if dims == 1
        return items
    else
        val_map = zeros(dims..., length(items))
        for i in 1:length(items)
            val_map[:, :, i] .= items[i]
        end
        return val_map
    end
end

function extract_dict(dict_item::Dict{Symbol, T}; key_library = [t_pars, t_conds]) where T <: Real
    #To extract the dictionary first we want to take a look at the first param to see if it is a parameter or a condition
    for key_names in key_library
        #println("Here")
        if dict_item.keys[1] ∈ key_names
            items = extract_dict(dict_item, key_names)
            return items
        end
    end
    println("if you made it here there was an error")
    throw("NoKeyLibraries")
end

function extract_dict(dict_item::Dict{Symbol, T}, dims...; dict_names = [t_pars, t_conds]) where T <: Real
    #Extract the dict_item
    vals = extract_dict(dict_item; dict_names = dict_names)
    dims = Int64.(dims) #Put this in to convert the dimensions  
    val_map = zeros(dims..., length(vals))
    for i in 1:length(vals)
        val_map[:,:,i] .= vals[i]
    end
    val_map
end

#This function converts the solution into a CPU based solution
function convert_to_cpu(model::RODESolution)
    prob_remade = remake(model.prob, u0 = model.prob.u0 |> Array)
    #println(prob_remade |> typeof)
    #lets go into remaking the problem
    alg = model.alg
    t = model.t
    u = map(u -> u |> Array, model.u) #This pa
    W = model.W |> Array
    return RetinalChaos.DiffEqBase.build_solution(prob_remade, alg, t, u, W=W)
end

"""
This function saves the solution. 
In order to save the solution correctly, ensure to run: 
    julia> sol = sol |> convert_to_cpu
    julia> save_solution(sol, save_path)

"""
function save_solution(sol_t::Vector{T}, sol_u::Array{Vector{T}, 1}, save_path::String; 
    name = "sol", partitions = 1, mode = :bson
    ) where T <: Real
    
    if mode == :bson
        if partitions == -1 #This means data is getting appended
            #check to see if a file_contents file exists
            file_contents = read_JSON(Dict{Symbol, Vector{String}}, "$(save_path)\\file_contents.json")
            push!(file_contents[:t], "$(name)_t.bson")
            push!(file_contents[:u], "$(name)_u.bson")
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            bson("$(save_path)\\$(name)_t.bson", Dict(:t => sol_t))
            bson("$(save_path)\\$(name)_u.bson", Dict(:u => sol_u))
        elseif partitions == 1
            file_contents = Dict(:t => ["$(name)_t.bson"], :u => ["$(name)_u.bson"])
            #write the details to a 
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            bson("$(save_path)\\$(name)_t.bson", Dict(:t => sol_t))
            bson("$(save_path)\\$(name)_u.bson", Dict(:u => sol_u))
        else
            file_contents = Dict(:t => ["$(name)1_t.bson"], :u => ["$(name)1_u.bson"])
            for i in 2:partitions
                push!(file_contents[:t], "$(name)$(i)_t.bson")
                push!(file_contents[:u], "$(name)$(i)_u.bson")
            end
            write_JSON(file_contents, "$(save_path)\\file_contents.json")
            partition_idxs = round.(Int64, LinRange(1, length(sol_t), partitions+1))
            for idx in 1:length(partition_idxs)-1
                if idx == 1
                    bson("$(save_path)\\$(name)$(idx)_t.bson", Dict(:t => sol_t[1:partition_idxs[idx+1]]))
                    bson("$(save_path)\\$(name)$(idx)_u.bson", Dict(:u => sol_u[1:partition_idxs[idx+1]]))
                else
                    bson("$(save_path)\\$(name)$(idx)_t.bson", Dict(:t => sol_t[partition_idxs[idx]+1:partition_idxs[idx+1]]))
                    bson("$(save_path)\\$(name)$(idx)_u.bson", Dict(:u => sol_u[partition_idxs[idx]+1:partition_idxs[idx+1]]))
                end
            end
        end
    else
        println("TODO implement JLD2")
    end
end

save_solution(sol, save_path::String; kwargs...) = save_solution(sol.t, sol.u; karwgs...)
save_solution(res::DiffEqCallbacks.SavedValues{Float64, Vector{Float64}}, save_path::String; kwargs...) = save_solution(res.t, res.saveval, save_path; kwargs...)

function load_solution(load_path)
    file_contents = read_JSON(Dict{Symbol, Vector{String}}, "$(load_path)\\file_contents.json")
    n_partitions = length(file_contents[:u])
    if n_partitions == 1 #there are no partitions in the file
        sol_t = BSON.load("$(load_path)\\$(file_contents[:t][1])")[:t]
        sol_u = BSON.load("$(load_path)\\$(file_contents[:u][1])")[:u]
    else
        sol_t = Float64[]
        sol_u = Vector{Vector{Float64}}()
        for i in 1:n_partitions
            sol_ti = BSON.load("$(load_path)\\$(file_contents[:t][i])")[:t]
            sol_ui = BSON.load("$(load_path)\\$(file_contents[:u][i])")[:u]
            push!(sol_t, sol_ti...)
            push!(sol_u, sol_ui...)
        end
    end

    p_dict = read_JSON(Dict{Symbol, Float32}, "$(load_path)\\params.json")
    p0 = p_dict |> extract_dict
    net = Network(p_dict[:nx], p_dict[:ny]; μ = p_dict[:μ])
    sol_prob = SDEProblem(net, noise, sol_u[1], (0f0 , p_dict[:t_warm]), p0)
    SciMLBase.build_solution(sol_prob, SOSRI(), sol_t, sol_u) #we can use this to build a solution without GPU
end

"""
This function runs the model using the indicated parameters
idxs: 
1 = Voltage
2 = n
3 = m
4 = h
5 = Ca
"""
function run_model(p_dict::Dict{Symbol,T}, u_dict::Dict{Symbol,T};
    tmax=120e3, xmax=64, ymax=64, warmup_tmax = 120e3, DEmodel = T_PDE_w_NA, verbose = true,
    idx = 5, #
    kwargs...
) where {T<:Real}

    if verbose
        print("[$(now())]: Setting up parameters, conditions, and network settings... ")
    end
    u0 = extract_dict(u_dict, t_conds, dims=(xmax, ymax))
    p = p_dict |> extract_dict
    if verbose
        println("Complete")
        print("[$(now())]: Warming up the model for $(round(warmup_tmax/1000))s... ")
    end
    prob = SDEProblem(DEmodel, noise, u0, (0.0, warmup_tmax), p)
    warmup = solve(prob, save_everystep=false, progress=true, progress_steps=1; kwargs...)
    if verbose
        println("Completed")
        print("[$(now())]: Simulating up the model for $(round(tmax/1000))s... ")
    end
    prob = SDEProblem(DEmodel, noise, warmup[end], (0.0, tmax), p)
    start_idx = idx*(xmax*ymax)-(xmax*ymax)+1
    end_idx = idx*(xmax*ymax)
    sol = solve(prob, progress=true, progress_steps=1, save_idxs=[start_idx:end_idx...]; kwargs...)
    if verbose
        println("Completed")
    end
    return sol
end 

function run_model(p_dict::Dict{Symbol,T}, u_dict::Dict{Symbol,T}, loc::String;
    plot_hists = false, 
    tmax=120e3, xmax=64, ymax=64, animate_dt = 60.0,
    kwargs...
) where {T<:Real}
    sol = run_model(p_dict, u_dict; tmax=tmax, xmax=xmax, ymax=ymax, kwargs...)
    if !isdir(loc) #If the directory doesn't exist, make it
        println("directory doesn't exist. Making it")
        mkdir(loc)
    end
    #animate_dt = 60.0
    #anim = @animate for t = 1.0:animate_dt:sol.t[end]
    #    println("[$(now())]: Animating simulation $(t) out of $(sol.t[end])...")
    #    frame_i = reshape(sol(t) |> Array, (nx, ny))
    #    heatmap(frame_i, ratio=:equal, grid=false, xaxis="", yaxis="", xlims=(0, nx), ylims=(0, ny), c=:curl, clims=(-90.0, 0.0))
    #end
    #gif(anim, "$(loc)/regular_animation.gif", fps=1000.0 / animate_dt)
    
    timestamps, data = timeseries_analysis(sol, loc)
    if plot_hists
        try
            hist_plot = plot_histograms(data, loc)
            return timestamps, data, hist_plot
        catch
            println("Something went wrong in plotting")
            return timestamps, data
        end
    else
        return timestamps, data
    end
end