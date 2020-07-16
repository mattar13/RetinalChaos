#Reading and writing JSON files
"""
This file writes either a named tuple or a dictionary into a JSON file
"""
function write_JSON(data, name_file)
    string_data = JSON2.write(data)
    open(name_file, "w") do f
        write(f, string_data)
    end
end

"""
This function reads JSON files into whatever type is_type is.
"""
function read_JSON(name_file::String; is_type = Dict{Symbol, Float64})
    nt = nothing
    open(name_file, "r") do f
        nt = JSON2.read(f, is_type)
    end
    nt
end

"""
Extract a parameter, condition dictionary
USAGE:
#p_dict -> [The parameter dictionary read directly from JSON]
> p0 = extract_dict(p_dict, BurstModel.params)

If your initial conditions are for a network, you can add the extra argument for dims

> In[1]: dims = (64,64);
> In[2]: p0 = extract_dict(p_dict, BurstModel.params, dims)
> Out[2]: 64x64xN AbstractArray:
"""
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}) = map(x -> Float64(dict_item[x]), pars)
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}, dims::Tuple) = cat(map(x -> fill(Float64(dict_item[x]), dims), pars)..., dims = length(dims)+1)

#get the index of the conditions in the list
u_find(cond::Symbol; list_u::Array{Symbol, 1} = tar_conds) = findall(c -> c == cond, list_u)[1]
#get the index of the parameter in the list
p_find(par::Symbol; list_p::Array{Symbol, 1}  = tar_pars) = findall(p -> p == par, list_p)[1]
#get the index of the conditions in the list
u_find(cond::Symbol; list_u = tar_conds) = findall(c -> c == cond, list_u)[1]
#get the index of the parameter in the list
p_find(par::Symbol; list_p = tar_pars) = findall(p -> p == par, list_p)[1]

"""
When using the Modeling Toolkit, the dictionary needs to be converted into an array of operations
"""
function extract_dict(dict_item::Dict{Symbol, Float64})
    par_set = nothing
    idx = 1
    for key in keys(dict_item)
        if par_set == nothing
            par_set = [Variable(key) => dict_item[key]]
        else
            push!(par_set, (Variable(key)=> dict_item[key]))
        end
    end
    par_set
end

"""
Multiple dispatch function that is called to extract a distribution dictionary
"""
function extract_dict(dict_item::Dict{Symbol, Tuple}, pars::Array{Symbol})
    d0 = Array{Distribution}([])
    for par in pars
        push!(d0, convert_dist(dict_item[par]))
    end
    d0
end

"""
This function looks in the path for the initial_condition file. If it finds it it returns it and a true success flag.
"""
function parse_ic(path, ic_path)
    success = false
    ic = nothing
    try 
        mkdir(path)
        println("[$(now())]: Path does not yet, exist, creating path")
    catch
        
        try 
            ic = jldopen(ic_path, "r") do file
                read(file, "ic")
            end
            println("[$(now())]: Previous solution found")
            success = true
        catch
            println("[$(now())]: Previous solution not found, warmup required")
        end
    end
    return success, ic
end

"""
This function creates a default set of distributions for the parameters and values entered
DEFAULT DISTRIBUTION
value > 0 ? truncated(Normal(value, 0.5), 0.0 )
"""
function create_distributions(pars, vals; std::Float64 = 0.01)
    dist_dict = Dict{Symbol, Tuple}()
    for (idx, par) in enumerate(pars)
        if vals[idx] > 0.0
            dist_dict[par] = (vals[idx], abs(vals[idx]*std), 0.0, Inf)
        elseif vals[idx] == 0.0
            dist_dict[par] = (vals[idx], abs(vals[idx]*std), -Inf, Inf)
        else
            dist_dict[par] = (vals[idx], abs(vals[idx]*std), -Inf, 0.0)
        end
    end
    dist_dict
end


function convert_dist(d_arr::NTuple{4,Union{Real, Nothing}})
    μ = d_arr[1]
    σ = d_arr[2]
    l = d_arr[3] != nothing ? d_arr[3] : -Inf
    u = d_arr[4] != nothing ? d_arr[4] : Inf
    if l == u
        Normal(μ, σ)
    else
        truncated(Normal(μ, σ ), l, u)
    end
end

"""
This function turns all conditions and parameters into a datasheet to write with excel
"""
function params_to_datasheet(timestamp, pars, conds)
    df_params = DataFrame()
    vals = (timestamp, conds..., pars...)
    keys = (:Date, BurstModel.syms..., BurstModel.params...)
    for (idx, key) in enumerate(keys)
        df_params[!, key] = [vals[idx]]
    end
    df_params
end

"""
These function are for reading and writing to excel spreadsheets
"""
function append_modeldata(filename, stats, params)
    try
        XLSX.openxlsx(filename, mode="rw") do xf
            #First we need to open the datasheets
            df_params_old = DataFrame(XLSX.readtable(filename, "Parameters")...)
            df_wavestats_old = DataFrame(XLSX.readtable(filename, "WaveStats")...)

            stats_sheet = xf[1]
            param_sheet = xf[2]

            df_params_new = vcat(df_params_old, params)
            df_wavestats_new = vcat(df_wavestats_old, stats)
            XLSX.writetable!(stats_sheet, collect(DataFrames.eachcol(df_wavestats_new)), DataFrames.names(df_wavestats_new))
            XLSX.writetable!(param_sheet, collect(DataFrames.eachcol(df_params_new)), DataFrames.names(df_params_new))
        end
    catch
        println("Excel sheet not yet made")
        XLSX.writetable(filename,
                WaveStats  = (collect(DataFrames.eachcol(stats)), DataFrames.names(stats)),
                Parameters = (collect(DataFrames.eachcol(params)), DataFrames.names(params))
            )
    end
end

function interleave(l1, l2)
    interleaved_list = Array{l1[1]|>typeof}([])
    @assert length(l1) == length(l2)
    for i = 1:length(l1)
        push!(interleaved_list, l1[i])
        push!(interleaved_list, l2[i])
    end
    interleaved_list
end

###############################################Opening and interacting with .abf files

"""
This function walks through the directory and locates any .abf file. 
The extension can be changed with the keyword argument extension
"""
function parse_abf(super_folder::String; extension::String = ".abf", verbose = false)
    file_list = []
    for (root, dirs, files) in walkdir(super_folder)
        for file in files
            if file[end-3:end] == extension
                path = joinpath(root, file)
                if verbose 
                    println(path) # path to files
                end
                push!(file_list, path)
            end
        end
    end
    file_list
end

function extract_abf(abf_path; verbose = false, v_offset = -25.0)
    if length(abf_path |> splitpath) > 1
        full_path = abf_path
    else
        full_path = joinpath(pwd(), abf_path)   
    end
    #extract the abf file by using pyABF
    exp_data = pyABF.ABF(full_path)
    #if the data is segmented into sweeps (which Jordans data is) concatenate all sweeps
    if length(exp_data.sweepList) > 1
        data = Float64[]
        time = Float64[]
        previous_time = 0.0
        for sweepNumber in exp_data.sweepList
            exp_data.setSweep(sweepNumber = sweepNumber, channel = 0);
            push!(data, exp_data.sweepY...);
            push!(time, (exp_data.sweepX.+previous_time)...);
            previous_time = time[end]
        end
        dt = time[2]*1000
    else
        exp_data.setSweep(sweepNumber = 0, channel = 0);
        data = Float64.(exp_data.sweepY);
        time = Float64.(exp_data.sweepX);
        dt = time[2]*1000
        
    end
    if verbose
        println("Data extracted from $full_path")
        println("Data from time stamp $(t[1]) s to $(t[end]+dt) s with dt = $dt ms")
        println("Data was acquired at $(1/dt/1000) Hz")
        println("$(length(t)) data points")
    end
    data, time, dt
end
