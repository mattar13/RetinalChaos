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
function read_JSON(name_file::String; is_type = NamedTuple{Symbol, Float64})
    nt = nothing
    open(name_file, "r") do f
        nt = JSON2.read(f, is_type)
    end
    nt
end

"""
Extract a parameter, condition dictionary
USAGE:
> p0 = p_dict |> x -> extract_dict(x, BurstModel.params)

OR

> p0 = extract_dict(p_dict, BurstModel.params)

"""
extract_dict(dict_item::Dict{Symbol, Float64}, pars::Array{Symbol}) = map(x -> Float64(dict_item[x]), pars)

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
        XLSX.openxlsx(filename, mode="w") do xf
            XLSX.writetable!(
                WaveStats  = (collect(DataFrames.eachcol(stats)), DataFrames.names(stats)),
                Parameters = (collect(DataFrames.eachcol(params)), DataFrames.names(params))
            )
        end
    end
end
