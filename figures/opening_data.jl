#%% Walk through each file analyzing using timescale analysis
spike_dur = 25
burst_dur = 1000
IBI_dur = 60e3

max_spike_duration = 100.0
max_spike_interval = 100.0
max_burst_duration = 10e5
max_burst_interval = 10e5
target_folder = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Patching\Jordans_Patch_Data\UsuableData"
paths = target_folder |> parseABF;

phys_baseline = Float64[]
phys_max = Float64[]
phys_min = Float64[]
phys_spike_durs = Float64[]
phys_isis = Float64[]
phys_burst_durs = Float64[]
phys_ibis = Float64[]
good_paths = [2, 3, 4, 5, 9] #I have further refined the paths 
for path in paths[good_paths]
     print("[$(now())] Opening $path... ")
     data_raw = readABF(path,
          channels=["Im_scaled"], stimulus_name="IN 2", time_unit=:ms, flatten_episodic=true
     )
     data_i = downsample(data_raw, 1 / 1.0) #downsample the data
     #These have to be done before highpass filtering the data
     push!(
          phys_baseline,
          sum(data_i.data_array[1, :, 1]) / length(data_i.data_array[1, :, 1])
     )
     push!(phys_min, minimum(data_i))
     push!(phys_max, maximum(data_i))

     data_i = ePhys.filter_data(data_raw, mode=:Highpass, freq_start=0.01)
     timestamps, data = timeseries_analysis(data_raw)

     #println(min_phys)

     push!(phys_spike_durs, data["SpikeDurs"]...)
     push!(phys_isis, data["ISIs"]...)
     push!(phys_burst_durs, data["BurstDurs"]...)
     push!(phys_ibis, data["IBIs"]...)

     println("Success")
end

#%%1) open the physiological data loaded from the single cell recordings I made
print("Opening physiological data... ")
dt = 1.0
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
dataŶ = readABF(target_file, channels=["Vm_prime4"], stimulus_name="Im_sec5", time_unit=:ms)
downsample!(dataŶ, 1 / dt)
tsPHYS, dataPHYS = timeseries_analysis(dataŶ)
println("Complete")

#%% Open all of the data for the wave models
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling"

#%% Extract the Isolated ===============================================#
print("Opening Isolated data")
isolated_path = "$(data_root)/isolated_model"
dataISO = load("$(isolated_path)/data.jld2")
tsISO = load("$(isolated_path)/timestamps.jld2")
tsISO = convert(Dict{String,Vector{Matrix{Float64}}}, tsISO)
#iso_xIdx = rand(findall(map(x -> size(x, 1) > 1, tsISO["Bursts"])))
iso_xIdx = 1684
println("Complete")

#%% Extract the No GABA ===============================================#
print("Opening No GABA data... ")
noGABA_path = "$(data_root)/no_GABA_model"
dataNG = load("$(noGABA_path)/data.jld2")
tsNG = load("$(noGABA_path)/timestamps.jld2")
tsNG = convert(Dict{String,Vector{Matrix{Float64}}}, tsNG)
#ng_xIdx = rand(findall(map(x -> size(x, 1) >= 2, tsNG["Bursts"])))
ng_xIdx = 3494
println("Complete")


#%% Extract the Reversal Potential Data ===============================================#
print("Opening Reversal data... ")
ECl55_path = "$(data_root)/ECl55_model"
dataEC = load("$(ECl55_path)/data.jld2")
tsEC = load("$(ECl55_path)/timestamps.jld2")
tsEC = convert(Dict{String,Vector{Matrix{Float64}}}, tsEC)
#ec_xIdx = rand(findall(map(x -> size(x, 1) >= 2, tsEC["Bursts"])))
ec_xIdx = 2699
println("Complete")

#%% Extract the Wave ===============================================#
print("Opening wave data... ")
wave_path = "$(data_root)/FastI_EDiffusion"
dataWAVE = load("$(wave_path)/data.jld2")
tsWAVE = load("$(wave_path)/timestamps.jld2")
tsWAVE = convert(Dict{String,Vector{Matrix{Float64}}}, tsWAVE)
#wave_xIdx = rand(findall(map(x -> size(x, 1) >= 2, tsWAVE["Bursts"])))
wave_xIdx = 1874
println("Complete")

#%% Extract the Fast diffusion ===============================================#
print("Opening wave data... ")
FASTe_path = wave_path = "$(data_root)/FastEdiffusion"
dataFASTe = load("$(FASTe_path)/data.jld2")
tsFASTe = load("$(FASTe_path)/timestamps.jld2")
tsFASTe = convert(Dict{String,Vector{Matrix{Float64}}}, tsFASTe)

print("Opening wave data... ")
FASTi_path = wave_path = "$(data_root)/FastIdiffusion"
dataFASTi = load("$(FASTi_path)/data.jld2")
tsFASTi = load("$(FASTi_path)/timestamps.jld2")
tsFASTi = convert(Dict{String,Vector{Matrix{Float64}}}, tsFASTi)

print("Opening wave data... ")
FASTei_path = wave_path = "$(data_root)/FastI_Ediffusion"
dataFASTei = load("$(FASTei_path)/data.jld2")
tsFASTei = load("$(FASTei_path)/timestamps.jld2")
tsFASTei = convert(Dict{String,Vector{Matrix{Float64}}}, tsFASTei)

#%% SLOW=====================================================================================#
SLOWe_path = wave_path = "$(data_root)/SlowEdiffusion"
dataSLOWe = load("$(SLOWe_path)/data.jld2")
tsSLOWe = load("$(SLOWe_path)/timestamps.jld2")
tsSLOWe = convert(Dict{String,Vector{Matrix{Float64}}}, tsSLOWe)

print("Opening wave data... ")
SLOWi_path = wave_path = "$(data_root)/SlowIdiffusion"
dataSLOWi = load("$(SLOWi_path)/data.jld2")
tsSLOWi = load("$(SLOWi_path)/timestamps.jld2")
tsSLOWi = convert(Dict{String,Vector{Matrix{Float64}}}, tsSLOWi)

print("Opening wave data... ")
SLOWei_path = wave_path = "$(data_root)/SlowI_Ediffusion"
dataSLOWei = load("$(SLOWei_path)/data.jld2")
tsSLOWei = load("$(SLOWei_path)/timestamps.jld2")
tsSLOWei = convert(Dict{String,Vector{Matrix{Float64}}}, tsSLOWei)
