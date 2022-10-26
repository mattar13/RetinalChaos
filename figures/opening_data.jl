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
#Filter out any spike duration under 5.0
phys_spike_durs = phys_spike_durs[phys_spike_durs.>5.0]
#%% Extract all of the histograms
sdur_hfit = fit(Histogram, phys_spike_durs, LinRange(0.0, 100.0, 50))
sdur_weights = sdur_hfit.weights / maximum(sdur_hfit.weights)
sdur_edges = collect(sdur_hfit.edges[1])[1:length(sdur_weights)]

isi_hfit = fit(Histogram, phys_isis, LinRange(0.0, 100.0, 50))
isi_weights = isi_hfit.weights / maximum(isi_hfit.weights)
isi_edges = collect(isi_hfit.edges[1])[1:length(isi_weights)]

bdur_hfit = fit(Histogram, phys_burst_durs, LinRange(0.0, 2000.0, 50))
bdur_weights = bdur_hfit.weights / maximum(bdur_hfit.weights)
bdur_edges = collect(bdur_hfit.edges[1])[1:length(bdur_weights)]

ibi_hfit = fit(Histogram, phys_ibis, LinRange(0.0, 120e3, 50))
ibi_weights = ibi_hfit.weights / maximum(ibi_hfit.weights)
ibi_edges = collect(ibi_hfit.edges[1])[1:length(ibi_weights)]

println("Return all the results")
n_recordings = length(paths)

#Average the stats
phys_baseline_avg = sum(phys_baseline) / n_recordings
phys_baseline_SEM = std(phys_baseline) / sqrt(n_recordings)

phys_min_avg = sum(phys_min) / n_recordings
phys_min_SEM = std(phys_min) / sqrt(n_recordings)

phys_max_avg = sum(phys_max) / n_recordings
phys_max_SEM = std(phys_max) / sqrt(n_recordings)

phys_spike_avg = sum(phys_spike_durs) / length(phys_spike_durs)
phys_spike_SEM = std(phys_spike_durs) / sqrt(length(phys_spike_durs))

phys_burst_avg = sum(phys_burst_durs) / length(phys_burst_durs)
phys_burst_SEM = std(phys_burst_durs) / sqrt(length(phys_burst_durs))

phys_ibis_avg = sum(phys_ibis) / length(phys_ibis)
phys_ibis_SEM = std(phys_ibis) / sqrt(length(phys_ibis))

#%%1) open the physiological data loaded from the single cell recordings I made
print("Opening physiological data... ")
dt = 0.5
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
dataŶ = readABF(target_file, channels=["Vm_prime4"], stimulus_name="Im_sec5", time_unit=:ms)
downsample!(dataŶ, 1 / dt)
tsPHYS, dataPHYS = timeseries_analysis(dataŶ)
println("Complete")

print("Calculating the loss... ")
spike_t_PHYS, spike_vt_PHYS = extract_spike_trace(tsPHYS, dataPHYS, idx=5, dt=dt, normalize=false)
burst_t_PHYS, burst_vt_PHYS = extract_burst_trace(tsPHYS, dataPHYS, idx=2, dt=dt, normalize=false)
IBI_t_PHYS, IBI_vt_PHYS = extract_IBI_trace(tsPHYS, dataPHYS, idx=2, dt=dt, normalize=false)
println("Complete")
#%% Open all of the data for the wave models
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"
#eXTRACT THE ISOLATED PATH 
isolated_path = "$(data_root)/isolated_model"
dataISO = load("$(isolated_path)/data.jld2")
tsISO = load("$(isolated_path)/timestamps.jld2")
tsISO = convert(Dict{String,Vector{Matrix{Float64}}}, tsISO)

iso_xIdx = rand(findall(map(x -> size(x, 1) > 1, tsISO["Bursts"])))
iso_xIdx = 2739
spike_t_ISO, spike_vt_ISO = extract_spike_trace(tsISO, dataISO, normalize=false, cell_n=iso_xIdx, idx=1)
burst_t_ISO, burst_vt_ISO = extract_burst_trace(tsISO, dataISO, normalize=false, cell_n=iso_xIdx, idx=1)
IBI_t_ISO, IBI_vt_ISO = extract_IBI_trace(tsISO, dataISO, normalize=false, cell_n=iso_xIdx, idx=1)

iso_baselines = calculate_threshold(dataISO["DataArray"], Z=1.0, dims=2)
iso_baseline_avg = sum(iso_baselines) / length(iso_baselines)
iso_baseline_SEM = std(iso_baselines) / sqrt(length(iso_baselines))

iso_min = minimum(dataISO["DataArray"], dims=2)
iso_min_avg = sum(iso_min) / length(iso_min)
iso_min_SEM = std(iso_min) / sqrt(length(iso_min))

iso_max = maximum(dataISO["DataArray"], dims=2)
iso_max_avg = sum(iso_max) / length(iso_max)
iso_max_SEM = std(iso_max) / sqrt(length(iso_max))

iso_sdur_hfit = fit(Histogram, dataISO["SpikeDurs"], LinRange(0.0, 100.0, 50))
iso_sdur_weights = iso_sdur_hfit.weights / maximum(iso_sdur_hfit.weights)
iso_sdur_edges = collect(iso_sdur_hfit.edges[1])[1:length(iso_sdur_weights)]

iso_isi_hfit = fit(Histogram, dataISO["ISIs"], LinRange(0.0, 100.0, 50))
iso_isi_weights = iso_isi_hfit.weights / maximum(iso_isi_hfit.weights)
iso_isi_edges = collect(iso_isi_hfit.edges[1])[1:length(iso_isi_weights)]

iso_bdur_hfit = fit(Histogram, dataISO["BurstDurs"], LinRange(0.0, 2000.0, 50))
iso_bdur_weights = iso_bdur_hfit.weights / maximum(iso_bdur_hfit.weights)
iso_bdur_edges = collect(iso_bdur_hfit.edges[1])[1:length(iso_bdur_weights)]

iso_ibi_hfit = fit(Histogram, dataISO["IBIs"], LinRange(0.0, 120000, 50))
iso_ibi_weights = iso_ibi_hfit.weights / maximum(iso_ibi_hfit.weights)
iso_ibi_edges = collect(iso_ibi_hfit.edges[1])[1:length(iso_ibi_weights)]

#%% We need to run this one again and for a longer period
print("Opening No GABA data... ")
noGABA_path = "$(data_root)/no_GABA_model"
dataNG = load("$(noGABA_path)/data.jld2")
tsNG = load("$(noGABA_path)/timestamps.jld2")
tsNG = convert(Dict{String,Vector{Matrix{Float64}}}, tsNG)
#ng_xIdx = rand(findall(tsNG["Bursts"] .!= nothing))
ng_xIdx = 2066
#ng_xIdx = rand(findall(map(x -> size(x, 1) >= 2, tsNG["Bursts"])))
spike_t_NG, spike_vt_NG = extract_spike_trace(tsNG, dataNG, normalize=false, cell_n=ng_xIdx, idx=1)
burst_t_NG, burst_vt_NG = extract_burst_trace(tsNG, dataNG, normalize=false, cell_n=ng_xIdx, idx=1)
IBI_t_NG, IBI_vt_NG = extract_IBI_trace(tsNG, dataNG, normalize=false, cell_n=ng_xIdx, idx=1)
println("Complete")

#Calculate the MSE between each of the data points
ng_baselines = calculate_threshold(dataNG["DataArray"], Z=1.0, dims=2)
ng_baseline_avg = sum(ng_baselines) / length(ng_baselines)
ng_baseline_SEM = std(ng_baselines) / sqrt(length(ng_baselines))

ng_min = minimum(dataNG["DataArray"], dims=2)
ng_min_avg = sum(ng_min) / length(ng_min)
ng_min_SEM = std(ng_min) / sqrt(length(ng_min))

ng_max = maximum(dataNG["DataArray"], dims=2)
ng_max_avg = sum(ng_max) / length(ng_max)
ng_max_SEM = std(ng_max) / sqrt(length(ng_max))

noGABA_sdur_hfit = fit(Histogram, dataNG["SpikeDurs"], LinRange(0.0, 100.0, 50))
noGABA_sdur_weights = noGABA_sdur_hfit.weights / maximum(noGABA_sdur_hfit.weights)
noGABA_sdur_edges = collect(noGABA_sdur_hfit.edges[1])[1:length(noGABA_sdur_weights)]

noGABA_isi_hfit = fit(Histogram, dataNG["ISIs"], LinRange(0.0, 100.0, 50))
noGABA_isi_weights = noGABA_isi_hfit.weights / maximum(noGABA_isi_hfit.weights)
noGABA_isi_edges = collect(noGABA_isi_hfit.edges[1])[1:length(noGABA_isi_weights)]

noGABA_bdur_hfit = fit(Histogram, dataNG["BurstDurs"], LinRange(0.0, 2000.0, 50))
noGABA_bdur_weights = noGABA_bdur_hfit.weights / maximum(noGABA_bdur_hfit.weights)
noGABA_bdur_edges = collect(noGABA_bdur_hfit.edges[1])[1:length(noGABA_bdur_weights)]

noGABA_ibi_hfit = fit(Histogram, dataNG["IBIs"], LinRange(0.0, 120000, 50))
noGABA_ibi_weights = noGABA_ibi_hfit.weights / maximum(noGABA_ibi_hfit.weights)
noGABA_ibi_edges = collect(noGABA_ibi_hfit.edges[1])[1:length(noGABA_ibi_weights)]

#%% Run the GABA reversal model
print("Opening No GABA data... ")
ECl55_path = "$(data_root)/ECl55_model"
dataEC = load("$(ECl55_path)/data.jld2")
tsEC = load("$(ECl55_path)/timestamps.jld2")
tsEC = convert(Dict{String,Vector{Matrix{Float64}}}, tsEC)

ec_xIdx = 3957
#ec_xIdx = rand(findall(map(x -> size(x, 1) >= 2, tsEC["Bursts"])))
spike_t_EC, spike_vt_EC = extract_spike_trace(tsEC, dataEC, normalize=false, cell_n=ec_xIdx, idx=1)
burst_t_EC, burst_vt_EC = extract_burst_trace(tsEC, dataEC, normalize=false, cell_n=ec_xIdx, idx=1)
IBI_t_EC, IBI_vt_EC = extract_IBI_trace(tsEC, dataEC, normalize=false, cell_n=ec_xIdx, idx=1)
println("Complete")

print("Calculating loss of No GABA data... ")
println("Complete")

#Calculate the MSE between each of the data points
ec_baselines = calculate_threshold(dataEC["DataArray"], Z=1.0, dims=2)
ec_baseline_avg = sum(ec_baselines) / length(ec_baselines)
ec_baseline_SEM = std(ec_baselines) / sqrt(length(ec_baselines))

ec_min = minimum(dataEC["DataArray"], dims=2)
ec_min_avg = sum(ec_min) / length(ec_min)
ec_min_SEM = std(ec_min) / sqrt(length(ec_min))

ec_max = maximum(dataEC["DataArray"], dims=2)
ec_max_avg = sum(ec_max) / length(ec_max)
ec_max_SEM = std(ec_max) / sqrt(length(ec_max))

ec_sdur_hfit = fit(Histogram, dataEC["SpikeDurs"], LinRange(0.0, 100.0, 50))
ec_sdur_weights = ec_sdur_hfit.weights / maximum(ec_sdur_hfit.weights)
ec_sdur_edges = collect(ec_sdur_hfit.edges[1])[1:length(ec_sdur_weights)]

ec_isi_hfit = fit(Histogram, dataEC["ISIs"], LinRange(0.0, 100.0, 50))
ec_isi_weights = ec_isi_hfit.weights / maximum(ec_isi_hfit.weights)
ec_isi_edges = collect(ec_isi_hfit.edges[1])[1:length(ec_isi_weights)]

ec_bdur_hfit = fit(Histogram, dataEC["BurstDurs"], LinRange(0.0, 2000.0, 50))
ec_bdur_weights = ec_bdur_hfit.weights / maximum(ec_bdur_hfit.weights)
ec_bdur_edges = collect(ec_bdur_hfit.edges[1])[1:length(ec_bdur_weights)]

ec_ibi_hfit = fit(Histogram, dataEC["IBIs"], LinRange(0.0, 120000, 50))
ec_ibi_weights = ec_ibi_hfit.weights / maximum(ec_ibi_hfit.weights)
ec_ibi_edges = collect(ec_ibi_hfit.edges[1])[1:length(ec_ibi_weights)]

dataEC["SpikeDurAvg"]
dataEC["SpikeDurSEM"]

dataEC["BurstDurAvg"]
dataEC["BurstDurSEM"]

dataEC["IBIAvg"]
dataEC["IBISEM"]
# Extract the wave model
print("Opening wave data... ")
wave_path = "$(data_root)/wave_model"
dataWAVE = load("$(wave_path)/data.jld2")
tsWAVE = load("$(wave_path)/timestamps.jld2")
tsWAVE = convert(Dict{String,Vector{Matrix{Float64}}}, tsWAVE)
#wave_xIdx = rand(findall(tsWAVE["Bursts"] .!= nothing))
ec_xIdx = rand(findall(map(x -> size(x, 1) >= 2, tsWAVE["Bursts"])))
wave_xIdx = 1801
spike_t_WAVE, spike_vt_WAVE = extract_spike_trace(tsWAVE, dataWAVE, normalize=false, cell_n=wave_xIdx, idx=1)
burst_t_WAVE, burst_vt_WAVE = extract_burst_trace(tsWAVE, dataWAVE, normalize=false, cell_n=wave_xIdx, idx=1)
IBI_t_WAVE, IBI_vt_WAVE = extract_IBI_trace(tsWAVE, dataWAVE, normalize=false, cell_n=wave_xIdx, idx=1)
println("Complete")

wave_baselines = calculate_threshold(dataWAVE["DataArray"], Z=1.0, dims=2)
wave_baseline_avg = sum(wave_baselines) / length(wave_baselines)
wave_baseline_SEM = std(wave_baselines) / sqrt(length(wave_baselines))

wave_min = minimum(dataWAVE["DataArray"], dims=2)
wave_min_avg = sum(wave_min) / length(wave_min)
wave_min_SEM = std(wave_min) / sqrt(length(wave_min))

wave_max = maximum(dataWAVE["DataArray"], dims=2)
wave_max_avg = sum(wave_max) / length(wave_max)
wave_max_SEM = std(wave_max) / sqrt(length(wave_max))

wave_sdur_hfit = fit(Histogram, dataWAVE["SpikeDurs"], LinRange(0.0, 100.0, 50))
wave_sdur_weights = wave_sdur_hfit.weights / maximum(wave_sdur_hfit.weights)
wave_sdur_edges = collect(wave_sdur_hfit.edges[1])[1:length(wave_sdur_weights)]

wave_isi_hfit = fit(Histogram, dataWAVE["ISIs"], LinRange(0.0, 100.0, 50))
wave_isi_weights = wave_isi_hfit.weights / maximum(wave_isi_hfit.weights)
wave_isi_edges = collect(wave_isi_hfit.edges[1])[1:length(wave_isi_weights)]

wave_bdur_hfit = fit(Histogram, dataWAVE["BurstDurs"], LinRange(0.0, 2000.0, 50))
wave_bdur_weights = wave_bdur_hfit.weights / maximum(wave_bdur_hfit.weights)
wave_bdur_edges = collect(wave_bdur_hfit.edges[1])[1:length(wave_bdur_weights)]

wave_ibi_hfit = fit(Histogram, dataWAVE["IBIs"], LinRange(0.0, 120e3, 50))
wave_ibi_weights = wave_ibi_hfit.weights / maximum(wave_ibi_hfit.weights)
wave_ibi_edges = collect(wave_ibi_hfit.edges[1])[1:length(wave_ibi_weights)]