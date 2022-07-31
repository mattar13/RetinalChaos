#%% Walk through each file analyzing using timescale analysis
max_spike_duration = 50.0
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
          channels=["Im_scaled"], time_unit=:ms,
          stimulus_name=nothing, flatten_episodic=true
     )
     data_i = downsample(data_raw, 1 / 1.0) #downsample the data
     data_i = highpass_filter(data_i, freq=0.01)
     timestamps, data = timeseries_analysis(data_i)

     push!(
          phys_baseline,
          sum(data_i.data_array[1, :, 1]) / length(data_i.data_array[1, :, 1])
     )

     push!(phys_min, minimum(data_i))
     push!(phys_max, maximum(data_i))
     #println(min_phys)

     push!(phys_spike_durs, data["SpikeDurs"]...)
     push!(phys_isis, data["ISIs"]...)
     push!(phys_burst_durs, data["BurstDurs"]...)
     push!(phys_ibis, data["IBIs"]...)

     println("Success")
end

#%% Extract all of the histograms
sdur_hfit = fit(Histogram, phys_spike_durs, LinRange(0.0, 50.0, 50))
sdur_weights = sdur_hfit.weights / maximum(sdur_hfit.weights)
sdur_edges = collect(sdur_hfit.edges[1])[1:length(sdur_weights)]

isi_hfit = fit(Histogram, phys_isis, LinRange(0.0, 50.0, 50))
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
phys_baseline_sem = std(phys_baseline) / sqrt(n_recordings)

phys_min_avg = sum(phys_min) / n_recordings
phys_min_sem = std(phys_min) / sqrt(n_recordings)

phys_max_avg = sum(phys_max) / n_recordings
phys_max_sem = std(phys_max) / sqrt(n_recordings)

phys_spike_avg = sum(phys_spike_durs) / length(phys_spike_durs)
phys_spike_sem = std(phys_spike_durs) / sqrt(length(phys_spike_durs))

phys_burst_avg = sum(phys_burst_durs) / length(phys_burst_durs)
phys_burst_sem = std(phys_burst_durs) / sqrt(length(phys_burst_durs))

phys_ibis_avg = sum(phys_ibis) / length(phys_ibis)
phys_ibis_sem = std(phys_ibis) / sqrt(length(phys_ibis))

#%%1) open the physiological data loaded from the single cell recordings I made
dt = 0.5
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
data = readABF(target_file, channels=["Vm_prime4"], stimulus_name=nothing, time_unit=:ms)
#data - 25.0
example_timestamps, example_data = timeseries_analysis(data)
ex_bursts = example_timestamps["Bursts"][1]

spike_t_PHYS, spike_vt_PHYS = extract_spike_trace(tsPHYS, dataPHYS, idx=3, dt=dt, normalize=false)
burst_t_PHYS, burst_vt_PHYS = extract_burst_trace(tsPHYS, dataPHYS, idx=2, dt=dt, normalize=false)
IBI_t_PHYS, IBI_vt_PHYS = extract_IBI_trace(tsPHYS, dataPHYS, idx=2, dt=dt, normalize=false)

#%% Open all of the data for the wave models
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"
region = (-100, 2500)
#eXTRACT THE ISOLATED PATH 
isolated_path = "$(data_root)/isolated_model"
isolated_data = load("$(isolated_path)/data.jld2")
isolated_timestamps = load("$(isolated_path)/timestamps.jld2")

iso_bursts = isolated_timestamps["Bursts"]
#iso_xIdx = rand(findall(iso_bursts .!= nothing))
iso_xIdx = 3213

spike_t_ISO, spike_vt_ISO = extract_spike_trace(isolated_timestamps, isolated_data, normalize=false, cell_n=iso_xIdx, idx=1)
burst_t_ISO, burst_vt_ISO = extract_burst_trace(isolated_timestamps, isolated_data, normalize=false, cell_n=iso_xIdx, idx=1)
IBI_t_ISO, IBI_vt_ISO = extract_IBI_trace(isolated_timestamps, isolated_data, normalize=false, cell_n=iso_xIdx, idx=1)

spikeMSE_ISO = SpikeLoss(example_timestamps, example_data, isolated_timestamps, isolated_data, normalize=true)
burstMSE_ISO = BurstLoss(example_timestamps, example_data, isolated_timestamps, isolated_data, normalize=true)
IBIMSE_ISO = IBILoss(example_timestamps, example_data, isolated_timestamps, isolated_data, normalize=true)

iso_sdur_hfit = fit(Histogram, isolated_data["SpikeDurs"], LinRange(0.0, 50.0, 50))
iso_sdur_weights = iso_sdur_hfit.weights / maximum(iso_sdur_hfit.weights)
iso_sdur_edges = collect(iso_sdur_hfit.edges[1])[1:length(iso_sdur_weights)]

iso_isi_hfit = fit(Histogram, isolated_data["ISIs"], LinRange(0.0, 50.0, 50))
iso_isi_weights = iso_isi_hfit.weights / maximum(iso_isi_hfit.weights)
iso_isi_edges = collect(iso_isi_hfit.edges[1])[1:length(iso_isi_weights)]

iso_bdur_hfit = fit(Histogram, isolated_data["BurstDurs"], LinRange(0.0, 2000.0, 50))
iso_bdur_weights = iso_bdur_hfit.weights / maximum(iso_bdur_hfit.weights)
iso_bdur_edges = collect(iso_bdur_hfit.edges[1])[1:length(iso_bdur_weights)]

iso_ibi_hfit = fit(Histogram, isolated_data["IBIs"], LinRange(0.0, 120000, 50))
iso_ibi_weights = iso_ibi_hfit.weights / maximum(iso_ibi_hfit.weights)
iso_ibi_edges = collect(iso_ibi_hfit.edges[1])[1:length(iso_ibi_weights)]

#We need to run this one again and for a longer period
noGABA_path = "$(data_root)/no_GABA_model"
noGABA_data = load("$(noGABA_path)/data.jld2")
noGABA_timestamps = load("$(noGABA_path)/timestamps.jld2")


ng_bursts = noGABA_timestamps["Bursts"]
ng_xIdx = rand(findall(ng_bursts .!= nothing))
ng_xIdx = 4065
ng_bursts[ng_xIdx]
noGABA_data["Time"]
noGABA_timestamps["Bursts"][ng_xIdx]

spike_t_NG, spike_vt_NG = extract_spike_trace(noGABA_timestamps, noGABA_data, normalize=false, cell_n=ng_xIdx, idx=1)
burst_t_NG, burst_vt_NG = extract_burst_trace(noGABA_timestamps, noGABA_data, normalize=false, cell_n=ng_xIdx, idx=1)
IBI_t_NG, IBI_vt_NG = extract_IBI_trace(noGABA_timestamps, noGABA_data, normalize=false, cell_n=ng_xIdx, idx=1)
spikeMSE_NG = SpikeLoss(example_timestamps, example_data, noGABA_timestamps, noGABA_data, normalize=true)
burstMSE_NG = BurstLoss(example_timestamps, example_data, noGABA_timestamps, noGABA_data, normalize=true)
IBIMSE_NG = IBILoss(example_timestamps, example_data, noGABA_timestamps, noGABA_data, normalize=true)
#Calculate the MSE between each of the data points
#spikeMSE_NG

noGABA_sdur_hfit = fit(Histogram, noGABA_data["SpikeDurs"], LinRange(0.0, 50.0, 50))
noGABA_sdur_weights = noGABA_sdur_hfit.weights / maximum(noGABA_sdur_hfit.weights)
noGABA_sdur_edges = collect(noGABA_sdur_hfit.edges[1])[1:length(noGABA_sdur_weights)]

noGABA_isi_hfit = fit(Histogram, noGABA_data["ISIs"], LinRange(0.0, 50.0, 50))
noGABA_isi_weights = noGABA_isi_hfit.weights / maximum(noGABA_isi_hfit.weights)
noGABA_isi_edges = collect(noGABA_isi_hfit.edges[1])[1:length(noGABA_isi_weights)]

noGABA_bdur_hfit = fit(Histogram, noGABA_data["BurstDurs"], LinRange(0.0, 2000.0, 50))
noGABA_bdur_weights = noGABA_bdur_hfit.weights / maximum(noGABA_bdur_hfit.weights)
noGABA_bdur_edges = collect(noGABA_bdur_hfit.edges[1])[1:length(noGABA_bdur_weights)]

noGABA_ibi_hfit = fit(Histogram, noGABA_data["IBIs"], LinRange(0.0, 120000, 50))
noGABA_ibi_weights = noGABA_ibi_hfit.weights / maximum(noGABA_ibi_hfit.weights)
noGABA_ibi_edges = collect(noGABA_ibi_hfit.edges[1])[1:length(noGABA_ibi_weights)]


# Extract the wave model
wave_path = "$(data_root)/wave_model"
wave_data = load("$(wave_path)/data.jld2")
wave_timestamps = load("$(wave_path)/timestamps.jld2")


wave_bursts = wave_timestamps["Bursts"]
wave_xIdx = rand(findall(wave_bursts .!= nothing))
wave_xIdx = 1838

spike_t_WAVE, spike_vt_WAVE = extract_spike_trace(wave_timestamps, wave_data, normalize=false, cell_n=wave_xIdx, idx=1)
burst_t_WAVE, burst_vt_WAVE = extract_burst_trace(wave_timestamps, wave_data, normalize=false, cell_n=wave_xIdx, idx=1)
IBI_t_WAVE, IBI_vt_WAVE = extract_IBI_trace(wave_timestamps, wave_data, normalize=false, cell_n=wave_xIdx, idx=1)
spikeMSE_WAVE = SpikeLoss(example_timestamps, example_data, wave_timestamps, wave_data, normalize=true)
burstMSE_WAVE = BurstLoss(example_timestamps, example_data, wave_timestamps, wave_data, normalize=true)
IBIMSE_WAVE = IBILoss(example_timestamps, example_data, wave_timestamps, wave_data, normalize=true)

wave_sdur_hfit = fit(Histogram, wave_data["SpikeDurs"], LinRange(0.0, 50.0, 50))
wave_sdur_weights = wave_sdur_hfit.weights / maximum(wave_sdur_hfit.weights)
wave_sdur_edges = collect(wave_sdur_hfit.edges[1])[1:length(wave_sdur_weights)]

wave_isi_hfit = fit(Histogram, wave_data["ISIs"], LinRange(0.0, 50.0, 50))
wave_isi_weights = wave_isi_hfit.weights / maximum(wave_isi_hfit.weights)
wave_isi_edges = collect(wave_isi_hfit.edges[1])[1:length(wave_isi_weights)]

wave_bdur_hfit = fit(Histogram, wave_data["BurstDurs"], LinRange(0.0, 2000.0, 50))
wave_bdur_weights = wave_bdur_hfit.weights / maximum(wave_bdur_hfit.weights)
wave_bdur_edges = collect(wave_bdur_hfit.edges[1])[1:length(wave_bdur_weights)]

wave_ibi_hfit = fit(Histogram, wave_data["IBIs"], LinRange(0.0, 120e3, 50))
wave_ibi_weights = wave_ibi_hfit.weights / maximum(wave_ibi_hfit.weights)
wave_ibi_edges = collect(wave_ibi_hfit.edges[1])[1:length(wave_ibi_weights)]



#%% Data testing 
#Find out what files are good
#=
#Good: 2,3,4,5, 9
#OK: 1, 15
#Bad: 6, 7, 8, 10, 11, 12, 13, 14, 16, 17
data_raw = readABF(paths[17],
     channels=["Im_scaled"], time_unit=:ms,
     stimulus_name=nothing, flatten_episodic=true
)
data_original = copy(data_raw)
data_original = downsample(data_raw, 1 / 1.0)
data_i = downsample(data_raw, 1 / 1.0) #downsample the data
data_i = highpass_filter(data_i, freq=0.01)
timestamps, data = timeseries_analysis(data_i)

t = data["Time"]
data_arr = data["DataArray"][1, :, 1]

fig, ax = plt.subplots(2, 3)
ax[1, 1].plot(t, data_original.data_array[1, :], c=:blue)
ax[1, 1].vlines(timestamps["Spikes"][1][:, 1], ymin=-90.0, ymax=10.0, color=:green)
ax[1, 1].vlines(timestamps["Spikes"][1][:, 2], ymin=-90.0, ymax=10.0, color=:red)

ax[2, 1].plot(t, data_arr, c=:black)
ax[2, 1].hlines(thresholds, xmin=0.0, xmax=t[end], color=:red)
ax[2, 1].vlines(timestamps["Spikes"][1][:, 1], ymin=-10.0, ymax=10.0, color=:green)
ax[2, 1].vlines(timestamps["Spikes"][1][:, 2], ymin=-10.0, ymax=10.0, color=:red)


sdur_hfit = fit(Histogram, data["SpikeDurs"], LinRange(0.0, 50.0, 50))
sdur_weights = sdur_hfit.weights / maximum(sdur_hfit.weights)
sdur_edges = collect(sdur_hfit.edges[1])[1:length(sdur_weights)]

isi_hfit = fit(Histogram, data["ISIs"], LinRange(0.0, 50.0, 50))
isi_weights = isi_hfit.weights / maximum(isi_hfit.weights)
isi_edges = collect(isi_hfit.edges[1])[1:length(isi_weights)]

bdur_hfit = fit(Histogram, data["BurstDurs"], LinRange(0.0, 2000.0, 50))
bdur_weights = bdur_hfit.weights / maximum(bdur_hfit.weights)
bdur_edges = collect(bdur_hfit.edges[1])[1:length(bdur_weights)]

ibi_hfit = fit(Histogram, data["IBIs"], LinRange(0.0, 120e3, 50))
ibi_weights = ibi_hfit.weights / maximum(ibi_hfit.weights)
ibi_edges = collect(ibi_hfit.edges[1])[1:length(ibi_weights)]


ax[1, 2].plot(sdur_edges, sdur_weights)
ax[2, 2].plot(isi_edges, isi_weights)
ax[1, 3].plot(bdur_edges, bdur_weights)
ax[2, 3].plot(ibi_edges, ibi_weights)
=#