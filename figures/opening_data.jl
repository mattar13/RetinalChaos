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
for path in paths
     print("[$(now())] Opening $path... ")
     data_raw = readABF(path,
          channels=["Im_scaled"],
          stimulus_name=nothing, flatten_episodic=true
     )
     data_i = data_raw
     #data_i = dwt_filter(data_raw, period_window=(11, 18))

     t = data_i.t * 1000 #Turn seconds into milliseconds
     vm_array = data_i.data_array[:, :, 1]
     thresholds = RetinalChaos.calculate_threshold(vm_array, Z=4.0, dims=2)
     println("Complete")
     print("[$(now())]: Extracting the spikes... ")
     spike_array = Matrix{Bool}(vm_array .> thresholds)
     spikes = RetinalChaos.get_timestamps(spike_array, t)
     spike_durs, isi = RetinalChaos.extract_interval(spikes,
          min_duration=1.0, min_interval=1.0,
          max_duration=max_spike_duration, max_interval=max_spike_interval
     )
     println("Complete")

     print("[$(now())]: Extracting the bursts... ")
     bursts = RetinalChaos.max_interval_algorithim(spikes)
     burst_durs, ibi = RetinalChaos.extract_interval(bursts[1], max_duration=max_burst_duration, max_interval=max_burst_interval, min_duration=1000.0)
     println("Complete")

     push!(
          phys_baseline,
          sum(data_i.data_array[1, :, 1]) / length(data_i.data_array[1, :, 1])
     )

     push!(phys_min, minimum(data_i))
     push!(phys_max, maximum(data_i))
     #println(min_phys)

     push!(phys_spike_durs, (spike_durs)...)
     push!(phys_isis, (isi)...)

     push!(phys_burst_durs, (burst_durs)...)
     push!(phys_ibis, (ibi)...)
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
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
data = readABF(target_file, channels=["Vm_prime4"], stimulus_name=nothing, time_unit=:ms)
data - 25.0
example_timestamps, example_data = timeseries_analysis(data.t, data.data_array[:, :, 1])
ex_bursts = timestamps["Bursts"][1]

t_phys_burst = ex_bursts[1, 1]-100:1.0:ex_bursts[1, 1]+2500
phys_burst_idxs = round.(Int64, t_phys_burst ./ data.dt)

t_phys_burst = t_phys_burst .- t_phys_burst[1]
vt_phys_burst = data.data_array[1, phys_burst_idxs, 1]

#%% Open all of the data for the wave models
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"
#eXTRACT THE ISOLATED PATH 
isolated_path = "$(data_root)/isolated_model"
isolated_data = load("$(isolated_path)/data.jld2")
isolated_timestamps = load("$(isolated_path)/timestamps.jld2")

iso_bursts = isolated_timestamps["Bursts"]
iso_xIdx = rand(findall(iso_bursts .!= nothing))
iso_burst = iso_bursts[iso_xIdx]
iso_burst_idx = round.(Int64, (iso_burst[1, 1]-100):1.0:(iso_burst[1, 1]+2500))

t_iso_burst = isolated_data["Time"][iso_burst_idx]
t_iso_burst .-= t_iso_burst[1]
vt_iso_burst = isolated_data["DataArray"][iso_xIdx, iso_burst_idx]

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

#eXTRACT THE NO GABA
noGABA_path = "$(data_root)/no_GABA_model"
noGABA_data = load("$(noGABA_path)/data.jld2")
noGABA_timestamps = load("$(noGABA_path)/timestamps.jld2")

ng_bursts = noGABA_timestamps["Spikes"]
ng_xIdx = rand(findall(ng_bursts .!= nothing))
ng_burst = ng_bursts[ng_xIdx]
ng_burst_idx = round.(Int64, (ng_burst[1, 1]-100):1.0:(ng_burst[1, 1]+2500))

t_ng_burst = noGABA_data["Time"][ng_burst_idx]
t_ng_burst .-= t_ng_burst[1]
vt_ng_burst = noGABA_data["DataArray"][ng_xIdx, ng_burst_idx]

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
wave_burst = wave_bursts[wave_xIdx]
wave_burst_idx = round.(Int64, (wave_burst[1, 1]-100):1.0:(wave_burst[1, 1]+2500))

t_wave_burst = wave_data["Time"][wave_burst_idx]
t_wave_burst .-= t_wave_burst[1]
vt_wave_burst = wave_data["DataArray"][wave_xIdx, wave_burst_idx]

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

#Good: 1,2,3
#=
data_raw = readABF(paths[3],
     channels=["Im_scaled"],
     stimulus_name=nothing, flatten_episodic=true
)

#A data filtering makes the analysis slightly easier

#Calculate the data for the spikes
data = dwt_filter(data_raw, period_window=(11, 18)) #This is the ideal setting


#data = data_raw
t = data.t * 1000 #Turn seconds into milliseconds
vm_array = data.data_array[:, :, 1]
thresholds = RetinalChaos.calculate_threshold(vm_array, Z=4.0, dims=2)
println("Complete")
print("[$(now())]: Extracting the spikes... ")
spike_array = Matrix{Bool}(vm_array .> thresholds)
spikes = RetinalChaos.get_timestamps(spike_array, t)
spike_durs, isi = RetinalChaos.extract_interval(spikes,
     min_duration=1.0, min_interval=1.0,
     max_duration=max_spike_duration, max_interval=max_spike_interval
)
println("Complete")

print("[$(now())]: Extracting the bursts... ")
bursts, spb = RetinalChaos.max_interval_algorithim(spikes)
burst_durs, ibi = RetinalChaos.extract_interval(bursts[1], max_duration=max_burst_duration, max_interval=max_burst_interval, min_duration=500.0, min_interval=10e3)
println("Complete")
blims = round.(Int64, bursts[1] / data.dt / 1000)

#
fig, ax = plt.subplots(3, 2)
ax[1, 1].plot(t, data_raw.data_array[1, :], c=:black)
ax[1, 1].vlines(spikes[1][:, 1], ymin=-90.0, ymax=10.0, color=:green)
ax[1, 1].vlines(spikes[1][:, 2], ymin=-90.0, ymax=10.0, color=:red)

ax[1, 2].plot(t, data.data_array[1, :], c=:blue)
ax[1, 2].hlines(thresholds, xmin=0.0, xmax=t[end], color=:red)

sdur_hfit = fit(Histogram, spike_durs, LinRange(0.0, 50.0, 50))
sdur_weights = sdur_hfit.weights / maximum(sdur_hfit.weights)
sdur_edges = collect(sdur_hfit.edges[1])[1:length(sdur_weights)]

isi_hfit = fit(Histogram, isi, LinRange(0.0, 50.0, 50))
isi_weights = isi_hfit.weights / maximum(isi_hfit.weights)
isi_edges = collect(isi_hfit.edges[1])[1:length(isi_weights)]

bdur_hfit = fit(Histogram, burst_durs, LinRange(0.0, 2000.0, 50))
bdur_weights = bdur_hfit.weights / maximum(bdur_hfit.weights)
bdur_edges = collect(bdur_hfit.edges[1])[1:length(bdur_weights)]

ibi_hfit = fit(Histogram, ibi, LinRange(0.0, 120e3, 50))
ibi_weights = ibi_hfit.weights / maximum(ibi_hfit.weights)
ibi_edges = collect(ibi_hfit.edges[1])[1:length(ibi_weights)]


ax[2, 1].plot(sdur_edges, sdur_weights)
ax[2, 2].plot(isi_edges, isi_weights)
ax[3, 1].plot(bdur_edges, bdur_weights)
ax[3, 2].plot(ibi_edges, ibi_weights)
=#