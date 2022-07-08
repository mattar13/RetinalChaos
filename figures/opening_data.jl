using Revise
using PhysAnalysis, ABFReader
using StatsBase, Statistics

#%% Run the data on the physiological data that Jordan has collected
println("[$(now())]: Opening Physiological data ...")
target_folder = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Patching\Jordans_Patch_Data\UsuableData"
paths = target_folder |> parseABF;

#Walk through each file analyzing using timescale analysis
max_spike_duration = 50.0
max_spike_interval = 100.0
max_burst_duration = 10e5
max_burst_interval = 10e5

phys_baseline = Float64[]
phys_max = Float64[]
phys_min = Float64[]
phys_spike_durs = Float64[]
phys_isis = Float64[]
phys_burst_durs = Float64[]
phys_ibis = Float64[]
for path in paths
     print("[$(now())] Opening $path... ")
     data_i = readABF(path,
          channels=["Im_scaled"],
          stimulus_name=nothing, flatten_episodic=true
     )
     t = data_i.t * 1000 #Turn seconds into milliseconds
     vm_array = data_i.data_array[:, :, 1]
     thresholds = RetinalChaos.calculate_threshold(vm_array, Z = 2.0, dims=2)
     println("Complete")
     print("[$(now())]: Extracting the spikes... ")
     spike_array = Matrix{Bool}(vm_array .> thresholds)
     spikes = RetinalChaos.get_timestamps(spike_array, t)
     spike_durs, isi = RetinalChaos.extract_interval(spikes, min_duration = 1.0, max_duration=max_spike_duration, max_interval=max_spike_interval)
     println("Complete")

     print("[$(now())]: Extracting the bursts... ")
     bursts = RetinalChaos.max_interval_algorithim(spikes)
     burst_durs, ibi = RetinalChaos.extract_interval(bursts[1], max_duration=max_burst_duration, max_interval=max_burst_interval, min_duration = 1000.0)
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

sdur_hfit = fit(Histogram, phys_spike_durs, LinRange(0.0, 50.0, 100))
sdur_weights = sdur_hfit.weights / maximum(sdur_hfit.weights)
sdur_edges = collect(sdur_hfit.edges[1])[1:length(sdur_weights)]

isi_hfit = fit(Histogram, phys_isis, LinRange(0.0, 100.0, 100))
isi_weights = isi_hfit.weights / maximum(isi_hfit.weights)
isi_edges = collect(isi_hfit.edges[1])[1:length(isi_weights)]

bdur_hfit = fit(Histogram, phys_burst_durs, LinRange(0.0, 2000.0, 100))
bdur_weights = bdur_hfit.weights / maximum(bdur_hfit.weights)
bdur_edges = collect(bdur_hfit.edges[1])[1:length(bdur_weights)]

ibi_hfit = fit(Histogram, phys_ibis, LinRange(0.0, 120e3, 100))
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

#%% Open all of the data for the wave models
data_root = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Modelling/figure_data"
wave_path = "$(data_root)/wave_model"
wave_data = load("$(wave_path)/data.jld2")

# Extract the histograms
wave_sdur_hfit = fit(Histogram, wave_data["SpikeDurs"], LinRange(0.0, 50.0, 100))
wave_sdur_weights = wave_sdur_hfit.weights / maximum(wave_sdur_hfit.weights)
wave_sdur_edges = collect(wave_sdur_hfit.edges[1])[1:length(wave_sdur_weights)]

wave_isi_hfit = fit(Histogram, wave_data["ISIs"], LinRange(0.0, 100.0, 100))
wave_isi_weights = wave_isi_hfit.weights / maximum(wave_isi_hfit.weights)
wave_isi_edges = collect(wave_isi_hfit.edges[1])[1:length(wave_isi_weights)]

wave_bdur_hfit = fit(Histogram, wave_data["BurstDurs"], LinRange(0.0, 2000.0, 100))
wave_bdur_weights = wave_bdur_hfit.weights / maximum(wave_bdur_hfit.weights)
wave_bdur_edges = collect(wave_bdur_hfit.edges[1])[1:length(wave_bdur_weights)]

wave_ibi_hfit = fit(Histogram, wave_data["IBIs"], LinRange(0.0, 120e3, 100))
wave_ibi_weights = wave_ibi_hfit.weights / maximum(wave_ibi_hfit.weights)
wave_ibi_edges = collect(wave_ibi_hfit.edges[1])[1:length(wave_ibi_weights)]

#eXTRACT THE ISOLATED PATH 
isolated_path = "$(data_root)/isolated_model"
isolated_data = load("$(isolated_path)/data.jld2")

iso_sdur_hfit = fit(Histogram, isolated_data["SpikeDurs"], LinRange(0.0, 50.0, 100))
iso_sdur_weights = iso_sdur_hfit.weights / maximum(iso_sdur_hfit.weights)
iso_sdur_edges = collect(iso_sdur_hfit.edges[1])[1:length(iso_sdur_weights)]

iso_isi_hfit = fit(Histogram, isolated_data["ISIs"], LinRange(0.0, 100.0, 100))
iso_isi_weights = iso_isi_hfit.weights / maximum(iso_isi_hfit.weights)
iso_isi_edges = collect(iso_isi_hfit.edges[1])[1:length(iso_isi_weights)]

iso_bdur_hfit = fit(Histogram, isolated_data["BurstDurs"], LinRange(0.0, 2000.0, 100))
iso_bdur_weights = iso_bdur_hfit.weights / maximum(iso_bdur_hfit.weights)
iso_bdur_edges = collect(iso_bdur_hfit.edges[1])[1:length(iso_bdur_weights)]

iso_ibi_hfit = fit(Histogram, isolated_data["IBIs"], LinRange(0.0, 120000, 100))
iso_ibi_weights = iso_ibi_hfit.weights / maximum(iso_ibi_hfit.weights)
iso_ibi_edges = collect(iso_ibi_hfit.edges[1])[1:length(iso_ibi_weights)]


#%% Open gradient data