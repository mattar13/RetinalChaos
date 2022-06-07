using Revise
using ABFReader
using StatsBase, Statistics
#%% Run the data on the physiological data that Jordan has collected
target_folder = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Patching\Jordans_Patch_Data\UsuableData"
paths = target_folder |> parse_abf;

phys_baseline = Float64[]
phys_max = Float64[]
phys_min = Float64[]
phys_spike_durs = Float64[]
phys_isis = Float64[]
phys_burst_durs = Float64[]
phys_ibis = Float64[]
for path in paths
     print("Opening $path ...")
     data_i = readABF(path,
          channels=["Im_scaled"],
          stimulus_name=nothing, flatten_episodic=true
     )
     spike_timestampsi = ABFReader.get_timestamps(data_i)
     burst_timestampsi = ABFReader.max_interval_algorithim(spike_timestampsi)
     println(burst_timestampsi)
     spike_durs, isi = ABFReader.extract_interval(spike_timestampsi)
     burst_durs, ibi = ABFReader.extract_interval(burst_timestampsi[1])
     println(burst_durs)
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

n_recordings = length(paths)
#Average the stats
phys_baseline_avg = sum(phys_baseline) / n_recordings
phys_baseline_sem = std(phys_baseline) / sqrt(n_recordings)

phys_min_avg = sum(phys_min) / n_recordings
phys_min_sem = std(phys_min) / sqrt(n_recordings)

phys_max_avg = sum(phys_max) / n_recordings
phys_max_sem = std(phys_max) / sqrt(n_recordings)

phys_spike_avg = sum(phys_spike_durs) / length(phys_spike_durs) * 1000
phys_spike_sem = std(phys_spike_durs) / sqrt(length(phys_spike_durs)) * 1000

phys_burst_avg = sum(phys_burst_durs) / length(phys_burst_durs)
phys_burst_sem = std(phys_burst_durs) / sqrt(length(phys_burst_durs))

phys_ibis_avg = sum(phys_ibis) / length(phys_ibis)
phys_ibis_sem = std(phys_ibis) / sqrt(length(phys_ibis))