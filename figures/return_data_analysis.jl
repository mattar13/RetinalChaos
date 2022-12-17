#Filter out any spike duration under 5.0
phys_spike_durs = phys_spike_durs[phys_spike_durs.>5.0]

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


#%% Print the results of all the traces
println("Data Values")
println("-----------------------------------------------")
println("""

ISO_n = $(length(iso_baselines))
ISO Baseline Avg: $(iso_baseline_avg|>x->round(x, digits = 2))±$(iso_baseline_SEM|>x->round(x, digits = 2))
ISO Minimum Avg: $(iso_min_avg|>x->round(x, digits = 2))±$(iso_min_SEM|>x->round(x, digits = 2))
ISO Maximum Avg: $(iso_max_avg|>x->round(x, digits = 2))±$(iso_max_SEM|>x->round(x, digits = 2))
ISO Spike Duration Avg: $(dataISO["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataISO["SpikeDurSEM"]|>x->round(x, digits = 2))
ISO Burst Duration Avg: $(dataISO["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataISO["BurstDurSEM"]|>x->round(x, digits = 2))
ISO IBI_Avg: $(dataISO["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataISO["IBISEM"]./1000|>x->round(x, digits = 2))

NO GABA n = $(length(ng_baselines))
NO GABA Baseline Avg: $(ng_baseline_avg|>x->round(x, digits = 2))±$(ng_baseline_SEM|>x->round(x, digits = 2))
NO GABA Minimum Avg: $(ng_min_avg|>x->round(x, digits = 2))±$(ng_min_SEM|>x->round(x, digits = 2))
NO GABA Maximum Avg: $(ng_max_avg|>x->round(x, digits = 2))±$(ng_max_SEM|>x->round(x, digits = 2))
NO GABA Spike Duration Avg: $(dataNG["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataNG["SpikeDurSEM"]|>x->round(x, digits = 2))
NO GABA Burst Duration Avg: $(dataNG["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataNG["BurstDurSEM"]|>x->round(x, digits = 2))
NO GABA IBI_Avg: $(dataNG["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataNG["IBISEM"]./1000|>x->round(x, digits = 2))

ECl -55 n = $(length(ec_baselines))
ECl -55 Baseline Avg: $(ec_baseline_avg|>x->round(x, digits = 2))±$(ec_baseline_SEM|>x->round(x, digits = 2))
ECl -55 Minimum Avg: $(ec_min_avg|>x->round(x, digits = 2))±$(ec_min_SEM|>x->round(x, digits = 2))
ECl -55 Maximum Avg: $(ec_max_avg|>x->round(x, digits = 2))±$(ec_max_SEM|>x->round(x, digits = 2))
ECl -55 Spike Duration Avg: $(dataEC["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataEC["SpikeDurSEM"]|>x->round(x, digits = 2))
ECl -55 Burst Duration Avg: $(dataEC["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataEC["BurstDurSEM"]|>x->round(x, digits = 2))
ECl -55 IBI_Avg: $(dataEC["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataEC["IBISEM"]./1000|>x->round(x, digits = 2))

WAVE n = $(length(wave_baselines))
WAVE Baseline Avg: $(wave_baseline_avg|>x->round(x, digits = 2))±$(wave_baseline_SEM|>x->round(x, digits = 2))
WAVE Minimum Avg: $(wave_min_avg|>x->round(x, digits = 2))±$(wave_min_SEM|>x->round(x, digits = 2))
WAVE Maximum Avg: $(wave_max_avg|>x->round(x, digits = 2))±$(wave_max_SEM|>x->round(x, digits = 2))
WAVE Spike Duration Avg: $(dataWAVE["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataWAVE["SpikeDurSEM"]|>x->round(x, digits = 2))
WAVE Burst Duration Avg: $(dataWAVE["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataWAVE["BurstDurSEM"]|>x->round(x, digits = 2))
WAVE IBI_Avg: $(dataWAVE["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataWAVE["IBISEM"]./1000|>x->round(x, digits = 2))

PHYS n = $(length(phys_baseline))
PHYS Baseline Avg: $(phys_baseline_avg|>x->round(x, digits = 2))±$(phys_baseline_SEM|>x->round(x, digits = 2))
PHYS Minimum Avg: $(phys_min_avg|>x->round(x, digits = 2))±$(phys_min_SEM|>x->round(x, digits = 2))
PHYS Maximum Avg: $(phys_max_avg|>x->round(x, digits = 2))±$(phys_max_SEM|>x->round(x, digits = 2))
PHYS Spike Duration Avg: $(phys_spike_avg|>x->round(x, digits = 2))±$(phys_spike_SEM |>x->round(x, digits = 2))
PHYS Burst Duration Avg: $(phys_burst_avg|>x->round(x, digits = 2))±$(phys_burst_SEM|>x->round(x, digits = 2))
PHYS IBI_Avg: $(phys_ibis_avg./1000|>x->round(x, digits = 2))±$(phys_ibis_SEM./1000|>x->round(x, digits = 2))
""")