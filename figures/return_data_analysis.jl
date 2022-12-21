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
println("Data Values")
println("-----------------------------------------------")

println("""
     PHYS n = $(length(phys_baseline))
     PHYS Baseline Avg: $(phys_baseline_avg|>x->round(x, digits = 2))±$(phys_baseline_SEM|>x->round(x, digits = 2))
     PHYS Minimum Avg: $(phys_min_avg|>x->round(x, digits = 2))±$(phys_min_SEM|>x->round(x, digits = 2))
     PHYS Maximum Avg: $(phys_max_avg|>x->round(x, digits = 2))±$(phys_max_SEM|>x->round(x, digits = 2))
     PHYS Spike Duration Avg: $(phys_spike_avg|>x->round(x, digits = 2))±$(phys_spike_SEM |>x->round(x, digits = 2))
     PHYS Burst Duration Avg: $(phys_burst_avg|>x->round(x, digits = 2))±$(phys_burst_SEM|>x->round(x, digits = 2))
     PHYS IBI_Avg: $(phys_ibis_avg./1000|>x->round(x, digits = 2))±$(phys_ibis_SEM./1000|>x->round(x, digits = 2))

""")

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


println("""
ISO_n = $(length(iso_baselines))
ISO Baseline Avg: $(iso_baseline_avg|>x->round(x, digits = 2))±$(iso_baseline_SEM|>x->round(x, digits = 2))
ISO Minimum Avg: $(iso_min_avg|>x->round(x, digits = 2))±$(iso_min_SEM|>x->round(x, digits = 2))
ISO Maximum Avg: $(iso_max_avg|>x->round(x, digits = 2))±$(iso_max_SEM|>x->round(x, digits = 2))
ISO Spike Duration Avg: $(dataISO["SpikeDurAvg"]|>x->round(x, digits = 2)) ± $(dataISO["SpikeDurSEM"]|>x->round(x, digits = 2))
ISO Burst Duration Avg: $(dataISO["BurstDurAvg"]|>x->round(x, digits = 2)) ± $(dataISO["BurstDurSEM"]|>x->round(x, digits = 2))
ISO IBI_Avg: $(dataISO["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataISO["IBISEM"]./1000|>x->round(x, digits = 2))
""")

wavesISO = DataFrame(XLSX.readtable("$(isolated_path)\\wave_data.xlsx", "Waves"))
eventsISO = DataFrame(XLSX.readtable("$(isolated_path)\\wave_data.xlsx", "WaveStats"))


println("""
ISO Wave/Burst ratio = $(round(sum(eventsISO.IsWave)/length(eventsISO.IsWave), digits = 3))
ISO Wave Duration = $(round(mean(wavesISO.WaveTime)/1000, digits = 2)) ± $(round(std(wavesISO.WaveTime)/sqrt(length(wavesISO.WaveTime))/1000, digits = 2))
ISO Wave Area = $(round(mean(wavesISO.WaveDistances), digits = 2)) ± $(round(std(wavesISO.WaveDistances)/sqrt(length(wavesISO.WaveTime)), digits = 2))
ISO Wave Velocity = $(round(mean(wavesISO.WaveVelocities)/1000, digits = 2)) ± $(round(std(wavesISO.WaveVelocities)/sqrt(length(wavesISO.WaveTime))/1000, digits = 2))
""")

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

println("""
     NO GABA n = $(length(ng_baselines))
     NO GABA Baseline Avg: $(ng_baseline_avg|>x->round(x, digits = 2))±$(ng_baseline_SEM|>x->round(x, digits = 2))
     NO GABA Minimum Avg: $(ng_min_avg|>x->round(x, digits = 2))±$(ng_min_SEM|>x->round(x, digits = 2))
     NO GABA Maximum Avg: $(ng_max_avg|>x->round(x, digits = 2))±$(ng_max_SEM|>x->round(x, digits = 2))
     NO GABA Spike Duration Avg: $(dataNG["SpikeDurAvg"]|>x->round(x, digits = 2)) ± $(dataNG["SpikeDurSEM"]|>x->round(x, digits = 2))
     NO GABA Burst Duration Avg: $(dataNG["BurstDurAvg"]|>x->round(x, digits = 2)) ± $(dataNG["BurstDurSEM"]|>x->round(x, digits = 2))
     NO GABA IBI_Avg: $(dataNG["IBIAvg"]./1000|>x->round(x, digits = 2)) ± $(dataNG["IBISEM"]./1000|>x->round(x, digits = 2))
""")

wavesNG = DataFrame(XLSX.readtable("$(noGABA_path)\\wave_data.xlsx", "Waves"))
eventsNG = DataFrame(XLSX.readtable("$(noGABA_path)\\wave_data.xlsx", "WaveStats"))

println("""
NO GABA Wave/Burst ratio = $(round(sum(eventsNG.IsWave)/length(eventsNG.IsWave), digits = 3))
NO GABA Wave Duration = $(round(mean(wavesNG.WaveTime)/1000, digits = 2)) ± $(round(std(wavesNG.WaveTime)/sqrt(length(wavesNG.WaveTime))/1000, digits = 2))
NO GABA Wave Area = $(round(mean(wavesNG.WaveDistances), digits = 2)) ± $(round(std(wavesNG.WaveDistances)/sqrt(length(wavesNG.WaveTime)), digits = 2))
NO GABA Wave Velocity = $(round(mean(wavesNG.WaveVelocities)/1000, digits = 2)) ± $(round(std(wavesNG.WaveVelocities)/sqrt(length(wavesNG.WaveTime))/1000, digits = 2))
""")

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

println("""
     ECl -55 n = $(length(ec_baselines))
     ECl -55 Baseline Avg: $(ec_baseline_avg|>x->round(x, digits = 2))±$(ec_baseline_SEM|>x->round(x, digits = 2))
     ECl -55 Minimum Avg: $(ec_min_avg|>x->round(x, digits = 2))±$(ec_min_SEM|>x->round(x, digits = 2))
     ECl -55 Maximum Avg: $(ec_max_avg|>x->round(x, digits = 2))±$(ec_max_SEM|>x->round(x, digits = 2))
     ECl -55 Spike Duration Avg: $(dataEC["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataEC["SpikeDurSEM"]|>x->round(x, digits = 2))
     ECl -55 Burst Duration Avg: $(dataEC["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataEC["BurstDurSEM"]|>x->round(x, digits = 2))
     ECl -55 IBI_Avg: $(dataEC["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataEC["IBISEM"]./1000|>x->round(x, digits = 2))
""")

wavesEC = DataFrame(XLSX.readtable("$(ECl55_path)\\wave_data.xlsx", "Waves"))
eventsEC = DataFrame(XLSX.readtable("$(ECl55_path)\\wave_data.xlsx", "WaveStats"))

println("""
ECl -55 Wave/Burst ratio = $(round(sum(eventsEC.IsWave)/length(eventsEC.IsWave), digits = 3))
ECl -55 Wave Duration = $(round(mean(wavesEC.WaveTime)/1000, digits = 2)) ± $(round(std(wavesEC.WaveTime)/sqrt(length(wavesEC.WaveTime))/1000, digits = 2))
ECl -55 Wave Area = $(round(mean(wavesEC.WaveDistances), digits = 2)) ± $(round(std(wavesEC.WaveDistances)/sqrt(length(wavesEC.WaveTime)), digits = 2))
ECl -55 Wave Velocity = $(round(mean(wavesEC.WaveVelocities)/1000, digits = 2)) ± $(round(std(wavesEC.WaveVelocities)/sqrt(length(wavesEC.WaveTime))/1000, digits = 2))
""")

wave_baselines = calculate_threshold(dataWAVE["DataArray"], Z=1.0, dims=2)
wave_baseline_avg = sum(wave_baselines) / length(wave_baselines)
wave_baseline_SEM = std(wave_baselines) / sqrt(length(wave_baselines))

wave_min = minimum(dataWAVE["DataArray"], dims=2)
wave_min_avg = sum(wave_min) / length(wave_min)
wave_min_SEM = std(wave_min) / sqrt(length(wave_min))

wave_max = maximum(dataWAVE["DataArray"], dims=2)
wave_max_avg = sum(wave_max) / length(wave_max)
wave_max_SEM = std(wave_max) / sqrt(length(wave_max))

println("""
     ECl -65 n = $(length(wave_baselines))
     ECl -65 Baseline Avg: $(wave_baseline_avg|>x->round(x, digits = 2))±$(wave_baseline_SEM|>x->round(x, digits = 2))
     ECl -65 Minimum Avg: $(wave_min_avg|>x->round(x, digits = 2))±$(wave_min_SEM|>x->round(x, digits = 2))
     ECl -65 Maximum Avg: $(wave_max_avg|>x->round(x, digits = 2))±$(wave_max_SEM|>x->round(x, digits = 2))
     ECl -65 Spike Duration Avg: $(dataWAVE["SpikeDurAvg"]|>x->round(x, digits = 2))±$(dataWAVE["SpikeDurSEM"]|>x->round(x, digits = 2))
     ECl -65 Burst Duration Avg: $(dataWAVE["BurstDurAvg"]|>x->round(x, digits = 2))±$(dataWAVE["BurstDurSEM"]|>x->round(x, digits = 2))
     ECl -65 IBI_Avg: $(dataWAVE["IBIAvg"]./1000|>x->round(x, digits = 2))±$(dataWAVE["IBISEM"]./1000|>x->round(x, digits = 2))
""")


wavesDE = DataFrame(XLSX.readtable("$(wave_path)\\wave_data.xlsx", "Waves"))
eventsDE = DataFrame(XLSX.readtable("$(wave_path)\\wave_data.xlsx", "WaveStats"))

println("""
ECl -65 Wave/Burst ratio = $(round(sum(eventsDE.IsWave)/length(eventsDE.IsWave), digits = 3))
ECl -65 Wave Duration = $(round(mean(wavesDE.WaveTime)/1000, digits = 2)) ± $(round(std(wavesDE.WaveTime)/sqrt(length(wavesDE.WaveTime))/1000, digits = 2))
ECl -65 Wave Area = $(round(mean(wavesDE.WaveDistances), digits = 2)) ± $(round(std(wavesDE.WaveDistances)/sqrt(length(wavesDE.WaveTime)), digits = 2))
ECl -65 Wave Velocity = $(round(mean(wavesDE.WaveVelocities)/1000, digits = 2)) ± $(round(std(wavesDE.WaveVelocities)/sqrt(length(wavesDE.WaveTime))/1000, digits = 2))
""")

