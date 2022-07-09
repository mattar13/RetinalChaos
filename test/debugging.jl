using Revise
using RetinalChaos
import RetinalChaos: calculate_threshold, get
import RetinalChaos: extract_equilibria, find_equilibria
include("../figures/figure_setup.jl");

max_spike_duration = 50.0
max_spike_interval = 100.0
max_burst_duration = 10e5
max_burst_interval = 10e5

#Lets look at each trace and figure out if there are bursts or not
target_folder = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Patching\Jordans_Patch_Data\UsuableData"
paths = target_folder |> parseABF;

path = paths[1]
print("[$(now())] Opening $path... ")
data_i = readABF(path,
     channels=["Im_scaled"],
     stimulus_name=nothing, flatten_episodic=true
)
t = data_i.t*1000 #Turn seconds into milliseconds
vm_array = data_i.data_array[:, :, 1]
thresholds = RetinalChaos.calculate_threshold(vm_array, Z=2, dims=2)

println("Complete")
print("[$(now())]: Extracting the spikes... ")
spike_array = Matrix{Bool}(vm_array .> thresholds)
spikes = RetinalChaos.get_timestamps(spike_array, t)
thresholds

spike_durs, isi = RetinalChaos.extract_interval(spikes, max_duration=max_spike_duration, max_interval=max_spike_interval)
println("Complete")

#%%
print("[$(now())]: Extracting the bursts... ")
bursts, spb = RetinalChaos.max_interval_algorithim(spikes)
burst_durs, ibi = RetinalChaos.extract_interval(bursts[1], max_duration=max_burst_duration, max_interval=max_burst_interval, min_duration=1000.0)
println("Complete")