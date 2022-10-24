#We need to make a new model for this
using Revise
using RetinalChaos
include("../figures/figure_setup.jl")
using Plots
#Step 1: Open the physiological patch-clamp data.
using ePhys
file_loc = "C:/Users/mtarc/OneDrive - The University of Akron/Data/Patching"
target_file = "$(file_loc)/2019_11_03_Patch/Animal_2/Cell_3/19n03042.abf"
data_raw = readABF(target_file, channels=["Vm_prime4", "Im_sec5"], stimulus_name="Im_sec5", time_unit=:ms)
#Parse through data for acceptable files#
data = deepcopy(data_raw)
ePhys.truncate_data!(data, t_begin=1.0, t_end=300e3)
data - 10.0

tmin = 0.0
dt = data.dt
tmax = data.t[end]
tstops = tmin:dt:tmax

# Run the model
import RetinalChaos.ODEModel #import the ODEModel
import RetinalChaos.SDEModel #import the ODEModel
import RetinalChaos.u0 #import the 
import RetinalChaos.parameters
reload_parameters()

parameters[I_app] = 0.0
parameters[g_Na] = 1.0
#Step 3: determine the timespan
#Step 4: Warmup the problem
warmupProb = SDEProblem(SDEModel, u0, (tmin, 300e3), parameters)
warmup = solve(warmupProb, callback=cbs, EM(), dt=dt, save_everystep=false)
probSDE = SDEProblem(SDEModel, warmup.u[end], (tmin, tmax), parameters)

I_PHYS = data.data_array[:, :, 2] #Enter the external current as the input to the model
affect!(integrator) = integrator.p[1] = I_PHYS[round(Int64, (integrator.t / dt) + 1)]
cbs = PresetTimeCallback(tstops, affect!)


#Step 5: Solve the problem
@time sol = solve(probSDE, callback=cbs, EM(), dt=dt, tstops=tstops); #So far the best method is SOSRI

p = plot(sol, idxs=[v, I_ext], layout=(2, 1))
plot!(p, data)

#%% Calculate the loss of the functions
tsPHYS, dataPHYS = timeseries_analysis(data)
tsSOL, dataSOL = timeseries_analysis(sol)
dataSOL

dataPHYS["SpikeDurAvg"]
dataPHYS["BurstDurAvg"]
dataPHYS["IBIAvg"]

dataSOL["SpikeDurAvg"]
dataSOL["BurstDurAvg"]
dataSOL["IBIAvg"]

spikeLOSS, burstLOSS, IBILOSS = TimescaleLoss(data, sol)