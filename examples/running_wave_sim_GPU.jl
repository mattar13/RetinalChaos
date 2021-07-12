using Revise
using RetinalChaos
using Dates, Plots, JLD2

#Configure the logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#Load the telegram client and env
dotenv("D:\\TelegramAccessEnv\\.env")
#Activate the GPU
RetinalChaos.CUDA.allowscalar(false)

#%% Load the needed files to run the model
save_file = "C:\\Users\\RennaLabSA1\\Documents\\ModellingData\\$(Date(Dates.now()))_test1\\"
p_dict = read_JSON("params\\params.json", is_type = Dict{Symbol, Float32}) 
u_dict = read_JSON("params\\conds.json", is_type = Dict{Symbol, Float32})
NetSol = load_model(save_file, p_dict, u_dict)
#%%
timestamps, data = timeseries_analysis(NetSol)

#%%Conduct the threshold analysis
thresholds = calculate_threshold(NetSol; dt = 500.0) #This takes really long
spikes = get_timestamps(NetSol, thresholds)
spike_durs, isi = extract_interval(spikes)
bursts, spb = max_interval_algorithim(spikes)
burst_durs, ibi = extract_interval(bursts)

JLD2.save("$(save_file)\\timestamps.jld2", 
    Dict(
        "Spikes" => spikes,
        "Bursts" => bursts
    )
)

JLD2.save("$(save_file)\\data.jld2", 
    Dict(
        "Thresholds" => thresholds,
        "SpikeDurs" => spike_durs, 
        "ISIs" => isi, 
        "BurstDurs" => burst_durs, 
        "IBIs" => ibi,
        "SpikesPerBurst" => spb
    )
)
BotNotify("{Wave} Finished running timeseries analysis")

#%%
timestamps = JLD2.load("$(save_file)\\timestamps.jld2")
data = JLD2.load("$(save_file)\\data.jld2")
#%%
run_loop = false #so far set this section to ignore
if run_loop
    #%% Run multiple samples
    n = 4
    for val in LinRange(1.0, 20.0, 10)
       
    end
end