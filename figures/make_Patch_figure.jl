#Load the example trace I collected
using Revise
using ABFReader
include("../../FunctionLibrary/FunctionLibrary.jl") #these are all of the functions that I have used to complete small tasks

target_file = raw"C:\Users\mtarc\OneDrive - The University of Akron\Data\Patching\2019_11_03_Patch\Animal_2\Cell_3\19n03042.abf"
data = readABF(target_file, channels=["Vm_prime4"], stimulus_name=nothing)

#%% Plotting data
width_inches = 6.0
height_inches = 6.0
fig1 = plt.figure("Good Trace", figsize=(width_inches, height_inches))
fig1.add_subplot(1, 1)
plot_experiment(fig1, data)