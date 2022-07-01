#Scientific reports figure sizes is 
#89mm (single column) to 183mm (double column) 3.5 to 7.5 inches

#Run this to set up the default parameters for plotting the figures
using Plots
using Plots.Measures
import PyCall as py
import PyCall.@py_str
import PyCall: @pyimport
import PyPlot as plt
import PyPlot: matplotlib, xlim, ylim, xlabel, ylabel
#Pyplot plot will be plt.plot
#normal plots will be plot

Plots.pyplot() #Switch the backend to pyplot
plt.pygui(true) #Make the GUI external to vscode
@pyimport matplotlib.gridspec as gspec #add the gridspec interface
@pyimport matplotlib.ticker as TICK #add the ticker interface
MultipleLocator = TICK.MultipleLocator #This is for formatting normal axis
LogLocator = TICK.LogLocator #This is for formatting the log axis]

#This functions allows us to get ranges of pyobjects in a Pythonic way
#In order to get slices we can wrap "getindex" 
import Base.getindex


slice(i, j) = py.pycall(py.pybuiltin("slice"), i, j)
#==============================These are the default parameters for plotting==============================#

print("Default plotting parameters loading... ")
rcParams = py.PyDict(matplotlib["rcParams"])

#Settings related the the DPI of the current plot interface
rcParams["figure.dpi"] = 60 #This is used to display the plot
rcParams["savefig.dpi"] = 600 #This is used to save the plot

rcParams["font.size"] = 18.0 #This controls the default font size
rcParams["font.family"] = "arial" #This controls the font family
rcParams["axes.spines.right"] = false #Make spines to the right invisible
rcParams["axes.spines.top"] = false #Make spines at the top invisible
rcParams["axes.linewidth"] = 1.0 #Make the spine width 1.0pts
rcParams["lines.linewidth"] = 0.7 #Set the lines.linewidth to 0.7pts
rcParams["xtick.major.size"] = 2.0 #Set the major x ticksize 2.0pts
rcParams["xtick.minor.size"] = 1.5
rcParams["xtick.major.pad"] = 2.0
rcParams["xtick.minor.pad"] = 2.0

rcParams["ytick.major.size"] = 2.0
rcParams["ytick.minor.size"] = 1.5
rcParams["ytick.major.pad"] = 2.0
rcParams["ytick.minor.pad"] = 2.0

rcParams["legend.frameon"] = false #This turns the legend frame off
rcParams["legend.labelspacing"] = 0.25 #this changes the spacing of the labels
rcParams["legend.borderpad"] = 0.2 #This changes the padding of the legend
#This changes the between the handle and label
rcParams["legend.handletextpad"] = -0.2
rcParams["legend.borderaxespad"] = 0.1
rcParams["legend.loc"] = "upper left"
rcParams["errorbar.capsize"] = 1.0 #Set the length of the errorbar cap
#set the background color for 
rcParams["figure.facecolor"] = (0.0, 0.0, 0.0, 0.0) #Make the figure background transparent white
rcParams["axes.facecolor"] = (0.0, 0.0, 0.0, 0.0) #Make the axes background transparent white

#These are the savefig params
rcParams["savefig.pad_inches"] = 0.0
println(" Completed")
