#Scientific reports figure sizes is 
#For guidance, Nature's standard figure sizes are 89 mm wide (single column) and 183 mm wide (double column). 
#The full depth of a Nature page is 247 mm. 
#Figures can also be a column-and-a-half where necessary (120–136 mm).
cm = 1 / 2.54

#Authors should check (using a reducing photocopier) that, 
#at the smallest possible size, lettering remains readable and lines are sufficiently (but not too) heavy to print 
#clearly. Line weights and strokes should be set between 0.25 and 1 pt at the final size 
#(lines thinner than 0.25 pt may vanish in print). Do not rasterize or outline these lines if possible.
lw_standard = 1.0
#Run this to set up the default parameters for plotting the figures
using Plots
using Plots.Measures
import PyCall as py
import PyCall.@py_str
import PyCall: @pyimport
import PyPlot as plt
import PyPlot: matplotlib, xlim, ylim, xlabel, ylabel, title
#Pyplot plot will be plt.plot
#normal plots will be plot

Plots.pyplot() #Switch the backend to pyplot
plt.pygui(true) #Make the GUI external to vscode
@pyimport matplotlib.colors as COLOR
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

rcParams["font.size"] = 12.0 #This controls the default font size
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
#rcParams["figure.facecolor"] = (0.0, 0.0, 0.0, 0.0) #Make the figure background transparent white
#rcParams["axes.facecolor"] = (0.0, 0.0, 0.0, 0.0) #Make the axes background transparent white

#These are the savefig params
rcParams["savefig.pad_inches"] = 0.0
println(" Completed")

function plot_histograms(data, loc::String; name="histogram_plot")

    if !isempty(data["SpikeDurs"])
        sdur_hfit = fit(Histogram, data["SpikeDurs"], LinRange(0.0, 50.0, 100))
        sdur_weights = sdur_hfit.weights / maximum(sdur_hfit.weights)
        sdur_edges = collect(sdur_hfit.edges[1])[1:length(sdur_weights)]
    else
        sdur_edges = sdur_weights = [0]
    end
    p1 = plot(sdur_edges, sdur_weights, xlabel="Spike Duration (ms)")

    if !isempty(data["ISIs"])
        isi_hfit = fit(Histogram, data["ISIs"], LinRange(0.0, 100.0, 100))
        isi_weights = isi_hfit.weights / maximum(isi_hfit.weights)
        isi_edges = collect(isi_hfit.edges[1])[1:length(isi_weights)]
    else
        isi_edges = isi_weights = [0]
    end
    p2 = plot(isi_edges, isi_weights, xlabel="Spike Interval (s)", xformatter=x -> x / 1000)

    if !isempty(data["BurstDurs"])
        bdur_hfit = fit(Histogram, data["BurstDurs"], LinRange(0.0, 2000.0, 100))
        bdur_weights = bdur_hfit.weights / maximum(bdur_hfit.weights)
        bdur_edges = collect(bdur_hfit.edges[1])[1:length(bdur_weights)]
    else
        bdur_edges = bdur_weights = [0]
    end
    p3 = plot(bdur_edges, bdur_weights, xlabel="Burst Duration (s)", xformatter=x -> x / 1000)

    if !isempty(data["IBIs"])
        ibi_hfit = fit(Histogram, data["IBIs"], LinRange(0.0, 120e3, 100))
        ibi_weights = ibi_hfit.weights / maximum(ibi_hfit.weights)
        ibi_edges = collect(ibi_hfit.edges[1])[1:length(ibi_weights)]
    else
        ibi_edges = ibi_weights = [0]
    end
    p4 = plot(ibi_edges, ibi_weights, xlabel="Interburst Interval (s)", xformatter=x -> x / 1000)

    if !isempty(data["Thresholds"])
        p5 = histogram(data["Thresholds"], yaxis=:log, xlabel="Voltage threshold")
    else
        p5 = plot(xlabel="Voltage threshold")
    end

    if !isempty(data["SpikesPerBurst"])
        p6 = histogram(data["SpikesPerBurst"], yaxis=:log, xlabel="Spikes per Burst")
    else
        p6 = plot(xlabel="Spike Per Burst")
    end
    hist_plot = plot(
        p1, p2, p3, p4, p5, p6,
        layout=grid(3, 2), ylabel="Counts",
        legend=false
    )
    savefig(hist_plot, "$(loc)/$(name).png")
    return hist_plot
end

function add_direction(ax, data_x, data_y;
    start_rng::Int64=1, end_rng::Union{Int64,Nothing}=nothing,
    color=:black,
    delta::Int64=100, nArrows::Int64=50,
    width=0.02, headwidth=7.0, headaxislength=5.0
) where {T}
    if isnothing(end_rng)
        rng = round.(Int64, LinRange(start_rng, length(data_x) - delta, nArrows))
    else
        rng = round.(Int64, LinRange(start_rng, end_rng, nArrows))
    end
    for i in rng
        x = data_x[i]
        y = data_y[i]
        dx = data_x[i+delta] - data_x[i]
        dy = data_y[i+delta] - data_y[i]
        #axAR.plot([x, x+dx], [y, y+dy], marker = "o", lw = 0.0, c = :blue)
        ax.quiver(x, y, dx, dy, color=color, angles="xy", pivot="mid",
            width=width, headwidth=headwidth, headaxislength=headaxislength)
    end
end
#%% If we want to run each script by itself use this
#include("making_figure1.jl")
#include("making_figure2.jl")
#include("making_figure3.jl")
