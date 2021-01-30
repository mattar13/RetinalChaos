#%% Running and analyzing the model using RetinalChaos.jl
using RetinalChaos

#Setup the fonts and stuff
font_title = font("Arial", 24)
font_axis = font("Arial", 12)
font_legend = font("Arial", 8)
pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)

#%% Making Figure 1
