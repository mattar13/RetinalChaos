using RetinalChaos
using Plots, Colors
using LaTeXStrings
using Plots.Measures
import RetinalChaos: M_INF, N_INF
#Define plotting attributes
font_title = Plots.font("Arial", 24)
font_axis = Plots.font("Arial", 12)
font_legend = Plots.font("Arial", 8)
pyplot(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
homedir = pwd()
pars_path = joinpath(homedir, "Settings")
figs_path = joinpath(homedir, "Examples")
#Set up default graph colors
v_color = :deepskyblue
n_color = :magenta
c_color = :green
a_color = :purple
b_color = :red
ach_color = :blue

#%% Setting up basic Settings
#Read JSON files for initial conditions and parameters
u0 = read_JSON(joinpath(pars_path,"conds.json")) |> extract_dict;
#Read JSON files for parameters and set the applied current to 10pA
p_dict = read_JSON(joinpath(pars_path, "params.json")) 
p_dict[:I_app] = 10.0
p = p_dict |> extract_dict
#Set the time span from 0s -> 60s
tspan = (0.0, 60e3)
#Create the problem
prob = ODEProblem(T_ode, u0, tspan, p)
println("Time it took to simulate 60s:")
@time sol = solve(prob); 

#%% Run all code here to generate Figure 1
begin
#Making Figure part A
xlims = (23.8e3, 24e3)
ellapsed_time = (xlims[end]-xlims[1])
xticks = (collect(xlims[1]:20.0:xlims[2]), collect(0:20.0:ellapsed_time))
###### Figure 1 A
fig1_Aa = plot(sol, vars = [:v, :n,], c = [v_color n_color], lw = 2.0, 
    ylabel = ["\$v_t \$(mV)" "\$n_t\$"], xlabel = ["" "time (ms)"], xlims = xlims,
    legend = :none, layout = grid(2, 1), 
)

vt_sol = map(t -> sol(t)[1], collect(xlims[1]:0.1:xlims[2]))
nt_sol = map(t -> sol(t)[2], collect(xlims[1]:0.1:xlims[2]))
fig1_Ab = plot(vt_sol, nt_sol, c = :black, label = "", 
    xlabel = "\$v_t\$ (mV)", ylabel = "\$n_t\$",
    title = "\$v_t\$ to \$n_t\$ Oscillation Diagram", titlefontsize = 12.0,
)

fig1_A = plot(fig1_Aa, fig1_Ab, layout = grid(1,2, widths = [0.7, 0.3]), grid = false);
title!(fig1_A[1], "A", title_location = :left)
###### Figure 1 B
v_rng = collect(-100:1.0:0)
C0 = p_dict[:C_0]; Ci = C0
δ = p_dict[:δ]; λ = p_dict[:λ]
#for calculating of I_ca
gCa = p_dict[:g_Ca]; ECa = p_dict[:E_Ca]
V1 = p_dict[:V1]; V2 = p_dict[:V2]
τc = p_dict[:τc]
f(v) = (C0+δ*(-gCa * M_INF(v, V1, V2) * (v - ECa)) - (λ*Ci))/τc
dc_rng = map(f, v_rng);


xlims = (22e3, 26e3) 
xticks = (collect(xlims[1]:500:xlims[2]), collect(0:0.5:(xlims[1]-xlims[2]/1000)))
fig1_Ba = plot(sol, vars = [:v, :c], c = [v_color c_color], lw = 2.0,
    ylabel = ["\$V_t\$ (mV)" "\$C_t\$ (mM)"], xlabel = ["" "time (s)"], 
    legend = :none, 
    layout = grid(2, 1), 
    xlims = xlims, 
    #xticks = collect(0:0.5:(xlims[2]-xlims[1]/1000))
)
fig1_Bb = plot(v_rng, dc_rng*1000, 
    c = c_color, linestyle = :dash, lw = 2.0, fill = (0, 0.5, c_color),
    title = "Instantaneous Calcium Change", titlefontsize = 12.0,
    xlabel = "\$V_t\$ (mV)",ylabel = "Inst. \$C_t\$ (mM/s)", label = "")

fig1_B = plot(fig1_Ba, fig1_Bb, 
    layout = grid(1, 2, widths=[0.7, 0.3]), size = (700, 350), grid = false)
title!(fig1_B[1], "B", title_pos = :left)

###### Figure 1 C
c_rng = collect(0.0:0.01:1.5)
ai = 0.0
α = p_dict[:α]
τa = p_dict[:τa]
fa(c) = ((α*c^4)/τa) - (ai/τa)
a_rng = map(fa, c_rng).*1000;

a_rng_2 = collect(0.0:0.01:1.5)
bi = 0.0
β = p_dict[:β]
τb = p_dict[:τb]
fb(a) = ((β*a^4)/τb) - (bi/τb)
b_rng = map(fb, a_rng_2).*1000;

xlims = (0.0, 40e3) 
xticks = (collect(xlims[1]:5e3:xlims[end]), collect(0:5.0:(xlims[end]-xlims[1]/1000)))

fig1_Ca = plot(
    sol, vars = [:c, :a, :b, :v], layout = grid(4,1), grid = false,
    ylabel = ["\$C_t\$ (mM)" "\$A_t\$ (mM)" "\$B_t\$ (mM)" "\$V_t\$ (mM)" ], xlabel = "",
    c = [c_color a_color b_color v_color], lw = 2.0, 
    legend = :none, 
    xlims = xlims, 
    #xticks = xticks, 
    ylims = [(0.0, 1.0) (0.0, 1.0) (0.0, 1.0) (-90.0, 0.0)],
)
fig1_Cb = plot(xlims = (0.0, 0.6), ylims = (0.0, 1.0), grid = false)
plot!(fig1_Cb, c_rng, a_rng, label = "", c = a_color, lw = 2.0, linestyle = :dash, 
    xlabel = "\$C_t\$ (mM)", ylabel = "Inst. \$A_t\$", 
    title = "\$C_t\$ and \$A_t\$", titlefontsize = 10.0, fill = (0, 0.2, a_color)
)
fig1_Cc = plot(xlims = (0.0, 0.6), ylims = (0.0, 1.0), grid = false)
plot!( fig1_Cc,
    a_rng_2, b_rng, label = "", 
    xlabel = "\$A_t\$", ylabel = "Inst. \$B_t\$",
    title = "\$A_t\$ and \$B_t\$", titlefontsize = 10.0,
    c = b_color, lw = 2.0, linestyle = :dash, 
    fill = (0, 0.2, b_color)
)
fig1_C = plot(fig1_Ca, fig1_Cb, fig1_Cc, layout = grid(1,3, widths = [0.6, 0.2, 0.2]), grid = false)
title!(fig1_C[1], "C", titlepos = :left);
###### Putting all parts together
fig1 = plot(fig1_A, fig1_B, fig1_C, layout = grid(3, 1, heights = [0.2, 0.3, 0.5]), size = (1000, 1000))
###### Saving
savefig(fig1, joinpath(figs_path, "Figure1_ModelDynamics.png"))