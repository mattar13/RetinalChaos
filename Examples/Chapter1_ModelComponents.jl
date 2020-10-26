### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 77bf8f30-f211-11ea-2516-3d678de83568
using PlutoUI

# ╔═╡ fe8ed2ce-f1e1-11ea-2c15-ad90425e6b01
using RetinalChaos

# ╔═╡ 6e849542-f2d5-11ea-0eba-632e96ea9ded
using DifferentialEquations

# ╔═╡ c52c2942-f207-11ea-33e2-1f597edebb9c
using Plots, Colors

# ╔═╡ 84194d80-fb8e-11ea-2d3f-a9a237c412ef
using Distributions, StatsPlots

# ╔═╡ d2bde4d0-f1e0-11ea-395a-21b3cc2800b1
md"# RetinalChaos

This is the tutorial for running the model in the paper ... Tarchick Et. al. In order to begin this tutorial the first thing we must do is import the package that this file came with."

# ╔═╡ 3826848e-f208-11ea-1b46-0d98f24a3da7
pyplot();

# ╔═╡ 55ce0450-f203-11ea-3d1b-d52ee88fffc3
import RetinalChaos: I_n, R_INF, Λ, Φ, ħ, ∇

# ╔═╡ d468e040-f1e1-11ea-152f-dff8d05c8594
md"
# **[1] Components of the Model**
- [1.1] $V_t$ Voltage
- [1.2] $N_t$ Potassium Gating
- [1.3] $C_t$ Calcium
- [1.4] $A_t$ cAMP breakdown)
- [1.5] $B_t$ TREK phosphorylation
- [1.6] $E_t$ Acetylcholine
- [1.7] $W_t$ Noise "

# ╔═╡ 0cdbf9e2-f1e1-11ea-3a06-a702d6cae6fb
md"## [1.1] Voltage ($V_t$)

The equation for the membrane voltage of a single cell is as follows. 

$$\begin{align}
 C_m\frac{dV}{dt} &= I_{leak}(V_t) + I_{Ca}(V_t) + I_K(V_t, N_t) + I_{TREK}(V_t, B_t) + I_{ACh}(V_t, E_t) + I_{app} + I_{noise}
\end{align}$$

Ohmic conductance relationships are used to describe the voltage. The cell membrane acts as both a capacitor and a conductor, so the voltage has to be described as a differential. These equations are similar the original equations used by the Morris-Lecar model (Morris1981). Sodium currents are not absent, however starburst cells are predominantly driven by calcium transients, so $I_{Na}$ is assumed to be low and does not have a significant effect on the  (Matzakos-karvouniari2019). The first term $C_m$ is the membrane capacitance. This is related to the total surface area of the cell and can easily be experimentally determined. 

### [1.1.a] Currents($I_n$) 

- We can use the function $I_n(g_n, R, V_t, E_n). 
To access more information about this function place your cursor in the next box and click Live Docs"

# ╔═╡ 6b1a70a0-f203-11ea-0026-b5faa986fecf
I_n

# ╔═╡ 6db280a0-f203-11ea-2cc1-0381b15e1410
md" ##### In the case of this model there are 5 currents. 
#### 1) Leaky current ($I_{leak}$): 
- This current represents a constituently active ionic channel. This channel is important for setting the baseline. Because this channel is always open, $R_{\infty} = 1.0$ 
- Leaky currents will repolarize the membrane if the membrane potential is over -70 and will depolarize the membrane if the membrane potential is under -70
    
$$\begin{align}
I_{leak} = -g_{leak} *(V_t - E_{leak})
\end{align}$$"

# ╔═╡ b82b71a0-f203-11ea-3dc9-93fc0632e420
md" 
($g_{leak}$): $(@bind g_leak1 NumberField(1:100, default = 2.0)) nS

($E_{leak}$) : $(@bind E_leak1 NumberField(-100:100, default = -70.0)) mV

We can show the instantaneous voltage for a range between -90 and 10
"

# ╔═╡ 9861a390-f207-11ea-3f45-6198c2f22003
begin
	f1 = v -> I_n(g_leak1, 1.0, v, E_leak1)
	plot(f1, -90, 10, label = "g:$(g_leak1)nS",lw = 2.0,
		xlabel = "Membrane Potential (mV)", xlims = (-90, 10),
		ylabel = "Leaky current (pA)", ylims = (-200, 200)
		)
end

# ╔═╡ c1e40650-f20b-11ea-34a8-e157d35f600e
md"### 2) Calcium Currents ($I_{Ca}$) 

$$I_{Ca}(V_t) = -g_{Ca} M_\infty(V_t, V_1, V_2) (V_t - E_{Ca})$$ 

- Reversal potentials are positive ($E_{Ca} = 50.0mV$) so channel is depolarizing.
- Voltage gated channels are gated by a modification of the $R_\infty$ function $M_\infty$
    - \$V_1$: is the slope
    - \$V_2$: as the half saturation value"

# ╔═╡ bda8ed50-f20e-11ea-2e23-075cedea76b7
R_INF

# ╔═╡ 0962f460-f210-11ea-21f1-0349a3e20635
md"""
First we will talk about the auxillary equation $R_{\infty}$

More information about $R_{\infty}$ exists in the Live Docs of the above function. This function uses 2 parameters as well as the variable $V_t$. The plotting example below shows the effect of varying 

V1: $(@bind V1p NumberField(-100:100, default = -20.0))

V2: $(@bind V2p NumberField(-100:100, default = 20.0))
"""

# ╔═╡ 8e05ba42-f210-11ea-0684-3972c2fffb56
begin
	f2 = v -> R_INF(v, V1p, V2p)
	plot(f2, -90.0, 10.0, 
		lw = 3.0, c = :green,
		xlabel = "Membrane Potential (mV)", xlims = (-90.0, 10.0), 
		ylabel = "Normalized Current", ylims = (0.0, 1.0)
	
	)
end

# ╔═╡ 7d839cb0-f20f-11ea-122d-cb56b26e2f77
md"($g_{Ca}$): $(@bind gCap NumberField(1:100, default = 10.0)) nS
	
($E_{Ca}$): $(@bind ECap NumberField(-100:100, default = 50.0)) mV"

# ╔═╡ 4df37fa2-f215-11ea-2a3f-f3c139376fad
begin 
	f3 = v -> I_n(gCap, f2(v), v, ECap)
	plot(f3, -90.0, 10.0,
		lw = 3.0, c = :blue,
		xlabel = "Membrane Potential (mV)", xlims = (-90.0, 10.0),
		ylabel = "Instant \$I_{Ca}\$ (pA)", ylims = (0.0, 800)
	)
end

# ╔═╡ 02fb4610-f2d1-11ea-0ea8-5d9a7cf72864
md"### The rest of the currents will be explained in their respective sections."

# ╔═╡ 9277f880-f2c0-11ea-3491-8f9f2ad9941d
md"
## [1.2] Fast Potassium Repolarization ($N_t$)

The next term describes the fast repolarization due to the opening of voltage gated potassium channels. This term is used to calculate the $I_K$ 

$\begin{align}
\tau_N\frac{dN}{dt} = \Lambda(V_t, V_3, V_4)(N_\infty(V_t, V_3, V_4) - N)
\end{align}$

- The decay time constant ($\tau_N$) represents the rate at which the Potassium gating decays. 
- The rate constant ($\Lambda$) is decsribed in the Live Docs below. To see the details of lambda put your cursor in the box below and click Live Docs. 
- The normalized activation ($N_\infty$) is just another refactorization of $R_\infty$

"

# ╔═╡ 75f593a0-f2cc-11ea-23d6-a71cc7db012a
md"
The two parameters for both $N_\infty$ and $\Lambda$ are:

V3: $(@bind V3p NumberField(0:100, default = -25.0))

V4: $(@bind V4p NumberField(0:100, default = 7.0))

and we can plot both $N_\infty$ and $\Lambda$
"

# ╔═╡ 0eb2c2d0-f23b-11ea-1b27-0541ac4d0b36
Λ

# ╔═╡ ec80e870-f2cd-11ea-1c3a-05beec6b0480
begin 
	NINF = v-> R_INF(v, V3p, V4p)
	LAM = v -> Λ(v, V3p, V4p)
	plt = plot(layout = grid(1,2))
	plot!(plt[1], NINF, -90.0, 10.0, 
		lw = 3.0, label = "\$N_\\infty\$", title = "\$N_\\infty\$",
		xlabel = "Membrane Potential (mV)", xlims = (-90.0, 10.0), 
		ylabel = "Normalized Current", ylims = (0.0, 1.0))
	plot!(plt[2], LAM, -90.0, 10.0, 
		lw = 3.0, label = "\$\\Lambda\$", title = "\$\\Lambda\$",
		xlabel = "Membrane Potential (mV)", xlims = (-90.0, 10.0), 
		ylabel = "Relaxation Time Constant", ylims = (0.0, 25.0)
	)
	plt
end

# ╔═╡ aace080e-f2d0-11ea-15cd-bb7c1926458f
md" 
##### In order to demonstrate how $N_t$ and $I_K$ work we have to run a small simulation
This is the first demonstration of actually running our model. Using the code below can act as a quick tutorial on how to run the model. Included in this section will be all the parameters mentioned so far as well as a few extras. All other parameters needed to run the model will be zeroed out. 
"

# ╔═╡ a79cda6e-f2d2-11ea-27d9-4df5215595fe
md" 
##### Parameters:

I_app: $(@bind Iapp NumberField(-100:100, default = 10.0)) pA

Cm: $(@bind C_m NumberField(-100:100, default = 22.0)) nF

($g_{leak}$): $(@bind g_leak NumberField(1:100, default = g_leak1)) nS

($E_{leak}$) : $(@bind E_leak NumberField(-100:10:100, default = E_leak1)) mV

V1: $(@bind V1 NumberField(-100:5:100, default = V1p))
V2: $(@bind V2 NumberField(-100:5:100, default = V2p))

($g_{Ca}$): $(@bind gCa NumberField(1:100, default = gCap)) nS

($E_{Ca}$): $(@bind ECa NumberField(-100:10:100, default = ECap)) mV

V3: $(@bind V3 NumberField(-100:5:100, default = V3p))
V4: $(@bind V4 NumberField(-100:5:100, default = V4p))

($g_K$): $(@bind gK NumberField(1:100, default = 10.0)) nS

($E_K$): $(@bind EK NumberField(-100:10:100, default = -90.0)) mV

#### Settings
tmax = $(@bind tmax1 NumberField(0:100:10e3, default = 300.0)) ms
"

# ╔═╡ 8510a620-f2d8-11ea-07bd-3146b9fd9fc5
begin
	#This block imports all the necessary conditions for the model. 
	p_dict1 = read_JSON("params.json"); 
	#Zero out all conductances to run this experiment
	p_dict1[:I_app] = Iapp
	p_dict1[:g_TREK] = 0.0; #Zero TREK params
	p_dict1[:g_ACh] = 0.0; #Zero Acetylcholine for now
	#Add in model parameters
	p_dict1[:C_m] = C_m
	p_dict1[:g_leak] = g_leak
	p_dict1[:E_leak] = E_leak
	p_dict1[:V1] = V1
	p_dict1[:V2] = V2
	p_dict1[:g_Ca] = gCa
	p_dict1[:E_Ca] = ECa
	p_dict1[:V3] = V3
	p_dict1[:V4] = V4
	p_dict1[:g_K] = gK
	p_dict1[:E_K] = EK
	#Import the initial conditions
	u0 = read_JSON("conds.json") |> extract_dict;
	p1 = p_dict1 |> extract_dict;
	#Finish
	prob1 = ODEProblem(T_ode, u0, (0.0, tmax1), p1)
	sol1 = solve(prob1);
	fig4_1 = plot(sol1, vars = [:v, :n], layout = grid(2,1), 
		c = [:deepskyblue :magenta], lw = 2.0,
		xlabel = ["" "Time (ms)"], 
		ylabel = ["\$V_t\$ (mV)" "\$N_t\$"], ylims = [(-90, 10) (0, 1)], 
	)
	fig4_2 = plot(sol1, vars = (:v, :n),  
		xlabel = "\$V_t\$ (mV)", xlims = (-90, 10.0), 
		ylabel = "\$N_t\$", ylims = (-0.01, 1.0)
	)
	plot(fig4_1, fig4_2, layout = grid(1,2))
end

# ╔═╡ 016e5160-f52a-11ea-11f1-a55e3757bfea
md" 
### 3) Potassium Currents ($I_K$).

##### We can use the value of N as a gating factor for Potassium currents. Therefore the value will fall in a range between [0,1]

$$I_K(V_t, N_t) = -g_K N_t (V_t-E_K)$$
"

# ╔═╡ 6587ed9e-f52a-11ea-13a4-319c3b549822
begin
	t_samples = 0.0:tmax1
	Vt1 = map(t -> sol1(t)[1], t_samples)
	Nt1 = map(t -> sol1(t)[2], t_samples)
	IK1 = -gK*Nt1.*(Vt1 .- EK)
	#ICA1 = -gCa*R_INF(Vt1, V1, V2).*(Vt1.-ECa)
	fig5_1 = plot(sol1, vars = [:v, :n], layout = grid(2,1), lw = 2.0,
		c = [:deepskyblue :magenta], 
		xlabel = ["" ""],
		ylabel = ["\$V_t\$ (mV)" "\$N_t\$"], ylims = [(-90, 10) (0, 1)]
	)
	fig5_2 = plot(t_samples, IK1, 
		ylabel = "\$I_K\$ (pA)", xlabel = "Time (ms)", label = "",
		lw = 2.0, linestyle = :dash, c = :green)
	plot(fig5_1, fig5_2, layout = grid(2,1, heights = (0.70, 0.30)))
end

# ╔═╡ d3d6af70-f3af-11ea-1227-f378f539174e
md"""
### [1.3] Calcium Influx ($C_t$)

$\begin{align}
   \tau_C\frac{dC}{dt} &= C_0 + \delta I_{Ca}(V_t) -\lambda  C_t
\end{align}$

Calcium concentration is dependent on the Calcium currents calculated above in the voltage equation. To prevent the equation from reaching negative concentrations of calcium $C_0$ represents a minimal calcium concentration similar to what was done in the Karvourniari model \cite{Matzakos-karvouniari2019}. Intracellular buffers such as **[Cl-]** and **Calmodulin** cause a decrease in calcium represented by $\lambda C_t$. Overtime calcium concentration will decay, and this is represented by $\tau_C$.

In order to demonstrate this we will simulate a external current injection like would be done using whole cell patch clamp.  
"""

# ╔═╡ 69e747e0-f3b0-11ea-3462-87019b1136b3
md" 
##### Parameters:

C_0: $(@bind C0 NumberField(0.0:0.01:1.0, default = 0.088))
δ: $(@bind δ NumberField(0.0:0.01:1.0, default = 0.010503))
λ: $(@bind λ NumberField(0.1:0.1:5.0, default = 2.702))
τc: $(@bind τc NumberField(1.0:500:5000, default = 2000))
##### Settings
tmax: $(@bind tmax2 NumberField(0:500:100e3, default = 4000)) ms
"

# ╔═╡ f3d20610-f3b1-11ea-3c0a-55e4b84f7abd
begin
	#This block imports all the necessary conditions for the model. 
	p_dict2 = read_JSON("params.json"); 
	#Zero out all conductances to run this experiment
	p_dict2[:I_app] = Iapp
	#p_dict2[:g_TREK] = 0.0; #Zero TREK params
	p_dict2[:g_ACh] = 0.0; #Zero Acetylcholine for now
	#Add in model parameters
	p_dict2[:C_0] = C0
	p_dict2[:τc] = τc
	p_dict2[:δ] = δ
	p_dict2[:λ] = λ
	#Remake the parameters
	p2 = p_dict2 |> extract_dict;
	
	#Make the problem
	prob2 = ODEProblem(T_ode, u0, (0.0, tmax2), p2)
	sol2 = solve(prob2);
	fig6_1 = plot(sol2, 
		vars = [:v, :c], c = [:deepskyblue :green :purple],
		label = ["" "" ""],
		ylabel = ["\$V_t\$ (mV)" "\$C_t\$ (μM)" "\$A_t\$"],
		ylims = [(-Inf, Inf) (-Inf, Inf) (0.0, 1.0)],
		xlabel = ["" "" "Time (ms)"],
		layout = grid(2,1), lw = 2.0
	)
end

# ╔═╡ f3ed2ce0-f52a-11ea-064a-bb3455fc392e
md"""
### [1.4] cAMP breakdown ($A_t$)

The sAHP intermediate is dependent on calcium. We utilize a logistic growth ODE to describe how as calcium concentration is tied to the evolution of second messenger processes. The model described by Ohadi indicates a potential mechanism for Calcium dependent breakdown of cAMP requires 4 Calcium molecules \cite{Ohadi2019}. For this reason $C_t$ is taken to the 4th power. 

$\begin{align}
\tau_A\frac{dA}{dt} &= \alpha  {C_t}^4  \left(1- A_t\right) - A_t
\end{align}$
"""

# ╔═╡ 5b3e1ad0-f535-11ea-19ba-ff8488165431
md" 
##### Parameters:

α: $(@bind α NumberField(0.0:125.0:1500.0, default = 625.0))
τA: $(@bind τA NumberField(0.0:500.0:15e3, default = 8300.0)) ms

##### Settings
tmax: $(@bind tmax3 NumberField(0:500:100e3, default = 40e3)) ms
"

# ╔═╡ 3887dae0-f69c-11ea-16fc-0725694a82c0
begin
	#This block imports all the necessary conditions for the model. 
	p_dict3 = read_JSON("params.json"); 
	#Zero out all conductances to run this experiment
	p_dict3[:I_app] = Iapp
	#p_dict2[:g_TREK] = 0.0; #Zero TREK params
	p_dict3[:g_ACh] = 0.0; #Zero Acetylcholine for now
	#Add in model parameters
	p_dict3[:α] = α
	p_dict3[:τa] = τA
	#Remake the parameters
	p3 = p_dict3 |> extract_dict;
	
	#Make the problem
	prob3 = ODEProblem(T_ode, u0, (0.0, tmax3), p3)
	sol3 = solve(prob3);
	fig7_1 = plot(sol3, 
		vars = [:v, :c, :a], c = [:deepskyblue :green :purple],
		label = ["" "" ""],
		ylabel = ["\$V_t\$ (mV)" "\$C_t\$ (μM)" "\$A_t\$"],
		ylims = [(-Inf, Inf) (-Inf, Inf) (0.0, 1.0)],
		xlabel = ["" "" "Time (ms)"],
		layout = grid(3,1), lw = 2.0
	)
	fig7_2 = plot(sol3, vars = (:c, :a), lw = 2.0, c = :purple,
		xlabel = "\$C_t\$ (μM)",
		ylabel = "\$A_t\$", ylims = (0.0, 1.0)		
	)
	plot(fig7_1, fig7_2, layout = grid(1,2))
end

# ╔═╡ 9cfc8200-f534-11ea-0941-3515f6eb0068
md"""
### [1.5] TREK Activation ($B_t$)

$\begin{align}
\tau_B\frac{dB}{dt} &= \beta  {A_t}^4  (1 - B_t) - B_t
\end{align}$

With high cAMP concentrations, the Ohadi model demonstrates that protein kineases are activated. Physiological imaging studies indicate that, much like cAMP, PKA concentrations also were inversely related to calcium transients \cite{Dunn2006}. The Ohadi model shows that PKA becomes active after binding with 2 cAMPs \cite{Ohadi2019}. Protein structure analysis of the TREK1 channel shows two serine residues (S-300 and S-333) need to be phosphorylated are needed to inactivate a single TREK pore. This analysis showed that PKA activation via Forskolin, adenyl cyclase and cAMP caused an inhibtion of TREK currents  However when mutations were made that removed one or both of the Serine residues, inhibition was no longer possible \cite{Murbartian2005}. 
"""

# ╔═╡ 9792ec10-f6a0-11ea-0a48-d7c842d29005
md" 
##### Parameters:

β: $(@bind β NumberField(0.0:10.0:500.0, default = 34.0))
τb: $(@bind τb NumberField(0.0:500.0:15e3, default = 8300.0)) ms

##### Settings
tmax: $(@bind tmax4 NumberField(0:500:100e3, default = 40e3)) ms
"

# ╔═╡ 9c15bfa0-f6a1-11ea-323a-276826c2a52b
begin
	#This block imports all the necessary conditions for the model. 
	p_dict4 = read_JSON("params.json"); 
	#Zero out all conductances to run this experiment
	p_dict4[:I_app] = Iapp
	#p_dict2[:g_TREK] = 0.0; #Zero TREK params
	p_dict4[:g_ACh] = 0.0; #Zero Acetylcholine for now
	#Add in model parameters
	p_dict4[:β] = β
	p_dict4[:τb] = τb
	#Remake the parameters
	p4 = p_dict4 |> extract_dict;
	
	#Make the problem
	prob4 = ODEProblem(T_ode, u0, (0.0, tmax4), p4)
	sol4 = solve(prob4);
	fig8_1 = plot(sol4, 
		vars = [:c, :a, :b], c = [:green :purple :red],
		label = ["" "" ""],
		ylabel = ["\$V_t\$ (mV)" "\$A_t\$" "\$B_t\$"],
		ylims = [(-Inf, Inf) (-Inf, Inf) (0.0, 1.0)],
		xlabel = ["" "" "Time (ms)"],
		layout = grid(3,1), lw = 2.0
	)
	fig8_2 = plot(sol4, vars = (:a, :b), lw = 2.0, c = :red,
		xlabel = "\$A_t\$",
		ylabel = "\$B_t\$", ylims = (0.0, 1.0)		
	)
	plot(fig8_1, fig8_2, layout = grid(1,2))
end

# ╔═╡ c4c3eb20-f534-11ea-21e2-c502cc5fef5c
md"""
### [1.6] Acetylcholine Release ($E_t$)

Modeling simple release of Acetylcholine utilizes this equation. 

$\begin{align}
    \tau_E\frac{dE}{dt} &=  \rho \Phi(V_t, V_s, V_0) - E_t
\end{align}$

The Acetylcholine receptors recieve acetylcholine and activate via the activation function: 

$\begin{align*}
\bar{H}(E_t, \kappa) &= \frac{E_t^2}{E_t^2 + \kappa^2}\\
\end{align*}$

##### Parameters:

Vs: $(@bind Vs NumberField(0.0:10.0:500.0, default = 0.2))
V0: $(@bind V0 NumberField(0.0:10.0:500.0, default = -40.0))

ρ: $(@bind ρ NumberField(0.0:1.0:50.0, default = 6.0))
τe: $(@bind τe NumberField(0.0:500.0:15e3, default = 540.0)) ms

κ: $(@bind κ NumberField(0.0:0.01:1.0, default = 0.1))

g_ACh: $(@bind gACh NumberField(0.0:0.1:10.0, default = 1.1))
E_ACh: $(@bind EACh NumberField(-40.0:40.0, default = 0.0))

##### Settings
tmax: $(@bind tmax5 NumberField(0:500:100e3, default = 4000)) ms
"""

# ╔═╡ 89c63ca0-f6c2-11ea-166d-9161ba68d706
Φ

# ╔═╡ 26629b30-f6c8-11ea-0086-71194b1b03eb
ħ

# ╔═╡ 7a1152c0-f6c4-11ea-22d8-254557f09f12
begin
	v_series = collect(-90:1.0:10.0)
	ρ_v = map(v -> Φ(v, Vs, V0), v_series)
	fig9_1a = plot(v_series, ρ_v, c = :blue, lw = 2.0,
		xlabel = "\$V_t\$ (mV)",
		ylabel = "Normalized ACh release", 
		title = "Voltage Dependent ACh release"
	)
	#This block imports all the necessary conditions for the model. 
	p_dict5 = read_JSON("params.json"); 
	#Zero out all conductances to run this experiment
	p_dict5[:I_app] = Iapp
	#p_dict2[:g_TREK] = 0.0; #Zero TREK params
	#p_dict5[:g_ACh] = 0.0; #Zero Acetylcholine for now
	#Add in model parameters
	p_dict5[:V0] = V0
	p_dict5[:k] = Vs
	p_dict5[:ρ] = ρ
	p_dict5[:τACh] = τe
	#Remake the parameters
	p5 = p_dict5 |> extract_dict;
	
	#Make the problem
	prob5 = ODEProblem(T_ode, u0, (0.0, tmax5), p5)
	sol5 = solve(prob5);
	fig9_2a = plot(sol5, 
		vars = [:v, :e], c = [:deepskyblue :blue],
		label = ["" ""],
		ylabel = ["\$V_t\$ (mV)" "\$E_t\$ (μM)"],
		xlabel = ["" "Time (ms)"],
		layout = grid(2,1), lw = 2.0
	)
	fig9a = plot(fig9_1a, fig9_2a, layout = grid(1,2))
	
	e_series = collect(0.0:0.05:1.0)
	h_e = map(e -> ħ(e, κ), e_series)
	fig9_1b = plot(e_series, h_e, 
		ylabel = "Normalized Current", ylims = (0.0, 1.0),
		xlabel = "Released ACh"
		
	)
	
	t_series5 = collect(0.0:0.01:tmax5)
	h_t = map(t -> -gACh*ħ(sol5(t)[6], κ)*(sol5(t)[1] - EACh), t_series5)
	fig9_2b = plot(t_series5, h_t, 
		ylabel = "IACh (pA)", xlabel = "Time (ms)"
	)
	
	fig9b = plot(fig9_1b, fig9_2b, layout = grid(1,2))
	plot(fig9a, fig9b, layout = grid(2,1))
end

# ╔═╡ 60671c80-f6c2-11ea-0513-9ff8bf60ae3f
md"""

#### 1) Crank-Nicholson Diffusion of Acetylcholine 

Unlike the previous equations, the diffusion of acetylcholine is a PDE. Once a cell releases Acetylcholine, $E_t$ will spill over into the next cell. This is described 
This equation describes the voltage dependent release of acetylcholine. It also describes the diffusion of Acetylcholine over a 2D grid $\nabla^2$ based on a diffusion coefficient (D). The diffusion of acetylcholine can be converted to the linear form by Crank-Nicholson diffusion scheme by $D \nabla^2 E_t$. 

$\begin{align}
    \tau_E\frac{dE}{dt} &=  D\nabla^2 E_t + \rho  \Phi(V_t, V_s, V_0) - E_t
\end{align}$

The 2D grid ($\nabla^2$) is the linear discretization of the PDE. This means that $D\nabla^2$ 2 dimensional tridiagonal diffusion coefficient matrix. 


In order to create a border, 2 is added to the top most and bottom most value to ensure that the acetylcholine does not diffuse out of the 2D grid which is padded with an infinite amount of zeros. 

D: $(@bind D NumberField(0.0:0.05:1.0, default = 0.1)) 
"""

# ╔═╡ 86967a30-fad8-11ea-131a-87b071185380
∇

# ╔═╡ ac2b4be0-fad8-11ea-08da-43e9440c5dc8
begin 
	dt = 5.0 #ms
	nx = 50; ny = 50
	u = zeros(nx,ny)
	local u
	#u[(nx/2) |> Int, (ny/2) |> Int] = 0.0
	du = similar(u)
	E_t = Float64[]
	t_rng = Float64[]
	anim = @animate for i = 0.0:dt:tmax5
		i % 1000.0 == 0 ? print(i) : nothing
		
		u[(nx/2) |> Int, (ny/2) |> Int] = sol5(i)[6]
		push!(E_t, sol5(i)[6])
		push!(t_rng, i)
		∇(du, u, D)
		u += du
		
		fig11_above = plot(t_rng, E_t, label = "", c = :blue, lw = 2.0,
			xlims = (0.0, tmax5), ylims = (0.0, 0.50)
		)
		
        fig11_1 = plot(collect(1:nx), collect(1:ny), u, 
			seriestype = :surface,
			c = :delta, clims = (0.01, 0.50),
            #ratio = :equal, 
			grid = false,
            xaxis = "", yaxis = "", xlims = (0, nx), ylims = (0, ny),
			zlims = (0.0, 0.3)
        )
		contourf!(fig11_1, u, c = :delta, clims = (0.01, 0.30),
			grid = false, xaxis = "", yaxis = "", 
			xlims = (0, nx), ylims = (0, ny)
		)
		
		u_vert = sum(u, dims = 1)
		fig11_below = plot(u',ylims = (0, 0.3), label = "")
		#u_horizon = sum(u, dims = 2)
		#fig11_side = plot(u_horizon)
		annotate!([10], [3], "t=$(i/1000)s", :white)
		plot(fig11_above, fig11_1, fig11_below, 
			layout = grid(3,1, heights = (0.2, 0.6, 0.2))
		)
	end
	gif(anim, fps = (tmax5/dt)*4)
end

# ╔═╡ 0901b4c2-f535-11ea-320b-d1b33a8f84a0
md"""
### [1.6] Ornstein-Uhlenbeck Stochastic activity ($W_t$)

In order to simulate noise, a Ornstein-Uhlenbeck process (OUP) is simulated (represented by $W_t$). This noise process is widely used in neuronal modelling (\cite{Tuckwell1988}\cite{Tuckwell1989} for reviews), and accounts for regularly times ionic diffusion events below the spiking threshold \cite{Lanska1994}. The equation for this process is represented below. 

$\begin{align}
   \tau_W\frac{dW}{dt} &= -W_t
\end{align}$

Similar to a simple white noise Wiener process, the OUP random walk includes an amplitude term ($\sigma$) corresponding to the mean of a Gaussian distribution. However, OUP also includes a term for correlation from one noise event to the next ($\tau_W$)  This allows the noise to be grouped into bursting events spaced out by the correlation term. Because the Gaussian process is a continuous function ranging from $-\infty$ to $\infty$, currents can be either positive or negative and therefore have a chance of either be hyperpolarizing or depolarizing. By adjusting the correlation time, noise bursts can occur at different time scales. 

σ: $(@bind σ NumberField(0.0:0.05:1.0, default = 0.1)) 
τw: $(@bind τw NumberField(1.0:500.0:100e3, default = 1e3))  

#### Settings
tmax: $(@bind tmax6 NumberField(1.0:500:100e3, default = 60e3)) 
"""

# ╔═╡ 4236e360-fb6a-11ea-0f67-1de40bb2f748
begin
	#This block imports all the necessary conditions for the model. 
	p_dict6 = read_JSON("params.json"); 
	#For this simulation, no parameters need to be zeroed out
	p_dict6[:σ] = σ
	p_dict6[:τw] = τw
	#Remake the parameters
	p6 = p_dict6 |> extract_dict;
	
	#Make the problem
	prob6 = SDEProblem(T_sde, u0, (0.0, tmax6), p6)
	@time sol6 = solve(prob6, SOSRI());
	fig11_1 = plot(sol6, vars = [:v, :W], label = ["" ""],
		c = [:deepskyblue :gray],lw = [2.0 2.0],
		xlabel = ["" "Time (ms)"], ylabel = ["\$V_t\$ (mV)" "\$W_t\$ (pA)"], 
		ylims = [(-90,10) (-15.0, 15.0)],
		layout = grid(2,1)
	)
	
	t_series2 = collect(0.0:1.0:sol6.t[end])
	Wt = map(t -> sol6(t)[7], t_series2)
	noise_fit = fit(Normal, Wt)
	mu = round(noise_fit.μ, digits = 2)
	sig = round(noise_fit.σ, digits = 2)
	
	fig11_2 = histogram(Wt, c = :gray, 
				xlabel = "Noise (pA)", 
				ylabel = "PDF", 
				label = "", 
				normalize = :pdf, 
				xlims = xlims
			)
	plot!(fig11_2, noise_fit, c = :black, 
    	label = "Fit (μ = $mu pA σ = $sig pA)",
		linestyle = :dash, lw = 4.0
	)
	
	
	plot(fig11_1, fig11_2, layout = grid(2,1, heights = (0.6, 0.3)))
end

# ╔═╡ cc2cc220-fb8c-11ea-1f3d-03247240956d


# ╔═╡ Cell order:
# ╠═77bf8f30-f211-11ea-2516-3d678de83568
# ╟─d2bde4d0-f1e0-11ea-395a-21b3cc2800b1
# ╠═fe8ed2ce-f1e1-11ea-2c15-ad90425e6b01
# ╠═6e849542-f2d5-11ea-0eba-632e96ea9ded
# ╠═c52c2942-f207-11ea-33e2-1f597edebb9c
# ╠═84194d80-fb8e-11ea-2d3f-a9a237c412ef
# ╠═3826848e-f208-11ea-1b46-0d98f24a3da7
# ╠═55ce0450-f203-11ea-3d1b-d52ee88fffc3
# ╟─d468e040-f1e1-11ea-152f-dff8d05c8594
# ╟─0cdbf9e2-f1e1-11ea-3a06-a702d6cae6fb
# ╟─6b1a70a0-f203-11ea-0026-b5faa986fecf
# ╟─6db280a0-f203-11ea-2cc1-0381b15e1410
# ╟─b82b71a0-f203-11ea-3dc9-93fc0632e420
# ╟─9861a390-f207-11ea-3f45-6198c2f22003
# ╟─c1e40650-f20b-11ea-34a8-e157d35f600e
# ╟─bda8ed50-f20e-11ea-2e23-075cedea76b7
# ╟─0962f460-f210-11ea-21f1-0349a3e20635
# ╟─8e05ba42-f210-11ea-0684-3972c2fffb56
# ╟─7d839cb0-f20f-11ea-122d-cb56b26e2f77
# ╟─4df37fa2-f215-11ea-2a3f-f3c139376fad
# ╟─02fb4610-f2d1-11ea-0ea8-5d9a7cf72864
# ╟─9277f880-f2c0-11ea-3491-8f9f2ad9941d
# ╟─75f593a0-f2cc-11ea-23d6-a71cc7db012a
# ╟─0eb2c2d0-f23b-11ea-1b27-0541ac4d0b36
# ╟─ec80e870-f2cd-11ea-1c3a-05beec6b0480
# ╟─aace080e-f2d0-11ea-15cd-bb7c1926458f
# ╟─a79cda6e-f2d2-11ea-27d9-4df5215595fe
# ╟─8510a620-f2d8-11ea-07bd-3146b9fd9fc5
# ╟─016e5160-f52a-11ea-11f1-a55e3757bfea
# ╟─6587ed9e-f52a-11ea-13a4-319c3b549822
# ╟─d3d6af70-f3af-11ea-1227-f378f539174e
# ╟─69e747e0-f3b0-11ea-3462-87019b1136b3
# ╟─f3d20610-f3b1-11ea-3c0a-55e4b84f7abd
# ╟─f3ed2ce0-f52a-11ea-064a-bb3455fc392e
# ╟─5b3e1ad0-f535-11ea-19ba-ff8488165431
# ╟─3887dae0-f69c-11ea-16fc-0725694a82c0
# ╟─9cfc8200-f534-11ea-0941-3515f6eb0068
# ╟─9792ec10-f6a0-11ea-0a48-d7c842d29005
# ╟─9c15bfa0-f6a1-11ea-323a-276826c2a52b
# ╟─c4c3eb20-f534-11ea-21e2-c502cc5fef5c
# ╟─89c63ca0-f6c2-11ea-166d-9161ba68d706
# ╟─26629b30-f6c8-11ea-0086-71194b1b03eb
# ╟─7a1152c0-f6c4-11ea-22d8-254557f09f12
# ╟─60671c80-f6c2-11ea-0513-9ff8bf60ae3f
# ╟─86967a30-fad8-11ea-131a-87b071185380
# ╠═ac2b4be0-fad8-11ea-08da-43e9440c5dc8
# ╟─0901b4c2-f535-11ea-320b-d1b33a8f84a0
# ╟─4236e360-fb6a-11ea-0f67-1de40bb2f748
# ╠═cc2cc220-fb8c-11ea-1f3d-03247240956d
