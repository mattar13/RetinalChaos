### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# ╔═╡ bc9b17c0-4674-44e7-af4d-e85ed0cf98e0
using Revise

# ╔═╡ fb0b35b0-9a27-11eb-1925-338bb9e3754f
using RetinalChaos

# ╔═╡ 0f1ae382-f251-49fa-8925-db18e832e228
md"
# Making figure 2
"

# ╔═╡ a2b3ce24-6874-4d54-bf98-c2bfe7649bf4
import RetinalChaos: Φ, diffuse, ħ

# ╔═╡ 4187d216-079d-46d7-a6e0-77d2e20f17e3
begin
	#Set up fonts 
	font_title = font("Arial", 24) #title font and size
	font_axis = font("Arial", 12) #axis font and size
	font_legend = font("Arial", 8) #legend font and size
	gr(titlefont=font_title, guidefont = font_axis, legendfont = font_legend)
end

# ╔═╡ d3b0c671-734f-4fff-a26c-0f6c2ded5d3a
begin
	#Set up the parameter files and save files
	param_root = joinpath(splitpath(pwd())[1:end-1]..., "params\\")
	params_file = joinpath(param_root, "params.json")
	conds_file = joinpath(param_root, "conds.json")
end

# ╔═╡ 15393dea-28b8-4234-8e03-35f7e0f8cf56
begin
	p = read_JSON(params_file) 
	p[:I_app] = 10.0
	#Set up the initial conditions
	u0 = read_JSON(conds_file);
	#Set up the Time span
	tspan = (0.0, 120e3)
	#set up the problem to solve
	prob = ODEProblem(T_ode, u0 |> extract_dict, tspan, p |> extract_dict);
end

# ╔═╡ 3d1a1189-dd11-4f12-add2-44ba7c0b26e2
#Solve the problem
@time sol = solve(prob, abstol = 2e-2, reltol = 2e-2, maxiters = 1e7); 

# ╔═╡ d26aa235-c403-4ce3-8934-a598c5566be5
begin 
	#Run analysis 
	dt = 0.1 #set the time differential
	v_thresh = calculate_threshold(sol, dt = dt)
	timestamps = get_timestamps(sol, dt = dt)
	burst_idxs, dur_list, spb_list, ibi_list = max_interval_algorithim(sol, dt = dt)
end

# ╔═╡ 74f7a4e8-cfd4-4b37-b835-cf42d21301d4
begin
	burst_lims = burst_idxs[2]
	A_dx = 500 #the tick interval is 20ms
	A_xlims = (burst_lims[1]-2000, burst_lims[2]+2000)
	A_xticks = A_xlims[1]:A_dx:A_xlims[2] #Do I actually need to collect here
	A_trng = A_xlims[1]:dt:A_xlims[2]
	
	#Pick 3 evenly spaced sections to plot
	frame_stops = LinRange(burst_lims[1], burst_lims[2]+500, 3)
	ach_stops = sol(frame_stops, idxs = 6)
	

end

# ╔═╡ 88b02dbb-2108-43f8-911c-949538105ac1
begin
	fig2_Aa1 = plot(A_trng, sol(A_trng, idxs = 1),
		xlabel = "", ylabel = "Vₜ (mV)", 
		lw = 3.0, c = v_color, grid = false,
		xlims = A_xlims, xticks = false
	)
	fig2_Aa2 = plot(A_trng, sol(A_trng, idxs = 6),
		xlabel = "Time (s)", ylabel = "Eₜ (μM)", 
		lw = 3.0, c = e_color, grid = false,
		xlims = A_xlims, xticks = A_xticks, 
		x_formatter = x -> (x-A_xlims[1])/1000
	)
	fig2_Aa = plot(fig2_Aa1, fig2_Aa2, layout = grid(2,1))	
	plot!(fig2_Aa[2], frame_stops, ach_stops, label = "Frame stops",
    	seriestype = :scatter, marker = :star, markersize = 10.0
	)
	
	fig2_Ab = plot(v -> Φ(v, p[:k], p[:V0]), -90, 10, 
        legend = false, xlabel = "Vₜ (mV)", ylabel = "Φ (Normalized ACh release)",
		c = e_color, lw = 3.0, grid = false, linestyle = :dash
    )
	
	fig2_A = plot(fig2_Aa, fig2_Ab, 
		layout = grid(1,2, widths = (0.75, 0.25)), margin = 0.0mm
	)
	title!(fig2_A[1], "A", titlepos = :left)
end

# ╔═╡ a9faf5c6-e40e-4b3d-bdc9-3ea0b3207c25
begin 
	#Run a short wave lattice simulation
	sim_rng = A_trng# burst_lims[1]:dt:burst_lims[2]
	nx, ny = (50, 50) #Set up the size of the sim
	#c1x, c1y = (round(Int, nx/2), round(Int, ny/2)) #Pick midpoints to place the cell
	c1x, c1y = (1, round(Int, nx/2)) #Pick midpoints to place the cell
	
	bp_model = Network(nx, ny, μ = 0.15, version = :gACh) #Make the network
	#Run over a period of time from A_xlims[1] -> A_xlims[2]
	lattice_c = zeros(nx, ny, length(sim_rng))
	println(size(lattice_c))
	

	
	for (idx,t) = enumerate(sim_rng)
		println(t)
		if idx == 1
			lattice_c[c1y, c1x, idx] = sol(t, idxs = 6)
		else
			lattice_c[:, :, idx] = diffuse(lattice_c[:,:,idx-1], p[:D], bp_model)
			lattice_c[c1y, c1x, idx] = sol(t, idxs = 6)
		end
	end
	
	
end	

# ╔═╡ 277d5329-62ae-461a-b094-bae821fd739f
begin
	#This is for coloring on B and C
	top_lim = 0.15
	n_samples = 6
	cell_samples = round.(Int64, LinRange(2, c1y, n_samples))
	mapping = LinRange(0.0, top_lim, n_samples)

	d_frame = 200
	anim = @animate for frame_i in 1:d_frame:size(lattice_c,3)-1
		println(frame_i)
		#println("Animating at time: $(sim_rng[frame_i])")

		plot(lattice_c[:,:, round(Int64, frame_i)],
			st = :contourf,
			ratio = :equal, grid = false,	
			xaxis = false, yaxis = false, 
			xticks = false, yticks = false,
			xlims = (0, nx), ylims = (0, ny),
			c = :thermal, clims = (0.0, top_lim), 		
		)
		
		for (color_i,cix) in enumerate(cell_samples)
			plot!([cix], [c1y],
				st = :scatter,
				markercolor = :jet,
				marker = :hex,
				marker_z = mapping[color_i],
				label = color_i == 1 ? "cell $color_i" : ""
			)
		end
		annotate!([25], [10], 
			"idx = $(frame_i)", :white
		)
		annotate!([25], [5], 
			"t = $(round((sim_rng[frame_i]-sim_rng[1])/1000, digits = 2)) ms", :white
		)
	end
	gif(anim, "ACh_Release.gif", fps = 2000/d_frame)
end

# ╔═╡ 21e5b5eb-d4ad-433f-b37e-2d7ee42ab0b3
begin
	fig2_B_map = plot(
		xlims = (0, nx), ylims = (0, ny), 
		layout = grid(1, length(frame_stops)), 
		colorbar_titlefontsize = 0.2,
		#margin = 1mm
	)
	
	for (idx, frame) in enumerate(frame_stops)
    	frame_idx = round(Int, (frame-A_xlims[1])/dt)
		plot!(fig2_B_map[idx], 
			lattice_c[:,:, round(Int64, frame_idx)], cbar = false,
			st = :contourf,
			ratio = :equal, grid = false,
			xaxis = false, yaxis = false, 
			xticks = false, yticks = false,
			c = :thermal, clims = (0.0, top_lim), 
			margins = 0.0mm
			)
		for (color_i,cix) in enumerate(cell_samples)
			plot!(fig2_B_map[idx],  [cix], [c1y],
				st = :scatter,
				markercolor = :jet,
				marker = :hex,
				marker_z = mapping[color_i],
				label = idx == 1 ? "cell $color_i" : ""
			)
		end
		annotate!(fig2_B_map[idx], [25], [5], 
			"t = $(round((sim_rng[frame_idx]-sim_rng[1])/1000, digits = 1)) ms", 				:white
				)
	end
	

	
	fig2_B_cbar = heatmap(
		repeat(collect(LinRange(0.0, top_lim, 16)), 1, 2)', 
		c = :thermal, cbar = false, 
		title = "Eₜ (μM)", titlefontsize = 0.2,
		yticks = false, 
		xticks = (1:16, round.(LinRange(0.0, top_lim, 16), digits = 2)),
		#size = (1000,100)
	)
	fig2_B = plot(fig2_B_map, fig2_B_cbar, 
		layout = grid(2,1, heights = (0.95, 0.05))
	)
	title!(fig2_B[1], "B", title_pos = :left)
end

# ╔═╡ eda807d8-a93d-45a4-be52-9f50174f091b
begin	
	fig2_Ca = plot()
	fig2_Cb = plot()
	for (color_i,cix) in enumerate(cell_samples)
		plot!(fig2_Ca, sim_rng, lattice_c[c1y, cix, :],
			c = :jet, line_z = mapping[color_i], 
			lw = 3.0, grid = false, cbar = false,
			label = "",
			xlabel = "Time (s)", ylabel = "Extracellular Eₜ (μM)", 
			xticks = sim_rng[1]:500:sim_rng[end], 
			x_formatter = x -> round(x -sim_rng[1])/1000
		)
		vc = -70.0
		I_ach(a) = -p[:g_ACh] * ħ(a, p[:k_d]) * (vc - p[:E_ACh])
		plot!(fig2_Cb, sim_rng, I_ach.(lattice_c[c1y, cix, :]),
			c = :jet, line_z = mapping[color_i], 
			lw = 3.0, grid = false, cbar = false,
			label = "",
			xlabel = "Time (s)", ylabel = "I_ACh (pA)", 
			xticks = sim_rng[1]:500:sim_rng[end], 
			x_formatter = x -> round(x -sim_rng[1])/1000		
		)
	end
	fig2_C = plot(fig2_Ca, fig2_Cb, layout = grid(1,2))
	#Make something with vclamp-esque properties
	title!(fig2_C, "C", title_pos = :left)
end

# ╔═╡ f3679820-7e5e-4bef-b6d6-412c3262b11f
fig2 = plot(fig2_A, fig2_B, fig2_C, 
	layout = grid(3,1, heights = (0.2, 0.5, 0.3)), 
	size = (800, 1000)
	)

# ╔═╡ ba39a93d-15e1-4f29-8f55-ddbcd2350aa5
savefig(fig2, "Fig2_Acteylcholine_Dynamics.png")

# ╔═╡ Cell order:
# ╟─0f1ae382-f251-49fa-8925-db18e832e228
# ╠═bc9b17c0-4674-44e7-af4d-e85ed0cf98e0
# ╠═fb0b35b0-9a27-11eb-1925-338bb9e3754f
# ╠═a2b3ce24-6874-4d54-bf98-c2bfe7649bf4
# ╠═4187d216-079d-46d7-a6e0-77d2e20f17e3
# ╠═d3b0c671-734f-4fff-a26c-0f6c2ded5d3a
# ╠═15393dea-28b8-4234-8e03-35f7e0f8cf56
# ╠═3d1a1189-dd11-4f12-add2-44ba7c0b26e2
# ╠═d26aa235-c403-4ce3-8934-a598c5566be5
# ╠═74f7a4e8-cfd4-4b37-b835-cf42d21301d4
# ╟─88b02dbb-2108-43f8-911c-949538105ac1
# ╟─a9faf5c6-e40e-4b3d-bdc9-3ea0b3207c25
# ╠═277d5329-62ae-461a-b094-bae821fd739f
# ╟─21e5b5eb-d4ad-433f-b37e-2d7ee42ab0b3
# ╟─eda807d8-a93d-45a4-be52-9f50174f091b
# ╟─f3679820-7e5e-4bef-b6d6-412c3262b11f
# ╠═ba39a93d-15e1-4f29-8f55-ddbcd2350aa5
