### A Pluto.jl notebook ###
# v0.14.9

using Markdown
using InteractiveUtils

# ╔═╡ 2756958e-e5b1-11eb-3da6-6f800a6d618c
using Revise

# ╔═╡ 09475831-f6a1-4989-8884-2369e29e16fc
using RetinalChaos

# ╔═╡ 63786ca7-4b1a-4805-be1e-eae87fedab19
using JLD2, StatsPlots

# ╔═╡ aa88f543-633a-476b-b6cb-441d020dd649
md"
For this one we need to open all of the simulation files I previously opened
"

# ╔═╡ 8c7722a0-22ff-449d-80e6-9639761f2915
#Lets find the file where all other files are located

# ╔═╡ c8d0590c-b9f5-4b2d-95eb-19ecaf816cee
begin
	directory = "E:\\Data\\Modelling\\mu_experiment"
	all_dirs = readdir(directory)
end

# ╔═╡ cc95c7d5-8190-4ad9-94f6-6d77e1afddb2
begin 
	#Lets try opening files
	tstamps_file = "$(directory)\\$(all_dirs[1])\\timestamps.jld2"
	data_file = "$(directory)\\$(all_dirs[1])\\data.jld2"
	dir1 = "$(directory)\\$(all_dirs[1])\\"
	timestamps = JLD2.load(tstamps_file)
	data = JLD2.load(data_file)
end

# ╔═╡ edd29aaa-3cce-4f37-a469-342b9798a8ff
begin
	p_dict = read_JSON("$(dir1)\\params.json")
	u_dict = read_JSON("$(dir1)\\conds.json")
	load_model("$(directory)\\$(all_dirs[1])", p_dict, u_dict)
end

# ╔═╡ 843793dc-da78-4446-9543-6711d486c915
tstamps_file

# ╔═╡ 8040d675-49f3-4ef2-91af-92e25edd05c6
histogram(data["BurstDurs"], yaxis = :log)

# ╔═╡ a6783568-a568-4b51-8813-36c799f51dbc
maximum(data["SpikeDurs"])

# ╔═╡ Cell order:
# ╠═2756958e-e5b1-11eb-3da6-6f800a6d618c
# ╠═09475831-f6a1-4989-8884-2369e29e16fc
# ╠═63786ca7-4b1a-4805-be1e-eae87fedab19
# ╟─aa88f543-633a-476b-b6cb-441d020dd649
# ╠═8c7722a0-22ff-449d-80e6-9639761f2915
# ╠═c8d0590c-b9f5-4b2d-95eb-19ecaf816cee
# ╠═cc95c7d5-8190-4ad9-94f6-6d77e1afddb2
# ╠═edd29aaa-3cce-4f37-a469-342b9798a8ff
# ╠═843793dc-da78-4446-9543-6711d486c915
# ╠═8040d675-49f3-4ef2-91af-92e25edd05c6
# ╠═a6783568-a568-4b51-8813-36c799f51dbc
