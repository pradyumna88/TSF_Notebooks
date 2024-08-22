### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ccca791d-b819-43aa-9c62-94df8657ac45
# ╠═╡ show_logs = false
begin
	using Revise
	using Pkg
	Pkg.activate("TSF_110")
	using QMCB
	using QMCB.Defs
	using QMCB.SimParameters
	using QMCB.Lattice
	using QMCB.Measurements
	using QMCB.MoveStruct
	using QMCB.WL_actions
	using QMCB.Results
	using QMCB.Plotter
	using Parameters
	using Dates
	using JLD2
	using Random
	using Parameters
	using Setfield
	using Distributed
	using StaticArrays
	using LinearAlgebra
	using Statistics
	using PlutoUI
	using Plots
end

# ╔═╡ db4e7db6-7f42-4da5-aec1-5e2ee05fd7ae
using DrWatson

# ╔═╡ 92fd972e-5a3e-11ef-27c3-cfcf217458f8
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2200px;
    	padding-left: max(50px, 1%);
    	padding-right: max(25px, 1%);
	}
</style>
"""

# ╔═╡ 71d10ec2-27fc-457a-9d65-4a436a57d714
mymap(f) = x -> map(f, x)

# ╔═╡ 2c4db510-1d7f-4d1a-82e4-936f994bfd46
pwd()

# ╔═╡ 13a16316-c6a5-4135-bbf2-b1e7b1a6caae
prefix_path = "/home/data/PD/Runs"

# ╔═╡ 47869b79-57ab-4c02-a1e7-aad9ffb0fa49
foldernames = prefix_path |> Results.foldernames_old

# ╔═╡ 2afcebe1-7ca7-4f01-af35-48a183cfa8ea
filepaths = prefix_path |> Results.filepaths_old

# ╔═╡ 9d3ce4a3-d7b9-4bfa-a8fa-462c210f1624
N1_selector = @bind N1 Select(filepaths |> Results.get_N1s_old);

# ╔═╡ bf1c51d2-93e2-4113-9e5d-ca020c8ed80a
T_selector = @bind T Select(filepaths |> x -> Results.get_Ts_old(x));

# ╔═╡ fb9562ac-3fac-459d-a93d-fc52e70ecbd2
n1_selector = @bind n1 Select(filepaths |> x -> Results.get_n1s_old(x));

# ╔═╡ 9719b0b1-5459-413e-a599-090a1c1a4766
phys_params = (N1, T, n1)

# ╔═╡ ea7b2b0a-0d8b-43a6-a6ee-e8fd5f4a6b30
run_version_selector = @bind run_version Select(filepaths |> x -> Results.get_run_version_old(x, phys_params));

# ╔═╡ 090048c4-19d9-4a5a-9bd4-6be6d34b8351
run_time = Results.filepaths_old(prefix_path) |> x -> Plotter.filter_datafilepaths_old(x, 5, (phys_params..., run_version))[end] |> x -> load(x, "run_time") |> x -> x/(60*60*24) |> x -> round(x, digits=2)

# ╔═╡ 376451b5-395e-44bb-a75e-6424f6f56375
current_meas = Results.filepaths_old(prefix_path) |> x -> Plotter.filter_datafilepaths_old(x, 5, (phys_params..., run_version))[end] |> x -> load(x, "current_meas") 

# ╔═╡ ea28abf7-61de-497b-8544-3f584adb35ac
ignore_selector = @bind ignore Select([:bare, :ignore_zeros, :ignore_non_trimer]);

# ╔═╡ 089c5d73-5dbb-495f-8766-a88ad25eb41f
md""" $N1_selector $T_selector $n1_selector $run_version_selector $ignore_selector"""

# ╔═╡ 9d187387-ac73-40e2-b594-e77eaf29ef09
begin
	filepath = Results.filepaths_old(prefix_path) |> x -> Plotter.filter_datafilepaths_old(x, 8, (phys_params...,run_version))[1]
	pd = plot(
		filepath |> 
		x -> Plotter.get_plotting_data(x, x -> nothing, x -> load(x, "superfluid_convergence") |> x -> getproperty(x, ignore) |> x -> (x.plus, x.minus)) |>
		x -> plot(x.xs, x.ys *phys_params[2]*100*3/(2*phys_params[3]) |> transpose, yerr = x.yers *phys_params[2]*100*3/(2*phys_params[3]) |> transpose, title = (run_time, current_meas/10^3)) |>
		x -> hline!(x, [100, 2*phys_params[2]*3*100/(pi*phys_params[3])]),
		filepath |> 
		x -> Plotter.get_plotting_data(x, x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).plus |> mymap(first) |> x -> x*4, 
			x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).plus |> mymap(last), yerr_flag=false) |>
		x -> bar(x.xs, x.ys ./ sum(x.ys), bar_width = .25),
		filepath |> 
		x -> Plotter.get_plotting_data(x, x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).minus |> mymap(first), 
			x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).minus |> mymap(last), yerr_flag=false) |>
		x -> bar(x.xs, x.ys ./ sum(x.ys), bar_width = .25, c = :red),
		filepath |> x -> load(x, "winding_numbers") |> x -> 2*@view(x[:,:,1]) + @view(x[:,:,2]) |> x -> reshape(x, :) |> x -> x[1:10:end] |> x -> plot(range(1, 100, 10000), x[Int64.(ceil.(length(x) *range(.01, 1, 10000)))]) |>
		x -> hline!(x, 4*(-2:2), label = false),
		# layout = (2, 2),size = (800, 800), 
		layout = (1, 4),size = (1800, 400), 
		margin = 9Plots.mm
	)
	md"""$pd"""
end

# ╔═╡ ef36fdc6-1900-4437-852c-da95c1869b87
# ╠═╡ disabled = true
#=╠═╡
plot(
	filepath |> x -> load(x, "winding_numbers") |> x -> 2*@view(x[:,:,1]) + @view(x[:,:,2]) |> x -> reshape(x, :) |> x -> x[1:10:end] |> plot,
	filepaths |> x -> Results.get_sim_params(x, phys_params, run_version) |> x -> x[sortperm(multi_ranker)] |> x -> x[1:20] |>
	x -> scatter(x |> mymap(first) |> mymap(x -> log2(x)), x |> mymap(last)),
	size = (800, 400)
);
  ╠═╡ =#

# ╔═╡ ce41903e-6bd2-4235-8895-a76b5592ce15
ranker(rev = false) = x -> x |> (x -> sortperm(x, rev=rev) .=> 1:length(x)) |> x -> sort(x, by = x -> x[1]) |> mymap(last)

# ╔═╡ afc91517-d6ee-44aa-b67e-1fd252192ca9
# scatter(
		# 	Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> 
		# 	x -> Plotter.world_lines(x, 2, 1),
		# 	Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> 
		# 	x -> Plotter.world_lines(x, 2, 2),
		# 	ms = 2, msw = 0, c = :red
		# ) |> 
		# x -> scatter!(x, 
		# 	Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> 
		# 	x -> Plotter.world_lines(x, 1, 1),
		# 	Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> 
		# 	x -> Plotter.world_lines(x, 1, 2),
		# 	ms = 2, msw = 0, c = :blue
		# ),

# ╔═╡ 794740ff-310a-40a1-b8af-ae9f0c3c46ff


# ╔═╡ Cell order:
# ╟─92fd972e-5a3e-11ef-27c3-cfcf217458f8
# ╟─db4e7db6-7f42-4da5-aec1-5e2ee05fd7ae
# ╠═ccca791d-b819-43aa-9c62-94df8657ac45
# ╟─71d10ec2-27fc-457a-9d65-4a436a57d714
# ╟─2c4db510-1d7f-4d1a-82e4-936f994bfd46
# ╟─13a16316-c6a5-4135-bbf2-b1e7b1a6caae
# ╠═47869b79-57ab-4c02-a1e7-aad9ffb0fa49
# ╠═2afcebe1-7ca7-4f01-af35-48a183cfa8ea
# ╠═9d3ce4a3-d7b9-4bfa-a8fa-462c210f1624
# ╠═bf1c51d2-93e2-4113-9e5d-ca020c8ed80a
# ╠═fb9562ac-3fac-459d-a93d-fc52e70ecbd2
# ╟─9719b0b1-5459-413e-a599-090a1c1a4766
# ╟─ea7b2b0a-0d8b-43a6-a6ee-e8fd5f4a6b30
# ╟─090048c4-19d9-4a5a-9bd4-6be6d34b8351
# ╟─376451b5-395e-44bb-a75e-6424f6f56375
# ╟─ea28abf7-61de-497b-8544-3f584adb35ac
# ╟─089c5d73-5dbb-495f-8766-a88ad25eb41f
# ╠═9d187387-ac73-40e2-b594-e77eaf29ef09
# ╟─ef36fdc6-1900-4437-852c-da95c1869b87
# ╟─ce41903e-6bd2-4235-8895-a76b5592ce15
# ╟─afc91517-d6ee-44aa-b67e-1fd252192ca9
# ╠═794740ff-310a-40a1-b8af-ae9f0c3c46ff
