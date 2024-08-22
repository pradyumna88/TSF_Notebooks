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
prefix_path = "../"

# ╔═╡ 47869b79-57ab-4c02-a1e7-aad9ffb0fa49
foldernames = prefix_path |> Results.foldernames

# ╔═╡ b0ef3f79-f6dc-462a-8ea9-f36dd9f1f583
filenames = prefix_path |> Results.filenames

# ╔═╡ 2afcebe1-7ca7-4f01-af35-48a183cfa8ea
filepaths = prefix_path |> Results.filepaths

# ╔═╡ 916bedc1-18e2-4a7b-b275-c202543600bf
phys_params_selector = @bind phys_params Select(filepaths |> Results.get_phys_params);

# ╔═╡ ea7b2b0a-0d8b-43a6-a6ee-e8fd5f4a6b30
run_version_selector = @bind run_version Select(filepaths |> x -> Results.get_run_version(x, phys_params));

# ╔═╡ 082d2074-fed5-4b9e-8dab-8c648cbbd25a
param_C_selector = @bind param_C Select(filepaths |> x -> Results.get_Cs(x, phys_params, run_version));

# ╔═╡ af611a51-a939-4ac3-949d-6b762eac94b9
param_Mbar_selector = @bind param_Mbar Select(filepaths |> x -> Results.get_Mbars(x, phys_params, run_version));

# ╔═╡ ea28abf7-61de-497b-8544-3f584adb35ac
ignore_selector = @bind ignore Select([:bare, :ignore_zeros, :ignore_non_trimer]);

# ╔═╡ 089c5d73-5dbb-495f-8766-a88ad25eb41f
md""" $phys_params_selector $run_version_selector $ignore_selector"""

# ╔═╡ ce41903e-6bd2-4235-8895-a76b5592ce15
ranker(rev = false) = x -> x |> (x -> sortperm(x, rev=rev) .=> 1:length(x)) |> x -> sort(x, by = x -> x[1]) |> mymap(last)

# ╔═╡ f9231ee4-9e0e-408d-a35d-de0aac7e74dd
ranker1 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5,6], (phys_params..., run_version)) |> mymap(x -> load(x, "superfluid")) |> mymap(x -> x.ignore_zeros.plus[1]) |> ranker(true)

# ╔═╡ 3fd6cf00-87a3-495f-88a0-f9e2abec657e
ranker2 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5,6], (phys_params..., run_version)) |> mymap(x -> load(x, "superfluid_winding_number_frequency")) |> mymap(x -> x.ignore_non_trimer.plus) |> mymap(length) |> ranker(true)

# ╔═╡ 1edf5403-af4f-4df4-9bde-c8543e9a96ef
ranker3 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5,6], (phys_params..., run_version)) |> mymap(x -> load(x, "superfluid_winding_number_frequency")) |> mymap(x -> x.ignore_non_trimer.plus) |>
mymap(mymap(x -> x[1] == 0 ? 0 : x[1] < 0 ? -x[2] : x[2])) |> mymap(x -> x/sum(abs.(x))) |> mymap(mean) |> mymap(abs) |> ranker()

# ╔═╡ f016d916-0ce8-4e59-80bb-1fdeff754245
ranker4 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5,6], (phys_params..., run_version)) |> mymap(x -> load(x, "superfluid_convergence").ignore_non_trimer.plus) |> mymap(mymap(first)) |> mymap(x -> x[end] - x[Int64(ceil(.8*length(x)))]) |> ranker()

# ╔═╡ efd358d3-e6c1-4abd-9a00-177f38583e67
multi_ranker = (ranker1, ranker2, ranker3) |> mymap(x -> 1 ./(x .+ 10)) |> sum |> ranker(true)

# ╔═╡ 3da8cbf1-c7bf-4adf-8cce-b1a2f54eebf9
# ╠═╡ show_logs = false
sim_params_params_selector = @bind sim_params Select(filepaths |> x -> Results.get_sim_params(x, phys_params, run_version) |> x -> x[sortperm(multi_ranker)]);

# ╔═╡ 090048c4-19d9-4a5a-9bd4-6be6d34b8351
run_time = Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> x -> load(x, "run_time") |> x -> x/(60*60*24) |> x -> round(x, digits=2)

# ╔═╡ 376451b5-395e-44bb-a75e-6424f6f56375
current_meas = Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> x -> load(x, "current_meas") 

# ╔═╡ 01e66054-1398-4878-bd36-95d7e7ad0923
 md""" $sim_params_params_selector """

# ╔═╡ 9d187387-ac73-40e2-b594-e77eaf29ef09
begin
	filepath = Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end]
	pd = plot(
		filepath |> 
		x -> Plotter.get_plotting_data(x, x -> nothing, x -> load(x, "superfluid_convergence") |> x -> getproperty(x, ignore) |> x -> (x.plus, x.minus)) |>
		x -> plot(x.xs, x.ys *phys_params[2]*100*3/(2*phys_params[3]) |> transpose, yerr = x.yers *phys_params[2]*100*3/(2*phys_params[3]) |> transpose, title = (run_time, current_meas/10^4)) |>
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
		layout = (1, 4),
		size = (1800, 400), 
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

# ╔═╡ c7abe550-993f-4684-a02d-883af41e8c1e
@bind N1 Select(filepaths |> Results.get_N1s)

# ╔═╡ f9f2b9e3-921e-43a5-abe8-3d374716e5cb
begin
	function get_stiffness_T()
		filepaths |> x -> Plotter.filter_datafilepaths(x, (2, 5, 6), (N1, .7, .2, "#")) |> x -> map(Results.get_parameter, x) .=> x
	end
	get_stiffness_T()
end

# ╔═╡ 2763c7fa-ef74-4417-a1a5-f833753dd783
param_Mbar_selector

# ╔═╡ b74db8bb-e9d9-4609-aa86-8a683ee8ae0b
filepaths |> x -> Plotter.filter_datafilepaths(x, 5, (phys_params..., param_Mbar, run_version)) |> x -> convert(Array{String}, x) |> x -> Plotter.get_plotting_data(x, 5, x -> load(x, "superfluid").ignore_zeros.plus) |>
x -> plot(x.xs, x.ys, yerr = x.yers)

# ╔═╡ 794740ff-310a-40a1-b8af-ae9f0c3c46ff


# ╔═╡ Cell order:
# ╟─92fd972e-5a3e-11ef-27c3-cfcf217458f8
# ╟─db4e7db6-7f42-4da5-aec1-5e2ee05fd7ae
# ╠═ccca791d-b819-43aa-9c62-94df8657ac45
# ╟─71d10ec2-27fc-457a-9d65-4a436a57d714
# ╠═2c4db510-1d7f-4d1a-82e4-936f994bfd46
# ╠═13a16316-c6a5-4135-bbf2-b1e7b1a6caae
# ╠═47869b79-57ab-4c02-a1e7-aad9ffb0fa49
# ╠═b0ef3f79-f6dc-462a-8ea9-f36dd9f1f583
# ╠═2afcebe1-7ca7-4f01-af35-48a183cfa8ea
# ╠═916bedc1-18e2-4a7b-b275-c202543600bf
# ╠═ea7b2b0a-0d8b-43a6-a6ee-e8fd5f4a6b30
# ╠═3da8cbf1-c7bf-4adf-8cce-b1a2f54eebf9
# ╠═082d2074-fed5-4b9e-8dab-8c648cbbd25a
# ╠═af611a51-a939-4ac3-949d-6b762eac94b9
# ╠═090048c4-19d9-4a5a-9bd4-6be6d34b8351
# ╠═376451b5-395e-44bb-a75e-6424f6f56375
# ╠═ea28abf7-61de-497b-8544-3f584adb35ac
# ╠═089c5d73-5dbb-495f-8766-a88ad25eb41f
# ╠═01e66054-1398-4878-bd36-95d7e7ad0923
# ╟─9d187387-ac73-40e2-b594-e77eaf29ef09
# ╟─ef36fdc6-1900-4437-852c-da95c1869b87
# ╟─f9231ee4-9e0e-408d-a35d-de0aac7e74dd
# ╟─3fd6cf00-87a3-495f-88a0-f9e2abec657e
# ╟─1edf5403-af4f-4df4-9bde-c8543e9a96ef
# ╟─f016d916-0ce8-4e59-80bb-1fdeff754245
# ╟─ce41903e-6bd2-4235-8895-a76b5592ce15
# ╟─efd358d3-e6c1-4abd-9a00-177f38583e67
# ╟─afc91517-d6ee-44aa-b67e-1fd252192ca9
# ╠═c7abe550-993f-4684-a02d-883af41e8c1e
# ╠═f9f2b9e3-921e-43a5-abe8-3d374716e5cb
# ╟─2763c7fa-ef74-4417-a1a5-f833753dd783
# ╠═b74db8bb-e9d9-4609-aa86-8a683ee8ae0b
# ╠═794740ff-310a-40a1-b8af-ae9f0c3c46ff
