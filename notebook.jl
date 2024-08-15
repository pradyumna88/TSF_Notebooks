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

# ╔═╡ 13a16316-c6a5-4135-bbf2-b1e7b1a6caae
prefix_path = "../../Cluster/curie/rough_runs_curie/"

# ╔═╡ 47869b79-57ab-4c02-a1e7-aad9ffb0fa49
foldernames = prefix_path |> Results.foldernames

# ╔═╡ b0ef3f79-f6dc-462a-8ea9-f36dd9f1f583
filenames = prefix_path |> Results.filenames

# ╔═╡ 2afcebe1-7ca7-4f01-af35-48a183cfa8ea
filepaths = prefix_path |> Results.filepaths

# ╔═╡ 916bedc1-18e2-4a7b-b275-c202543600bf
@bind phys_params Select(filepaths |> Results.get_phys_params)

# ╔═╡ 7d10fc71-459b-4504-8cd2-14d01d1bb0ef
@bind run_version Select(filepaths |> x -> Results.get_run_version(x, phys_params))

# ╔═╡ 9e313fc0-a9f9-4473-a369-fc9ff6589cb8
@bind sim_params Select(filepaths |> x -> Results.get_sim_params(x, phys_params, run_version))

# ╔═╡ 786911f2-719b-48d6-9a86-51f296ea5642
@bind param_C Select(filepaths |> x -> Results.get_Cs(x, phys_params, run_version))

# ╔═╡ 63d19ba1-0160-498c-8d64-ba120d27e7ef
@bind param_Mbar Select(filepaths |> x -> Results.get_Mbars(x, phys_params, run_version))

# ╔═╡ b65486c7-abca-43e9-8e47-e2dc70c260c3
filepaths |> filter(x -> Results.get_parameter(x)[[1:4; 7; 6]] == (phys_params..., run_version, param_Mbar)) |> mymap(x -> load(x, "counts", typemap = Dict("Main.Results.Counts" => Counts))) |> rand 

# ╔═╡ Cell order:
# ╟─92fd972e-5a3e-11ef-27c3-cfcf217458f8
# ╟─db4e7db6-7f42-4da5-aec1-5e2ee05fd7ae
# ╠═ccca791d-b819-43aa-9c62-94df8657ac45
# ╟─71d10ec2-27fc-457a-9d65-4a436a57d714
# ╟─13a16316-c6a5-4135-bbf2-b1e7b1a6caae
# ╠═47869b79-57ab-4c02-a1e7-aad9ffb0fa49
# ╠═b0ef3f79-f6dc-462a-8ea9-f36dd9f1f583
# ╠═2afcebe1-7ca7-4f01-af35-48a183cfa8ea
# ╠═916bedc1-18e2-4a7b-b275-c202543600bf
# ╠═7d10fc71-459b-4504-8cd2-14d01d1bb0ef
# ╠═9e313fc0-a9f9-4473-a369-fc9ff6589cb8
# ╠═786911f2-719b-48d6-9a86-51f296ea5642
# ╠═63d19ba1-0160-498c-8d64-ba120d27e7ef
# ╠═b65486c7-abca-43e9-8e47-e2dc70c260c3
