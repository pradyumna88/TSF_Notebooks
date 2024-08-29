### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try
            Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value
        catch
            b -> missing
        end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 7fee701b-2d8e-4b84-a588-f2d751be0b6f
# NOTE: html to use full screen
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

# ╔═╡ 34a1d51b-b7cc-49f7-be9f-cb5128599f81
# NOTE: Packages
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

# ╔═╡ e71d3946-8f41-4d7d-bdd4-920571683a9c
# NOTE: prefix path and filepaths
begin
    prefix_path = "NewRuns/"
    filepaths = prefix_path |> Results.filepaths
end;

# ╔═╡ e2cb5371-dabe-4ab9-98fa-97fd1dd3ab55
# NOTE: phys_params
begin
    phys_params_selector = @bind phys_params Select(filepaths |> Results.get_phys_params)
    ignore_selector = @bind ignore Select([:bare, :ignore_zeros, :ignore_non_trimer])
end;

# ╔═╡ 2d9fd99a-d590-406f-9894-363fa5b025b9
# NOTE: run_version selector
begin
    run_version_selector = @bind run_version Select(filepaths |> x -> Results.get_run_version(x, phys_params))
end;

# ╔═╡ f3033ccd-aedb-4196-b08a-c6d20f74dd80
# NOTE: rankers and sim_params_selector
begin
    ranker(rev=false) = x -> x |> (x -> sortperm(x, rev=rev) .=> 1:length(x)) |> x -> sort(x, by=x -> x[1]) |> mymap(last)
    ranker1 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5, 6], (phys_params..., run_version)) |>
                                mymap(x -> load(x, "superfluid")) |> mymap(x -> x.ignore_zeros.plus[1]) |> ranker(true)
    ranker2 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5, 6], (phys_params..., run_version)) |>
                                mymap(x -> load(x, "superfluid_winding_number_frequency")) |> mymap(x -> x.ignore_non_trimer.plus) |> mymap(length) |> ranker(true)
    ranker3 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5, 6], (phys_params..., run_version)) |>
                                mymap(x -> load(x, "superfluid_winding_number_frequency")) |> mymap(x -> x.ignore_non_trimer.plus) |>
                                mymap(mymap(x -> x[1] == 0 ? 0 : x[1] < 0 ? -x[2] : x[2])) |> mymap(x -> x / sum(abs.(x))) |> mymap(mean) |> mymap(abs) |> ranker()
    ranker4 = filepaths |> x -> Plotter.filter_datafilepaths(x, [5, 6], (phys_params..., run_version)) |>
                                mymap(x -> load(x, "superfluid_convergence").ignore_non_trimer.plus) |> mymap(mymap(first)) |> mymap(x -> x[end] - x[Int64(ceil(0.8 * length(x)))]) |> ranker()
    multi_ranker = (ranker1,) |> mymap(x -> 1 ./ (x .+ 10)) |> sum |> ranker(true)
    sim_params_list = filepaths |> x -> Results.get_sim_params(x, phys_params, run_version) |> x -> x[sortperm(multi_ranker)]
    sim_params_selector = @bind sim_params Select(sim_params_list)
end;

# ╔═╡ 10d11908-cfdc-4b39-95ea-999bd28a5566
# NOTE: run_time and current_meas
begin
    run_time = Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> x -> load(x, "run_time") |> x -> x / (60 * 60 * 24) |> x -> round(x, digits=2)
    current_meas = Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end] |> x -> load(x, "current_meas")
end;


# ╔═╡ 256e8e23-feeb-42eb-b2ec-ad798803d6f4
# NOTE: show the selectors
md""" $phys_params_selector $run_version_selector $ignore_selector"""


# ╔═╡ 582dbc1f-0623-40e8-a891-639803bd3a19
md""" $sim_params_selector """


# ╔═╡ 11848192-9459-41da-a69a-935f8121b230
# NOTE: plot convergence, histograms
begin
    filepath = Results.filepaths(prefix_path) |> x -> Plotter.filter_datafilepaths(x, 8, (phys_params..., sim_params..., run_version))[end]
    plot(
        filepath |>
        x -> Plotter.get_plotting_data(x, x -> nothing, x -> load(x, "superfluid_convergence") |> x -> getproperty(x, ignore) |> x -> (x.plus, x.minus)) |>
             x -> plot(x.xs, x.ys * phys_params[2] * 100 * 3 / (2 * phys_params[3]) |> transpose, yerr=x.yers * phys_params[2] * 100 * 3 / (2 * phys_params[3]) |> transpose, title=(run_time, current_meas / 10^4)) |>
                  x -> hline!(x, [100, 2 * phys_params[2] * 3 * 100 / (pi * phys_params[3])]),
        filepath |>
        x -> Plotter.get_plotting_data(x, x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).plus |> mymap(first) |> x -> x * 4,
            x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).plus |> mymap(last), yerr_flag=false) |>
             x -> bar(x.xs, x.ys ./ sum(x.ys), bar_width=0.25),
        filepath |>
        x -> Plotter.get_plotting_data(x, x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).minus |> mymap(first),
            x -> load(x, "superfluid_winding_number_frequency") |> x -> getproperty(x, ignore).minus |> mymap(last), yerr_flag=false) |>
             x -> bar(x.xs, x.ys ./ sum(x.ys), bar_width=0.25, c=:red),
        # filepath |> x -> load(x, "winding_numbers") |> x -> 2 * @view(x[:, :, 1]) + @view(x[:, :, 2]) |> x -> reshape(x, :) |> x -> x[1:10:end] |> x -> plot(range(1, 100, 10000), x[Int64.(ceil.(length(x) * range(0.01, 1, 10000)))]) |>
        #                                                                                                                                                 x -> hline!(x, 4 * (-2:2), label=false),
        layout=(2, 2), size=(1000, 800), margin=9Plots.mm
    )
end


# ╔═╡ ec229c11-74bf-4d70-945f-d7fe13b62e4e
begin
    param_C_selector = @bind param_C Select(filepaths |> x -> Results.get_Cs(x, phys_params, run_version))
    param_Mbar_selector = @bind param_Mbar Select(filepaths |> x -> Results.get_Mbars(x, phys_params, run_version))
end;


# ╔═╡ b31fa4ab-1fca-40fd-8d41-d98235abf512
begin
    md""" $param_Mbar_selector """
end

# ╔═╡ cbc03d9e-42ad-45f1-8b7d-e3acee7f9b8b
begin
    filepaths |> x -> Plotter.filter_datafilepaths(x, 5, (phys_params..., param_Mbar, run_version)) |> (x -> convert(Array{String}, x)) |>
                      (x -> Plotter.get_plotting_data(x, 5, x -> load(x, "superfluid").ignore_zeros.plus)) |>
                      (x -> plot(x.xs, x.ys |> x -> x .* (100 * 3 * 2 * phys_params[2] / (pi * phys_params[3])),
                          yerr=x.yers |> x -> x .* (100 * 3 * 2 * phys_params[2] / (pi * phys_params[3]))))
end

# ╔═╡ Cell order:
# ╟─7fee701b-2d8e-4b84-a588-f2d751be0b6f
# ╟─34a1d51b-b7cc-49f7-be9f-cb5128599f81
# ╟─e71d3946-8f41-4d7d-bdd4-920571683a9c
# ╟─e2cb5371-dabe-4ab9-98fa-97fd1dd3ab55
# ╟─f3033ccd-aedb-4196-b08a-c6d20f74dd80
# ╟─2d9fd99a-d590-406f-9894-363fa5b025b9
# ╟─10d11908-cfdc-4b39-95ea-999bd28a5566
# ╟─256e8e23-feeb-42eb-b2ec-ad798803d6f4
# ╟─582dbc1f-0623-40e8-a891-639803bd3a19
# ╟─11848192-9459-41da-a69a-935f8121b230
# ╟─ec229c11-74bf-4d70-945f-d7fe13b62e4e
# ╟─b31fa4ab-1fca-40fd-8d41-d98235abf512
# ╟─cbc03d9e-42ad-45f1-8b7d-e3acee7f9b8b
