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
end;

# ╔═╡ a6a04ae7-1925-4cbe-86a3-fe6ccf39d657
# NOTE: make cells wider
html"""
<style>
	main {
		margin: 0 auto;
		max-width: 2200px;
    	padding-left: max(50px, 1%);
    	padding-right: max(25px, 1%);
	}
</style>
<style>
  pluto-helpbox {
    display: none;
  }
</style>
"""

# ╔═╡ 45865586-8c02-4698-ae29-41fd3981197e
# NOTE: Packages required
begin
  using Revise
  using Pkg
  Pkg.activate("TSF_110")
  using DrWatson
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
  using LsqFit
  using PrettyTables
end;



# ╔═╡ ff8cbd1e-19a2-4524-9a01-8c195126df35
# NOTE: refresh button for filepaths
begin
  refresh = ""
  refresh_button = @bind refresh Button("refresh")
end;

# ╔═╡ 77394bb8-9b65-4153-9bb3-9c8c846a3815
# NOTE: refresh button to trigger re-evaluation of cell
# filepaths_old/new : array containing path to the results
# parameters_old/new : array of pairs from parameters to filepaths
begin
  refresh
  filepaths_old = Results.filepaths_old("Runs")
  filepaths_new = Results.filepaths_new("NewRuns")
  parameters_old = Results.get_parameter_old.(filepaths_old) .=> filepaths_old
  parameters_new = Results.get_parameter.(filepaths_new) .=> filepaths_new
end;

# ╔═╡ c3f7756d-65ef-41ac-8b6a-a61cf950b33f
# NOTE:
begin
  Results.filepaths_new("NewRuns/") |> mymap(Results.get_parameter)
end

# ╔═╡ d96bcb93-62cc-467b-8884-180757f88a30
# NOTE: myfunctions: print, filter, partition
begin
  function myprint(arr)
    arr |> x -> join(x, '\n') |> println
  end
  function myfilter(str)
    filter(x -> occursin(str, x))
  end
  function mypartition(f, arr)
    part_keys = arr |> mymap(f) |> unique
    try
      sort!(part_keys)
    catch
    end
    part_keys |> mymap(
      x -> filter(y -> f(y) == x, arr)
    )
  end
  function mypretty_table(header="header")
    x -> x |> keys |> collect |> (y -> pretty_table(x, header=[header]))
  end
  myload(x) = load(x, "data")
  mysplit(x) = y -> split(y, x)
  myread(x) = read(x, String)
  myreplace(x) = y -> replace(y, x)
end;

# ╔═╡ f479f3d8-fe68-4ac9-974f-c19e4f017685
# NOTE: N1/T/n1/runversion_old/new_lists : corresponding arrays unique and sorted
begin
  N1_old_list = parameters_old |> mymap(first) |> mymap(x -> x[1]) |> unique |> sort
  T_old_list = parameters_old |> mymap(first) |> mymap(x -> x[2]) |> unique |> sort
  n1_old_list = parameters_old |> mymap(first) |> mymap(x -> x[3]) |> unique |> sort
  runversion_old_list = parameters_old |> mymap(first) |> mymap(x -> x[4]) |> unique |> sort
  N1_new_list = parameters_new |> mymap(first) |> mymap(x -> x[1]) |> unique |> sort
  T_new_list = parameters_new |> mymap(first) |> mymap(x -> x[2]) |> unique |> sort
  n1_new_list = parameters_new |> mymap(first) |> mymap(x -> x[3]) |> unique |> sort
  CMbar_new_list = parameters_new |> mymap(first) |> mymap(x -> x[5:6]) |> unique |> sort
  runversion_new_list = parameters_new |> mymap(first) |> mymap(x -> x[7]) |> unique |> sort
end;

# ╔═╡ 72caad3c-1798-4755-84ea-224c9b9f412b
# NOTE: save vardicts/new_old_compare/(N1/T/n1/runversion).jld2 
begin
  N1_default, T_default, n1_default, runversion_default = parameters_old |> first |> first
  save("vardicts/new_old_compare/N1.jld2", Dict("data" => N1_default))
  save("vardicts/new_old_compare/T.jld2", Dict("data" => T_default))
  save("vardicts/new_old_compare/n1.jld2", Dict("data" => n1_default))
  save("vardicts/new_old_compare/runversion.jld2", Dict("data" => runversion_default))
end;

# ╔═╡ 2779d4d1-bb0c-4f32-bf93-daaf03fad042
# NOTE: load n1 and make n1_selector with loaded n1 as default if available
begin
  n1 = myload("vardicts/new_old_compare/n1.jld2")
  n1_selector = @bind n1 parameters_old |> mymap(first) |> mymap(x -> x[3]) |>
                         unique |> sort |> (x -> Select(x, default=n1 in x ? n1 : x[1]))
end;

# ╔═╡ b004d4f3-f9e4-46fa-8aa3-772e7cc6ac09
# NOTE: load T and make T_selector based on n1 with loaded T as default if available
begin
  T = myload("vardicts/new_old_compare/T.jld2")
  T_selector = @bind T parameters_old |> mymap(first) |> filter(x -> x[3] == n1) |>
                       mymap(x -> x[2]) |> unique |> sort |> (x -> Select(x, default=T in x ? T : x[1]))
end;

# ╔═╡ dd9165f4-8af8-4b81-b96e-ecdfde94ddb1
# NOTE: load N1 and make N1_selector based on n1 and T with loaded N1 as default if available
begin
  N1 = myload("vardicts/new_old_compare/N1.jld2")
  N1_selector = @bind N1 parameters_old |> mymap(first) |> filter(x -> x[2:3] == (T, n1)) |>
                         mymap(x -> x[1]) |> unique |> sort |> (x -> Select(x, default=N1 in x ? N1 : x[1]))
end;

# ╔═╡ fd352f67-0bb8-49f3-b0d9-4d5aea31f2ca
# NOTE: load runversion and make runversion_selector based on n1, T and N1 with loaded runversion as default
# if available
begin
  runversion = myload("vardicts/new_old_compare/runversion.jld2")
  runversion_selector = @bind runversion parameters_old |> mymap(first) |>
                                         filter(x -> x[1:3] == (N1, T, n1)) |> mymap(x -> x[end]) |>
                                         unique |> sort |>
                                         (x -> Select(x, default=runversion in x ? runversion : x[1]))
end;

# ╔═╡ 5e537877-b1a0-4b4f-b321-61a9bb56d85e
# NOTE: save N1, T, n1, and runversion to be used as defaults
begin
  save("vardicts/new_old_compare/N1.jld2", Dict("data" => N1))
  save("vardicts/new_old_compare/T.jld2", Dict("data" => T))
  save("vardicts/new_old_compare/n1.jld2", Dict("data" => n1))
  save("vardicts/new_old_compare/runversion.jld2", Dict("data" => runversion))
end;

# ╔═╡ b1648906-dbec-4d4b-82eb-f1aa6d4c9f38
# NOTE: factor: 100*(3/n1). also runtime_old and measpercent_old 
# rhoconv_old : rho convergence for old data, both plus and minus
# windingplus/minus_old : histograms for old data
begin
  factor = 100 * 3 / n1
  runtime_old = parameters_old |> filter(x -> x[1] == (N1, T, n1, runversion)) |> first |>
                (x -> load(x[2], "run_time")) |> (x -> x / (60 * 60 * 24)) |> (x -> round(x, digits=2))
  measpercent_old = parameters_old |> filter(x -> x[1] == (N1, T, n1, runversion)) |> first |>
                    (x -> load(x[2], "current_meas"))
  rhoconv_old = parameters_old |> filter(x -> x[1] == (N1, T, n1, runversion)) |> first |>
                (x -> load(x[2], "superfluid_convergence")) |> (x -> (x.bare.plus, x.bare.minus)) |>
                stack |> mymap(x -> x .* (T / 2) .* factor) |>
                (x -> plot(x |> mymap(first), yerr=x |> mymap(last), label=["ρ₊" "ρ₋"])) |>
                (x -> hline!(x, [100, factor * 2T / π], label=false)) |>
                (x -> plot(x, xlabel="measure%", ylabel="ρ%", legendfontsize=10))
  windingplus_old = parameters_old |> filter(x -> x[1] == (N1, T, n1, runversion)) |> first |>
                    (x -> load(x[2], "superfluid_winding_number_frequency").bare.plus) |>
                    (x -> bar(
                      x |> mymap(first) |> (x -> x .* 4),
                      x |> mymap(last) |> (x -> x ./ sum(x)),
                      bar_width=0.25, xlabel="winding numbers", ylabel="probability",
                      label="W₊", legendfontsize=10
                    ))
  windingminus_old = parameters_old |> filter(x -> x[1] == (N1, T, n1, runversion)) |> first |>
                     (x -> load(x[2], "superfluid_winding_number_frequency").bare.minus) |>
                     (x -> bar(
                       x |> mymap(first),
                       x |> mymap(last) |> (x -> x ./ sum(x)),
                       bar_width=0.25, xlabel="winding numbers", ylabel="probability",
                       label="W₋", legendfontsize=10
                     ))
end;

# ╔═╡ 7ca0e64d-0ed2-4250-b959-cd1138d3fc5f
# NOTE: CMbars_list : for given N1, T, n1 and runversion; empty if no new runs 
# CMbars_selector only if CMbars_list is non mepty, CMbar is (0,0) if no new runs
# make runsim_checkbox if no new runs
# make CMbarRange_textfield if no new runs
begin
  CMbars_list = try
    parameters_new |> filter(x -> x[1][[1, 2, 3, end]] == (N1, T, n1, runversion)) |>
    (x -> sort(x, by=x -> load(x[2], "superfluid").bare.plus[1], rev=true)) |>
    mymap(x -> x[1][5:6])
  catch
    []
  end
  if !isempty(CMbars_list)
    CMbar = CMbars_list |> first
    CMbars_selector = @bind CMbar Select(CMbars_list)
  else
    CMbar = (0, 0)
    runsim_cb = false
    runsim_checkbox = @bind runsim_cb CheckBox(default=false)
    CMbarRange = ""
    CMbarRange_textfield = @bind CMbarRange TextField(default="0, 9, 10, 1, 2, 2")
  end
end;

# ╔═╡ 88321f27-b09c-4531-87fc-7100016b3d5c
# NOTE: runtime_new, measpercent_new, rhoconv_new and windingplus/minus_new analogous to old if new runs exist
begin
  if CMbar != (0, 0)
    runtime_new = parameters_new |> filter(x -> x[1] == (N1, T, n1, 0.1, CMbar..., runversion)) |> first |>
                  (x -> load(x[2], "run_time")) |> (x -> x / (60 * 60 * 24)) |> (x -> round(x, digits=2))
    measpercent_new = parameters_new |> filter(x -> x[1] == (N1, T, n1, 0.1, CMbar..., runversion)) |> first |>
                      (x -> load(x[2], "current_meas"))
    lastmeas = "Logs/log" |> readlines |> mymap(mysplit('\t')) |>
               filter(x -> x[1] == (join([N1, T, n1, 0.1, runversion], '_') |> myreplace("0." => "."))) |>
               filter(x -> x[3] == string(CMbar)) |> first |> (x -> x[5])
    rhoconv_new = parameters_new |> filter(x -> x[1] == (N1, T, n1, 0.1, CMbar..., runversion)) |>
                  first |> (x -> load(x[2], "superfluid_convergence")) |> (x -> (x.bare.plus, x.bare.minus)) |>
                  stack |> mymap(x -> x .* (T / 2) .* factor) |>
                  (x -> plot(x |> mymap(first), yerr=x |> mymap(last), label=["ρ₊" "ρ₋"])) |>
                  (x -> hline!(x, [100, factor * 2T / π], label=false)) |>
                  (x -> plot(x, xlabel="measure%", ylabel="ρ%", legendfontsize=10))
    windingplus_new = parameters_new |> filter(x -> x[1] == (N1, T, n1, 0.1, CMbar..., runversion)) |>
                      first |> (x -> load(x[2], "superfluid_winding_number_frequency").bare.plus) |>
                      (x -> bar(
                        x |> mymap(first) |> (x -> x .* 4),
                        x |> mymap(last) |> (x -> x ./ sum(x)),
                        bar_width=0.25, xlabel="winding numbers", ylabel="probability",
                        label="W₊", legendfontsize=10
                      ))
    windingminus_new = parameters_new |> filter(x -> x[1] == (N1, T, n1, 0.1, CMbar..., runversion)) |>
                       first |> (x -> load(x[2], "superfluid_winding_number_frequency").bare.minus) |>
                       (x -> bar(
                         x |> mymap(first),
                         x |> mymap(last) |> (x -> x ./ sum(x)),
                         bar_width=0.25, xlabel="winding numbers", ylabel="probability",
                         label="W₋", legendfontsize=10
                       ))
  end
end;

# ╔═╡ 0577ae28-e3fc-4a5b-a1a1-74ecbab78929
# NOTE: nbsp (non breaking space) and display n1, T, N1, runversion and Cmbars if available
begin
  if isempty(CMbars_list)
    nbsp = html"&nbsp"
    md"""
    n1 = $n1_selector $nbsp T = $T_selector $nbsp N1 = $N1_selector $nbsp $runversion_selector
    $nbsp $runsim_checkbox $nbsp $CMbarRange_textfield $nbsp $refresh_button
    """
  else
    nbsp = html"&nbsp"
    md"""
    n1 = $n1_selector $nbsp T = $T_selector $nbsp N1 = $N1_selector $nbsp $runversion_selector 
    $nbsp $CMbars_selector $nbsp $refresh_button
    """
  end
end

# ╔═╡ 9498041f-546b-4a3a-8f19-030c88f58c51
# NOTE: combine rhoconv, windingplus/minus plots and new ones if available
begin
  old_data = plot(
    rhoconv_old, windingplus_old, windingminus_old,
    plot_title="OLD DATA; $runtime_old days, $measpercent_old measurements",
    layout=(1, 3), size=(2000, 500), margin=10Plots.mm,
  )
  if CMbar != (0, 0)
    new_data = plot(
      rhoconv_new, windingplus_new, windingminus_new,
      plot_title="new DATA; $runtime_new days, $measpercent_new measurements; last meas at $lastmeas",
      layout=(1, 3), size=(2000, 500), margin=10Plots.mm,
    )
  end
  if CMbar != (0, 0)
    plot(old_data, new_data, layout=(2, 1), size=(2000, 1000))
  else
    old_data
  end
end


# ╔═╡ afca303a-2f81-4d58-b4a7-87712a248c8d
# NOTE:
begin
  if isempty(CMbars_list)
    if runsim_cb
      Cstart, Cstop, Clength, Mstart, Mstop, Mlength = CMbarRange |> x -> split(x, ',')
      @async run(
        `ssh curie "/usr/people/snirgaz/prady/data/PD/Simulations/NewRuns/Slurm/run_scripts/run.sh $N1 $T $n1 $runversion $Cstart $Cstop $Clength $Mstart $Mstop $Mlength"`
      )
    end
  end
end

# ╔═╡ 4623dd78-75a9-49f3-937f-cbaa1d2b8f89
# NOTE: 
begin
  parameters_new |> rand |> last |> load |> keys |> collect |>
  (x -> pretty_table(x; header=["result keys"]))
end

# ╔═╡ Cell order:
# ╟─a6a04ae7-1925-4cbe-86a3-fe6ccf39d657
# ╟─45865586-8c02-4698-ae29-41fd3981197e
# ╟─ff8cbd1e-19a2-4524-9a01-8c195126df35
# ╟─77394bb8-9b65-4153-9bb3-9c8c846a3815
# ╠═c3f7756d-65ef-41ac-8b6a-a61cf950b33f
# ╟─d96bcb93-62cc-467b-8884-180757f88a30
# ╟─f479f3d8-fe68-4ac9-974f-c19e4f017685
# ╟─72caad3c-1798-4755-84ea-224c9b9f412b
# ╟─2779d4d1-bb0c-4f32-bf93-daaf03fad042
# ╟─b004d4f3-f9e4-46fa-8aa3-772e7cc6ac09
# ╟─dd9165f4-8af8-4b81-b96e-ecdfde94ddb1
# ╟─fd352f67-0bb8-49f3-b0d9-4d5aea31f2ca
# ╟─5e537877-b1a0-4b4f-b321-61a9bb56d85e
# ╟─b1648906-dbec-4d4b-82eb-f1aa6d4c9f38
# ╟─7ca0e64d-0ed2-4250-b959-cd1138d3fc5f
# ╟─88321f27-b09c-4531-87fc-7100016b3d5c
# ╟─0577ae28-e3fc-4a5b-a1a1-74ecbab78929
# ╟─9498041f-546b-4a3a-8f19-030c88f58c51
# ╟─afca303a-2f81-4d58-b4a7-87712a248c8d
# ╟─4623dd78-75a9-49f3-937f-cbaa1d2b8f89
