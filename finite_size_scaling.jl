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

# ╔═╡ 2e1fc825-7fad-47f7-967e-dae00c2595fa
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
"""


# ╔═╡ b011d064-a39f-4ec1-bacf-983d6cc0f116
# NOTE: include relevant packages
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
  using LsqFit
end;

# ╔═╡ 1e6efcc4-860b-4a27-9f6a-8b3ff26d230a
# NOTE: custom functions: myprint, myfilter, mypartition
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
  mymap2(f) = mymap(x -> x[1] => f(x[2]))
end;

# ╔═╡ def69dd5-35c6-490c-a273-818edbba3fb0
# NOTE: 
# filepaths : array of strings 
# parameters : 4 tuple of params => pair of filepaths
# rhoplus : 4 tuple => 2 tuple of val and err
# rhominus : 4 tuple => 2 tuple of val and err
begin
  filepaths = Results.filepaths_old("Runs")
  parameters = Results.get_parameter_old.(filepaths) .=> filepaths
  rhoplus = parameters |> mymap(
    x -> x[1] => load(x[2], "superfluid").bare.plus .* 100 .* 3 .* x[1][2] ./ (2 .* x[1][3])
  )
  rhominus = parameters |> mymap(
    x -> x[1] => load(x[2], "superfluid").bare.minus .* 100 .* 3 .* x[1][2] ./ (2 .* x[1][3])
  )
end;

# ╔═╡ f947b6b7-4797-48e2-b602-3125004a2e71
# NOTE:
# partitioned_vals : function to partiion rhoplus/minus by n1, then T, then N1
# run_avg_vals : function to partition, then combine all runs, then average over runs, 
# fss_vals : funmction to finite size scale upto N1 particles, then flatten the partition
#        pair N1, T, n1 => val, err
begin
  partitioned_vals(x) = mypartition(x -> x[1][3], x) |>
                        mymap(x -> mypartition(x -> x[1][2], x)) |>
                        mymap(mymap(x -> mypartition(x -> x[1][1], x)))
  run_avg_vals(x) = x |> partitioned_vals |> mymap(mymap(mymap(
                      x -> x[1][1][1:end-1] => x |> mymap(last)
                    ))) |>
                    mymap(mymap(mymap(
                      x -> x[1] => x[2] |> (x -> tuple(x...)) ∘ (x -> mean(x, dims=1)) ∘ transpose ∘ stack
                    )))
  fss_vals(x) = x |> run_avg_vals |>
                mymap(mymap(
                  x -> x[1][1][2:end] => x |> mymap(x -> x[1][1] => x[2]) |> fss_all
                )) |>
                mymap(mymap(
                  x -> x[2] |> mymap(y -> (y[1], x[1]...) => y[2])
                )) |> x -> vcat(x...) |> x -> vcat(x...)
end;

# ╔═╡ b4b6a281-1d26-4844-a139-b81a08459c10
# NOTE:
# fss : function to finite size scale taking a N1 => 2 tuple of val, err
# fss_all : function to fss upto N1 particles
begin
  function fss(ps)
    if length(ps) == 1
      return ps[1][2]
    end
    ps = sort(ps, by=x -> x[1])
    N1s = ps |> mymap(first)
    ys = ps |> mymap(last) |> mymap(first)
    yers = ps |> mymap(last) |> mymap(last)
    yhighs = ys + yers
    ylows = ys - yers
    model(x, p) = p[1] .* x .+ p[2]
    y = curve_fit(model, 1 ./ N1s, ys, rand(2)).param[2] |> x -> max(x, 0)
    yh = curve_fit(model, 1 ./ N1s, yhighs, rand(2)).param[2] |> x -> max(x, 0)
    yl = curve_fit(model, 1 ./ N1s, ylows, rand(2)).param[2] |> x -> max(x, 0)
    return (y, abs.(yh - yl) / 2) |> x -> round.(x, digits=2)
  end
  function fss_all(ps)
    ps = sort(ps, by=x -> x[1])
    N1s = ps |> mymap(first)
    new_ps = []
    for i in eachindex(N1s)
      push!(new_ps, N1s[i] => fss(ps[1:i]))
    end
    return new_ps
  end
  # tmp[5][4][2] |> fss_all
end;

# ╔═╡ 6b8e4605-e19b-46e0-91d4-50b9605ac148
# NOTE: rhoplus_vals and rhominus_vals -> after application of fss_vals
begin
  rhoplus_run_avg = run_avg_vals(rhoplus) |> x -> vcat(x...) |> x -> vcat(x...)
  rhominus_run_avg = run_avg_vals(rhominus) |> x -> vcat(x...) |> x -> vcat(x...)
  rhoplus_vals = fss_vals(rhoplus)
  rhominus_vals = fss_vals(rhominus)
end

# ╔═╡ feeb078d-e923-42c9-a200-8b003a7e368a
# NOTE: n1 selector with default value of 0.71
begin
  n1 = 0.71
  n1selector = @bind n1 (rhoplus_vals |> mymap(x -> x[1][3]) |> unique |> sort |>
                         x -> Select(x, default=n1))
end;

# ╔═╡ c56c0b07-e968-406b-960b-083940016d6b
# NOTE: N1 selector with default value of 6
begin
  N1 = 4
  N1selector = @bind N1 (rhoplus_vals |> mymap(x -> x[1][1]) |> unique |> sort |>
                         x -> Select(x, default=N1))
end;



# ╔═╡ a6c7fc1a-587c-4551-9eee-c24e3b89063a
md"""
N1: $N1selector

n1: $n1selector
"""

# ╔═╡ 739253b7-fd27-4152-92a7-4314dc9c3616
# NOTE:
begin
  fig1 = plot(title="finite size scaled")
  rhoplus_vals |> filter(x -> x[1][[1, 3]] == (N1, n1)) |> mymap(x -> (x[1][2], x[2]...)) |>
  stack |> transpose |> (x -> plot!(
    fig1, x[:, 1], x[:, 2], yerr=x[:, 3], xlim=[0, 2], label="ρ₊"
  ))
  rhominus_vals |> filter(x -> x[1][[1, 3]] == (N1, n1)) |> mymap(x -> (x[1][2], x[2]...)) |>
  stack |> transpose |> (x -> plot!(
    fig1, x[:, 1], x[:, 2], yerr=x[:, 3], xlim=[0, 2], label="ρ₋"
  ))
  plot!(fig1, range(0, 0.8, 10), x -> 2x * 100 * 3 / (n1 * π), label=false)

  fig2plus = plot()
  rhoplus_run_avg |> filter(x -> x[1][1] <= N1 && x[1][3] == n1) |> mymap(x -> (x[1][1:2]..., x[2]...)) |>
  (x -> mypartition(x -> x[1], x)) |> mymap(stack) |> mymap(transpose) |> x -> foldl(
    (fg, x) -> plot!(
      fg, x[:, 2], x[:, 3], yerr=x[:, 4], xlim=[0, 2], label="N1 = $(x[1, 1])", title="ρ₊"
    ), x, init=fig2plus
  )
  plot!(fig2plus, range(0, 0.8, 10), x -> 2x * 100 * 3 / (n1 * π), label=false)

  fig2minus = plot()
  rhominus_run_avg |> filter(x -> x[1][1] <= N1 && x[1][3] == n1) |> mymap(x -> (x[1][1:2]..., x[2]...)) |>
  (x -> mypartition(x -> x[1], x)) |> mymap(stack) |> mymap(transpose) |> x -> foldl(
    (fg, x) -> plot!(
      fg, x[:, 2], x[:, 3], yerr=x[:, 4], xlim=[0, 2], label="N1 = $(x[1, 1])", title="ρ₋"
    ), x, init=fig2minus
  )
  plot!(fig2minus, range(0, 0.8, 10), x -> 2x * 100 * 3 / (n1 * π), label=false)

  fig3 = plot()
  rhoplus_run_avg |> filter(x -> x[1][[1, 3]] == (N1, n1)) |> mymap(x -> (x[1][2], x[2]...)) |>
  stack |> transpose |> (x -> plot!(
    fig3, x[:, 1], x[:, 2], yerr=x[:, 3], xlim=[0, 2], label="ρ₊", title="N1 = $N1"
  ))
  rhominus_run_avg |> filter(x -> x[1][[1, 3]] == (N1, n1)) |> mymap(x -> (x[1][2], x[2]...)) |>
  stack |> transpose |> (x -> plot!(
    fig3, x[:, 1], x[:, 2], yerr=x[:, 3], xlim=[0, 2], label="ρ₋"
  ))
  plot!(fig3, range(0, 0.8, 10), x -> 2x * 100 * 3 / (n1 * π), label=false)

  plot(fig1, fig3, fig2plus, fig2minus, layout=(2, 2), size=(1000, 800), xlabel="T", margin=5Plots.mm)




end

# ╔═╡ Cell order:
# ╟─2e1fc825-7fad-47f7-967e-dae00c2595fa
# ╟─b011d064-a39f-4ec1-bacf-983d6cc0f116
# ╟─1e6efcc4-860b-4a27-9f6a-8b3ff26d230a
# ╟─def69dd5-35c6-490c-a273-818edbba3fb0
# ╟─f947b6b7-4797-48e2-b602-3125004a2e71
# ╟─b4b6a281-1d26-4844-a139-b81a08459c10
# ╟─6b8e4605-e19b-46e0-91d4-50b9605ac148
# ╟─feeb078d-e923-42c9-a200-8b003a7e368a
# ╟─c56c0b07-e968-406b-960b-083940016d6b
# ╟─a6c7fc1a-587c-4551-9eee-c24e3b89063a
# ╟─739253b7-fd27-4152-92a7-4314dc9c3616
