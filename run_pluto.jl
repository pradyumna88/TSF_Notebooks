using Pluto

if length(ARGS) == 2
  notebook = ARGS[2] |> x -> split(x, '/', limit=2)[end] |>
                             x -> x == "jl" ? ARGS[2] : readdir(ARGS[2], join=true)
else
  notebook = empty(["a"])
end


Pluto.run(port=parse(Int64, ARGS[1]),
  launch_browser=false,
  require_secret_for_open_links=false,
  require_secret_for_access=false,
  auto_reload_from_file=true,
  run_notebook_on_load=true,
  notebook=notebook)
