# FEM release gate

This directory is intentionally outside `test/`. Routine package tests and
coverage do not discover it.

Run the gate with Julia 1.12 after installing Gmsh in the active environment and
setting `GETDP_EXECUTABLE` to a GetDP 3.5.0 executable:

```sh
GETDP_EXECUTABLE=/absolute/path/to/getdp \
  julia --project=. -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path=pwd()); Pkg.add("Gmsh"); include("integration/fem/runtests.jl")'
```

The manual GitHub workflow downloads the official Linux archive and checks
SHA-256 `d3c28fa18f20d6147b4c7367d4dd802e9f7ddb58c608688bbb71919dbca8041d`
before running the 19 assertions with coverage disabled.
