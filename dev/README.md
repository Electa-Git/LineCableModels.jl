# Manual development scripts

Use the active Julia environment when including these scripts from an IDE or
REPL. The scripts do not select or install a project on the caller's behalf.

| Script | Purpose |
| --- | --- |
| `run_pscad.jl` | Runs a two-wire PSCAD toy case using the private `local-pscad.toml` beside the script. |
| `run_quasi_full.jl` | Runs the disposable quasi-full Gmsh/GetDP investigation in a fresh local directory. |
| `run_three_bare_wires_baseline.jl` | Captures the three-bare-wire analytical baseline. |
| `verify_getdp_artifact.jl` | Verifies the packaged GetDP executable and source assets. |

For PSCAD, edit the ignored `dev/local-pscad.toml` with your station settings,
then run `include("dev/run_pscad.jl")` from the repository root. A sanitized
configuration template is available at `examples/pscad/remote.example.toml`.
Each execution creates a fresh run in the configured exchange directory and
leaves `result`, `Zresult`, and `Yresult` available for inspection.

The `plotting/` directory has its own project and manual GL/WGL gallery scripts.
Activate `dev/plotting` explicitly when running those examples.
