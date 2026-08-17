# PSCAD harvester

The maintained harvester is Julia code. Install its project-local Julia app
once from the repository root and ensure `~/.julia/bin` is on `PATH`:

```powershell
julia -e 'using Pkg; Pkg.Apps.develop(path="gauntlet")'
```

The installed command is `linecablebenchmark`. Gauntlet owns planning,
completeness checks, hashing, capture, normalization, and reporting.

Static commands need neither PSCAD nor Python:

```powershell
linecablebenchmark pending dataset_manifest.json
linecablebenchmark verify dataset_manifest.json runs/pendingcases
linecablebenchmark ingest . staging --amendments runs/pendingcases
```

Live PSCAD automation is the one external runtime boundary. Point PythonCall at
the Python installation containing the external `mhi.pscad` package, then
instantiate its isolated environment once. The CLI starts that environment and
loads the Gauntlet extension only for a live command:

```powershell
$env:JULIA_CONDAPKG_BACKEND = "Null"
$env:JULIA_PYTHONCALL_EXE = "C:\path\to\python.exe"
julia --project=gauntlet/harvest/pscad/live -e 'using Pkg; Pkg.instantiate()'
linecablebenchmark inspect SOURCE.pscx
linecablebenchmark case SOURCE.pscx --output runs
```

`benchmark.toml` defines the frequency sweep, soil resistivities, reduction
states, and formulation axes. `sources.toml` pins external donors and hashes.
No project-owned Python source is retained.

Each successful run records the normalized `.cli` or `.tli` input, ordinary
and detailed outputs, a `.clo` or `.tlo` fit when emitted, hashes, PSCAD
version, source definition, mutations, and elapsed time. Terminal PSCAD solver
rejections remain evidence. Gauntlet's Julia ingester deduplicates these runs,
normalizes units and ordering, and stages the immutable dataset archives.

The source campaign used for `pscad-raw-v1` and `pscad-normalized-v1` is
inventoried in the repository's `GAUNTLET.md`. Raw artifact publication still
requires a human licensing review.
