# PSCAD harvester boundary

This directory preserves the narrow Python boundary that launches PSCAD's
Line Constants Program and captures content-addressed evidence. It is not used
by Gauntlet's parser, smoke suite, documentation, or full-corpus validation.
Those paths are implemented in Julia and operate without Python or PSCAD.

The harvester requires Windows, a licensed PSCAD installation, `mhi.pscad`,
`mhi.xml`, and PyYAML. It compiles selected frequency-dependent line or cable
definitions; it does not run an EMTDC time-domain simulation.

From this directory in PowerShell:

```powershell
.\linecablebenchmark.cmd --help
.\linecablebenchmark.cmd harvest --verify-only --no-extract
.\linecablebenchmark.cmd inspect --file SOURCE.pscx
.\linecablebenchmark.cmd case --file SOURCE.pscx --dry-run
.\linecablebenchmark.cmd batch
```

`benchmark.yml` defines the frequency sweep, soil resistivities, reduction
states, and formulation axes. `official_sources.yml` pins the external donor
sources and hashes. Generated donor trees, copied projects, runs, caches,
diagnostics, and dataset manifests deliberately do not live here.

Each successful run records the normalized `.cli` or `.tli` input, ordinary
and detailed outputs, a `.clo` or `.tlo` fit when emitted, hashes, PSCAD
version, source definition, mutations, and elapsed time. Terminal PSCAD solver
rejections remain evidence. Gauntlet's Julia ingester deduplicates these runs,
normalizes units and ordering, and stages the immutable corpus archives.

The source campaign used for `pscad-raw-v1` and `pscad-normalized-v1` is
inventoried in the repository's `GAUNTLET.md`. Raw artifact publication still
requires a human licensing review.
