# # Gauntlet
#
# Gauntlet runs cases manually and retains numerical results. This page is its
# only report: a compact comparison of the **applicable default formulations**
# for each case and backend, read from completed calculation files.
# Documentation generation never starts PSCAD, FEM, analytical calculations or
# uncertainty propagation. Full-catalogue and UQ results remain stored for
# analysis through the ordinary result, observation and plotting APIs.
#
# ## Recorded default comparisons
#
# `ε` is the element-wise RMS relative difference over the complete stored
# frequency vector, calculated by [`compare`](@ref LineCableModels.Engine.compare).
# Each table entry is the largest matrix-entry difference, expressed as a
# percentage. An infinite value means a nonzero difference against a reference
# entry whose RMS is zero. No numerical floor hides small entries.
#
# Only selections explicitly recorded as all-`:default` are included. Defaults
# are contextual and may implement different approximations in each backend.
# FEM and PSCAD are each compared with the analytical (`coaxial`) calculation;
# they are not compared against one another here. The analytical baseline only
# supplies the relative-difference denominator: no backend is ground truth or
# an approved CI reference. Inputs, terminal order, basis and frequencies must
# match; there is no implicit interpolation or conversion. No detailed HTML
# pages, plots, input-object dumps or numerical-file copies are published here.
#
# <!-- GAUNTLET_REPORT -->
#
# ## Run and select stored results
#
# ```bash
# lcm gauntlet run --directory /path/to/campaign \
#   --backends coaxial,fem,pscad --formulas default
# lcm gauntlet status --directory /path/to/campaign
# lcm gauntlet resume --directory /path/to/campaign
#
# LINECABLEMODELS_GAUNTLET_RESULTS=/path/to/campaign \
#   julia --project=docs docs/make.jl
# ```
#
# Omit `--cases` to run the indexed catalogue, or pass comma-separated case IDs.
# `--formulas catalogue` retains the broader formula sweep without expanding
# this summary. Formula grids, explicit lossy selections and uncertainty runs
# are documented in the
# [Gauntlet CLI guide](https://github.com/Electa-Git/LineCableModels.jl/blob/main/test/gauntlet/README.md).
# FEM Monte Carlo remains disabled.
#
# Multiple campaign directories can be selected with the platform path-list
# separator (`:` on Unix, `;` on Windows). Select at most one completed default
# per case/backend; duplicate identical files are included once. Interrupted
# campaigns contribute their completed defaults and an unfinished count.
# No latest-run selection, case reconstruction or report bundle is required.
#
# ## Numerical references for CI
#
# Stored Gauntlet results are not automatically approved CI references. The
# separate [numerical-reference gate](https://github.com/Electa-Git/LineCableModels.jl/tree/main/test/numerical)
# requires reviewed results, explicit tolerances and artifact bindings. It never
# starts a Gauntlet campaign or refreshes a reference. Approval remains separate
# from collecting results or displaying this summary.
#
# ## Benchmark data API
#
# ```@docs
# LineCableModels.Engine.RMSError
# LineCableModels.Engine.LineParametersBenchmark
# LineCableModels.Engine.compare
# LineCableModels.Engine.absolute_error
# LineCableModels.Engine.relative_error
# ```
