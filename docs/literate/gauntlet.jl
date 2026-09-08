# # Gauntlet
#
# Gauntlet runs cases manually and retains numerical results. This page is its
# only report: a compact comparison of the **applicable baseline formulations**
# for each case, read from explicitly configured, persisted benchmarks.
# Documentation generation never starts PSCAD, FEM, analytical calculations or
# uncertainty propagation. Full-catalogue and UQ results remain stored for
# analysis through the ordinary result, observation and plotting APIs.
#
# ## Recorded baseline comparisons
#
# Both errors are calculated element-wise by
# [`compare`](@ref LineCableModels.Engine.compare), with **A = the benchmark's
# reference** and **B = its candidate**. Neither role is inferred from a backend.
#
# ```math
# \mathrm{NRMSE}=100\sqrt{\frac{\sum_k|B_k-A_k|^2}{\sum_k|A_k|^2}},
# \qquad
# \mathrm{RMS}_{\mathrm{pointwise}}=100\sqrt{\frac1N\sum_k\left|\frac{B_k-A_k}{A_k}\right|^2}.
# ```
#
# **One row is one benchmark.** Reference and candidate columns identify the
# backend and requested formulation selection. Z and Y stay side by side;
# each cell contains **maximum error in percent (response, excitation)**.
# The maximum is across matrix entries, not a whole-matrix error. The two
# normalizations can attain their maxima at different entries. Expand the
# benchmark identities below each full-band table to see the terminal order.
#
# **The full-band summary comes first.** Frequency slices follow in separate
# sections, with the same column layout and their actual stored ranges.
# No interpolation or additional simulations are performed.
#
# <!-- GAUNTLET_REPORT -->
#
# ## Comparison conventions
#
# No denominator floor or solver-data modification is applied. When **both**
# traces are below the recorded quantity-specific numerical-zero tolerance,
# their error is zero. Defaults are 1e-10 Ω/m
# for R, 1e-15 H/m for L, 1e-12 S/m for G and 1e-16 F/m for C. The Z and Y
# thresholds follow R + 2πfL and G + 2πfC, respectively. Callers may override
# them in `compare(...; atol=(G=..., C=...))`. Unsupported comparisons and
# empty bands return `missing` with a reason. An infinite error otherwise
# means a nonzero difference against an exact zero normalization. In the
# pointwise metric, exact 0/0 contributes zero and nonzero/0 contributes infinity;
# no samples are omitted. This differs from normalizing by the complete reference
# trace's RMS. Absolute RMS arrays and all individual relative errors remain in
# the stored comparison objects.
#
# The deterministic summary displays **Z and Y only**. Shunt conductance
# `G = real(Y)` remains available for an explicitly requested loss study through
# `compare(reference, candidate, G; ...)`; it is not an additional default KPI.
# Saved G comparisons are retained, but do not expand the Z/Y summary.
# UQ moment benchmarks, when selected, have separate mean and standard-deviation
# sections; their R/L/C/G observables are not mixed into deterministic Z/Y tables.
#
# Two profiles are included: all-`:default`, and `:Ametani2004` for both
# insulation and semicon with all other slots at `:default`. Defaults
# are contextual and may implement different approximations in each backend.
# PSCAD's scalar export matches the selected radial dielectric admittance at
# its base frequency (50 Hz in this campaign), before the native loss-tangent
# cap of 10. That fit is not an arbitrary broadband constitutive law. The owned
# engine and FEM evaluate the retained dielectric constituents at each frequency;
# homogeneous radial equivalence does not assert full-geometry field equivalence.
# In this campaign, each benchmark explicitly selects PSCAD or FEM as reference
# and coaxial as candidate. This is a choice in its definition, not a Gauntlet
# rule: same-backend and other cross-backend pairs are equally valid.
# No backend is ground truth or
# an approved CI reference. Inputs, terminal order, basis and frequencies must
# match; there is no implicit interpolation or conversion. No detailed HTML
# pages, plots, input-object dumps or numerical-file copies are published here.
#
# ## Run and select stored results
#
# ```bash
# lcm gauntlet run --directory /path/to/campaign \
#   --backends coaxial,fem,pscad --formulas default --frequency-range 0.1,1e6
# lcm gauntlet status --directory /path/to/campaign
# lcm gauntlet resume --directory /path/to/campaign
# lcm gauntlet run --directory /path/to/lossy-campaign \
#   --backends coaxial,fem,pscad --formulas default --dielectric Ametani2004 \
#   --frequency-range 0.1,1e6
#
# lcm gauntlet compare --definition /path/to/benchmarks.toml \
#   --output /path/to/comparisons
# LINECABLEMODELS_GAUNTLET_RESULTS=/path/to/comparisons \
#   julia --project=docs docs/make.jl
# ```
#
# Omit `--cases` to run the indexed catalogue, or pass comma-separated case IDs.
# `--frequency-range` sets the same 101 logarithmic samples for every case;
# without it, each case retains its own extent, with a 0.1 Hz minimum.
# `--formulas catalogue` retains the broader formula sweep without expanding
# this summary. Formula grids, explicit lossy selections and uncertainty runs
# are documented in the
# [Gauntlet CLI guide](https://github.com/Electa-Git/LineCableModels.jl/blob/main/test/gauntlet/README.md).
# FEM Monte Carlo remains disabled.
#
# A benchmark definition names exactly one reference and one candidate artifact
# (paths and SHA-256 checksums), plus quantities, bands and normalizations. See the
# CLI guide for the TOML format. Comparing saved files never reruns their solvers.
# Missing operands are errors; they never trigger a replacement reference.
# Multiple comparison directories can be selected with the platform path-list
# separator (`:` on Unix, `;` on Windows); identical benchmark records appear once.
# This completed-only campaign retains 86 calculations and 46 benchmark pairs:
# 38 PSCAD-reference and 8 FEM-reference. The 32 unfinished FEM selections remain
# excluded. No simulations are resumed by this page. Batch elapsed-at-completion
# values are not advertised as per-selection cold or warmed execution timings.
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
