"""
    Gauntlet

Case ingestion, declared benchmark execution, comparisons and immutable artifacts.
Computations use LineCableModels.compute and retain their declared inputs.
"""
module Gauntlet
using LineCableModels
using LineCableModels.Engine
using LineCableModels: PSCAD, AbstractCoreResult, AbstractFormulation, AbstractGrid,
    AbstractParametricResult, AbstractUncertaintyResult, Formulation, Grid, Gridspace,
    LineParametersProblem, ParametricResult, PhaseDomain, build, description, details,
    formula, nominal, quantity
using LineCableModels.Engine: Engine
import LineCableModels.Grammar
import LineCableModels.ImportExport
import LineCableModels.Units
import Pkg
using BenchmarkTools: BenchmarkTools
using LinearAlgebra: BLAS
using Statistics: Statistics, median
import LineCableModels: compute, validate
export ARTIFACT_ROOT, ARTIFACTS_TOML, SNAPSHOT_SCHEMA_VERSION,
       WORK_ROOT,
       AbstractCaseVariation, CaseDefinition, CaseParameter, CompositeVariation,
       ExactOverrides, LoadedCase, NoVariation, ParameterGrids,
       RelativeStandardUncertainty,
       BenchmarkCalculation, MomentBenchmark, MomentResult,
       BenchmarkDefinition,
       UQ_MONTE_CARLO_TRIALS,
       artifact_name, bind_published_artifact,
       collection_archive_name, collection_release, collection_stage,
       benchmark_local, benchmark_stage,
       benchmark_definition, calculation_record,
       cleanup_work,
       case_definition, case_index, case_parameter, compose_variations,
       correlation_record,
       gauntlet_instrumented, finalize_staging,
       formulation_record,
       extract_moments, moment_comparison_passes, moment_error_summary,
       parameter_manifest,
       numerical_input_sha256, implementation_record, repository_revision,
       load_case, performance_comparison, package_collection, prepare_staging, release_tag,
       read_collection,
       run_benchmark, uq_inner_formulation, uq_moment_tolerances,
       variation_record

include("artifacts.jl")
include("cases.jl")
include("records.jl")
include("uq_benchmarks.jl")
include("performance.jl")
include("benchmarks.jl")
include("campaigns.jl")
include("reporting.jl")
end
