# Reproduce accepted earth-only matrices with the production closure and record
# warmed backend cost. No internal/insulation contribution or FEM solve is used.
# julia --project=. test/gauntlet/unified_earth_validation.jl [output-directory]
using LineCableModels, LinearAlgebra, TOML, JSON3, SHA
const E=LineCableModels.Engine
const ROOT=normpath(joinpath(@__DIR__, "../.."))
const FIXTURE=joinpath(ROOT, "test/fixtures/reference/unified_earth_return.toml")
const SOURCE_FILES=("src/engine/earthkernels.jl", "src/engine/earthreturn.jl",
    "src/engine/integration.jl", "src/engine/spectralsampling.jl")
function source_hashes()
    Dict(file=>bytes2hex(sha256(read(joinpath(ROOT, file)))) for file in SOURCE_FILES)
end
const SOURCE_HASHES=source_hashes()

function validate_earth(directory)
    fixture=TOML.parsefile(FIXTURE)
    rows=NamedTuple[]
    for n in (2, 3), method in (:quad, :trapz, :cim)

        expected=only(filter(fixture["cases"]) do row
            row["model"]=="proposed_full_current"&&length(row["positions_m"])==n&&row["frequency"]==1e6
        end)
        positions=expected["positions_m"]
        geometry=E.EarthReturnGeometry(first.(positions), last.(positions), fill(expected["radius_m"], n))
        workspace=E.EarthReturnWorkspace(geometry)
        s=2pi*1e6im
        sigma=[0.0, inv(expected["earth_resistivity_ohm_m"])]
        epsilon=fill(fixture["analytic_epsilon"], 2)
        mu=fill(fixture["analytic_mu"], 2)
        state=(jω = s, Γ = zero(s), sigma, epsilon, mu,
            gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
        integration=E.computation_options(
            E.SpectralIntegral, (method, options = (rtol = method===:quad ? 1e-9 : 1e-6,)))
        E.unified_earth!(workspace, state, integration)
        times=[@elapsed(E.unified_earth!(workspace, state, integration)) for _ in 1:3]
        bytes=@allocated E.unified_earth!(workspace, state, integration)
        matrices=map(((:Ze, 1e-14), (:Pe, 1e-14), (:Ye, 1e-10))) do (kind, floor)
            wanted=reshape(complex.(expected[string(kind) * "_real"], expected[string(kind) * "_imag"]), n, n)
            actual=getproperty(workspace, kind)
            tolerance=method===:quad ? 1e-7 : 1e-5
            ratios=[abs(component(actual[i]-wanted[i]))/(floor+tolerance*abs(component(wanted[i])))
                    for component in (real, imag) for i in eachindex(actual)]
            @assert maximum(ratios)<=1 (n, method, kind, maximum(ratios))
            (quantity = kind, real = vec(real.(actual)),
                imag = vec(imag.(actual)), maximum_tolerance_ratio = maximum(ratios))
        end
        report=workspace.scratch.report
        row=(conductors = n, frequency_Hz = 1e6, method,
            minimum_warmed_seconds = minimum(times),
            allocated_bytes = bytes, integrals = report.integrals[], kernel_evaluations = report.evaluations[],
            matrix_refinements = report.refinements[], cutoff = report.cutoff[], matrices)
        push!(rows, row)
        println((n, method, seconds = row.minimum_warmed_seconds,
            bytes, evaluations = row.kernel_evaluations))
        flush(stdout)
    end
    mkpath(directory)
    path=joinpath(directory, "production-validation.json")
    sources=source_hashes()
    @assert sources==SOURCE_HASHES "Sources changed during the benchmark; rerun with a stable working tree"
    open(path, "w") do io
        JSON3.pretty(io,
            (julia_version = string(VERSION), blas_threads = BLAS.get_num_threads(),
                fixture_sha256 = bytes2hex(sha256(read(FIXTURE))), sources, rows))
    end
    println(path)
end
validate_earth(isempty(ARGS) ?
               joinpath(ROOT, ".linecablemodels/qa/unified-earth-production") : only(ARGS))
