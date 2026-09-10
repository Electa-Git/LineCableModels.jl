# Separate warmed construction from image reuse on the frozen earth-only case.
# julia --project=. test/gauntlet/spectral_backend_cost.jl [output-directory]
using LineCableModels, LinearAlgebra, TOML, JSON3, SHA
const E=LineCableModels.Engine
const ROOT=normpath(joinpath(@__DIR__, "../.."))
const FILES=("src/engine/integration.jl", "src/engine/spectralsampling.jl",
    "src/engine/compleximages.jl", "src/engine/earthreturn.jl", "src/engine/earthkernels.jl",
    "src/engine/input.jl", "src/engine/Engine.jl")
hashes()=Dict(file=>bytes2hex(sha256(read(joinpath(ROOT, file)))) for file in FILES)
const INITIAL_HASHES=hashes()

function measure_backends(directory)
    fixture=TOML.parsefile(joinpath(ROOT, "test/fixtures/reference/unified_earth_return.toml"))
    rows=NamedTuple[]
    for n in (2, 3), method in (:quad, :trapz, :cim)
        expected=only(filter(row->row["model"]=="proposed_full_current"&&
            length(row["positions_m"])==n&&row["frequency"]==1e6, fixture["cases"]))
        geometry=E.EarthReturnGeometry(first.(expected["positions_m"]), last.(expected["positions_m"]),
            fill(expected["radius_m"], n))
        workspace=E.EarthReturnWorkspace(geometry)
        s=2pi*1e6im
        sigma=[0.0, inv(expected["earth_resistivity_ohm_m"])]
        epsilon=fill(fixture["analytic_epsilon"], 2); mu=fill(fixture["analytic_mu"], 2)
        state=(jω = s, Γ = zero(s), sigma, epsilon, mu,
            gamma_medium_squared = s .* mu .* (sigma .+ s .* epsilon))
        controls=E.computation_options(E.SpectralIntegral, (method,
            options = (rtol = method===:quad ? 1e-9 : 1e-6,)))
        E.unified_earth!(workspace, state, controls)
        for mode in (:construction, :reuse)
            records=map(1:3) do _
                mode===:construction && empty!(workspace.scratch.numerical.cim.fits)
                timed=@timed E.unified_earth!(workspace, state, controls)
                cache=workspace.scratch.numerical.cim
                (seconds = timed.time, bytes = timed.bytes,
                    statistics = map(x->x[], cache.statistics),
                    images = [length(fit.images) for fit in cache.fits],
                    evaluations = workspace.scratch.report.evaluations[])
            end
            measured=records[argmin(map(row->row.seconds, records))]
            for (kind, floor) in ((:Ze, 1e-14), (:Pe, 1e-14), (:Ye, 1e-10))
                actual=getproperty(workspace, kind)
                wanted=reshape(complex.(expected[string(kind)*"_real"], expected[string(kind)*"_imag"]), n, n)
                for component in (real, imag), i in eachindex(actual)
                    @assert isapprox(component(actual[i]), component(wanted[i]); rtol = 1e-5, atol = floor) (n, method, mode, kind, i)
                end
            end
            row=(conductors = n, frequency_Hz = 1e6, method, mode, measured...)
            push!(rows, row)
            println(row); flush(stdout)
        end
    end
    @assert hashes()==INITIAL_HASHES "Source changed during measurement"
    mkpath(directory)
    open(joinpath(directory, "spectral-backend-cost.json"), "w") do io
        JSON3.pretty(io, (julia_version = string(VERSION), blas_threads = BLAS.get_num_threads(),
            sources = INITIAL_HASHES, rows))
    end
end
measure_backends(isempty(ARGS) ? joinpath(ROOT, ".linecablemodels/qa/spectral-backend-cost") : only(ARGS))
