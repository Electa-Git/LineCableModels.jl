# Current construction recipes. These are inputs and routing controls, never
# numerical snapshots. Scientific expectations belong beside the relevant test.
module CurrentScenarios
    export conductor_material, coaxial_design, three_phase_system, two_wire_system,
        line_parameters_problem, channel_value, two_conductor_results, cable_monte_carlo_result
    using LineCableModels
    using Statistics: mean

    conductor_material() = Material(kind=:conductor, rho=2e-8, eps_r=1.0,
        mu_r=1.0, T0=20.0, alpha=0.004)

    function coaxial_design(; scale=1.0, name="concentric-control")
        dielectric = Material(kind=:insulator, rho=1e8, eps_r=3.0, mu_r=1.0)
        return build(CableDesign, name, Stack(
            terminal(:core, Region(:metal, Disk(0.005scale), conductor_material())),
            Region(:dielectric, Shell(0.005scale), dielectric),
            terminal(:sheath, Region(:return, Shell(0.001scale), conductor_material())),
            Region(:cover, Shell(0.001scale), dielectric)))
    end

    function three_phase_system(; line_length=600.0, spacing=0.08)
        designs = [coaxial_design(; name="phase-$phase") for phase in 1:3]
        return build(LineCableSystem, designs,
            [(-spacing, -1.0), (0.0, -1.0-spacing), (spacing, -1.0)];
            system_id="current-three-phase", line_length,
            connections=[Dict(:core=>phase, :sheath=>0) for phase in 1:3])
    end

    function two_wire_system()
        designs = [build(CableDesign, "wire-$i", Stack(
            terminal(:core, Region(:metal, Disk(0.005), conductor_material())),
            Region(:cover, Shell(0.001), Material(kind=:insulator, rho=1e8, eps_r=3.0))))
            for i in 1:2]
        return build(LineCableSystem, designs, [(0.0,-1.0), (0.2,-1.3)];
            system_id="current-two-wire", line_length=600.0,
            connections=[Dict(:core=>i) for i in 1:2])
    end

    # Passing the system is essential when a consumer checks its terminal identity.
    function line_parameters_problem(system=three_phase_system(); frequencies=[50.0])
        return LineParametersProblem(system; temperature=20.0,
            earth_props=homogeneous(rho=100.0, eps_r=10.0, mu_r=1.0), frequencies)
    end

    # Deliberately asymmetric metadata: this tests routing, not passive physics.
    # Every quantity, coordinate and frequency has a distinct declared value.
    channel_value(::Val{:R}, i,j,k) = 1e-4*(1+2i+3j+k)
    channel_value(::Val{:L}, i,j,k) = 1e-7*(2+i+4j+2k)
    channel_value(::Val{:G}, i,j,k) = 1e-9*(3+3i+j+3k)
    channel_value(::Val{:C}, i,j,k) = 1e-10*(4+4i+2j+k)
    function two_conductor_results(; frequencies=[10.0,100.0,1000.0])
        f=collect(frequencies)
        channels=NamedTuple{(:R,:L,:G,:C)}(Tuple(
            [channel_value(Val(q),i,j,k) for i in 1:2,j in 1:2,k in eachindex(f)]
            for q in (:R,:L,:G,:C)))
        omega=reshape(2pi.*f,1,1,:)
        return LineParameters(channels.R .+ im.*omega.*channels.L,
            channels.G .+ im.*omega.*channels.C, f; basis=:pul,
            details=ComputationDetails(;coordinates=["west","east"],))
    end

    # Independent arrays and histogram storage for each observable. This object
    # tests publication/retention only; actual Monte Carlo validation runs compute.
    function cable_monte_carlo_result()
        samples=(R=reshape([2.0,3.0,5.0,8.0],1,:),
            L=reshape([11.0,13.0,17.0,19.0].*1e-6,1,:),
            C=reshape([23.0,29.0,31.0,37.0].*1e-10,1,:),
            G=reshape([41.0,43.0,47.0,53.0].*1e-9,1,:))
        statistics=map(x->[SampleSummary(vec(x))],samples)
        histograms=map(x->[HistogramDensity(vec(x);bins=2)],samples)
        representation=CableConstants(mean(samples.R),mean(samples.L),
            mean(samples.C),mean(samples.G))
        representation=LineCableModels.materialize(representation,statistics)
        formulation=MonteCarlo(Formulation();trials=4,seed=2027,
            return_samples=true,return_histograms=true)
        return MonteCarloResult(formulation,[representation],[statistics],
            [samples],[histograms],UInt64(2027),UInt64[2039],[4])
    end
end
