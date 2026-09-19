# Current construction recipes. These are inputs and routing controls, never
# numerical snapshots. Scientific expectations belong beside the relevant test.
module CurrentScenarios
    export conductor_material, coaxial_design, three_phase_system, two_wire_system,
        line_parameters_problem, channel_value, two_conductor_results, cable_monte_carlo_result,
        three_bare_wires_problem, three_bare_wires_layouts
    using LineCableModels
    using LineCableModels.DocStringExtensions: TYPEDSIGNATURES
    using Statistics: mean

    "Signed wire heights \\[m\\] for the five three-wire interface placements."
    const three_bare_wires_layouts = (
        all_air=(1.0, 1.0, 1.0), all_earth=(-1.0, -1.0, -1.0),
        air_1=(1.0, -1.0, -1.0), air_2=(-1.0, 1.0, -1.0),
        air_3=(-1.0, -1.0, 1.0))

    """
    $(TYPEDSIGNATURES)

    Construct three identical bare copper wires with separate left-to-right
    terminals. This physical fixture selects no analytical or FEM formulation.

    # Keywords

    - `heights`, `horizontal`: Three signed vertical and horizontal coordinates \\[m\\].
    - `radius`: Wire radius \\[m\\].
    - `copper`: Conductor material; defaults to the materials library's copper.
    - `rho`: Soil resistivity \\[Ω·m\\].
    - `eps_r`, `mu_r`: Soil relative permittivity and permeability \\[dimensionless\\].
    - `temperature`: Operating temperature \\[°C\\].
    - `line_length`: Line length \\[m\\].
    - `frequencies`: Analysis frequencies \\[Hz\\].
    - `name`: Physical system identifier.

    # Returns

    - A `LineParametersProblem` with no prescribed longitudinal propagation constant.
    """
    function three_bare_wires_problem(; heights=(1.0, 1.0, 1.0),
            horizontal=(0.0, 1.0, 2.0), radius=0.0425,
            copper=Material(MaterialsLibrary(add_defaults=true), :copper),
            rho=0.1, eps_r=1.0, mu_r=1.0, temperature=20.0,
            line_length=1.0, frequencies=10.0 .^ (-1:7), name="three_bare_wires")
        length(heights) == length(horizontal) == 3 ||
            throw(ArgumentError("three wire coordinates are required"))
        design = build(CableDesign, name, Stack(
            Group(:core, Region(:core_metal, Disk(radius), copper))))
        system = build(LineCableSystem, fill(design, 3),
            [Pose2(horizontal[i], heights[i]) for i in 1:3];
            connections=[Dict(:core=>i) for i in 1:3], system_id=name, line_length)
        return LineParametersProblem(system; temperature, frequencies,
            earth_props=homogeneous(; rho, eps_r, mu_r))
    end

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
