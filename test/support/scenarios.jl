# Current construction recipes. These are inputs and routing controls, never
# numerical snapshots. Scientific expectations belong beside the relevant test.
module CurrentScenarios
    export conductor_material, coaxial_design, three_phase_system, two_wire_system,
        line_parameters_problem, channel_value, two_conductor_results, cable_monte_carlo_result,
        three_bare_wires_problem, three_bare_wires_layouts, known_case_problem,
        preservation_corpus,
        internal_shunt_test_design, screened_cable_design, armored_cable_design,
        three_core_design, screened_three_core_design, offset_conduit_design,
        single_cable_system, trefoil_system, bipole_system,
        exported_text, reduction_corpus
    using LineCableModels
    using LineCableModels.DocStringExtensions: TYPEDSIGNATURES
    using Logging: NullLogger, with_logger
    using Statistics: mean

    include("internal_shunt.jl")

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

    - `heights`, `horizontal`: three signed vertical and horizontal coordinates \\[m\\].
    - `radius`: wire radius \\[m\\].
    - `copper`: conductor material. Defaults to the materials library's copper.
    - `rho`: soil resistivity \\[Ω·m\\].
    - `eps_r`, `mu_r`: soil relative permittivity and permeability \\[dimensionless\\].
    - `temperature`: operating temperature \\[°C\\].
    - `line_length`: line length \\[m\\].
    - `frequencies`: analysis frequencies \\[Hz\\].
    - `name`: physical system identifier.

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

    """
    $(TYPEDSIGNATURES)

    Construct a case of the FEM investigation, as its script in
    `LinecableModels-fem-investigation/test/manual/cases/<case>.jl` builds it: the same
    geometry, materials, placements, connections, earth and temperature.
    `:two_bare_wires` and `:two_insulated_wires` are copper wires 1 m apart and 1 m deep,
    bare or with 1 mm of insulation. `:solid_1000mm2_single` is an aluminum core with XLPE.
    `:cable_525kv_land_no_armour_dc_bipole` is a solid-core land cable bipole.
    `:cable_320kv_no_armour_dc_bipole` is a stranded-core cable bipole.
    `:cable_380kv_no_armour_ac_flat` is three stranded-core cables in flat formation.
    `:cable_132kv_cigre_tb880_case0_630cu_trefoil` is three solid-core cables in trefoil.
    `:cable_132kv_630mm2_flathor` is three stranded-core cables in horizontal flat formation.
    `:cable_18kv_1000mm2_trefoil` is three stranded-core 18 kV cables in trefoil.
    `:cable_18kv_1000mm2_trefoil_homogenized` is three homogenized 18 kV cables in trefoil.
    `:cable_220kv_milliken_1x2500_252_trefoil` is three Milliken-core cables in trefoil.
    `:cable_30kv_na2xs2y_630mm2_trefoil` is three stranded-core 30 kV cables in trefoil.
    `:cable_320kv_armoured_dc_bipole` is an armored stranded-core cable bipole.
    `:cable_380kv_2000mm2_flatver` is three stranded-core cables in vertical flat formation.
    `:cable_380kv_armoured_ac_flat` is three armored stranded-core cables in flat formation.
    `:cable_525kv_1600mm2_bipole` is two stranded-core 525 kV cables in a bipole.
    `:cable_525kv_land_no_armour_ac_flat` is three solid-core land cables in flat formation.
    `:cable_525kv_subsea_armoured_ac_flat` is three armored subsea cables in flat formation.
    `:cable_525kv_subsea_armoured_dc_bipole` is an armored subsea cable bipole.
    `:cable_640kv_2000mm2_bipole` is two stranded-core 640 kV cables in a bipole.
    This physical fixture selects no formula.

    # Keywords

    - `frequencies`: analysis frequencies \\[Hz\\]. The case scripts use 101 points from
      0.1 Hz to 10 MHz.
    """
    known_case_problem(case::Symbol; kwargs...) = known_case_problem(Val(case); kwargs...)

    include("known_cases/cable_132kv_630mm2_flathor.jl")
    include("known_cases/cable_132kv_cigre_tb880_case0_630cu_trefoil.jl")
    include("known_cases/cable_18kv_1000mm2_trefoil.jl")
    include("known_cases/cable_18kv_1000mm2_trefoil_homogenized.jl")
    include("known_cases/cable_220kv_milliken_1x2500_252_trefoil.jl")
    include("known_cases/cable_30kv_na2xs2y_630mm2_trefoil.jl")
    include("known_cases/cable_320kv_armoured_dc_bipole.jl")
    include("known_cases/cable_320kv_no_armour_dc_bipole.jl")
    include("known_cases/cable_380kv_2000mm2_flatver.jl")
    include("known_cases/cable_380kv_armoured_ac_flat.jl")
    include("known_cases/cable_380kv_no_armour_ac_flat.jl")
    include("known_cases/cable_525kv_1600mm2_bipole.jl")
    include("known_cases/cable_525kv_land_no_armour_ac_flat.jl")
    include("known_cases/cable_525kv_land_no_armour_dc_bipole.jl")
    include("known_cases/cable_525kv_subsea_armoured_ac_flat.jl")
    include("known_cases/cable_525kv_subsea_armoured_dc_bipole.jl")
    include("known_cases/cable_640kv_2000mm2_bipole.jl")
    include("known_cases/solid_1000mm2_single.jl")
    include("known_cases/two_bare_wires.jl")
    include("known_cases/two_insulated_wires.jl")

    conductor_material() = Material(kind=:conductor, rho=2e-8, eps_r=1.0,
        mu_r=1.0, T0=20.0, alpha=0.004)

    function coaxial_design(; scale=1.0, name="concentric-control", nominal_data=nothing)
        dielectric = Material(kind=:insulator, rho=1e8, eps_r=3.0, mu_r=1.0)
        return build(CableDesign, name, Stack(
            terminal(:core, Region(:metal, Disk(0.005scale), conductor_material())),
            Region(:dielectric, Shell(0.005scale), dielectric),
            terminal(:sheath, Region(:return, Shell(0.001scale), conductor_material())),
            Region(:cover, Shell(0.001scale), dielectric)); nominal_data)
    end

    function three_phase_system(; line_length=600.0, spacing=0.08, nominal_data=nothing)
        designs = [coaxial_design(; name="phase-$phase", nominal_data) for phase in 1:3]
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
    # tests publication and retention only. Actual Monte Carlo validation runs compute.
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

    # Counters of the equivalence check: `materialized()` runs once per point a
    # parametric scenario realizes, and `nominal_data` reaches its cable designs.
    const NO_COUNTERS = (materialized=Returns(nothing), nominal_data=nothing)

    # Parametric study: earth resistivity and spacing, combined as `combine` asks, crossed
    # with two insulation admittance formulations.
    study_space(frequencies, combine, counters) = Gridspace{LineParametersProblem}(
        (rho, spacing) -> (counters.materialized();
            LineParametersProblem(three_phase_system(; spacing, counters.nominal_data);
                temperature=20.0, earth_props=homogeneous(; rho, eps_r=10.0, mu_r=1.0),
                frequencies)), (Grid((100.0, 1000.0)), Grid((0.08, 0.1))); combine)
    study_formulations() = Formulation(insulation_admittance=Grid((:default, :lossy)))

    # One uncertain cable spacing [m]. Its evaluation requires the Measurements extension.
    uncertain_space(frequencies, counters) = Gridspace{LineParametersProblem}(
        spacing -> (counters.materialized(); line_parameters_problem(
            three_phase_system(; spacing, counters.nominal_data); frequencies)),
        (Grid(0.08, AbsoluteError(0.002)),))

    """
    $(TYPEDSIGNATURES)

    The preservation corpus has one scenario per computation shape. The allocation
    ceilings, the timing comparison and the equivalence check share it.
    Each scenario is the tuple of positional arguments of `compute`, with `n` analysis
    frequencies between 10 Hz and 100 kHz. The CableConstants problem has its own
    single frequency. LinearError and MonteCarlo require the Measurements extension.
    `trefoil` is the only scenario with stranded cores, wire screens, a tape and
    semiconducting screens. New scenarios are appended, never inserted: the allocation counts depend on what a
    process computed before. `counters` instruments the parametric scenarios for the
    equivalence check and leaves their results unchanged.
    """
    function preservation_corpus(n::Integer; counters=NO_COUNTERS)
        frequencies = collect(10.0 .^ range(1, 5; length=n))
        coaxial = line_parameters_problem(three_phase_system(); frequencies)
        # On an earth of three soil layers, the default earth formulas consume an
        # equivalent earth.
        layered = LineParametersProblem(three_phase_system(); temperature=20.0,
            earth_props=build(EarthModel, (layer(rho=100.0, eps_r=10.0, thickness=0.5),
                layer(rho=500.0, eps_r=20.0, thickness=2.0), layer(rho=50.0, eps_r=15.0))),
            frequencies)
        before_fd = formula(:default; equivalent_earth=formula(:bottommost; order=:before))
        # The formulation of the known-cases test, `test/integration/known_cases.jl`.
        known_cases = Formulation(earth_impedance=:unified, earth_admittance=:unified,
            shunt_model=:equivalent,
            options=(reduce_bundle=false, kron_reduction=false, ideal_transposition=false))
        return (
            coaxial=(coaxial, Formulation()),
            overhead=(three_bare_wires_problem(;
                heights=three_bare_wires_layouts.all_air, frequencies), Formulation()),
            buried=(three_bare_wires_problem(;
                heights=three_bare_wires_layouts.all_earth, frequencies), Formulation()),
            modal=(coaxial, Formulation(), ModalAnalysisFormulation(:default)),
            cable_constants=(CableConstantsProblem(coaxial_design()), CableConstantsFormulation()),
            product=(ParametricProblem(study_space(frequencies, :product, counters)),
                Combinatorial(study_formulations())),
            zip=(ParametricProblem(study_space(frequencies, :zip, counters)),
                Combinatorial(study_formulations())),
            linear_error=(ParametricProblem(uncertain_space(frequencies, counters)),
                LinearError(Formulation())),
            monte_carlo=(ParametricProblem(uncertain_space(frequencies, counters)),
                MonteCarlo(Formulation(); trials=4, seed=2027)),
            three_layer=(layered, Formulation()),
            three_layer_before_fd=(layered, Formulation(earth_properties=formula(:alipio2014),
                earth_impedance=before_fd, earth_admittance=before_fd)),
            trefoil=(line_parameters_problem(trefoil_system(); frequencies), Formulation()),
            cable_132kv_630mm2_flathor=(
                known_case_problem(Val(:cable_132kv_630mm2_flathor); frequencies), known_cases),
            cable_132kv_cigre_tb880_case0_630cu_trefoil=(
                known_case_problem(Val(:cable_132kv_cigre_tb880_case0_630cu_trefoil); frequencies), known_cases),
            cable_18kv_1000mm2_trefoil=(
                known_case_problem(Val(:cable_18kv_1000mm2_trefoil); frequencies), known_cases),
            cable_18kv_1000mm2_trefoil_homogenized=(
                known_case_problem(Val(:cable_18kv_1000mm2_trefoil_homogenized); frequencies), known_cases),
            cable_220kv_milliken_1x2500_252_trefoil=(
                known_case_problem(Val(:cable_220kv_milliken_1x2500_252_trefoil); frequencies), known_cases),
            cable_30kv_na2xs2y_630mm2_trefoil=(
                known_case_problem(Val(:cable_30kv_na2xs2y_630mm2_trefoil); frequencies), known_cases),
            cable_320kv_armoured_dc_bipole=(
                known_case_problem(Val(:cable_320kv_armoured_dc_bipole); frequencies), known_cases),
            cable_320kv_no_armour_dc_bipole=(
                known_case_problem(Val(:cable_320kv_no_armour_dc_bipole); frequencies), known_cases),
            cable_380kv_2000mm2_flatver=(
                known_case_problem(Val(:cable_380kv_2000mm2_flatver); frequencies), known_cases),
            cable_380kv_armoured_ac_flat=(
                known_case_problem(Val(:cable_380kv_armoured_ac_flat); frequencies), known_cases),
            cable_380kv_no_armour_ac_flat=(
                known_case_problem(Val(:cable_380kv_no_armour_ac_flat); frequencies), known_cases),
            cable_525kv_1600mm2_bipole=(
                known_case_problem(Val(:cable_525kv_1600mm2_bipole); frequencies), known_cases),
            cable_525kv_land_no_armour_ac_flat=(
                known_case_problem(Val(:cable_525kv_land_no_armour_ac_flat); frequencies), known_cases),
            cable_525kv_land_no_armour_dc_bipole=(
                known_case_problem(Val(:cable_525kv_land_no_armour_dc_bipole); frequencies), known_cases),
            cable_525kv_subsea_armoured_ac_flat=(
                known_case_problem(Val(:cable_525kv_subsea_armoured_ac_flat); frequencies), known_cases),
            cable_525kv_subsea_armoured_dc_bipole=(
                known_case_problem(Val(:cable_525kv_subsea_armoured_dc_bipole); frequencies), known_cases),
            cable_640kv_2000mm2_bipole=(
                known_case_problem(Val(:cable_640kv_2000mm2_bipole); frequencies), known_cases),
            solid_1000mm2_single=(
                known_case_problem(Val(:solid_1000mm2_single); frequencies), known_cases),
            two_bare_wires=(
                known_case_problem(Val(:two_bare_wires); frequencies), known_cases),
            two_insulated_wires=(
                known_case_problem(Val(:two_insulated_wires); frequencies), known_cases))
    end

    """
    $(TYPEDSIGNATURES)

    Construct the 18/30 kV single-core cable of tutorial 2 (`examples/tutorial2.jl`). A
    compacted stranded aluminum core sits under semiconducting tapes and screens and PE
    insulation. One terminal holds 49 copper screen wires and a helical tape of the library
    material `tape_material`, copper in the tutorial. An aluminum laminate and two PE
    jackets follow.
    """
    function screened_cable_design(; tape_material::Symbol = :copper)
        materials = MaterialsLibrary(add_defaults = true)
        aluminum = Material(materials, :aluminum)
        copper = Material(materials, :copper)
        polyacrylate = Material(materials, :polyacrylate)
        pe = Material(materials, :pe)
        d_core = 38.1e-3
        t_sc_in = 0.6e-3
        t_ins = 8e-3
        t_sc_out = 0.3e-3
        d_ws = 0.95e-3
        t_sct = 0.3e-3
        stranded_core = stranded(aluminum; shape = Disk(4.7e-3 / 2),
            lay = LayRatio(15.0, 13.5, 12.5, 11.0), compact = true,
            boundary = Disk(d_core / 2))
        screen_wire_locus = d_core / 2 + t_sct + t_sc_in + t_ins + t_sc_out + t_sct +
                            d_ws / 2
        return @cable "18kV_1000mm2" begin
            @terminal :core begin
                stranded_core
                screen(polyacrylate; t = t_sct)
                screen(Material(materials, :semicon1); t = t_sc_in)
                insulation(pe; t = t_ins)
                screen(Material(materials, :semicon2); t = t_sc_out)
                screen(polyacrylate; t = t_sct)
            end
            @terminal :sheath begin
                wires(copper; shape = Disk(d_ws / 2), n = 49, r = screen_wire_locus,
                    lay = LayRatio(10), tag = :screen_wire)
                tape(Material(materials, tape_material); section = Rectangle(10e-3, 0.1e-3),
                    n = 1, lay = LayRatio(10), tag = Symbol(tape_material, "_", "tape"))
            end
            screen(polyacrylate; t = 0.3e-3)
            @terminal :jacket begin
                sheath(aluminum; t = 0.15e-3)
            end
            jacket(pe; t = 0.05e-3)
            jacket(pe; t = 2.4e-3)
        end
    end

    """
    $(TYPEDSIGNATURES)

    Construct the 525 kV single-core cable of tutorial 3 (`examples/tutorial3.jl`). A
    compacted stranded copper core sits under semiconducting screens, PE insulation and a
    tape. A lead sheath, a PE sheath and PP bedding follow, then 68 steel armor wires and a
    PP jacket.
    """
    function armored_cable_design()
        materials = MaterialsLibrary(add_defaults = true)
        pe = Material(materials, :pe)
        pp = Material(materials, :pp)
        core = stranded(Material(materials, :copper); shape = Disk(3.6649e-3 / 2),
            lay = LayRatio(11.0), compact = true, boundary = Disk(0.0463 / 2))
        return @cable "525kV_1600mm2" begin
            @terminal :core begin
                core
                screen(Material(materials, :semicon1); t = 2e-3)
                insulation(pe; t = 26e-3)
                screen(Material(materials, :semicon2); t = 1.8e-3)
                screen(Material(materials, :polyacrylate); t = 0.3e-3)
            end
            @terminal :sheath begin
                sheath(Material(materials, :lead); t = 3.3e-3)
            end
            jacket(pe; t = 3e-3)
            bedding(pp; t = 3e-3)
            @terminal :armor begin
                armor(Material(materials, :steel); shape = Disk(5.827e-3 / 2), n = 68,
                    lay = LayRatio(10), tag = :armor_wire)
            end
            jacket(pp; t = 10e-3)
        end
    end

    """
    $(TYPEDSIGNATURES)

    Construct three insulated cores of 19 stranded copper wires on one ring, without a
    common screen, as `test/unit/datamodel/assembly_lowering.jl` builds them.
    """
    function three_core_design()
        copper = Material(:conductor, 1.72e-8, 1.0, 1.0)
        wire = Disk(0.5e-3)
        member = terminal(:core,
            stranded(copper; shape = wire, boundary = Disk(5wire.r), lay = LayRatio(12)),
            insulation(Material(:insulator, 1e14, 2.3, 1.0); t = 0.5e-3))
        return build(CableDesign, "three-core", assembly(member;
            pattern = Ring(3; r = 0.01, φ0 = 0.27), names = (:phase_1, :phase_2, :phase_3)))
    end

    """
    $(TYPEDSIGNATURES)

    Construct a three-core cable whose cores each have their own copper screen. Each core
    is a copper conductor under XLPE insulation and a copper sheath. The cores lie in a
    trefoil, in a polypropylene filler under a common 2.5 mm PE jacket. `flatten` rejects
    the design, because it cannot reduce groups that share a container yet.
    """
    function screened_three_core_design()
        materials = MaterialsLibrary(add_defaults = true)
        copper = Material(materials, :copper)
        xlpe = Material(materials, :xlpe)
        function screened_core(prefix)
            return Stack(
                terminal(Symbol(prefix, :_core), core(copper; r = 4.0e-3),
                    insulation(xlpe; t = 3.4e-3)),
                terminal(Symbol(prefix, :_screen), sheath(copper; t = 0.1e-3)))
        end
        # The screened cores have a 7.5 mm outer radius. Their axes lie 8.7 mm from the
        # cable axis, so neighboring screens are 0.07 mm apart.
        spacing = 8.7e-3
        return @cable "screened-three-core" begin
            @pipe shape=Disk(16.3e-3) fill=Material(materials, :pp) wall=jacket(
                Material(materials, :pe); t = 2.5e-3) begin
                @at screened_core(:a) (0.0, spacing)
                @at screened_core(:b) (-spacing * cos(pi / 6), -spacing * sin(pi / 6))
                @at screened_core(:c) (spacing * cos(pi / 6), -spacing * sin(pi / 6))
            end
        end
    end

    """
    $(TYPEDSIGNATURES)

    Construct one screened cable resting on the bottom of an insulating conduit. The cable
    is a copper conductor under XLPE insulation, a copper sheath and a PE jacket. Air fills
    the PE conduit. `flatten` rejects the design, because it cannot reduce a group off the
    center of its container yet.
    """
    function offset_conduit_design()
        materials = MaterialsLibrary(add_defaults = true)
        copper = Material(materials, :copper)
        pe = Material(materials, :pe)
        cable = Stack(
            terminal(:core, core(copper; r = 4.0e-3),
                insulation(Material(materials, :xlpe); t = 3.4e-3)),
            terminal(:screen, sheath(copper; t = 0.1e-3)),
            jacket(pe; t = 1.0e-3))
        # The cable's 8.5 mm outer radius touches the 12 mm bore of the conduit.
        return @cable "offset-conduit" begin
            @pipe shape=Disk(12.0e-3) fill=Material(materials, :air) wall=jacket(
                pe; t = 2.0e-3) begin
                @at cable (0.0, -3.5e-3)
            end
        end
    end

    # One cable at a depth of 1 m. Each terminal is a phase of its own.
    function single_cable_system(design)
        connections = Dict(terminal => phase
            for (phase, terminal) in enumerate(design.terminal_order))
        return build(LineCableSystem, [design], [(0.0, -1.0)];
            system_id = design.cable_id, line_length = 1000.0, connections = [connections])
    end

    # The trefoil of `screened_cable_design()` cables that tutorial 2 builds.
    function trefoil_system()
        design = screened_cable_design()
        placements = @trefoil design spacing=70e-3 center=(0.0, -1.0) core=(1, 2, 3) sheath=0 jacket=0
        return build(LineCableSystem, placements; system_id = "18kV_1000mm2_trefoil",
            line_length = 1000.0)
    end

    # The bipole of `armored_cable_design()` cables that tutorial 3 builds.
    function bipole_system()
        design = armored_cable_design()
        return build(LineCableSystem, [
                (@at design (-0.5, -1.0) connections = (core = 1, sheath = 0, armor = 0)),
                (@at design (0.5, -1.0) connections = (core = 2, sheath = 0, armor = 0))];
            system_id = "525kV_1600mm2_bipole", line_length = 1000.0)
    end

    """
    $(TYPEDSIGNATURES)

    Write `system` on `earth` in `format` to a temporary file and return its text in pieces.
    A piece ends after each newline and between adjacent XML elements, and the pieces join
    back to the file. The equivalence check then reports the element that changed. A
    `temperature` \\[°C\\] or a `formulation` other than `nothing` reaches the exporter as
    its keyword.
    """
    function exported_text(format::Symbol, system, earth, temperature = nothing,
            formulation = nothing)
        extension = (atp = ".xml", pscad = ".pscx", tralin = ".f05")[format]
        keywords = merge(temperature === nothing ? (;) : (; temperature),
            formulation === nothing ? (;) : (; formulation))
        return mktempdir() do directory
            path = with_logger(NullLogger()) do
                export_data(format, system, earth;
                    file_name = joinpath(directory, "export" * extension), keywords...)
            end
            String.(split(read(path, String), r"(?<=\n)|(?<=>)(?=<)"))
        end
    end

    """
    $(TYPEDSIGNATURES)

    The reduction corpus contains observables of the reduction from a cable design to
    coaxial layers. The preservation corpus lacks them. Only the equivalence check runs this
    corpus. Each scenario starts with a function, and its positional arguments follow.

    - Export files. ATP, PSCAD and TRALIN write the three-phase system of the preservation
      corpus, the 18 kV trefoil, the 525 kV bipole and the three-core assembly. PSCAD also
      writes the trefoil at 90 °C, and the three formats write it at 90 °C under the lossy
      laws.
    - Computations. LineParameters at two frequencies and CableConstants run with the
      default formulations. The cables are the open wire screen, the same screen with a
      strip, the 18 kV cable, the 525 kV cable and the three-core assembly. Each cable runs
      as built and after `homogenize`. The open wire screen also runs with the boundary
      shunt model at a coarse resolution.
    - Temperature. The 18 kV cable with an aluminum tape among its copper screen wires, a
      terminal of mixed conductor materials, at 90 °C in both computations.
    - Rejections. The three-core cable with screened cores and the cable off the center of
      its conduit have no coaxial reduction. The corpus records the messages of `flatten`,
      of both computations and of the three exports.
    """
    function reduction_corpus()
        earth = homogeneous(rho = 100.0, eps_r = 10.0, mu_r = 1.0)
        frequencies = collect(10.0 .^ range(1, 5; length = 2))
        boundary = formula(:boundary;
            options = (resolution = (wire = 16, order = 8, quadrature = 64, modes = 128),))
        designs = (wire_screen = internal_shunt_test_design(; tapes = false),
            wire_strip = internal_shunt_test_design(), screened = screened_cable_design(),
            armored = armored_cable_design(), three_core = three_core_design())
        systems = (three_phase = three_phase_system(), trefoil = trefoil_system(),
            bipole = bipole_system(), three_core = single_cable_system(designs.three_core))
        scenarios = Pair{Symbol, Tuple}[]
        for (name, system) in pairs(systems), format in (:atp, :pscad, :tralin)
            push!(scenarios,
                Symbol(format, :_, name) => (exported_text, format, system, earth))
        end
        push!(scenarios,
            :pscad_trefoil_90 => (exported_text, :pscad, systems.trefoil, earth, 90.0))
        lossy = Formulation(insulation_admittance = :lossy, semicon_admittance = :lossy)
        for format in (:atp, :pscad, :tralin)
            push!(scenarios, Symbol(format, :_trefoil_lossy_90) =>
                (exported_text, format, systems.trefoil, earth, 90.0, lossy))
        end
        for (name, design) in pairs(designs), (suffix, built) in
                (("", design), ("_homogenized", homogenize(design)))
            problem = line_parameters_problem(single_cable_system(built); frequencies)
            push!(scenarios,
                Symbol(name, :_line_parameters, suffix) =>
                    (compute, problem, Formulation()),
                Symbol(name, :_cable_constants, suffix) =>
                    (compute, CableConstantsProblem(built), CableConstantsFormulation()))
            name === :wire_screen && push!(scenarios,
                Symbol(name, :_line_parameters_boundary, suffix) =>
                    (compute, problem, Formulation(shunt_model = boundary)),
                Symbol(name, :_cable_constants_boundary, suffix) => (compute,
                    CableConstantsProblem(built),
                    CableConstantsFormulation(shunt_model = boundary)))
        end
        mixed = screened_cable_design(tape_material = :aluminum)
        push!(scenarios,
            :mixed_screen_line_parameters_90 => (compute, LineParametersProblem(
                single_cable_system(mixed); temperature = 90.0, earth_props = earth,
                frequencies), Formulation()),
            :mixed_screen_cable_constants_90 => (compute,
                CableConstantsProblem(mixed; temperature = 90.0), CableConstantsFormulation()))
        rejected = (screened_three_core = screened_three_core_design(),
            offset_conduit = offset_conduit_design())
        for (name, design) in pairs(rejected)
            system = single_cable_system(design)
            push!(scenarios,
                Symbol(name, :_flatten) => (LineCableModels.DataModel.flatten, design),
                Symbol(name, :_line_parameters) =>
                    (compute, line_parameters_problem(system; frequencies), Formulation()),
                Symbol(name, :_cable_constants) => (compute, CableConstantsProblem(design),
                    CableConstantsFormulation()))
            for format in (:atp, :pscad, :tralin)
                push!(scenarios,
                    Symbol(format, :_, name) => (exported_text, format, system, earth))
            end
        end
        return NamedTuple(scenarios)
    end
end
