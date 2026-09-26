# Run in an environment containing LineCableModels and GLMakie:
#     include("test/manual/calculations/run_two_bare_wires.jl")
# The active project and LOAD_PATH belong to the caller and remain unchanged.
using LineCableModels, GLMakie

# Two bare copper wires: 42.5 mm radius, 1 m spacing, buried 1 m deep.
# Geometry, temperature and static earth inputs match the two_bare_wires case.
materials = MaterialsLibrary(add_defaults=true)
copper = Material(materials, :copper)
design = build(CableDesign, "two_bare_wires",
    Stack(Group(:core, Region(:core_metal, Disk(0.0425), copper))))
system = build(LineCableSystem, [design, design],
    [Pose2(0.0, -1.0), Pose2(1.0, -1.0)];
    connections=[Dict(:core => 1), Dict(:core => 2)],
    system_id="two_bare_wires", line_length=1.0)
earth = homogeneous(rho=0.1, eps_r=1.0, mu_r=1.0)
frequency = 10.0 .^ range(-1, 7; length=17) # Hz
problem = LineParametersProblem(system;
    temperature=20.0, earth_props=earth, frequencies=frequency)

# Both calculations use Unified earth return. Compare constant earth properties
# with the implemented Longmire–Smith dispersion model, keeping static inputs
# fixed. The legend identifies the differing soil law automatically.
formulations = [Formulation(
    earth_impedance=:unified, earth_admittance=:unified,
    earth_properties=soil_law,
    options=(reduce_bundle=false, kron_reduction=false, ideal_transposition=false))
    for soil_law in (:constant, :longmire1975)]
results = @time compute(problem, formulations)
plots = LineCableModels.plot(results;
    ydata=(R, X, G, B), length_unit=:base, backend=:gl)
nothing
