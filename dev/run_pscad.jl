# Include from an IDE/REPL or run with julia --project=. dev/run_pscad.jl.
# Uses the active environment. Each execution starts a fresh PSCAD calculation.
using LineCableModels

# Edit this filename to use another station configuration.
config_path = joinpath(@__DIR__, "local-pscad.toml")
station = PSCAD.RemoteConfig(config_path)

# Two insulated copper wires, buried 1 m deep and spaced 1 m apart.
copper = Material(:conductor, 1.72e-8, 1.0)
dielectric = Material(:insulator, 1e14, 2.3)
design = build(CableDesign, "toy-insulated-wire", terminal(:core,
    solid(copper, Disk(0.004)), insulation(dielectric; t=0.002)))
system = build(LineCableSystem, [design, design], [Pose2(0, -1), Pose2(1, -1)];
    connections=[Dict(:core => 1), Dict(:core => 2)])
problem = LineParametersProblem(system; earth_props=homogeneous(rho=100.0),
    frequencies=10.0 .^ range(log10(50.0), log10(1000.0); length=101))
selected = Formulation(:pscad; earth_impedance=:wedepohl1973)

@time result = compute(problem, selected;
    options=(remote=station, verbosity=(default=0, PSCAD=1)))
Zresult = Z(result) # Series impedance in ohm/m: terminal × terminal × frequency.
Yresult = Y(result) # Shunt admittance in S/m, with the same ordering.
println("PSCAD run: ", details(result).data.execution.source_run)
