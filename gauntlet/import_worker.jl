using LineCableModels

length(ARGS) == 3 || error("usage: import_worker.jl SOURCE OUTPUT CASE_ID")
source, output, raw_id = abspath(ARGS[1]), abspath(ARGS[2]), ARGS[3]
isfile(source) || throw(ArgumentError("case source does not exist: $source"))

module ImportedGauntletCase end
problem = Base.include(ImportedGauntletCase, source)
problem isa LineCableModels.Engine.LineParametersProblem || throw(ArgumentError(
    "$source must evaluate to one concrete LineParametersProblem; got $(typeof(problem))",
))

system = problem.system
normalized_system = LineCableModels.build(
    LineCableModels.LineCableSystem,
    system.designs,
    system.positions;
    connections = system.connections,
    environment = system.environment,
    system_id = raw_id,
    line_length = system.line_length
)
normalized = LineCableModels.Engine.LineParametersProblem(
    normalized_system;
    temperature = problem.temperature,
    earth_props = problem.earth_props,
    frequencies = problem.frequencies,
    Γ = problem.Γ
)
LineCableModels.export_data(:json, normalized; file_name = output)
