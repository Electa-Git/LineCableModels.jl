@testmodule TestFixtures begin
    include("scenarios.jl")
    using .CurrentScenarios
end

@testsnippet CableSystemFixture begin
    cable_system=TestFixtures.three_phase_system()
    problem_atp=TestFixtures.line_parameters_problem(cable_system)
    earth_props=problem_atp.earth_props
    freqs=problem_atp.frequencies
    num_phases=LineCableModels.nphases(cable_system)
end
