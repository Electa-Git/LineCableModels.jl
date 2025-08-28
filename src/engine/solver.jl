function compute!(problem::LineParametersProblem,
    formulation::CoaxialFormulation)

    # Flatten the system into 1D arrays for computation
    workspace = init_workspace(problem, formulation)

    T = eltype(workspace.freq)
    ComplexT = Complex{T}

    n_frequencies = workspace.n_frequencies
    n_phases = workspace.n_phases
    Z = zeros(ComplexT, n_phases, n_phases, n_frequencies)
    Y = zeros(ComplexT, n_phases, n_phases, n_frequencies)

    @show workspace
    @show typeof(workspace)

    @info "Starting line parameters computation"
    line_params = nothing

    @info "Line parameters computation completed successfully"
    return workspace, line_params
end