# The shunt stage completes the boundary coefficients after flattening, for each
# formulation. The solved-domain inventory lives only for one shunt stage.
const ShuntDomainReport = NamedTuple{
    (:design, :terminals, :requested, :effective, :reason, :message),
    Tuple{Int, UnitRange{Int}, Symbol, Symbol, Symbol, String}}

function internal_shunt_response(selected::ShuntModelFormulation, design,
        geometry, T, methods, solutions, design_index)
    domains=internal_shunt_domains(design, geometry, T; design_index)
    return internal_shunt_response(selected, domains, methods, solutions)
end

function internal_shunt_response(
        selected::Formula{:equivalent},
        design::CableDesign, geometry, ::Type{T}, methods, solutions, design_index) where {T}
    requested = formula_id(selected)
    reports = ShuntDomainReport[(design_index, indices, requested,
        :equivalent, :selected,
        "Equivalent annular layers; no boundary extraction or audit")
        for indices in geometry.assembly_ranges]
    return (blocks = InternalShuntBlock{T}[],
        details =
        (requested, effective = :equivalent, solves = 0, domains = reports,
            diagnostics = InternalShuntDiagnostic[]))
end

function internal_shunt_response(selected::Formula{:boundary},
        domains::Vector{InternalShuntDomain{T}}, methods, solutions) where {T}
    blocks = InternalShuntBlock{T}[]
    diagnostics = InternalShuntDiagnostic[]
    reports = ShuntDomainReport[]
    for domain in domains
        try
            _shunt_lossless(methods) || throw(BoundarySolveError(:unsupported,
                (; design = domain.design, terminals = domain.terminals),
                "boundary shunt requires the built-in lossless dielectric laws"))
            # A solve depends on its physical inputs, the terminal of each wire and tape, and
            # its numerical controls.
            values = _shunt_values(domain, methods)
            topology = (length(domain.terminals), length(domain.left), length(domain.right),
                [wire.terminal for wire in domain.wires],
                [tape.terminal for tape in domain.tapes])
            controls = selected.options.data
            previous = findfirst(
                solution -> isequal(solution.controls, controls) &&
                            isequal(solution.topology, topology) &&
                            same_physical_state(solution.values, values),
                solutions)
            if previous === nothing
                # The admitted lossless laws are independent of frequency and
                # temperature. Their dielectric descriptors retain UQ sources.
                result = internal_shunt_response(selected, values, domain;
                    level = controls.resolution, integration = controls.integration,
                    audit = controls.audit)
                C = Matrix{T}(result.C)
                P = lu(C) \ Matrix{T}(I, size(C, 1), size(C, 1))
                diagnostic = result.diagnostic
                push!(solutions, (; values, topology, controls, C, P, diagnostic))
            else
                source = solutions[previous]
                C, P, diagnostic = source.C, source.P, source.diagnostic
            end
            push!(blocks, InternalShuntBlock(domain.assembly, domain.terminals, C, P))
            push!(diagnostics, diagnostic)
            push!(reports, (
                domain.design, domain.terminals, :boundary, :boundary, :resolved, ""))
        catch exception
            exception isa BoundarySolveError || rethrow()
            selected.parameters.fallback === :equivalent || rethrow()
            @warn "Boundary shunt replaced by the equivalent annular layer" design=domain.design terminals=domain.terminals reason=exception.category
            push!(reports,
                (domain.design, domain.terminals, :boundary,
                    :equivalent, exception.category, sprint(showerror, exception)))
        end
    end
    solved = Base.IdSet{Matrix{T}}()
    foreach(block -> push!(solved, block.C), blocks)
    effective = isempty(reports) ? :equivalent :
                all(r -> r.effective === :boundary, reports) ? :boundary :
                all(r -> r.effective === :equivalent, reports) ? :equivalent : :mixed
    return (blocks,
        details = (requested = :boundary, effective,
            solves = length(solved), domains = reports, diagnostics))
end
