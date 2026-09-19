# Boundary coefficients are completed while constructing coaxial blueprints.
# The solved-domain inventory lives only for this flattening call.
const ShuntDomainReport = NamedTuple{
    (:design, :terminals, :requested, :effective, :reason, :message),
    Tuple{Int, UnitRange{Int}, Symbol, Symbol, Symbol, String}}

# The local model owns its blueprint dependency closure. Other shunt models
# retain the supplied dielectric selections unless they declare a narrower set.
function blueprint_dependencies(::ShuntModelFormulation, methods)
    methods[(:shunt_model, :insulation_admittance, :semicon_admittance)]
end
blueprint_dependencies(::Formula{:coaxial}, methods) = methods[(:shunt_model,)]

function internal_shunt_response(selected::ShuntModelFormulation, design,
        geometry, T, methods, solutions, design_index)
    domains=internal_shunt_domains(design, geometry, T; design_index)
    return internal_shunt_response(selected, domains, methods, solutions)
end

function internal_shunt_response(
        selected::Formula{:coaxial},
        design::CableDesign, geometry, ::Type{T}, methods, solutions, design_index) where {T}
    requested = formula_id(selected)
    reports = ShuntDomainReport[(design_index, indices, requested,
        :coaxial, :selected, "Coaxial annuli; no boundary extraction or audit")
        for indices in geometry.assembly_ranges]
    return (blocks = InternalShuntBlock{T}[],
        details =
        (requested, effective = :coaxial, solves = 0, domains = reports,
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
            previous = findfirst(
                value -> isequal(value.methods, methods) &&
                         _shunt_domain_equal(value.domain, domain),
                solutions)
            if previous === nothing
                # The admitted lossless laws are independent of frequency and
                # temperature. Their dielectric descriptors retain UQ sources.
                result = internal_shunt_response(
                    selected, _shunt_values(domain, methods), domain;
                    level = selected.options.data.resolution,
                    integration = selected.options.data.integration, audit = selected.options.data.audit)
                C = Matrix{T}(result.C)
                P = lu(C) \ Matrix{T}(I, size(C, 1), size(C, 1))
                diagnostic = result.diagnostic
                push!(solutions, (; domain, methods, C, P, diagnostic))
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
            selected.parameters.fallback === :coaxial || rethrow()
            @warn "Boundary shunt replaced by coaxial annuli" design=domain.design terminals=domain.terminals reason=exception.category
            push!(reports,
                (domain.design, domain.terminals, :boundary,
                    :coaxial, exception.category, sprint(showerror, exception)))
        end
    end
    solved = Base.IdSet{Matrix{T}}()
    foreach(block -> push!(solved, block.C), blocks)
    effective = isempty(reports) ? :coaxial :
                all(r -> r.effective === :boundary, reports) ? :boundary :
                all(r -> r.effective === :coaxial, reports) ? :coaxial : :mixed
    return (blocks,
        details = (requested = :boundary, effective,
            solves = length(solved), domains = reports, diagnostics))
end
