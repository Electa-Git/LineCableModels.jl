# Boundary coefficients are completed while constructing coaxial blueprints.
# The solved-domain inventory lives only for this flattening call.
const ShuntDomainReport = NamedTuple{
    (:design, :terminals, :requested, :effective, :reason, :message),
    Tuple{Int, UnitRange{Int}, Symbol, Symbol, Symbol, String}}

"""
$(TYPEDSIGNATURES)

Construct cable blueprints for each selected formulation. Identical local
selections share the completed blueprints; equivalent lossless domains share
their coefficient matrices. Earth-return choices do not enter this calculation.
The returned outer vector follows formulation order, and each inner vector
follows design order. Sharing is confined to this construction call.
"""
function flatten(engine::LineCableModelsCoaxial, designs::AbstractVector,
        ::Type{T}, formulations::AbstractVector{<:AbstractFormulation}) where {T <: Real}
    solutions = NamedTuple[]
    selections = NamedTuple[]
    blueprints = Vector{CableBlueprint{T}}[]
    for formulation in formulations
        methods = formulation.methods
        selected = formula_id(methods.shunt_model) === :boundary ?
                   methods[(
            :shunt_model, :insulation_admittance, :semicon_admittance)] :
                   methods[(:shunt_model,)]
        previous = findfirst(value -> isequal(value, selected), selections)
        current = previous === nothing ?
                  CableBlueprint{T}[flatten(engine, design, T, selected, solutions, index)
                                    for (index, design) in pairs(designs)] :
                  blueprints[previous]
        push!(selections, selected)
        push!(blueprints, current)
    end
    return blueprints
end

function internal_shunt_response(
        selected::Union{ShuntModel.Formula{:default}, ShuntModel.Formula{:coaxial}},
        domains::Vector{InternalShuntDomain{T}}, methods, solutions) where {T}
    requested = formula_id(selected)
    reports = ShuntDomainReport[(d.design, d.terminals, requested,
                                    :coaxial, :selected, "") for d in domains]
    return (blocks = InternalShuntBlock{T}[],
        details =
        (requested, effective = :coaxial, solves = 0, domains = reports,
            diagnostics = InternalShuntDiagnostic[]))
end

function internal_shunt_response(selected::ShuntModel.Formula{:boundary},
        domains::Vector{InternalShuntDomain{T}}, methods, solutions) where {T}
    blocks = InternalShuntBlock{T}[]
    diagnostics = InternalShuntDiagnostic[]
    reports = ShuntDomainReport[]
    for domain in domains
        try
            _shunt_lossless(methods) || throw(BoundarySolveError(:unsupported,
                (; design = domain.design, terminals = domain.terminals),
                "boundary shunt requires the built-in lossless dielectric laws without overrides"))
            previous = findfirst(
                value -> isequal(value.methods, methods) &&
                         _shunt_domain_equal(value.domain, domain),
                solutions)
            if previous === nothing
                # The admitted lossless laws are independent of frequency and
                # temperature. Their dielectric descriptors retain UQ sources.
                result = internal_shunt_response(_shunt_values(domain, methods), domain;
                    level = selected.options.resolution,
                    integration = selected.options.integration, audit = selected.options.audit)
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
