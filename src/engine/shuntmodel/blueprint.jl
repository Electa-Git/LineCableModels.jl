# The shunt stage completes the internal potential coefficients after flattening, for each
# formulation. Its records of the solved groups live only for one stage.
const ShuntDomainReport = NamedTuple{
    (:design, :terminals, :requested, :reason, :message),
    Tuple{Int, UnitRange{Int}, Symbol, Symbol, String}}

function dielectric_groups(selected::Formula{:equivalent}, design::CableDesign,
        blueprint::CableBlueprint{T}, formulas, design_index) where {T}
    reports = ShuntDomainReport[(design_index, indices, formula_id(selected), :selected,
        "Equivalent annular layers; no boundary extraction")
        for indices in blueprint.assembly_ranges]
    return (groups = DielectricGroup{T}[], reports)
end

function dielectric_groups(selected::Formula{:boundary}, design::CableDesign,
        blueprint::CableBlueprint{T}, formulas, design_index) where {T}
    extraction = _shunt_extract(design, blueprint, T, design_index, formulas)
    validate(blueprint, selected, extraction, design)
    (; groups) = extraction
    isempty(groups) || validate(formulas, selected, design)
    reports = ShuntDomainReport[(group.design, group.terminals, :boundary, :resolved,
        (first(group.terminals) == first(group.assembly) ?
         "Equivalent inner core. " : "Closed inner anchor. ") *
        "Outer equipotential $(group.reference) [m].") for group in groups]
    for assembly in blueprint.assembly_ranges
        any(group -> group.assembly == assembly, groups) && continue
        push!(reports, (design_index, assembly, :boundary, :radial,
            "No discrete outer arrangement; exact radial dielectric solution. " *
            "The inner core retains the equivalent-core assumption."))
    end
    return (; groups, reports)
end

# Every conductor outside the inner core is a closed concentric tube or lies in a group.
# Otherwise the reason is the one that extraction gave for its dielectric group.
function validate(blueprint::CableBlueprint, ::Formula{:boundary}, extraction,
        design::CableDesign)
    for assembly in blueprint.assembly_ranges
        x, y = blueprint.positions[blueprint.conductors[first(assembly)].assembly]
        for index in (first(assembly)+1):last(assembly)
            _shunt_closed(design, index, x, y, blueprint.conductors[index].r_ex) !== nothing && continue
            any(group -> group.assembly == assembly &&
                first(group.terminals) < index <= last(group.terminals), extraction.groups) && continue
            rejection = findfirst(candidate -> index in candidate.terminals, extraction.rejected)
            rejection === nothing || throw(extraction.rejected[rejection].reason)
            throw(BoundarySolveError(:unsupported,
                (; cable = design.cable_id, terminal = blueprint.conductors[index].terminal),
                "boundary shunt cannot represent this conductor: it is neither a closed " *
                "concentric tube nor a wire or tape in an annular dielectric fill"))
        end
    end
    return blueprint
end

function validate(formulas::NamedTuple, ::Formula{:boundary}, design::CableDesign)
    _shunt_lossless(formulas) || throw(BoundarySolveError(:unsupported,
        (; cable = design.cable_id),
        "boundary shunt requires the built-in lossless dielectric laws"))
    return formulas
end
