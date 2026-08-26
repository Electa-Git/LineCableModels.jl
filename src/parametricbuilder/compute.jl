struct CableConstantsMaterializer end

function (::CableConstantsMaterializer)(design, separation, earth_resistivity)
    Engine.CableConstantsProblem(design; separation, earth_resistivity)
end

function _cable_constants_space(
        design,
        separation,
        earth_resistivity;
        combine::Symbol
)
    values = (design, separation, earth_resistivity)
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{Engine.CableConstantsProblem}(
        CableConstantsMaterializer(),
        grids;
        combine
    )
end

function Engine.CableConstantsProblem(
        design::Gridspace{<:DataModel.CableDesign};
        separation = nothing,
        earth_resistivity = 100.0,
        combine::Symbol = :product
)
    combine in (:product, :zip) ||
        throw(ArgumentError("combine must be :product or :zip"))
    return _cable_constants_space(
        design, separation, earth_resistivity; combine)
end

function Engine._cable_constants_problem(
        design::DataModel.CableDesign,
        separation,
        earth_resistivity
)
    values = (design, separation, earth_resistivity)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        return invoke(
            Engine._cable_constants_problem,
            Tuple{DataModel.CableDesign, Union{Nothing, Real}, Real},
            design,
            separation,
            earth_resistivity
        )
    return _cable_constants_space(
        design, separation, earth_resistivity; combine = :product)
end

_append_result(::Nothing, value) = typeof(value)[value]

_computation_owner(formulation) = Base.typename(typeof(formulation)).wrapper

function _append_result(values::Vector{T}, value) where {T}
    value isa T && (push!(values, value); return values)
    W = typejoin(T, typeof(value))
    W === Any && throw(ArgumentError(
        "Gridspace results do not share one core result grammar",
    ))
    widened = Vector{W}(undef, length(values) + 1)
    copyto!(widened, values)
    widened[end] = value
    return widened
end

function _traverse(problem::ParametricProblem, formulation, result_type)
    length(problem.space) == 0 && throw(ArgumentError(
        "higher-order problem space must contain at least one core problem",
    ))
    values = nothing
    retained = nothing
    for point in points(problem.space)
        core_problem = materialize(point)
        core_result = compute(
            core_problem,
            formulation.inner;
            options = problem.options
        )
        values = _append_result(values, core_result)
        if formulation.options.retain_details
            record = computation_details(
                Val(_computation_owner(formulation.inner)),
                core_result
            )
            retained = _append_result(retained, record)
        end
    end
    typed_values = values === nothing ? AbstractProblemResult[] : values
    retained_details = formulation.options.retain_details ?
                       (points = something(retained, NamedTuple[]),) : (;)
    return result_type(formulation, typed_values, retained_details)
end

function compute(problem::ParametricProblem, formulation::Combinatorial)
    _traverse(problem, formulation, ParametricResult)
end
