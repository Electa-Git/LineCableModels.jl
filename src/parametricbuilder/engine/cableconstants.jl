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

function Engine.cable_constants_problem(
        design::DataModel.CableDesign,
        separation,
        earth_resistivity
)
    values = (design, separation, earth_resistivity)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) ||
        return invoke(
            Engine.cable_constants_problem,
            Tuple{DataModel.CableDesign, Union{Nothing, Real}, Real},
            design,
            separation,
            earth_resistivity
        )
    return _cable_constants_space(
        design, separation, earth_resistivity; combine = :product)
end
