function _lift_pose(build, values; combine::Symbol)
    combine in (:product, :zip) || throw(ArgumentError(
        "combine must be :product or :zip"
    ))
    grids = map(values) do value
        value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
    end
    return Gridspace{DataModel.Pose2}(build, grids; combine)
end

"Construct one physical cable pose, or a finite space of poses."
function at(; x = 0, y = 0, φ = 0, combine::Symbol = :product)
    values = (x, y, φ)
    any(value -> value isa Union{AbstractGrid, Gridspace}, values) &&
        return _lift_pose(DataModel.Pose2, values; combine)
    return DataModel.Pose2(values...)
end

function _trefoil_poses(x, y, spacing)
    spacing > zero(spacing) || throw(DomainError(
        spacing,
        "trefoil spacing must be positive"
    ))
    values = DataModel.trefoil_formation(x, y, spacing / 2)
    return DataModel.Pose2[
        DataModel.Pose2(values[1], values[2], 0),
        DataModel.Pose2(values[3], values[4], 0),
        DataModel.Pose2(values[5], values[6], 0)
    ]
end

"Return three physical poses in an equilateral trefoil formation."
function trefoil(; x = 0, y, spacing, combine::Symbol = :product)
    values = (x, y, spacing)
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{Vector{DataModel.Pose2}}(_trefoil_poses, grids; combine)
    end
    return _trefoil_poses(values...)
end

function _flat_poses(x, y, spacing, n, vertical)
    n isa Integer && !(n isa Bool) && n > 0 || throw(ArgumentError(
        "flat formation cardinality must be a positive integer"
    ))
    spacing > zero(spacing) || throw(DomainError(
        spacing,
        "flat formation spacing must be positive"
    ))
    return DataModel.Pose2[DataModel.Pose2(
                               vertical ? x : x + index * spacing,
                               vertical ? y - index * spacing : y,
                               0
                           ) for index in 0:(Int(n) - 1)]
end

function _flat(kind::Val, x, y, spacing, n; combine::Symbol)
    values = (x, y, spacing, n)
    build = kind === Val(:vertical) ?
            ((a, b, c, d) -> _flat_poses(a, b, c, d, true)) :
            ((a, b, c, d) -> _flat_poses(a, b, c, d, false))
    if any(value -> value isa Union{AbstractGrid, Gridspace}, values)
        grids = map(values) do value
            value isa Union{AbstractGrid, Gridspace} ? value : Grid((value,))
        end
        return Gridspace{Vector{DataModel.Pose2}}(build, grids; combine)
    end
    return build(values...)
end

"Return `n` physical poses in a horizontal flat formation."
function hflat(; x = 0, y = 0, spacing, n = 3, combine::Symbol = :product)
    _flat(Val(:horizontal), x, y, spacing, n; combine)
end

"Return `n` physical poses in a vertical flat formation."
function vflat(; x = 0, y = 0, spacing, n = 3, combine::Symbol = :product)
    _flat(Val(:vertical), x, y, spacing, n; combine)
end
