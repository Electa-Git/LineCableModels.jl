@inline _owner_name(::Type{T}) where {T} = String(nameof(T))
@inline _field(value, name::Symbol) = getproperty(value, name)

function _real(owner::Type, field::Symbol, value)
    value isa Real || throw(ArgumentError(
        "[$(_owner_name(owner))] $field must be real; got $(typeof(value))",
    ))
    return value
end

function apply(rule::Finite, value, owner::Type)
    candidate = _real(owner, rule.field, _field(value, rule.field))
    isfinite(candidate) || throw(DomainError(
        candidate,
        "[$(_owner_name(owner))] $(rule.field) must be finite"
    ))
end

function apply(rule::Nonnegative, value, owner::Type)
    candidate = _real(owner, rule.field, _field(value, rule.field))
    candidate >= zero(candidate) || throw(DomainError(
        candidate,
        "[$(_owner_name(owner))] $(rule.field) must be nonnegative"
    ))
end

function apply(rule::Positive, value, owner::Type)
    candidate = _real(owner, rule.field, _field(value, rule.field))
    candidate > zero(candidate) || throw(DomainError(
        candidate,
        "[$(_owner_name(owner))] $(rule.field) must be positive"
    ))
end

function apply(rule::IntegerField, value, owner::Type)
    candidate = _field(value, rule.field)
    candidate isa Integer || throw(ArgumentError(
        "[$(_owner_name(owner))] $(rule.field) must be an integer; got $(typeof(candidate))",
    ))
end

function _ordered(value, owner::Type, left::Symbol, right::Symbol, predicate, relation)
    lhs = _real(owner, left, _field(value, left))
    rhs = _real(owner, right, _field(value, right))
    predicate(lhs, rhs) || throw(DomainError(
        (lhs, rhs),
        "[$(_owner_name(owner))] $left must be $relation $right"
    ))
end

function apply(rule::Less, value, owner::Type)
    _ordered(value, owner, rule.left, rule.right, <, "less than")
end
function apply(rule::LessEqual, value, owner::Type)
    _ordered(value, owner, rule.left, rule.right, <=, "less than or equal to")
end
function apply(rule::Greater, value, owner::Type)
    _ordered(value, owner, rule.left, rule.right, >, "greater than")
end
function apply(rule::GreaterEqual, value, owner::Type)
    _ordered(value, owner, rule.left, rule.right, >=, "greater than or equal to")
end

function apply(rule::IsA{T}, value, owner::Type) where {T}
    candidate = _field(value, rule.field)
    candidate isa T || throw(ArgumentError(
        "[$(_owner_name(owner))] $(rule.field) must be $T; got $(typeof(candidate))",
    ))
end

function apply(rule::OneOf, value, owner::Type)
    candidate = _field(value, rule.field)
    candidate in rule.values || throw(ArgumentError(
        "[$(_owner_name(owner))] $(rule.field) must be one of $(collect(rule.values)); " *
        "got $candidate",
    ))
end

function apply(rule::PhysicalFillLimit, value, owner::Type)
    count = _field(value, rule.count)
    count isa Integer || return nothing
    geometry = map(field -> _field(value, field), rule.geometry)
    all(candidate -> candidate isa Real, geometry) || return nothing
    limit = maxfill(owner, geometry...)
    count <= limit || throw(DomainError(
        count,
        "[$(_owner_name(owner))] $(rule.count) exceeds packing limit $limit"
    ))
end

"Apply `rules(owner)` to `value` and return `value` unchanged."
function check(owner::Type, value)
    foreach(rule -> apply(rule, value, owner), rules(owner))
    return value
end
