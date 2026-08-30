TextDisplay.@showfields Finite "Finite" rule -> (field = rule.field,)
TextDisplay.@showfields Nonnegative "Nonnegative" rule -> (field = rule.field,)
TextDisplay.@showfields Positive "Positive" rule -> (field = rule.field,)
TextDisplay.@showfields IntegerField "IntegerField" rule -> (field = rule.field,)
TextDisplay.@showfields Less "Less" rule -> (left = rule.left, right = rule.right)
TextDisplay.@showfields LessEqual "LessEqual" rule -> (left = rule.left, right = rule.right)
TextDisplay.@showfields Greater "Greater" rule -> (left = rule.left, right = rule.right)
TextDisplay.@showfields GreaterEqual "GreaterEqual" rule -> (left = rule.left, right = rule.right)

TextDisplay.name(::Type{<:IsA}) = "IsA"
Base.summary(io::IO, ::IsA) = print(io, "Type validation rule")
function Base.show(io::IO, rule::IsA{T}) where {T}
    print(io, "IsA(field=:", rule.field, ", type=", nameof(T), ")")
end
Base.show(io::IO, ::MIME"text/plain", rule::IsA) = show(io, rule)

TextDisplay.name(::Type{<:OneOf}) = "OneOf"
Base.summary(io::IO, ::OneOf) = print(io, "Allowed-values validation rule")
function Base.show(io::IO, rule::OneOf)
    print(io, "OneOf(field=:", rule.field, ", values=", length(rule.values), ")")
end
function Base.show(io::IO, ::MIME"text/plain", rule::OneOf)
    get(io, :compact, false) && return show(io, rule)
    return TextDisplay.tree(io, "OneOf :$(rule.field)", Tuple(
        (label = sprint(show, value; context = :compact => true), noun = "values")
        for value in rule.values
    ); noun = "values")
end

TextDisplay.name(::Type{<:OwnerRule}) = "OwnerRule"
Base.summary(io::IO, rule::OwnerRule) = print(io, "Owner validation rule :", rule.name)
Base.show(io::IO, rule::OwnerRule) = print(io, "OwnerRule(:", rule.name, ")")
Base.show(io::IO, ::MIME"text/plain", rule::OwnerRule) = show(io, rule)
