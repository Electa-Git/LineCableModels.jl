"""
    LineCableModels.InputValidation

Check materialized inputs before they enter another construction or numerical
operation. Concrete input types implement [`validate`](@ref) beside their type
definitions. A successful method returns its argument unchanged; rejected
inputs raise an actionable native exception.
"""
module InputValidation

export validate

using DocStringExtensions: TYPEDSIGNATURES
using RequiredInterfaces: @required
import ..Grammar: AbstractProblemDefinition

"""
$(TYPEDSIGNATURES)

Check that a materialized input satisfies the invariants of its concrete type
and return that same input unchanged.

# Errors

- Throws a native exception identifying the invalid field, value, and required
  condition.
- Throws `RequiredInterfaces.NotImplementedError` when a concrete problem type
  does not implement input validation.
"""
function validate end

@required AbstractProblemDefinition begin
    validate(::AbstractProblemDefinition)
end

end # module InputValidation
