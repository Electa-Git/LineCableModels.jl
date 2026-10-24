"""
    LineCableModels.Commons

Define the grammar that the other modules extend with methods.
"""
module Commons

export build, resolve, validate

"""
    validate(subject, context...)

Reject an invalid input before its consumer uses it. Return `subject`, the first positional
argument, unchanged.

`subject` is user input or an input derived from it. The optional context says what the
subject is checked for: its consumer, such as the problem, formulation, physical model or
result being built, or, for the result of a law or record, the function that produced it.
Data that the check needs follows the context. A method can warn on a value that it can
compute but that is questionable, such as one outside a formula's range of validity.

# Errors

- Throws a native exception that identifies the invalid field, value, and required
  condition.
"""
function validate end

"""
    resolve(context, declaration)

Return the geometry that `declaration` takes in `context`. The outline of a part gives its
geometry to the layer or part declared around it, and a placement locates a cross-section in
its frame.
"""
function resolve end

"""
    build(Target, declarations...)

Construct one completed `Target` from its declarations, such as the design of a cable from
its declared parts.
"""
function build end

end
