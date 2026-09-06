"""
    CssEditor

Describe a browser-local CSS editor; this descriptor never binds Julia state.
Unsupported CSS properties remain read-only regardless of component overrides.
"""
struct CssEditor
    "Editor family: length, number, choice, color, text, or readonly."
    kind::Symbol
    "Allowed choices, including CSS keywords."
    choices::Vector{String}
    "Allowed units for a simple numeric value."
    units::Vector{String}
    "Optional lower bound for a simple numeric value."
    minimum::Union{Nothing,Float64}
    "Optional upper bound for a simple numeric value."
    maximum::Union{Nothing,Float64}
    "Positive spinner increment."
    step::Float64
end

"""
    CssEditor(kind; choices=[], units=[], minimum=nothing, maximum=nothing, step=0.1)

Construct editor hints for authored CSS declarations.

# Arguments
- `kind`: One of `:length`, `:number`, `:choice`, `:color`, `:text`, `:readonly`.

# Keywords
- `choices`, `units`: Allowed CSS keywords and unit suffixes.
- `minimum`, `maximum`: Optional bounds on simple numeric values.
- `step`: Positive, finite spinner increment; defaults to `0.1`.

# Returns
- A `CssEditor` descriptor.

# Errors
- `ArgumentError` for an unknown kind, invalid bounds, or invalid increment.
"""
function CssEditor(kind::Symbol; choices=String[], units=String[],
        minimum=nothing, maximum=nothing, step::Real=0.1)
    kind in (:length, :number, :choice, :color, :text, :readonly) ||
        throw(ArgumentError("unknown CSS editor kind: $kind"))
    isfinite(step) && step > 0 || throw(ArgumentError("CSS editor step must be positive and finite"))
    for value in (minimum, maximum)
        isnothing(value) || isfinite(value) || throw(ArgumentError("CSS editor bounds must be finite"))
    end
    isnothing(minimum) || isnothing(maximum) || minimum <= maximum ||
        throw(ArgumentError("CSS editor minimum exceeds maximum"))
    return CssEditor(kind, string.(collect(choices)), string.(collect(units)),
        isnothing(minimum) ? nothing : Float64(minimum),
        isnothing(maximum) ? nothing : Float64(maximum), Float64(step))
end

"""
    css_editors(component)

Return property-name → `CssEditor` overrides for an owned component. Specialize
this method to constrain editors without duplicating stylesheet declarations.
The default is an empty dictionary. Code parameters and bindings are unaffected.
"""
css_editors(::Any) = Dict{String,CssEditor}()

function css_editor_payload(editor::CssEditor)
    return Dict("kind" => string(editor.kind), "choices" => editor.choices,
        "units" => editor.units, "minimum" => editor.minimum,
        "maximum" => editor.maximum, "step" => editor.step)
end

const CSS_EDITORS = let catalogue = Dict{String,CssEditor}()
    for name in split("width height min-width max-width min-height max-height gap row-gap column-gap padding padding-top padding-right padding-bottom padding-left border-width border-top-width border-right-width border-bottom-width border-left-width border-radius font-size flex-basis")
        catalogue[name] = CssEditor(:length; units=["px", "rem", "em", "%", "vw", "vh", "dvh", "ch"], minimum=0)
    end
    for name in split("margin margin-top margin-right margin-bottom margin-left top right bottom left letter-spacing word-spacing")
        catalogue[name] = CssEditor(:length; units=["px", "rem", "em", "%", "vw", "vh", "dvh", "ch"])
    end
    catalogue["opacity"] = CssEditor(:number; minimum=0, maximum=1, step=0.05)
    catalogue["line-height"] = CssEditor(:length; units=["", "px", "rem", "em", "%"], minimum=0)
    for name in ("flex-grow", "flex-shrink")
        catalogue[name] = CssEditor(:number; minimum=0)
    end
    catalogue["font-weight"] = CssEditor(:choice; choices=["normal", "bold", "400", "500", "600", "700"])
    for name in split("color background-color border-color border-top-color border-right-color border-bottom-color border-left-color outline-color fill stroke")
        catalogue[name] = CssEditor(:color)
    end
    for (name, choices) in (
        "display" => ["block", "inline", "inline-block", "flex", "inline-flex", "grid", "inline-grid", "none"],
        "flex-direction" => ["row", "column", "row-reverse", "column-reverse"],
        "flex-wrap" => ["nowrap", "wrap", "wrap-reverse"],
        "align-items" => ["normal", "stretch", "start", "end", "center", "baseline"],
        "justify-content" => ["normal", "start", "end", "center", "space-between", "space-around", "space-evenly"],
        "text-align" => ["start", "end", "left", "right", "center", "justify"],
        "white-space" => ["normal", "nowrap", "pre", "pre-wrap", "pre-line", "break-spaces"],
        "visibility" => ["visible", "hidden"],
        "box-sizing" => ["border-box", "content-box"],
        "border-style" => ["none", "solid", "dashed", "dotted", "double"],
    )
        catalogue[name] = CssEditor(:choice; choices)
    end
    for name in ("overflow", "overflow-x", "overflow-y")
        catalogue[name] = CssEditor(:choice; choices=["visible", "hidden", "clip", "auto", "scroll"])
    end
    for name in ("grid-template-columns", "grid-template-rows", "grid-auto-columns", "grid-auto-rows", "border", "outline")
        catalogue[name] = CssEditor(:text)
    end
    catalogue
end
