abstract type AbstractFormControl end

const FORM_NAME_PATTERN = r"^[a-z][a-z0-9_]*$"

function form_name(value, kind)
    name = Symbol(value)
    occursin(FORM_NAME_PATTERN, string(name)) || throw(ArgumentError(
        "$kind must use lowercase snake_case and start with a letter: $value"
    ))
    return name
end

control_name(control::AbstractFormControl) = control.name
control_id(control::AbstractFormControl) = control.id
control_required(control::AbstractFormControl) = control.required
reset_control!(control::AbstractFormControl) = (control.value[] = control.initial; control)

function normalized_options(options)
    resolved = Pair{String,String}[]
    for option in options
        push!(resolved, option isa Pair ? string(first(option)) => string(last(option)) :
            string(option) => string(option))
    end
    isempty(resolved) && throw(ArgumentError("a selection control requires an option"))
    allunique(first.(resolved)) || throw(ArgumentError("option values must be unique"))
    return resolved
end

function assert_selected(selected, options; multiple=false)
    allowed = Set(first.(options))
    values = multiple ? string.(collect(selected)) : [string(selected)]
    all(in(allowed), values) || throw(ArgumentError("selected value is not an option"))
    return multiple ? values : only(values)
end

function native_attributes(control; input_type=nothing)
    attributes = Dict{Symbol,Any}(
        :id => control.id,
        :name => string(control.name),
        :required => control.required,
        :disabled => control.disabled,
    )
    isnothing(input_type) || (attributes[:type] = input_type)
    return attributes
end

"""Create a single-line text control with observable session-local state."""
struct TextInput <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable DOM identifier."
    id::String
    "Current text."
    value::Observable{String}
    "Value restored by form reset."
    initial::String
    "Placeholder shown while empty."
    placeholder::String
    "Whether browser validation requires a value."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function TextInput(name; value::AbstractString="", id::AbstractString=string(name),
        placeholder::AbstractString="", required::Bool=false, disabled::Bool=false)
    initial = string(value)
    return TextInput(form_name(name, "text-input name"), string(id), Observable(initial),
        initial, string(placeholder), required, disabled)
end

function Bonito.jsrender(session::Session, control::TextInput)
    node = DOM.input(; native_attributes(control; input_type="text")...,
        value=control.value, placeholder=control.placeholder,
        class="lc-control-input lc-form-control",
        oninput=js"event => $(control.value).notify(event.currentTarget.value)")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

"""Create a multiline text control with observable session-local state."""
struct TextAreaInput <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable DOM identifier."
    id::String
    "Current text."
    value::Observable{String}
    "Value restored by form reset."
    initial::String
    "Placeholder shown while empty."
    placeholder::String
    "Visible text rows."
    rows::Int
    "Whether browser validation requires a value."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function TextAreaInput(name; value::AbstractString="", id::AbstractString=string(name),
        placeholder::AbstractString="", rows::Integer=4, required::Bool=false,
        disabled::Bool=false)
    rows > 0 || throw(ArgumentError("text-area rows must be positive"))
    initial = string(value)
    return TextAreaInput(form_name(name, "text-area name"), string(id), Observable(initial),
        initial, string(placeholder), Int(rows), required, disabled)
end

function Bonito.jsrender(session::Session, control::TextAreaInput)
    node = DOM.textarea(; native_attributes(control)..., rows=control.rows,
        value=control.value,
        placeholder=control.placeholder, class="lc-control-textarea lc-form-control",
        oninput=js"event => $(control.value).notify(event.currentTarget.value)")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

"""
    SecretInput(name; placeholder="", autocomplete="current-password", required=false)

Create a masked browser-local input. Its value is deliberately absent from the
Julia object, reactive bindings, and X-ray metadata. A containing [`Form`](@ref)
may transmit it only when the user explicitly submits that form.
"""
struct SecretInput <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable DOM identifier."
    id::String
    "Placeholder shown while empty."
    placeholder::String
    "Browser autocomplete policy."
    autocomplete::String
    "Whether browser validation requires a value."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
    "Whether the reveal control is available."
    revealable::Bool
end

function SecretInput(name; id::AbstractString=string(name), placeholder::AbstractString="",
        autocomplete::AbstractString="current-password", required::Bool=false,
        disabled::Bool=false, revealable::Bool=true)
    return SecretInput(form_name(name, "secret-input name"), string(id), string(placeholder),
        string(autocomplete), required, disabled, revealable)
end

reset_control!(::SecretInput) = nothing

function secret_icon()
    return SVG.svg(
        SVG.path(; d="M2.5 12s3.5-6 9.5-6 9.5 6-3.5 6-9.5 6S2.5 12 2.5 12Z"),
        SVG.circle(; cx="12", cy="12", r="2.5"); viewBox="0 0 24 24", fill="none",
        stroke="currentColor", var"stroke-width"="1.7", var"stroke-linecap"="round",
        var"stroke-linejoin"="round", var"aria-hidden"="true")
end

function Bonito.jsrender(session::Session, control::SecretInput)
    input = DOM.input(; native_attributes(control; input_type="password")...,
        placeholder=control.placeholder, autocomplete=control.autocomplete,
        class="lc-control-input lc-form-control lc-secret-input")
    reveal = control.revealable ? DOM.button(secret_icon(); type="button",
        class="lc-secret-reveal", title="Show secret", var"aria-label"="Show secret") : nothing
    behavior = control.revealable ? DOM.script(js"""
    (() => {
        const field = $(input);
        const button = $(reveal);
        button.addEventListener('pointerdown', event => event.preventDefault());
        button.addEventListener('click', () => {
            const visible = field.type === 'text';
            field.type = visible ? 'password' : 'text';
            button.title = visible ? 'Show secret' : 'Hide secret';
            button.setAttribute('aria-label', button.title);
        });
    })();
    """) : nothing
    node = DOM.div(input, reveal, behavior; class="lc-secret-control")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

"""Create a mutually exclusive set of labelled radio choices."""
struct RadioGroup <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable group identifier."
    id::String
    "Value-label option pairs."
    options::Vector{Pair{String,String}}
    "Currently selected value."
    value::Observable{String}
    "Value restored by form reset."
    initial::String
    "Layout direction."
    orientation::Symbol
    "Whether browser validation requires a choice."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function RadioGroup(name, options; selected=nothing, id::AbstractString=string(name),
        orientation::Symbol=:vertical, required::Bool=false, disabled::Bool=false)
    orientation in (:horizontal, :vertical) || throw(ArgumentError(
        "radio orientation must be :horizontal or :vertical"))
    resolved = normalized_options(options)
    initial = assert_selected(isnothing(selected) ? first(first(resolved)) : selected, resolved)
    return RadioGroup(form_name(name, "radio-group name"), string(id), resolved,
        Observable(initial), initial, orientation, required, disabled)
end

control_id(control::RadioGroup) = "$(control.id)-1"

function render_choice_group(session, control, root_class, role)
    orientation = hasfield(typeof(control), :orientation) ? control.orientation : :horizontal
    nodes = map(enumerate(control.options)) do (index, option)
        option_value, label = option
        input_id = "$(control.id)-$index"
        checked = map(session, control.value) do current
            current == option_value
        end
        input = DOM.input(; id=input_id, name=string(control.name), value=option_value,
            type="radio", checked, required=control.required && index == 1,
            disabled=control.disabled,
            onchange=js"event => event.currentTarget.checked && $(control.value).notify(event.currentTarget.value)")
        DOM.label(input, DOM.span(label); var"for"=input_id, class="$root_class-choice")
    end
    node = DOM.div(nodes...; id=control.id,
        class="$root_class is-$orientation", role)
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

Bonito.jsrender(session::Session, control::RadioGroup) =
    render_choice_group(session, control, "lc-radio-group", "radiogroup")

"""Create a compact mutually exclusive segmented selector."""
struct SegmentedControl <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable group identifier."
    id::String
    "Value-label option pairs."
    options::Vector{Pair{String,String}}
    "Currently selected value."
    value::Observable{String}
    "Value restored by form reset."
    initial::String
    "Whether browser validation requires a choice."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function SegmentedControl(name, options; selected=nothing, id::AbstractString=string(name),
        required::Bool=false, disabled::Bool=false)
    resolved = normalized_options(options)
    initial = assert_selected(isnothing(selected) ? first(first(resolved)) : selected, resolved)
    return SegmentedControl(form_name(name, "segmented-control name"), string(id), resolved,
        Observable(initial), initial, required, disabled)
end

control_id(control::SegmentedControl) = "$(control.id)-1"
Bonito.jsrender(session::Session, control::SegmentedControl) =
    render_choice_group(session, control, "lc-segmented-control", "radiogroup")

"""Create a native single-selection control suitable for long option lists."""
struct ComboBox <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable DOM identifier."
    id::String
    "Value-label option pairs."
    options::Vector{Pair{String,String}}
    "Currently selected value."
    value::Observable{String}
    "Value restored by form reset."
    initial::String
    "Whether browser validation requires a choice."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function ComboBox(name, options; selected=nothing, id::AbstractString=string(name),
        required::Bool=false, disabled::Bool=false)
    resolved = normalized_options(options)
    initial = assert_selected(isnothing(selected) ? first(first(resolved)) : selected, resolved)
    return ComboBox(form_name(name, "combo-box name"), string(id), resolved,
        Observable(initial), initial, required, disabled)
end

function Bonito.jsrender(session::Session, control::ComboBox)
    options = map(control.options) do option
        option_value, label = option
        selected = map(session, control.value) do current
            current == option_value
        end
        DOM.option(label; value=option_value, selected)
    end
    node = DOM.select(options...; native_attributes(control)...,
        class="lc-control-select lc-form-control",
        onchange=js"event => $(control.value).notify(event.currentTarget.value)")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

"""Create a native multiple-selection control with observable selected values."""
struct MultiSelect <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable DOM identifier."
    id::String
    "Value-label option pairs."
    options::Vector{Pair{String,String}}
    "Currently selected values."
    value::Observable{Vector{String}}
    "Values restored by form reset."
    initial::Vector{String}
    "Visible option rows."
    size::Int
    "Whether browser validation requires a choice."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function MultiSelect(name, options; selected=String[], id::AbstractString=string(name),
        size::Integer=4, required::Bool=false, disabled::Bool=false)
    size > 1 || throw(ArgumentError("multi-select size must exceed one"))
    resolved = normalized_options(options)
    initial = assert_selected(selected, resolved; multiple=true)
    return MultiSelect(form_name(name, "multi-select name"), string(id), resolved,
        Observable(initial), copy(initial), Int(size), required, disabled)
end

function Bonito.jsrender(session::Session, control::MultiSelect)
    options = map(control.options) do option
        option_value, label = option
        selected = map(session, control.value) do current
            option_value in current
        end
        DOM.option(label; value=option_value, selected)
    end
    wire = Observable(String(JSON3.write(control.value[])))
    on(session, wire) do payload
        control.value[] = JSON3.read(payload, Vector{String})
        return nothing
    end
    node = DOM.select(options...; native_attributes(control)..., multiple=true,
        size=control.size, class="lc-control-select lc-form-control lc-multi-select",
        onchange=js"event => $(wire).notify(JSON.stringify(Array.from(event.currentTarget.selectedOptions, option => option.value)))")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

"""Create a bounded numeric input with an adjacent engineering-unit label."""
struct UnitNumberInput <: AbstractFormControl
    "Form payload key."
    name::Symbol
    "Stable DOM identifier."
    id::String
    "Current numeric value."
    value::Observable{Float64}
    "Value restored by form reset."
    initial::Float64
    "Smallest accepted value."
    minimum::Float64
    "Largest accepted value."
    maximum::Float64
    "Native stepping interval."
    step::Float64
    "Displayed engineering unit."
    unit::String
    "Whether browser validation requires a value."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function UnitNumberInput(name; value::Real=0, minimum::Real=-Inf,
        maximum::Real=Inf, step::Real=1, unit::AbstractString="",
        id::AbstractString=string(name), required::Bool=false, disabled::Bool=false)
    minimum <= value <= maximum || throw(ArgumentError("numeric value is outside bounds"))
    step > 0 || throw(ArgumentError("numeric step must be positive"))
    initial = Float64(value)
    return UnitNumberInput(form_name(name, "number-input name"), string(id),
        Observable(initial), initial, Float64(minimum), Float64(maximum),
        Float64(step), string(unit), required, disabled)
end

finite_attribute(value) = isfinite(value) ? value : nothing

function Bonito.jsrender(session::Session, control::UnitNumberInput)
    input = DOM.input(; native_attributes(control; input_type="number")...,
        value=control.value, min=finite_attribute(control.minimum),
        max=finite_attribute(control.maximum), step=control.step,
        class="lc-control-input lc-form-control lc-unit-number-input",
        oninput=js"event => { const value = event.currentTarget.valueAsNumber; if (Number.isFinite(value)) $(control.value).notify(value); }")
    node = DOM.div(input, DOM.span(control.unit; class="lc-unit-number-unit");
        class="lc-unit-number")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

"""Create paired lower and upper numeric entries with shared bounds and unit."""
struct RangeInput <: AbstractFormControl
    "Form payload key prefix."
    name::Symbol
    "Stable group identifier."
    id::String
    "Current lower bound."
    lower::Observable{Float64}
    "Current upper bound."
    upper::Observable{Float64}
    "Lower bound restored by form reset."
    initial_lower::Float64
    "Upper bound restored by form reset."
    initial_upper::Float64
    "Smallest accepted value."
    minimum::Float64
    "Largest accepted value."
    maximum::Float64
    "Native stepping interval."
    step::Float64
    "Displayed engineering unit."
    unit::String
    "Whether browser validation requires both values."
    required::Bool
    "Whether interaction is disabled."
    disabled::Bool
end

function RangeInput(name; lower::Real, upper::Real, minimum::Real=-Inf,
        maximum::Real=Inf, step::Real=1, unit::AbstractString="",
        id::AbstractString=string(name), required::Bool=false, disabled::Bool=false)
    minimum <= lower <= upper <= maximum || throw(ArgumentError(
        "range values must satisfy minimum ≤ lower ≤ upper ≤ maximum"))
    step > 0 || throw(ArgumentError("range step must be positive"))
    return RangeInput(form_name(name, "range-input name"), string(id),
        Observable(Float64(lower)), Observable(Float64(upper)),
        Float64(lower), Float64(upper), Float64(minimum), Float64(maximum),
        Float64(step), string(unit), required, disabled)
end

control_id(control::RangeInput) = "$(control.id)-lower"
reset_control!(control::RangeInput) = begin
    control.lower[] = control.initial_lower
    control.upper[] = control.initial_upper
    control
end

function Bonito.jsrender(session::Session, control::RangeInput)
    lower = DOM.input(; id="$(control.id)-lower", name="$(control.name)_lower",
        type="number", value=control.lower, min=finite_attribute(control.minimum),
        max=finite_attribute(control.maximum), step=control.step,
        required=control.required, disabled=control.disabled,
        class="lc-control-input lc-form-control",
        oninput=js"event => { const value = event.currentTarget.valueAsNumber; if (Number.isFinite(value)) $(control.lower).notify(value); }")
    upper = DOM.input(; id="$(control.id)-upper", name="$(control.name)_upper",
        type="number", value=control.upper, min=finite_attribute(control.minimum),
        max=finite_attribute(control.maximum), step=control.step,
        required=control.required, disabled=control.disabled,
        class="lc-control-input lc-form-control",
        oninput=js"event => { const value = event.currentTarget.valueAsNumber; if (Number.isFinite(value)) $(control.upper).notify(value); }")
    node = DOM.div(
        DOM.label(DOM.span("Lower"), lower; var"for"="$(control.id)-lower"),
        DOM.span("to"; class="lc-range-separator"),
        DOM.label(DOM.span("Upper"), upper; var"for"="$(control.id)-upper"),
        DOM.span(control.unit; class="lc-unit-number-unit");
        id=control.id, class="lc-range-input", role="group")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, control))
end

"""
    Field(label, control; hint=nothing, error=nothing)

Compose a label, optional hint, and live validation message around one control.
"""
struct Field{C<:AbstractFormControl}
    "Visible field label."
    label::String
    "Owned form control."
    control::C
    "Optional explanatory text."
    hint::Union{Nothing,String}
    "Live validation message."
    error::Observable{Union{Nothing,String}}
    "Whether the label displays the required marker."
    required::Bool
end

function Field(label::AbstractString, control::C;
        hint::Union{Nothing,AbstractString}=nothing,
        error::Union{Nothing,AbstractString}=nothing,
        required::Bool=control_required(control)) where {C<:AbstractFormControl}
    return Field{C}(string(label), control,
        isnothing(hint) ? nothing : string(hint),
        Observable{Union{Nothing,String}}(isnothing(error) ? nothing : string(error)),
        required)
end

field_error!(field::Field, error::Nothing) = (field.error[] = nothing; field)
field_error!(field::Field, error::AbstractString) = (field.error[] = string(error); field)
reset_control!(field::Field) = (reset_control!(field.control); field_error!(field, nothing); field)

function Bonito.jsrender(session::Session, field::Field)
    message_class = map(session, field.error) do message
        isnothing(message) ? "lc-field-message" : "lc-field-message is-visible"
    end
    message = map(session, field.error) do value
        something(value, "")
    end
    node = DOM.div(
        DOM.label(field.label,
            field.required ? DOM.span(" required"; class="lc-field-required") : nothing;
            var"for"=control_id(field.control), class="lc-field-label"),
        field.control,
        isnothing(field.hint) ? nothing : DOM.div(field.hint; class="lc-field-hint"),
        DOM.div(message; class=message_class, role="status");
        class="lc-field", var"data-field-name"=string(control_name(field.control)))
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, field))
end

"""
    Form(fields...; name=:form, validator=nothing, submit_label="Apply",
         reset_label="Reset", cancel_label=nothing)

Compose typed fields into a native form. Browser validation runs first; the
plain payload is transferred to Julia only on explicit submission. `validator`
may return `nothing` or a dictionary keyed by field name.
"""
struct Form{F,V}
    "Stable form name."
    name::Symbol
    "Ordered fields."
    fields::F
    "Optional payload validation callback."
    validator::V
    "Latest successfully submitted payload."
    submitted::Observable{Union{Nothing,Dict{String,Any}}}
    "Form lifecycle state."
    status::Observable{Symbol}
    "Current field errors."
    errors::Observable{Dict{String,String}}
    "Submit button label."
    submit_label::String
    "Reset button label."
    reset_label::String
    "Optional cancel button label."
    cancel_label::Union{Nothing,String}
    "Number of cancel requests."
    cancelled::Observable{Int}
end

function Form(fields::Field...; name=:form, validator=nothing,
        submit_label::AbstractString="Apply", reset_label::AbstractString="Reset",
        cancel_label::Union{Nothing,AbstractString}=nothing)
    names = [control_name(field.control) for field in fields]
    allunique(names) || throw(ArgumentError("form field names must be unique"))
    isnothing(validator) || applicable(validator, Dict{String,Any}()) || throw(
        ArgumentError("form validator must accept a dictionary payload"))
    return Form(form_name(name, "form name"), fields, validator,
        Observable{Union{Nothing,Dict{String,Any}}}(nothing), Observable(:pristine),
        Observable(Dict{String,String}()), string(submit_label), string(reset_label),
        isnothing(cancel_label) ? nothing : string(cancel_label), Observable(0))
end

function normalize_form_errors(errors)
    isnothing(errors) && return Dict{String,String}()
    errors isa AbstractDict || throw(ArgumentError(
        "form validator must return nothing or a dictionary"))
    return Dict(string(key) => string(value) for (key, value) in pairs(errors))
end

function assign_form_errors!(form::Form, errors)
    resolved = normalize_form_errors(errors)
    form.errors[] = resolved
    for field in form.fields
        field_error!(field, get(resolved, string(control_name(field.control)), nothing))
    end
    return resolved
end

function accept_form_payload!(form::Form, payload::Dict{String,Any})
    errors = try
        isnothing(form.validator) ? Dict{String,String}() :
            normalize_form_errors(form.validator(payload))
    catch error
        form.status[] = :error
        assign_form_errors!(form, Dict("_form" => sprint(showerror, error)))
        return false
    end
    assign_form_errors!(form, errors)
    if isempty(errors)
        form.submitted[] = payload
        form.status[] = :submitted
        return true
    end
    form.status[] = :invalid
    return false
end

"""Return the latest accepted plain payload, or `nothing` before submission."""
submit_payload(form::Form) = form.submitted[]

"""Restore every field to its declared initial value and clear validation state."""
function reset_form!(form::Form)
    foreach(reset_control!, form.fields)
    assign_form_errors!(form, nothing)
    form.submitted[] = nothing
    form.status[] = :pristine
    return form
end

function Bonito.jsrender(session::Session, form::Form)
    submission_wire = Observable("{}")
    dirty_wire = Observable(false)
    reset_wire = Observable(0)
    cancel_wire = Observable(0)
    on(session, submission_wire) do encoded
        accept_form_payload!(form, JSON3.read(encoded, Dict{String,Any}))
        return nothing
    end
    on(session, dirty_wire) do dirty
        dirty && form.status[] in (:pristine, :submitted) && (form.status[] = :dirty)
        return nothing
    end
    on(session, reset_wire) do _
        reset_form!(form)
        return nothing
    end
    on(session, cancel_wire) do _
        form.cancelled[] += 1
        return nothing
    end
    cancel_button = isnothing(form.cancel_label) ? nothing : DOM.button(
        form.cancel_label; type="button", class="lc-button lc-button-secondary",
        onclick=js"event => $(cancel_wire).notify($(cancel_wire).value + 1)")
    reset_button = DOM.button(form.reset_label; type="reset",
        class="lc-button lc-button-secondary")
    submit_button = DOM.button(form.submit_label; type="submit",
        class="lc-button lc-button-primary")
    form_error = map(session, form.errors) do errors
        get(errors, "_form", "")
    end
    form_error_class = map(session, form_error) do message
        isempty(message) ? "lc-form-error" : "lc-form-error is-visible"
    end
    node = DOM.form(
        DOM.div(form.fields...; class="lc-form-fields"),
        DOM.div(form_error; class=form_error_class, role="alert"),
        DOM.div(cancel_button, reset_button, submit_button; class="lc-form-actions");
        id=string(form.name), class="lc-form", novalidate=true,
        oninput=js"event => $(dirty_wire).notify(true)",
        onreset=js"event => $(reset_wire).notify($(reset_wire).value + 1)",
        onsubmit=js"""
        event => {
            event.preventDefault();
            const form = event.currentTarget;
            if (!form.checkValidity()) {
                form.reportValidity();
                return;
            }
            const payload = {};
            for (const [key, value] of new FormData(form).entries()) {
                if (Object.prototype.hasOwnProperty.call(payload, key)) {
                    payload[key] = Array.isArray(payload[key]) ? [...payload[key], value] : [payload[key], value];
                } else {
                    payload[key] = value;
                }
            }
            $(submission_wire).notify(JSON.stringify(payload));
        }
        """)
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, form))
end

function control_inspection(control; parameters=[])
    bindings = hasfield(typeof(control), :value) ? [
        ComponentXRay.BindingInspection(:value, control.value; notes="session-local field state")
    ] : ComponentXRay.BindingInspection[]
    return ComponentXRay.ComponentInspection(control;
        name=string(nameof(typeof(control))), source=toolkit_source(@__FILE__, @__LINE__),
        parameters=[ComponentXRay.PropertyInspection(:name, control.name),
            ComponentXRay.PropertyInspection(:required, control.required),
            ComponentXRay.PropertyInspection(:disabled, control.disabled), parameters...],
        bindings, css_scopes=[".lc-form-control"])
end

ComponentXRay.inspection(control::TextInput) = control_inspection(control;
    parameters=[ComponentXRay.PropertyInspection(:placeholder, control.placeholder)])
ComponentXRay.inspection(control::TextAreaInput) = control_inspection(control;
    parameters=[ComponentXRay.PropertyInspection(:rows, control.rows)])

function ComponentXRay.inspection(control::SecretInput)
    return ComponentXRay.ComponentInspection(control; name="SecretInput",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:name, control.name),
            ComponentXRay.PropertyInspection(:autocomplete, control.autocomplete),
            ComponentXRay.PropertyInspection(:required, control.required),
            ComponentXRay.PropertyInspection(:disabled, control.disabled),
            ComponentXRay.PropertyInspection(:revealable, control.revealable)],
        actions=[ComponentXRay.ActionInspection(:click, "toggle visibility", nothing)],
        css_scopes=[".lc-secret-control", ".lc-secret-input", ".lc-secret-reveal"],
        notes=["Secret value is intentionally absent from Julia bindings and diagnostics."])
end

ComponentXRay.inspection(control::RadioGroup) = control_inspection(control;
    parameters=[ComponentXRay.PropertyInspection(:options, length(control.options)),
        ComponentXRay.PropertyInspection(:orientation, control.orientation)])
ComponentXRay.inspection(control::SegmentedControl) = control_inspection(control;
    parameters=[ComponentXRay.PropertyInspection(:options, length(control.options))])
ComponentXRay.inspection(control::ComboBox) = control_inspection(control;
    parameters=[ComponentXRay.PropertyInspection(:options, length(control.options))])
ComponentXRay.inspection(control::MultiSelect) = control_inspection(control;
    parameters=[ComponentXRay.PropertyInspection(:options, length(control.options)),
        ComponentXRay.PropertyInspection(:size, control.size)])
ComponentXRay.inspection(control::UnitNumberInput) = control_inspection(control;
    parameters=[ComponentXRay.PropertyInspection(:bounds, (control.minimum, control.maximum)),
        ComponentXRay.PropertyInspection(:step, control.step),
        ComponentXRay.PropertyInspection(:unit, control.unit)])

function ComponentXRay.inspection(control::RangeInput)
    return ComponentXRay.ComponentInspection(control; name="RangeInput",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:name, control.name),
            ComponentXRay.PropertyInspection(:bounds, (control.minimum, control.maximum)),
            ComponentXRay.PropertyInspection(:step, control.step),
            ComponentXRay.PropertyInspection(:unit, control.unit),
            ComponentXRay.PropertyInspection(:required, control.required),
            ComponentXRay.PropertyInspection(:disabled, control.disabled)],
        bindings=[ComponentXRay.BindingInspection(:lower, control.lower),
            ComponentXRay.BindingInspection(:upper, control.upper)],
        css_scopes=[".lc-range-input", ".lc-form-control"])
end

function ComponentXRay.inspection(field::Field)
    return ComponentXRay.ComponentInspection(field; name="Field",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:label, field.label),
            ComponentXRay.PropertyInspection(:hint, something(field.hint, "")),
            ComponentXRay.PropertyInspection(:required, field.required),
            ComponentXRay.PropertyInspection(:control, ComponentXRay.short_type(field.control))],
        bindings=[ComponentXRay.BindingInspection(:error, field.error)],
        css_scopes=[".lc-field", ".lc-field-label", ".lc-field-hint", ".lc-field-message"])
end

function ComponentXRay.inspection(form::Form)
    return ComponentXRay.ComponentInspection(form; name="Form",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:name, form.name),
            ComponentXRay.PropertyInspection(:fields, length(form.fields)),
            ComponentXRay.PropertyInspection(:validator, isnothing(form.validator) ? :none : :custom)],
        bindings=[ComponentXRay.BindingInspection(:status, form.status),
            ComponentXRay.BindingInspection(:errors, form.errors),
            ComponentXRay.BindingInspection(:cancelled, form.cancelled)],
        actions=[ComponentXRay.ActionInspection(:submit, form.submit_label, form.validator),
            ComponentXRay.ActionInspection(:reset, form.reset_label, reset_form!)],
        css_scopes=[".lc-form", ".lc-form-fields", ".lc-form-actions"])
end
