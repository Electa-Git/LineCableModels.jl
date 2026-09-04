module Toolkit

using Bonito
using JSON3
using Observables
using UUIDs
using ..ComponentXRay

const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TOOLKIT_STYLE_PATHS = (
    joinpath(TOOLKIT_ROOT, "assets", "forms.css"),
    joinpath(TOOLKIT_ROOT, "assets", "overlays.css"),
    joinpath(TOOLKIT_ROOT, "assets", "data-views.css"),
)
foreach(include_dependency, TOOLKIT_STYLE_PATHS)
const TOOLKIT_STYLES = join((read(path, String) for path in TOOLKIT_STYLE_PATHS), '\n')

toolkit_source(file, line) = ComponentXRay.source_reference(@__MODULE__, file, line)

include("Forms.jl")
include("Overlays.jl")
include("DataViews.jl")

export ComboBox,
    ConfirmDialog,
    DataTable,
    Dialog,
    DialogAction,
    DialogEvent,
    Disclosure,
    Field,
    Form,
    FormDialog,
    InlineNotice,
    MessageDialog,
    MultiSelect,
    PropertyGrid,
    PropertyItem,
    RadioGroup,
    RangeInput,
    SecretInput,
    SegmentedControl,
    TableColumn,
    TextAreaInput,
    TextInput,
    ToastCenter,
    ToastEntry,
    UnitNumberInput,
    ViewportFrame,
    clear_toasts!,
    close_dialog!,
    dismiss_notice!,
    dismiss_toast!,
    field_error!,
    open_dialog!,
    push_toast!,
    reset_form!,
    set_viewport_state!,
    submit_payload

end
