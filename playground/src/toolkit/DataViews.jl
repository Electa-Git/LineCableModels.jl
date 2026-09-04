"""Create a native expandable content section with observable open state."""
struct Disclosure{C}
    "Visible summary label."
    label::String
    "Renderable persistent content."
    content::C
    "Whether the section is expanded."
    open::Observable{Bool}
    "Whether interaction is disabled."
    disabled::Bool
end

Disclosure(label::AbstractString, content; open::Bool=false, disabled::Bool=false) =
    Disclosure(string(label), content, Observable(open), disabled)

function Bonito.jsrender(session::Session, disclosure::Disclosure)
    changed = Observable(disclosure.open[])
    on(session, changed) do value
        disclosure.disabled || (disclosure.open[] = value)
        return nothing
    end
    class = disclosure.disabled ? "lc-disclosure is-disabled" : "lc-disclosure"
    node = DOM.details(
        DOM.summary(disclosure.label),
        DOM.div(disclosure.content; class="lc-disclosure-body");
        class, open=disclosure.open,
        ontoggle=js"event => $(changed).notify(event.currentTarget.open)")
    synchronize = DOM.script(js"""
    (() => {
        const element = $(node);
        const state = $(disclosure.open);
        state.on(value => { if (element.open !== value) element.open = value; });
        if ($(disclosure.disabled)) {
            element.querySelector('summary').addEventListener('click', event => event.preventDefault());
        }
    })();
    """)
    return Bonito.jsrender(session, ComponentXRay.instrument(session,
        DOM.div(node, synchronize; style="display: contents;"), disclosure))
end

"""Describe one labelled value in a [`PropertyGrid`](@ref)."""
struct PropertyItem{V}
    "Visible property label."
    label::String
    "Renderable or observable property value."
    value::V
    "Optional engineering unit."
    unit::Union{Nothing,String}
    "Optional explanatory text."
    description::Union{Nothing,String}
    "Semantic value state."
    state::Symbol
end

function PropertyItem(label::AbstractString, value; unit=nothing, description=nothing,
        state::Symbol=:normal)
    state in (:normal, :success, :warning, :danger, :muted) || throw(ArgumentError(
        "property state must be :normal, :success, :warning, :danger, or :muted"))
    return PropertyItem(string(label), value,
        isnothing(unit) ? nothing : string(unit),
        isnothing(description) ? nothing : string(description), state)
end

"""Lay out labelled engineering values as a compact description list."""
struct PropertyGrid{I<:Tuple}
    "Ordered property items."
    items::I
    "Accessible grid label."
    label::String
    "Number of responsive columns."
    columns::Int
end

function PropertyGrid(items::PropertyItem...; label::AbstractString="Properties",
        columns::Integer=1)
    1 <= columns <= 4 || throw(ArgumentError("property-grid columns must be 1 through 4"))
    return PropertyGrid(items, string(label), Int(columns))
end

function Bonito.jsrender(session::Session, item::PropertyItem)
    node = DOM.div(
        DOM.dt(item.label),
        DOM.dd(DOM.span(item.value; class="lc-property-value"),
            isnothing(item.unit) ? nothing : DOM.span(item.unit; class="lc-property-unit")),
        isnothing(item.description) ? nothing :
            DOM.div(item.description; class="lc-property-description");
        class="lc-property-item is-$(item.state)")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, item))
end

function Bonito.jsrender(session::Session, grid::PropertyGrid)
    node = DOM.dl(grid.items...; class="lc-property-grid columns-$(grid.columns)",
        var"aria-label"=grid.label)
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, grid))
end

"""Describe one column in a [`DataTable`](@ref)."""
struct TableColumn{F}
    "Row lookup key."
    key::Symbol
    "Visible column label."
    label::String
    "Cell formatter."
    format::F
    "Horizontal alignment."
    align::Symbol
    "Whether users may sort the column."
    sortable::Bool
end

function TableColumn(key, label::AbstractString; format=identity,
        align::Symbol=:left, sortable::Bool=true)
    align in (:left, :center, :right) || throw(ArgumentError(
        "table-column alignment must be :left, :center, or :right"))
    return TableColumn(form_name(key, "table-column key"), string(label),
        format, align, sortable)
end

function table_value(row, key::Symbol)
    if row isa AbstractDict
        haskey(row, key) && return row[key]
        haskey(row, string(key)) && return row[string(key)]
        return missing
    end
    hasproperty(row, key) && return getproperty(row, key)
    return missing
end

function formatted_value(column::TableColumn, row)
    value = table_value(row, column.key)
    return string(column.format(value))
end

"""
    DataTable(columns...; rows=[], row_key=nothing, label="Data table")

Create a reactive compact table with client-triggered filtering, stable row
selection, and sortable declared columns. Rows may be named tuples, structs, or
dictionaries. The component never mutates source rows.
"""
struct DataTable{C,K}
    "Declared table columns."
    columns::C
    "Source rows."
    rows::Observable{Vector{Any}}
    "Function returning stable row identity."
    row_key::K
    "Selected row identity."
    selected::Observable{Union{Nothing,String}}
    "Current case-insensitive filter."
    filter::Observable{String}
    "Current sort column."
    sort_by::Observable{Union{Nothing,Symbol}}
    "Current sort direction."
    sort_direction::Observable{Symbol}
    "Accessible table label."
    label::String
    "Whether the filter control is rendered."
    filterable::Bool
end

function DataTable(columns::TableColumn...; rows=Any[], row_key=nothing,
        label::AbstractString="Data table", filterable::Bool=true)
    isempty(columns) && throw(ArgumentError("data table requires a column"))
    allunique(column.key for column in columns) || throw(ArgumentError(
        "data-table column keys must be unique"))
    keyer = isnothing(row_key) ? ((row, index) -> string(index)) : row_key
    return DataTable(columns, Observable(Any[rows...]), keyer,
        Observable{Union{Nothing,String}}(nothing), Observable(""),
        Observable{Union{Nothing,Symbol}}(nothing), Observable(:ascending),
        string(label), filterable)
end

function table_rows(table::DataTable)
    rows = collect(enumerate(table.rows[]))
    query = lowercase(strip(table.filter[]))
    if !isempty(query)
        rows = filter(rows) do (_, row)
            any(column -> occursin(query, lowercase(formatted_value(column, row))), table.columns)
        end
    end
    key = table.sort_by[]
    if !isnothing(key)
        sort!(rows; by=pair -> begin
            value = table_value(last(pair), key)
            ismissing(value) ? "" : value
        end, rev=table.sort_direction[] == :descending)
    end
    return rows
end

function table_payload(table::DataTable)
    rows = [Dict{String,Any}(
        "key" => string(table.row_key(row, index)),
        "cells" => [formatted_value(column, row) for column in table.columns],
        "align" => [string(column.align) for column in table.columns],
    ) for (index, row) in table_rows(table)]
    allunique(row["key"] for row in rows) || throw(ArgumentError(
        "data-table row identities must be unique"))
    return String(JSON3.write(rows))
end

function Bonito.jsrender(session::Session, table::DataTable)
    rows_payload = Observable(table_payload(table))
    sort_wire = Observable("")
    refresh = () -> begin
        rows_payload[] = table_payload(table)
    end
    for observable in (table.rows, table.filter, table.sort_by, table.sort_direction)
        on(session, observable) do _
            refresh()
            return nothing
        end
    end
    on(session, sort_wire) do requested
        key = Symbol(requested)
        if table.sort_by[] == key
            table.sort_direction[] = table.sort_direction[] == :ascending ?
                :descending : :ascending
        else
            table.sort_by[] = key
            table.sort_direction[] = :ascending
        end
        return nothing
    end
    filter_control = table.filterable ? DOM.input(; type="search",
        class="lc-control-input lc-form-control lc-data-filter",
        placeholder="Filter rows", var"aria-label"="Filter table rows",
        value=table.filter,
        oninput=js"event => $(table.filter).notify(event.currentTarget.value)") : nothing
    headers = map(table.columns) do column
        if column.sortable
            DOM.th(DOM.button(column.label; type="button", class="lc-data-sort",
                onclick=js"""
                event => {
                    $(sort_wire).notify($(string(column.key)));
                }
                """); scope="col", class="align-$(column.align)")
        else
            DOM.th(column.label; scope="col", class="align-$(column.align)")
        end
    end
    body = DOM.tbody()
    node = DOM.div(
        table.filterable ? DOM.div(filter_control; class="lc-data-tools") : nothing,
        DOM.div(DOM.table(DOM.thead(DOM.tr(headers...)), body;
            var"aria-label"=table.label); class="lc-data-table-scroll");
        class="lc-data-table")
    behavior = DOM.script(js"""
    (() => {
        const body = $(body);
        const payload = $(rows_payload);
        const selected = $(table.selected);

        const synchronizeSelection = key => {
            for (const row of body.rows) {
                const active = row.dataset.rowKey === key;
                row.classList.toggle('is-selected', active);
                row.setAttribute('aria-selected', String(active));
            }
        };

        const renderRows = encoded => {
            const fragment = document.createDocumentFragment();
            for (const record of JSON.parse(encoded)) {
                const row = document.createElement('tr');
                row.className = 'lc-data-row';
                row.dataset.rowKey = record.key;
                row.tabIndex = 0;
                row.setAttribute('aria-selected', 'false');
                record.cells.forEach((value, index) => {
                    const cell = document.createElement('td');
                    cell.className = 'align-' + record.align[index];
                    cell.textContent = value;
                    row.appendChild(cell);
                });
                const choose = () => selected.notify(record.key);
                row.addEventListener('click', choose);
                row.addEventListener('keydown', event => {
                    if (event.key !== 'Enter' && event.key !== ' ') return;
                    event.preventDefault();
                    choose();
                });
                fragment.appendChild(row);
            }
            body.replaceChildren(fragment);
            synchronizeSelection(selected.value);
        };

        payload.on(renderRows);
        selected.on(synchronizeSelection);
        renderRows(payload.value);
    })();
    """)
    root = DOM.div(node, behavior; style="display: contents;")
    return Bonito.jsrender(session, ComponentXRay.instrument(session, root, table))
end

"""
    ViewportFrame(title, content; state=:ready, message="", tools=nothing, footer=nothing)

Wrap persistent content with a compact title bar and non-destructive loading,
empty, and error overlays. State changes never replace or remount `content`.
"""
struct ViewportFrame{C,T,F}
    "Visible viewport title."
    title::String
    "Persistent viewport content."
    content::C
    "Optional title-bar tools."
    tools::T
    "Optional footer content."
    footer::F
    "Current viewport state."
    state::Observable{Symbol}
    "Current overlay message."
    message::Observable{String}
end

function ViewportFrame(title::AbstractString, content; state::Symbol=:ready,
        message::AbstractString="", tools=nothing, footer=nothing)
    state in (:ready, :loading, :empty, :error) || throw(ArgumentError(
        "viewport state must be :ready, :loading, :empty, or :error"))
    return ViewportFrame(string(title), content, tools, footer,
        Observable(state), Observable(string(message)))
end

"""Change a viewport overlay state without replacing its mounted content."""
function set_viewport_state!(viewport::ViewportFrame, state::Symbol; message=nothing)
    state in (:ready, :loading, :empty, :error) || throw(ArgumentError(
        "viewport state must be :ready, :loading, :empty, or :error"))
    viewport.state[] = state
    isnothing(message) || (viewport.message[] = string(message))
    return viewport
end

function Bonito.jsrender(session::Session, viewport::ViewportFrame)
    class = map(session, viewport.state) do state
        "lc-viewport-frame is-$state"
    end
    label = map(session, viewport.state, viewport.message) do state, message
        !isempty(message) && return message
        state == :loading && return "Loading view"
        state == :empty && return "No data available"
        state == :error && return "View unavailable"
        return ""
    end
    node = DOM.section(
        DOM.header(DOM.strong(viewport.title), viewport.tools;
            class="lc-viewport-header"),
        DOM.div(DOM.div(viewport.content; class="lc-viewport-content"),
            DOM.div(DOM.span(; class="lc-viewport-spinner"), DOM.span(label);
                class="lc-viewport-state", role="status");
            class="lc-viewport-body"),
        isnothing(viewport.footer) ? nothing :
            DOM.footer(viewport.footer; class="lc-viewport-footer");
        class)
    return Bonito.jsrender(session, ComponentXRay.instrument(session, node, viewport))
end

function ComponentXRay.inspection(disclosure::Disclosure)
    return ComponentXRay.ComponentInspection(disclosure; name="Disclosure",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:label, disclosure.label),
            ComponentXRay.PropertyInspection(:disabled, disclosure.disabled)],
        bindings=[ComponentXRay.BindingInspection(:open, disclosure.open)],
        css_scopes=[".lc-disclosure", ".lc-disclosure-body"])
end

function ComponentXRay.inspection(item::PropertyItem)
    return ComponentXRay.ComponentInspection(item; name="PropertyItem",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:label, item.label),
            ComponentXRay.PropertyInspection(:unit, something(item.unit, "")),
            ComponentXRay.PropertyInspection(:state, item.state)],
        css_scopes=[".lc-property-item"])
end

function ComponentXRay.inspection(grid::PropertyGrid)
    return ComponentXRay.ComponentInspection(grid; name="PropertyGrid",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:items, length(grid.items)),
            ComponentXRay.PropertyInspection(:columns, grid.columns)],
        css_scopes=[".lc-property-grid"])
end

function ComponentXRay.inspection(table::DataTable)
    return ComponentXRay.ComponentInspection(table; name="DataTable",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:columns, length(table.columns)),
            ComponentXRay.PropertyInspection(:filterable, table.filterable)],
        bindings=[ComponentXRay.BindingInspection(:rows, table.rows),
            ComponentXRay.BindingInspection(:selected, table.selected),
            ComponentXRay.BindingInspection(:filter, table.filter),
            ComponentXRay.BindingInspection(:sort_by, table.sort_by)],
        css_scopes=[".lc-data-table", ".lc-data-table-scroll"])
end

function ComponentXRay.inspection(viewport::ViewportFrame)
    return ComponentXRay.ComponentInspection(viewport; name="ViewportFrame",
        source=toolkit_source(@__FILE__, @__LINE__), parameters=[
            ComponentXRay.PropertyInspection(:title, viewport.title),
            ComponentXRay.PropertyInspection(:tools, !isnothing(viewport.tools)),
            ComponentXRay.PropertyInspection(:footer, !isnothing(viewport.footer))],
        bindings=[ComponentXRay.BindingInspection(:state, viewport.state),
            ComponentXRay.BindingInspection(:message, viewport.message)],
        css_scopes=[".lc-viewport-frame", ".lc-viewport-header",
            ".lc-viewport-content", ".lc-viewport-state"])
end
