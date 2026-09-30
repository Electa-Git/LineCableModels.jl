# Live axis actions and native toolbar ownership. These operations retain no
# scientific source and subscribe only to the native resources they own.
const _TOOLBAR_STANDARD_CONTROLS = (:reset, :xlog, :ylog, :export_svg, :legend)

function _panel_axes(p, panel)
    panel===nothing && return p.axes
    panel=_panel_identity(panel)
    haskey(p.plot_state.panel_data, panel) ||
        throw(ArgumentError("panel $(repr(panel)) is absent from this figure"))
    return [p.plot_state.panel_data[panel].axis]
end

function _numeric_axis(axis, dimension)
    index=dimension===:x ? 1 : 2
    return !(axis.aspect[] isa DataAspect) &&
           getproperty(axis, Symbol(:dim, index, :_conversion))[]===nothing
end

function LineCableModels.axisscale!(p::LineCableModels.UIPlot, dimension::Symbol, scale; panel = nothing)
    dimension in (:x, :y) || throw(ArgumentError("axis dimension must be :x or :y"))
    axes=_panel_axes(p, panel)
    entries=filter(
        entry -> entry.axis in axes && _numeric_axis(entry.axis, dimension),
        p.plot_state.axis_bindings)
    isempty(entries) &&
        throw(ArgumentError("no eligible numeric $dimension axis in the selected panels"))
    _set_axis!(entries, dimension, scale)
    return p
end

function LineCableModels.resetview!(p::LineCableModels.UIPlot; panel = nothing, x::Bool = true, y::Bool = true)
    axes=_panel_axes(p, panel)
    for entry in p.plot_state.axis_bindings
        entry.axis in axes || continue
        entry.reset(; xauto = x, yauto = y)
    end
    return p
end

function _detach!(object)
    gc=GridLayoutBase.gridcontent(object)
    gc===nothing || GridLayoutBase.remove_from_gridlayout!(gc)
    return nothing
end

function _delete_subtree!(grid::GridLayout)
    for child in copy(grid.content)
        object=child.content
        object isa GridLayout ? _delete_subtree!(object) : delete!(object)
    end
    _detach!(grid)
    return nothing
end

function _belongs_to_slot(object, slot)
    current=object
    while current!==nothing
        current===slot && return true
        gc=GridLayoutBase.gridcontent(current)
        current=gc===nothing ? nothing : gc.parent
    end
    return false
end

function _repack_toolbar!(p)
    state=p.plot_state
    for (index, key) in enumerate(state.widget_order)
        state.shell.toolbar[1, index]=state.widget_bindings[key].slot
        colsize!(state.shell.toolbar, index, Auto(true))
    end
    GridLayoutBase.trim!(state.shell.toolbar)
    return nothing
end

function _expected_callback_error(error)
    error isa Union{ArgumentError, DomainError, BoundsError,
        DimensionMismatch, SystemError, Base.IOError}
end

function _widget!(builder, p, key; event = nothing, callback = nothing,
        success = nothing, standard = false)
    state=p.plot_state
    state.controls_enabled ||
        throw(ArgumentError("this figure was constructed with controls=false"))
    key isa Symbol || throw(ArgumentError("widget keys must be Symbols"))
    haskey(p.controls, key) && throw(ArgumentError("widget key $key is already registered"))
    standard || key ∉ _TOOLBAR_STANDARD_CONTROLS ||
        throw(ArgumentError("widget key $key is reserved"))
    (event===nothing)==(callback===nothing) ||
        throw(ArgumentError("event and callback must be supplied together"))
    slot=GridLayout(; tellwidth = true, tellheight = true)
    state.shell.toolbar[1, length(state.widget_order) + 1]=slot
    subscriptions=Any[]
    try
        control=builder(p, slot[1, 1])
        native_block=!(control isa GridLayout) && hasproperty(control,:blockscene)
        (native_block || control isa GridLayout) ||
            throw(ArgumentError("widget builders must return a native block or GridLayout"))
        _belongs_to_slot(control, slot) ||
            throw(ArgumentError("widget content must belong to its allocated slot"))
        if event!==nothing
            observable=event(control)
            owner=native_block ? control : nothing
            if owner===nothing
                owner=findfirst(p.figure.content) do block
                    _belongs_to_slot(block, slot) &&
                        any(propertynames(typeof(block))) do name
                            getproperty(block, name)===observable
                        end
                end
                owner=owner===nothing ? nothing : p.figure.content[owner]
            end
            scene=owner===nothing ? p.figure.scene : owner.blockscene
            push!(subscriptions, on(scene, observable) do value
                try
                    callback(p, value)
                    success===nothing || (p.status[]=string(success))
                catch error
                    _expected_callback_error(error) || rethrow()
                    p.status[]=sprint(showerror, error)
                end
                return nothing
            end)
        end
        p.controls[key]=control
        state.widget_bindings[key]=(; slot, subscriptions, standard)
        push!(state.widget_order, key)
        _repack_toolbar!(p)
        return control
    catch
        foreach(off, subscriptions)
        _delete_subtree!(slot)
        GridLayoutBase.trim!(state.shell.toolbar)
        rethrow()
    end
end

function LineCableModels.addwidget!(builder, p::LineCableModels.UIPlot, key::Symbol;
        event = nothing, callback = nothing, success = nothing)
    return _update_figure_layout!(p) do
        _widget!(builder, p, key; event, callback, success)
    end
end

function LineCableModels.removewidget!(p::LineCableModels.UIPlot, key::Symbol)
    key ∉ _TOOLBAR_STANDARD_CONTROLS ||
        throw(ArgumentError("standard widget $key cannot be removed"))
    haskey(p.plot_state.widget_bindings, key) ||
        throw(ArgumentError("unknown widget $key"))
    return _update_figure_layout!(p) do
        owned=pop!(p.plot_state.widget_bindings, key)
        foreach(off, owned.subscriptions)
        _delete_subtree!(owned.slot)
        delete!(p.controls, key)
        filter!(!=(key), p.plot_state.widget_order)
        _repack_toolbar!(p)
        return p
    end
end

function _controls!(p, xsetters, ysetters)
    p.plot_state.controls_enabled || return p
    if !isempty(p.axes)
        _widget!(
            (p, slot) -> Button(slot; label = _toolbar_icon(_TOOLBAR_REFRESH_ICON),
                width = _TOOLBAR_BUTTON_SIZE, height = _TOOLBAR_BUTTON_SIZE, buttoncolor = _TOOLBAR_BUTTON_BACKGROUND),
            p,
            :reset;
            event = button -> button.clicks, callback = (
                p, _) -> LineCableModels.resetview!(p),
            success = "Axis limits reset", standard = true)
    end
    if Base.get_extension(LineCableModels, :LineCableModelsCairoMakieExt)!==nothing
        _widget!(
            (p, slot) -> Button(slot; label = _toolbar_icon(_TOOLBAR_SAVE_ICON),
                width = _TOOLBAR_BUTTON_SIZE, height = _TOOLBAR_BUTTON_SIZE, buttoncolor = _TOOLBAR_BUTTON_BACKGROUND),
            p,
            :export_svg;
            event = button -> button.clicks, callback = (p, _) -> begin
                output=LineCableModels.export_svg(p)
                p.status[]="Saved SVG to $output"
            end, standard = true)
    end
    for (dimension, entries) in ((:x, xsetters), (:y, ysetters))
        isempty(entries) && continue
        changing=Ref(false)
        caption=Ref{Any}(nothing)
        key=Symbol(dimension, :log)
        toggle=_widget!(p, key; event = control -> control.active, standard = true,
            callback = (p, enabled) -> begin
                changing[] && return nothing
                changing[]=true
                try
                    eligible=filter(
                        entry -> getproperty(entry.axis, Symbol(dimension, :scale))[] in (identity, log10) ||
                                 getproperty(entry.axis, Symbol(dimension, :scale))[] isa Makie.ReversibleScale{_SignedLog10},
                        entries)
                    _set_axis!(eligible, dimension, enabled ? :log10 : :linear)
                    p.status[]=enabled ? "$dimension-axis logarithmic view" :
                               "$dimension-axis linear view"
                catch
                    p.controls[key].active[]=!enabled
                    rethrow()
                finally
                    changing[]=false
                end
            end) do p, slot
            group=GridLayout(slot; halign = :left, tellwidth = true)
            result=Toggle(group[1, 1]; active = false)
            caption[]=Label(group[1, 2], "log $dimension")
            colgap!(group, 4)
            result
        end
        for entry in entries
            subscription=on(
                toggle.blockscene, getproperty(entry.axis, Symbol(dimension, :scale));
                update = true) do _
                kinds=unique(getproperty(item.axis, Symbol(dimension, :scale))[]===log10 ?
                             "log" :
                             getproperty(item.axis, Symbol(dimension, :scale))[]===identity ?
                             "linear" :
                             getproperty(item.axis, Symbol(dimension, :scale))[] isa Makie.ReversibleScale{_SignedLog10} ?
                             "signed log" : "custom" for item in entries)
                mode=length(kinds)==1 ? only(kinds) : "mixed"
                caption[].text[]=mode=="linear" ? "log $dimension" : "$mode $dimension"
                if !changing[]
                    changing[]=true
                    try
                        active=all(kind -> kind in ("log", "signed log"), kinds)
                        toggle.active[]==active || (toggle.active[]=active)
                    finally
                        changing[]=false
                    end
                end
                nothing
            end
            push!(p.plot_state.widget_bindings[key].subscriptions, subscription)
        end
    end
    p.legend===nothing || (p.controls[:legend]=p.legend)
    return p
end
