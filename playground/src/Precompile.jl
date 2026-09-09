using PrecompileTools: @compile_workload

# Compile representative presentation code while building the package image,
# never when registering applications or accepting a browser request. These
# sessions have neither a network connection nor an asset server. In particular
# this does not prepare scientific profiles, connect to NATS, or create uploads.
function precompile_ui_workload()
    for factory in (form_toolkit_widget, control_panel_widget,
            () -> TemplateWorkbench.app(; xray=true))
        session = Bonito.Session(Bonito.NoConnection(); asset_server=Bonito.NoServer())
        try
            app = factory()
            dom = Bonito.rendered_dom(session, app, Bonito.HTTP.Request("GET", "/"))
            Bonito.jsrender(session, dom)
            session.init_error[] === nothing || error("UI precompile rendering failed")
        finally
            close(session)
        end
    end
    return nothing
end

@compile_workload begin
    precompile_ui_workload()
end
