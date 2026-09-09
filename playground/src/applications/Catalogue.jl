"""
    ApplicationCatalogue

Declare application metadata without importing Bonito, a broker or numerical
packages. Both static publishing and the runtime registry consume these rows.
"""
module ApplicationCatalogue

scientific_requirements() = (
    (role="parameters", profiles=("line-parameters",), required=true),
    (role="power-flow", profiles=("power-flow",), required=true),
    (role="terminal", profiles=("julia-terminal",), required=false),
)

"""
    entries()

Return the immutable catalogue. A UI entry identifies an approved local driver;
it is never supplied by the browser. Profiles are declared here, not prepared.
"""
entries() = (
    (id="ichqp-showcase", title="Frequency-dependent cable parameters", kind=:presentation,
     entrypoint="/presentations/showcase.html", version=v"1.0.0", visibility=:public,
     description="ICHQP2026: cable models, line parameters and an OHL/UGC application case.",
     ui=:showcase, entry_surface=:published, requirements=scientific_requirements()),
    (id="cable-study", title="CableStudy", kind=:workbench,
     entrypoint="/workbenches/cable-study", version=v"1.0.0", visibility=:public,
     description="Cable construction, frequency-dependent line parameters and OHL/UGC impedance sensitivity.",
     ui=:cable_study, entry_surface=:ui, requirements=scientific_requirements()),
    (id="template-workbench", title="Workbench foundation", kind=:workbench,
     entrypoint="/workbenches/template", version=v"1.0.0", visibility=:developer,
     description="The reusable engineering shell with local demonstration state and X-ray.",
     ui=:template, entry_surface=:ui, requirements=()),
    (id="toolkit-gallery", title="Interactive toolkit gallery", kind=:workbench,
     entrypoint="/widgets/index.html", version=v"1.0.0", visibility=:developer,
     description="The original registered widget factories inside one isolated UI host.",
     ui=:gallery, entry_surface=:published, requirements=()),
    (id="starter-deck", title="Starter deck", kind=:presentation,
     entrypoint="/presentations/starter.html", version=v"1.0.0", visibility=:developer,
     description="Copyable Markdown, native incremental lists and a persistent live view.",
     ui=:starter, entry_surface=:published, requirements=()),
    (id="hostile-deck", title="Hostile specimen", kind=:presentation,
     entrypoint="/presentations/specimen.html", version=v"1.0.0", visibility=:developer,
     description="Layout, interaction, math annotations and live-view acceptance cases.",
     ui=:specimen, entry_surface=:published, requirements=()),
)

"""
    public_entries()

Return browser-safe catalogue rows; omit the trusted UI driver implementation.
"""
public_entries() = map(entries()) do entry
    (; entry.id, entry.title, kind=String(entry.kind), entry.entrypoint,
        version=string(entry.version), visibility=String(entry.visibility),
        entry.description, entry_surface=String(entry.entry_surface), entry.requirements)
end

end
