const GEOGRAPHIC_MAP_STYLES_PATH = normpath(
    joinpath(@__DIR__, "..", "assets", "geographic-map.css")
)
const OPENLAYERS_STYLES_PATH = normpath(
    joinpath(@__DIR__, "..", "assets", "vendor", "openlayers.css")
)
const GEOGRAPHIC_MAP_SCRIPT_PATH = normpath(
    joinpath(@__DIR__, "..", "assets", "vendor", "geographic-map.bundle.js")
)
include_dependency(GEOGRAPHIC_MAP_STYLES_PATH)
include_dependency(OPENLAYERS_STYLES_PATH)
const GEOGRAPHIC_MAP_STYLES = read(GEOGRAPHIC_MAP_STYLES_PATH, String)
const OPENLAYERS_STYLES = read(OPENLAYERS_STYLES_PATH, String)
const GEOGRAPHIC_MAP_SCRIPT = Bonito.Asset(GEOGRAPHIC_MAP_SCRIPT_PATH)

const KML_MAX_BYTES = 16 * 1024 * 1024
const KMZ_MAX_ENTRIES = 64
const KMZ_MAX_EXPANDED_BYTES = 32 * 1024 * 1024
const KMZ_MAX_COMPRESSION_RATIO = 250.0

"A saved, user-defined point attached to one imported geographic feature."
struct MapReference
    uid::UUID
    source_feature::String
    longitude::Float64
    latitude::Float64
    identifier::String
    tag::String
end

"A KML or KMZ layer supplied through the reusable transactional uploader."
struct KMLUploadSource{Field<:UploadField}
    upload::Field
end

"A persistent browser map with line-snapped, editable reference points."
struct GeographicMap{Source<:KMLUploadSource}
    source::Source
    tags::Vector{Pair{String,String}}
    basemap::Symbol
    references::Observable{Vector{MapReference}}
    reference_payload::Observable{String}
end

function normalized_map_tags(tags)
    normalized = Pair{String,String}[]
    for option in tags
        option isa Pair || throw(ArgumentError(
            "map tags must be value => label pairs"
        ))
        value = strip(string(first(option)))
        label = strip(string(last(option)))
        isempty(value) && throw(ArgumentError("map tag values cannot be empty"))
        isempty(label) && throw(ArgumentError("map tag labels cannot be empty"))
        push!(normalized, value => label)
    end
    isempty(normalized) && throw(ArgumentError("at least one map tag is required"))
    values = first.(normalized)
    allunique(values) || throw(ArgumentError("map tag values must be unique"))
    return normalized
end

reference_record(reference::MapReference) = (
    uid=string(reference.uid),
    source_feature=reference.source_feature,
    longitude=reference.longitude,
    latitude=reference.latitude,
    identifier=reference.identifier,
    tag=reference.tag,
)

reference_payload(references) = JSON3.write(reference_record.(references))

function references_from_payload(payload::AbstractString)
    document = JSON3.read(payload)
    document isa JSON3.Array || throw(ArgumentError(
        "map reference payload must be an array"
    ))
    references = MapReference[]
    identifiers = Set{String}()
    for record in document
        uid = UUID(String(record.uid))
        longitude = Float64(record.longitude)
        latitude = Float64(record.latitude)
        isfinite(longitude) && -180 <= longitude <= 180 || throw(ArgumentError(
            "map reference longitude must be between -180 and 180 degrees"
        ))
        isfinite(latitude) && -90 <= latitude <= 90 || throw(ArgumentError(
            "map reference latitude must be between -90 and 90 degrees"
        ))
        identifier = strip(String(record.identifier))
        isempty(identifier) && throw(ArgumentError("map reference ID cannot be empty"))
        canonical_identifier = lowercase(identifier)
        canonical_identifier in identifiers && throw(ArgumentError(
            "map reference IDs must be unique"
        ))
        push!(identifiers, canonical_identifier)
        push!(references, MapReference(
            uid,
            String(record.source_feature),
            longitude,
            latitude,
            identifier,
            String(record.tag)
        ))
    end
    return references
end

function GeographicMap(
        source::KMLUploadSource;
        tags,
        basemap::Symbol=:osm,
        references=MapReference[]
    )
    basemap in (:none, :osm) || throw(ArgumentError(
        "geographic-map basemap must be :none or :osm"
    ))
    normalized_tags = normalized_map_tags(tags)
    initial_references = MapReference[references...]
    allowed_tags = Set(first.(normalized_tags))
    all(reference -> reference.tag in allowed_tags, initial_references) ||
        throw(ArgumentError("initial map references contain an unknown tag"))

    reference_observable = Observable(initial_references)
    payload_observable = Observable(reference_payload(initial_references))
    synchronizing = Ref(false)
    on(payload_observable) do payload
        synchronizing[] && return nothing
        parsed = try
            references_from_payload(string(payload))
        catch error
            @warn "Rejected invalid browser map-reference payload" exception=(
                error,
                catch_backtrace()
            )
            return nothing
        end
        all(reference -> reference.tag in allowed_tags, parsed) || begin
            @warn "Rejected browser map-reference payload containing an unknown tag"
            return nothing
        end
        synchronizing[] = true
        try
            reference_observable[] = parsed
        finally
            synchronizing[] = false
        end
        return nothing
    end
    on(reference_observable) do current
        synchronizing[] && return nothing
        all(reference -> reference.tag in allowed_tags, current) || throw(ArgumentError(
            "map references contain an unknown tag"
        ))
        synchronizing[] = true
        try
            payload_observable[] = reference_payload(current)
        finally
            synchronizing[] = false
        end
        return nothing
    end
    return GeographicMap(
        source,
        normalized_tags,
        basemap,
        reference_observable,
        payload_observable
    )
end

function safe_archive_member(name::AbstractString)
    normalized = replace(string(name), '\\' => '/')
    startswith(normalized, '/') && return false
    parts = split(normalized, '/'; keepempty=false)
    return !isempty(parts) && all(part -> part != "..", parts)
end

function validate_kml_bytes(bytes::AbstractVector{UInt8})
    length(bytes) <= KML_MAX_BYTES || throw(ArgumentError(
        "KML document exceeds the $(KML_MAX_BYTES)-byte expanded limit"
    ))
    isvalid(String, bytes) || throw(ArgumentError("KML document is not valid UTF-8"))
    source = String(copy(bytes))
    occursin(r"(?is)<!DOCTYPE|<!ENTITY", source) && throw(ArgumentError(
        "KML documents containing DTD or entity declarations are not accepted"
    ))
    EzXML.parsexml(source; noerror=true, nowarning=true)
    occursin(r"(?is)<(?:[A-Za-z_][A-Za-z0-9_.-]*:)?kml(?:\s|>)", source) ||
        throw(ArgumentError("input does not contain a KML root element"))
    occursin(r"(?is)<(?:[A-Za-z_][A-Za-z0-9_.-]*:)?Placemark(?:\s|>)", source) ||
        throw(ArgumentError("KML document contains no placemarks"))
    has_line = occursin(
        r"(?is)<(?:[A-Za-z_][A-Za-z0-9_.-]*:)?LineString(?:\s|>)",
        source
    )
    has_point = occursin(
        r"(?is)<(?:[A-Za-z_][A-Za-z0-9_.-]*:)?Point(?:\s|>)",
        source
    )
    (has_line || has_point) || throw(ArgumentError(
        "KML document contains no supported Point or LineString placemarks"
    ))
    return nothing
end

function kmz_kml_entry(files)
    kml_files = filter(file -> endswith(lowercase(file.name), ".kml"), files)
    root_files = filter(file -> lowercase(replace(file.name, '\\' => '/')) == "doc.kml", kml_files)
    length(root_files) == 1 && return only(root_files)
    isempty(root_files) && length(kml_files) == 1 && return only(kml_files)
    isempty(kml_files) && throw(ArgumentError(
        "KMZ archive does not contain a KML document"
    ))
    throw(ArgumentError(
        "KMZ archive has multiple KML documents and no unambiguous root doc.kml"
    ))
end

function validate_kmz(path::AbstractString)
    reader = ZipFile.Reader(path)
    try
        length(reader.files) <= KMZ_MAX_ENTRIES || throw(ArgumentError(
            "KMZ archive contains more than $KMZ_MAX_ENTRIES entries"
        ))
        total_expanded = UInt64(0)
        for file in reader.files
            safe_archive_member(file.name) || throw(ArgumentError(
                "KMZ archive contains an unsafe member path"
            ))
            total_expanded += file.uncompressedsize
            total_expanded <= KMZ_MAX_EXPANDED_BYTES || throw(ArgumentError(
                "KMZ archive exceeds the $(KMZ_MAX_EXPANDED_BYTES)-byte expanded limit"
            ))
            if file.compressedsize > 0
                ratio = Float64(file.uncompressedsize) / Float64(file.compressedsize)
                ratio <= KMZ_MAX_COMPRESSION_RATIO || throw(ArgumentError(
                    "KMZ archive contains an excessively compressed entry"
                ))
            end
        end
        kml_file = kmz_kml_entry(reader.files)
        validate_kml_bytes(read(kml_file))
    finally
        close(reader)
    end
    return nothing
end

"Validate the bounded KML/KMZ subset consumed by `GeographicMap`."
function validate_kml_upload(path::AbstractString)
    zipped = open(path, "r") do io
        signature = read(io, min(4, filesize(path)))
        return length(signature) >= 2 && signature[1:2] == UInt8[0x50, 0x4b]
    end
    zipped && return validate_kmz(path)
    validate_kml_bytes(read(path))
    return nothing
end

function map_add_icon()
    return SVG.svg(
        SVG.path(; d="M12 5v14"),
        SVG.path(; d="M5 12h14");
        viewBox="0 0 24 24",
        width="16",
        height="16",
        fill="none",
        stroke="currentColor",
        var"stroke-width"="2",
        var"stroke-linecap"="round",
        var"aria-hidden"="true"
    )
end

function Bonito.jsrender(session::Session, component::GeographicMap)
    upload = component.source.upload
    library = DOM.script(src=GEOGRAPHIC_MAP_SCRIPT)
    configuration = JSON3.write((
        basemap=string(component.basemap),
        tags=[(; value=first(tag), label=last(tag)) for tag in component.tags],
    ))
    map_canvas = DOM.div(
        tabindex="0",
        var"aria-label"="Interactive geographic map",
        var"data-map-role"="canvas",
        class="lc-map-canvas"
    )
    editor = DOM.section(
        DOM.h3("Reference"; var"data-map-role"="editor-title"),
        DOM.label(
            DOM.span("ID"),
            DOM.input(
                type="text",
                autocomplete="off",
                var"data-map-role"="identifier"
            );
            class="lc-map-field"
        ),
        DOM.label(
            DOM.span("Tag"),
            DOM.select(
                (
                    DOM.option(last(tag); value=first(tag))
                    for tag in component.tags
                )...;
                var"data-map-role"="tag"
            );
            class="lc-map-field"
        ),
        DOM.p(""; var"data-map-role"="editor-error", class="lc-map-editor-error"),
        DOM.div(
            DOM.button(
                "Cancel";
                type="button",
                var"data-map-role"="cancel"
            ),
            DOM.button(
                "Delete";
                type="button",
                var"data-map-role"="delete"
            ),
            DOM.button(
                "Save";
                type="button",
                var"data-map-role"="save"
            );
            class="lc-map-editor-actions"
        );
        hidden=true,
        var"data-map-role"="editor",
        class="lc-map-editor"
    )
    root = DOM.div(
        DOM.style(OPENLAYERS_STYLES),
        DOM.style(GEOGRAPHIC_MAP_STYLES),
        library,
        upload,
        DOM.div(
            DOM.div(
                map_canvas,
                DOM.div(
                    DOM.button(
                        map_add_icon(),
                        DOM.span("Add reference");
                        type="button",
                        var"aria-pressed"="false",
                        var"data-map-role"="add-reference",
                        class="lc-map-place-reference"
                    ),
                    DOM.output(
                        "Choose a KML or KMZ file";
                        var"data-kind"="empty",
                        var"data-map-role"="status",
                        class="lc-map-status"
                    );
                    class="lc-map-command-bar"
                ),
                DOM.div(
                    "Load route geometry, enable Add reference, then click a line.";
                    var"data-map-role"="empty",
                    class="lc-map-empty"
                );
                class="lc-map-stage"
            ),
            DOM.aside(
                DOM.header(
                    DOM.h2("References"),
                    DOM.output("0"; var"data-map-role"="reference-count");
                    class="lc-map-inspector-header"
                ),
                editor,
                DOM.div(
                    var"data-map-role"="reference-list",
                    class="lc-map-reference-list"
                );
                class="lc-map-inspector"
            );
            class="lc-map-workspace"
        );
        class="lc-map-component"
    )

    Bonito.onload(session, root, js"""
        (element) => {
            const library = $(library);
            const referencePayload = $(component.reference_payload);
            const configuration = JSON.parse($configuration);
            let mounted = false;

            const start = () => {
                if (mounted || !globalThis.LCMGeographicMap) return;
                mounted = true;
                try {
                    element._lcmGeographicMap = globalThis.LCMGeographicMap.mountGeographicMap(
                        element,
                        configuration,
                        referencePayload,
                    );
                } catch (error) {
                    const status = element.querySelector('[data-map-role="status"]');
                    status.dataset.kind = "error";
                    status.textContent = error instanceof Error
                        ? error.message
                        : "Unable to initialize geographic map";
                }
            };

            if (globalThis.LCMGeographicMap) start();
            else {
                library.addEventListener("load", start, {once: true});
                library.addEventListener("error", () => {
                    const status = element.querySelector('[data-map-role="status"]');
                    status.dataset.kind = "error";
                    status.textContent = "Unable to load the geographic map library";
                }, {once: true});
            }

            window.addEventListener("pagehide", () => {
                element._lcmGeographicMap?.destroy();
            }, {once: true});
        }
    """)
    return Bonito.jsrender(session, root)
end
