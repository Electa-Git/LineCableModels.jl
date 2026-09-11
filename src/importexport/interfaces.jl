"""
$(TYPEDSIGNATURES)

Export LineCableModels data in the format selected by `backend`.

# Arguments

- `backend`: Format selector.
- `args`: Inputs required by the selected format.

# Keywords

- Format-specific output options.

# Returns

- The output path or value defined by the selected format.

# Methods

$(METHODLIST)
"""
function export_data(backend::Symbol, args...; kwargs...)
    return export_data(Val(backend), args...; kwargs...)
end

function export_data(
        ::Val{:xlsx},
        line_parameters::LineParameters;
        file_name::Union{String, Nothing} = nothing,
        cable_system::Union{LineCableSystem, Nothing} = nothing
)::String
    artifact = ReportBuilder.report(
        ReportBuilder.XLSXReportDefinition(; file_name, cable_system),
        line_parameters
    )
    return artifact.output
end

"""
$(TYPEDSIGNATURES)

Import data in the format selected by `backend`.

# Arguments

- `backend`: Format selector.
- `args`: Inputs required by the selected format.

# Keywords

- Format-specific input options.

# Returns

- Materialised objects defined by the selected format.

# Methods

$(METHODLIST)
"""
function import_data(backend::Symbol, args...; kwargs...)
    return import_data(Val(backend), args...; kwargs...)
end
