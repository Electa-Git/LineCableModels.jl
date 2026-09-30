"""
$(TYPEDSIGNATURES)

Export LineCableModels data in the format selected by `format`.

# Arguments

- `format`: Format selector.
- `args`: Inputs required by the selected format.

# Keywords

- Format-specific output options.

# Returns

- The output path or value defined by the selected format.

# Methods

$(METHODLIST)
"""
function export_data(format::Symbol, args...; kwargs...)
    return export_data(Val(format), args...; kwargs...)
end

function export_data(
        ::Val{:xlsx},
        line_parameters::Union{LineParameters,Grammar.ObservedResult,AbstractVector{<:Grammar.ObservedResult}};
        file_name::Union{String, Nothing} = nothing,
        cable_system::Union{LineCableSystem, Nothing} = nothing,
        overwrite::Bool = false
)
    artifact = ReportBuilder.report(
        ReportBuilder.XLSXReportDefinition(; file_name, cable_system, overwrite),
        line_parameters
    )
    return artifact.output
end

"""
$(TYPEDSIGNATURES)

Import data in the format selected by `format`.

# Arguments

- `format`: Format selector.
- `args`: Inputs required by the selected format.

# Keywords

- Format-specific input options.

# Returns

- Materialized objects defined by the selected format.

# Methods

$(METHODLIST)
"""
function import_data(format::Symbol, args...; kwargs...)
    return import_data(Val(format), args...; kwargs...)
end
