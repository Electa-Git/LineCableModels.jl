# import DataFrames: DataFrame
# using DataFrames

# ---------------------------------------------------------------------------
# String conversion preserves uncertainty text in Excel.
# ---------------------------------------------------------------------------
stringify(x) = string(x)
stringify(::Missing) = ""
stringify(x::Real) = @sprintf("%.12g", float(x))

function df_to_strings(df::DataFrame)
    DataFrame((name => stringify.(df[!, name]) for name in names(df))...; copycols = false)
end

# Read the units dictionary from DataFrame metadata.
_get_units(df::DataFrame) =
    try
        DataFrames.metadata(df, "units", style = :note)
    catch
        try
            DataFrames.metadata(df, "units")
        catch
            nothing
        end
    end

# ---------------------------------------------------------------------------
# Reuse or rename Sheet1 for the first write to avoid an empty sheet.
# ---------------------------------------------------------------------------
function _write_sheet!(xf, sheetname::String, df::DataFrame; use_first_sheet::Bool)
    units = _get_units(df)
    df_str = df_to_strings(df)

    ws = nothing
    if use_first_sheet
        ws = try
            xf["Sheet1"]
        catch
            nothing
        end
        ws = ws === nothing ? XLSX.addsheet!(xf, sheetname) : ws
        # Write to the existing sheet if this XLSX version cannot rename it.
        try
            XLSX.rename!(ws, sheetname)
        catch
        end
    else
        ws = XLSX.addsheet!(xf, sheetname)
    end

    # First output row.
    start_row = 1

    # Optional units block from DataFrame metadata.
    if units isa AbstractDict
        for name in names(df)
            ws[start_row, 1] = String(name)
            u = get(units, name, get(units, Symbol(name), ""))
            ws[start_row, 2] = String(u)
            start_row += 1
        end
        start_row += 1
    end

    # XLSX requires a CellRef rather than a String for anchor_cell.
    XLSX.writetable!(
        ws,
        Tables.columntable(df_str);
        anchor_cell = XLSX.CellRef(start_row, 1)
    )
    return nothing
end

# ---------------------------------------------------------------------------
# Main export
# ---------------------------------------------------------------------------
function export_data(
        ::Val{:xlsx},
        line_params::LineParameters;
        file_name::Union{String, Nothing} = nothing,
        cable_system::Union{LineCableSystem, Nothing} = nothing
)::String

    # ---- Resolve final file_name (exactly as requested) --------------------
    if isnothing(file_name)
        if isnothing(cable_system)
            file_name = joinpath(@__DIR__, "ZY_export.xlsx")
        else
            file_name = joinpath(@__DIR__, "$(cable_system.system_id)_ZY_export.xlsx")
        end
    else
        requested = isabspath(file_name) ? file_name : joinpath(@__DIR__, file_name)
        if isnothing(cable_system)
            file_name = requested
        else
            dir = dirname(requested)
            base = basename(requested)
            file_name = joinpath(dir, "$(cable_system.system_id)_$base")
        end
    end

    # Build the frequency-indexed tables once.
    df_z, df_y = DataFrame(line_params)

    # Matrix dimensions.
    nzx, nzy = size(df_z)
    nyx, nyy = size(df_y)

    # Detect modal parameters represented by diagonal matrices.
    is_diagonal(matrix) = isapprox(
        matrix, Diagonal(diag(matrix)); rtol = 1.0e-8, atol = 1.0e-8
    )
    Z_isdiag = is_diagonal(observe(line_params, Z, :, :, 1))
    Y_isdiag = is_diagonal(observe(line_params, Y, :, :, 1))

    if Z_isdiag
        @warn "Z is diagonal within the selected tolerance. Exporting Z[i,i] and omitting off-diagonal elements."
    end
    if Y_isdiag
        @warn "Y is diagonal within the selected tolerance. Exporting Y[i,i] and omitting off-diagonal elements."
    end

    # ---- Write XLSX --------------------------------------------------------
    first_sheet = true
    XLSX.openxlsx(file_name, mode = "w") do xf
        # Z sheets
        if Z_isdiag
            for i in 1:min(nzx, nzy)
                _write_sheet!(xf, "Z($i,$i)", df_z[i, i]; use_first_sheet = first_sheet)
                first_sheet = false
            end
        else
            for i in 1:nzx, j in 1:nzy

                _write_sheet!(xf, "Z($i,$j)", df_z[i, j]; use_first_sheet = first_sheet)
                first_sheet = false
            end
        end

        # Y sheets
        if Y_isdiag
            for i in 1:min(nyx, nyy)
                _write_sheet!(xf, "Y($i,$i)", df_y[i, i]; use_first_sheet = first_sheet)
                first_sheet = false
            end
        else
            for i in 1:nyx, j in 1:nyy

                _write_sheet!(xf, "Y($i,$j)", df_y[i, j]; use_first_sheet = first_sheet)
                first_sheet = false
            end
        end
    end

    return file_name
end
