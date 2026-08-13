# import DataFrames: DataFrame
# using DataFrames

# ---------------------------------------------------------------------------
# Stringification helpers (keeps uncertainties visible in Excel)
# ---------------------------------------------------------------------------
stringify(x) = string(x)  # fallback (rarely reached)
stringify(::Missing) = ""
stringify(x::Real) = @sprintf("%.12g", float(x))
stringify(x::Measurements.Measurement) =
	@sprintf("%.12g ± %.6g", Measurements.value(x), Measurements.uncertainty(x))

function df_to_strings(df::DataFrame)
	DataFrame((name => stringify.(df[!, name]) for name in names(df))...; copycols = false)
end

# helper to fetch the units dict from df.metadata
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
# XLSX sheet writer: reuses/renames Sheet1 for the first write to avoid blanks
# ---------------------------------------------------------------------------
function _write_sheet!(xf, sheetname::String, df::DataFrame; use_first_sheet::Bool)
	units = _get_units(df)
	df_str = df_to_strings(df)

	ws = nothing
	if use_first_sheet
		ws = try
			xf["Sheet1"]            # reuse default first sheet
		catch
			nothing
		end
		ws = ws === nothing ? XLSX.addsheet!(xf, sheetname) : ws
		# If rename! exists, great; if not, we still write so Sheet1 isn't blank.
		try
			XLSX.rename!(ws, sheetname)
		catch
		end
	else
		ws = XLSX.addsheet!(xf, sheetname)
	end

	# Start row for writing
	start_row = 1

	# Optional UNITS block (Column | Unit) from DataFrame metadata
	if units isa AbstractDict
		for name in names(df)
			ws[start_row, 1] = String(name)
			u = get(units, name, get(units, Symbol(name), ""))
			ws[start_row, 2] = String(u)
			start_row += 1
		end
		start_row += 1   # spacer line
	end

	# IMPORTANT: anchor_cell must be a CellRef, not a String
	XLSX.writetable!(
		ws,
		Tables.columntable(df_str);
		anchor_cell = XLSX.CellRef(start_row, 1),
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
	cable_system::Union{LineCableSystem, Nothing} = nothing,
)::Union{String, Nothing}

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

	# ---- Build the DataFrames once (uses LP.f internally) ------------------
	df_z, df_y = DataFrame(line_params)  # each is Matrix{DataFrame}

	# Shapes
	nzx, nzy = size(df_z)
	nyx, nyy = size(df_y)

	# Diagonal-only logic (modal parameters)
	Z_isdiag = isdiag_approx(line_params.Z[:, :, 1])
	Y_isdiag = isdiag_approx(line_params.Y[:, :, 1])

	if Z_isdiag
		@warn "Z appears modal/diagonal (isdiag_approx=true). Exporting ONLY diagonal elements Z[i,i]; off-diagonals are intentionally omitted."
	end
	if Y_isdiag
		@warn "Y appears modal/diagonal (isdiag_approx=true). Exporting ONLY diagonal elements Y[i,i]; off-diagonals are intentionally omitted."
	end

	# ---- Write XLSX --------------------------------------------------------
	try
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
	catch err
		# If anything explodes (e.g., filesystem perms), return nothing.
		# Let the caller decide whether to rethrow.
		@error "Failed to export XLSX: $(err)"
		return nothing
	end
end
