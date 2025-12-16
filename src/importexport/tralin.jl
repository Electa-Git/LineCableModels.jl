const _TRALIN_COMP = ("CORE", "SHEATH", "ARMOUR")


function export_data(::Val{:tralin},
	cable_system::LineCableSystem,
	earth_props::EarthModel;
	freq = f₀,
	file_name::Union{String, Nothing} = nothing,
)::Union{String, Nothing}

	# -- helpers ---------------------------------------------------------------
	_freqs(x) = x isa AbstractVector ? collect(x) : [x]
	_fmt(x) = string(round(Float64(to_nominal(x)); digits = 6))
	_maybe(x) = (x === nothing) ? "" : _fmt(x)

	# Resolve output file name (prefix "tr_"; mirror  XML semantics)
	if isnothing(file_name)
		file_name = joinpath(@__DIR__, "tr_$(cable_system.system_id).f05")
	else
		dir = dirname(file_name)
		fname = basename(file_name)
		# Ensure filename has "tr_" prefix, but preserve user's name
		prefixed_fname = startswith(fname, "tr_") ? fname : "tr_" * fname
		# Rejoin with original path, handling relative vs absolute
		file_name =
			isabspath(file_name) ? joinpath(dir, prefixed_fname) :
			joinpath(@__DIR__, dir, prefixed_fname)
	end

	num_phases = length(cable_system.cables)
	freqs = map(f -> to_nominal(f), _freqs(freq))

	# -- build TRALIN lines ----------------------------------------------------
	lines = String[]

	push!(lines, "TRALIN")
	push!(lines, "TEXT,MODULE,LineCableModels run")
	push!(lines, "OPTIONS")
	push!(lines, "UNITS,METRIC")
	push!(lines, "RUN-IDENTIFICATION,$(cable_system.system_id)")
	push!(lines, "SEQUENCE,ON")
	push!(lines, "MULTILAYER,ON")
	push!(lines, "CONDUCTANCE,ON")
	push!(lines, "!KEEP_CIRCUIT_MODE")

	push!(lines, "PARAMETERS")
	push!(lines, "BASE-VALUES")
	push!(lines, "ACCURACY,1e-7")
	push!(lines, "BESSEL")
	push!(lines, "TERMS,300")
	for f in freqs
		push!(lines, "FREQUENCY,$(_fmt(f))")
	end
	push!(lines, "INTEGRATION,AUTO-ADJUST,9")
	push!(lines, "STEP,1e-6")
	push!(lines, "UPPER-LIMIT,5.")
	push!(lines, "SERIES-TERMS,300")

	nlayers = length(earth_props.layers)

	if nlayers == 2
		# [AIR, SOIL] => uniform semi-infinite earth
		soil  = earth_props.layers[end]
		rho   = _fmt(getfield(soil, :base_rho_g))
		mu_r  = hasfield(typeof(soil), :mu_r) ? _fmt(getfield(soil, :mu_r)) : "1"
		eps_r = hasfield(typeof(soil), :eps_r) ? _fmt(getfield(soil, :eps_r)) : "1"
		push!(lines, "SOIL-TYPE")
		push!(lines, "UNIFORM,$rho,$mu_r,$eps_r")
	else
		# [AIR, TOP, (CENTRAL...), BOTTOM] => HORIZONTAL
		push!(lines, "SOIL-TYPE")
		push!(lines, "HORIZONTAL")

		# AIR: no thickness -> explicit empty field `,,`
		push!(lines, "    LAYER,AIR,1e+18,,1,1")

		n_earth = nlayers - 1
		names =
			n_earth == 1 ? ["TOP"] :
			n_earth == 2 ? ["TOP", "BOTTOM"] :
			vcat("TOP", fill("CENTRAL", n_earth - 2), "BOTTOM")

		for (eidx, (lname, layer)) in enumerate(zip(names, earth_props.layers[2:end]))
			rho   = _fmt(getfield(layer, :base_rho_g))
			mu_r  = hasfield(typeof(layer), :base_mur_g) ? _fmt(getfield(layer, :base_mur_g)) : "1"
			eps_r = hasfield(typeof(layer), :base_epsr_g) ? _fmt(getfield(layer, :base_epsr_g)) : "1"

			if eidx == n_earth
				# BOTTOM: no thickness -> explicit empty field `,,`
				push!(lines, "    LAYER,$lname,$rho,,$mu_r,$eps_r")
			else
				# TOP/CENTRAL: include thickness if available; otherwise leave it empty to keep the slot
				thk =
					(
						hasfield(typeof(layer), :t) &&
						getfield(layer, :t) !== nothing
					) ?
					_fmt(getfield(layer, :t)) : ""
				push!(lines, "    LAYER,$lname,$rho,$thk,$mu_r,$eps_r")
			end
		end
	end

	push!(lines, "SYSTEM")

	for (pidx, cable) in enumerate(cable_system.cables)
		# Phase group position
		push!(lines, "GROUP,PH-$(pidx),$(_fmt(cable.horz)),$(_fmt(cable.vert))")

		comps_vec = cable.design_data.components  # assumed Vector in your corrected model
		ncomp = length(comps_vec)
		if ncomp > 3
			throw(
				ArgumentError(
					"TRALIN supports at most 3 concentric components (CORE/SHEATH/ARMOR); got $ncomp for cable index $pidx.",
				),
			)
		end
		# Outer radius for CABLE line
		outer_R = to_nominal(comps_vec[end].insulator_group.radius_ext)
		push!(lines, "CABLE,CA-$(pidx),$(_fmt(outer_R))")

		# Strict connection vector
		conn = getfield(cable, :conn)
		if !(conn isa AbstractVector)
			throw(
				ArgumentError(
					"cable.conn must be a Vector of Int mappings (0 or 1..$num_phases) for cable index $pidx.",
				),
			)
		end
		if length(conn) < ncomp
			throw(
				ArgumentError(
					"cable.conn length $(length(conn)) < number of components $ncomp for cable index $pidx.",
				),
			)
		end

		# Emit COMPONENT lines (same syntax for CORE/SHEATH/ARMOR)
		for i in 1:ncomp
			label = _TRALIN_COMP[i]
			comp = comps_vec[i]
			comp_id = String(getfield(comp, :id))  # <-- component name from your datamodel

			conn_val = Int(conn[i])  # 0 or 1..N phases

			cond_group = comp.conductor_group
			ins_group  = comp.insulator_group
			cond_props = comp.conductor_props
			ins_props  = comp.insulator_props

			rin  = _fmt(cond_group.radius_in)
			rex  = _fmt(cond_group.radius_ext)
			rho  = _fmt(cond_props.rho/ρ₀) # values in TRALIN are normalized to match the annealed copper
			muC  = _fmt(cond_props.mu_r)
			epsI = _fmt(ins_props.eps_r)  # coating εr

			# COMPONENT,<id-string>,<conn-int>,<Rout>,<Rin>,<rho>,<mu_r>,0,<eps>
			push!(lines, "$label,$comp_id,$conn_val,$rex,$rin,$rho,$muC,0,$epsI")
		end
	end

	push!(lines, "ENDPROGRAM")

	try
		open(file_name, "w") do fid
			for ln in lines
				write(fid, ln);
				write(fid, '\n')
			end
		end
		@info "TRALIN file saved to: $(display_path(file_name))"
		return file_name
	catch e
		@error "Failed to write TRALIN file '$(display_path(file_name))'" exception =
			(e, catch_backtrace())
		return nothing
	end
end


# --- internal utility: slice a block between an anchor and the next page header ---
# Finds the first line that contains `anchor` and returns the lines up to (but not including)
# the next "TRALIN package - PAGE" header. Throws if not found.
function _block_after_anchor(fileLines::Vector{String}, anchor::AbstractString)
	start_idx = findfirst(l -> occursin(anchor, l), fileLines)
	start_idx === nothing && throw(ArgumentError("Anchor not found: $anchor"))

	# page header appears after each page break; we stop before it
	page_hdr = "TRALIN package - PAGE"
	stop_idx = findnext(l -> occursin(page_hdr, l), fileLines, start_idx + 1)
	stop_idx === nothing && (stop_idx = length(fileLines) + 1)

	# drop the anchor line itself and the terminating page header (if any)
	return fileLines[(start_idx+1):(stop_idx-1)]
end

function _infer_tralin_order(file_or_lines)::Int
	fileLines =
		file_or_lines isa AbstractString ? readlines(String(file_or_lines)) : file_or_lines

	block = _block_after_anchor(
		fileLines,
		"CHARACTERISTICS OF ALL CONDUCTORS",
	)

	# Table rows look like:
	#    1      1     1     1     1     core      0.00000   0.01885  ...
	# Columns (first 5 numbers): CONDUCTOR, GROUP, CABLE, COAX, PHASE
	# We capture the 5th integer (PHASE) and keep nonzero uniques.
	phase_set = Set{Int}()
	row_re = r"^\s*\d+\s+\d+\s+\d+\s+\d+\s+(\d+)\s+\S+"

	for ln in block
		m = match(row_re, ln)
		if m !== nothing
			ph = parse(Int, m.captures[1])
			if ph != 0
				push!(phase_set, ph)
			end
		end
	end

	isempty(phase_set) && throw(
		ArgumentError(
			"Could not infer phase count from the 'CHARACTERISTICS OF ALL CONDUCTORS' table.",
		),
	)
	return length(phase_set)
end

# --- public: extract the frequency vector from the "FREQUENCY OF HARMONIC CURRENT" section ---
"""
	extract_tralin_frequencies(file_or_lines) -> Vector{Float64}

Parses the list of operating frequencies from the `FREQUENCY OF HARMONIC CURRENT:` section
up to the next page header. Returns a `Vector{Float64}` in \\[Hz\\].

Accepts either a filename (`AbstractString`) or a preloaded `Vector{String}` with file lines.
"""
function _extract_tralin_frequencies(file_or_lines)::Vector{Float64}
	fileLines =
		file_or_lines isa AbstractString ? readlines(String(file_or_lines)) : file_or_lines

	block = _block_after_anchor(
		fileLines,
		"FREQUENCY OF HARMONIC CURRENT:",
	)

	# Data lines look like:
	#       1       1.00
	#       6      0.215E+04
	# We capture the second column as a float (supports E-notation).
	freqs = Float64[]
	row_re = r"^\s*\d+\s+([+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?)\s*$"

	for ln in block
		m = match(row_re, ln)
		if m !== nothing
			push!(freqs, parse(Float64, m.captures[1]))
		end
	end

	isempty(freqs) && throw(
		ArgumentError("No frequency lines found under 'FREQUENCY OF HARMONIC CURRENT:'."),
	)

	return freqs
end

"""
	parse_tralin_file(filename)

Parse a TRALIN file and extract impedance, admittance, and potential coefficient matrices
for multiple frequency samples.
"""
function parse_tralin_file(filename)
	fileLines = readlines(filename)

	ord = _infer_tralin_order(fileLines)
	freqs = _extract_tralin_frequencies(fileLines)

	# Get all occurrences of "GROUND WIRES ELIMINATED"
	limited_str = "GROUND WIRES ELIMINATED"
	all_idx = findall(row -> occursin(limited_str, row), fileLines)

	# Initialize arrays to store matrices for all frequency samples
	Z_matrices = Vector{Matrix{ComplexF64}}(undef, length(all_idx))
	Y_matrices = Vector{Matrix{ComplexF64}}(undef, length(all_idx))
	P_matrices = Vector{Matrix{ComplexF64}}(undef, length(all_idx))

	# Loop through each frequency block
	for (k, start_idx) in enumerate(all_idx)

		# Slice the file from the current "GROUND WIRES ELIMINATED" position to end
		block_lines = fileLines[start_idx:end]

		# Extract matrices for this frequency sample, ensuring output is ComplexF64
		Z_matrices[k] = Complex{Float64}.(
			extract_tralin_variable(
				block_lines,
				ord,
				"SERIES IMPEDANCES - (ohms/kilometer)",
				"SHUNT ADMITTANCES (microsiemens/kilometer)",
			),
		)
		Y_matrices[k] = Complex{Float64}.(
			extract_tralin_variable(
				block_lines,
				ord,
				"SHUNT ADMITTANCES (microsiemens/kilometer)",
				"SERIES ADMITTANCES (siemens.kilometer)",
			),
		)
		P_matrices[k] = Complex{Float64}.(
			extract_tralin_variable(
				block_lines,
				ord,
				"POTENTIAL COEFFICIENTS (meghoms.kilometer)",
				"SERIES IMPEDANCES - (ohms/kilometer)",
			),
		)
	end

	# Convert lists of matrices into 3D arrays for each matrix type
	Z_stack = reshape(hcat(Z_matrices...), ord, ord, length(Z_matrices))
	Y_stack = reshape(hcat(Y_matrices...), ord, ord, length(Y_matrices))
	P_stack = reshape(hcat(P_matrices...), ord, ord, length(P_matrices))

	Z_stack = Z_stack ./ 1000
	Y_stack = Y_stack .* 1e-6 ./ 1000
	P_stack = P_stack .* 1e6 .* 1000

	return freqs, Z_stack, Y_stack, P_stack
end

"""
	extract_tralin_variable(fileLines, order, str_init, str_final)

Extracts matrix data between specified headers in `fileLines`, handling complex formatting.
"""
function extract_tralin_variable(fileLines, order, str_init, str_final)
	# Locate header and footer lines
	variable_init = findfirst(line -> occursin(str_init, line), fileLines)
	variable_final = findfirst(line -> occursin(str_final, line), fileLines)

	if isnothing(variable_init) || isnothing(variable_final)
		println("Could not locate start or end of the block.")
		return zeros(ComplexF64, order, order)
	end

	# Parse the relevant lines into a list of complex numbers
	variable_list_number = []
	for line in fileLines[(variable_init+15):(variable_final-1)]
		numbers = take_complex_list(line)
		if !isempty(numbers)
			push!(variable_list_number, numbers)
		end
	end

	# Process, clean, and arrange data into matrix form
	variable_list_number = clean_variable_list(variable_list_number, order)

	# Initialize matrix and fill, with padding if necessary
	matrix = zeros(ComplexF64, order, order)
	for (i, row) in enumerate(variable_list_number)
		matrix[i, 1:length(row)] = row
	end

	# Make symmetric by filling lower triangle
	matrix += tril(matrix, -1)'

	return matrix
end


"""
	take_complex_list(s)

Parses a string to identify real and complex numbers, with conditional scaling for scientific notation.
"""
function take_complex_list(s)
	numbers = []

	# Match the first real number (decimal, integer, or scientific notation)
	first_real_pattern = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)"
	first_real_match = match(first_real_pattern, s)
	if !isnothing(first_real_match)
		real_part_str = strip(first_real_match.match)
		real_value =
			occursin(r"[Ee]", real_part_str) ? parse(Float64, real_part_str) :
			parse(Float64, real_part_str) * 1
		push!(numbers, real_value)
	end

	# Match complex numbers (handles scientific notation or regular float, allowing extra whitespace before 'j')
	complex_pattern =
		r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)\s*\+\s*j\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)"
	for m in eachmatch(complex_pattern, s)
		real_part_str, imag_part_str = m.captures
		real_value =
			occursin(r"[Ee]", real_part_str) ? parse(Float64, real_part_str) :
			parse(Float64, real_part_str) * 1
		imag_value =
			occursin(r"[Ee]", imag_part_str) ? parse(Float64, imag_part_str) :
			parse(Float64, imag_part_str) * 1
		push!(numbers, Complex(real_value, imag_value))
	end

	return numbers
end


"""
	clean_variable_list(variable_list_number, order)

Cleans and arranges extracted list into a proper matrix format.
"""
function clean_variable_list(data, order)
	# Remove entries that lack values, filter short lists
	filter!(lst -> length(lst) > 1, data)

	# Trim row label elements and only keep the actual data
	data = [lst[2:end] for lst in data]

	# Apply padding to each row as needed to align with specified order
	data_padded = [vcat(lst, fill(0.0 + 0.0im, order - length(lst))) for lst in data]

	# Ensure `data_padded` has `order` rows; add extra rows of zeros if required
	if length(data_padded) < order
		for _ in 1:(order-length(data_padded))
			push!(data_padded, fill(0.0 + 0.0im, order))
		end
	end

	return data_padded
end

# -- Direct TRALIN constructor
function LineParameters(::Val{:tralin}, file_name::AbstractString)
	f, Z_tralin, Y_tralin, _ = parse_tralin_file(file_name)

	# Normalize types (ComplexF64 / Float64 by default; tweak if you need Measurements etc.)
	Z = ComplexF64.(Z_tralin)
	Y = ComplexF64.(Y_tralin)
	fv = Float64.(f)

	return LineParameters(SeriesImpedance(Z), ShuntAdmittance(Y), fv)
end

# -- Format-auto convenience (add branches as you implement other parsers)
function LineParameters(file_name::AbstractString; format::Symbol = :auto)
	fmt =
		format === :auto ? (endswith(lowercase(file_name), ".f09") ? :tralin : :unknown) :
		format
	if fmt === :tralin
		return LineParameters(Val(:tralin), file_name)
	else
		throw(
			ArgumentError("Unknown/unsupported format for '$file_name' (format=$format)."),
		)
	end
end

# helpful fallback for unknown symbols (better than a MethodError)
LineParameters(::Val{fmt}, args...; kwargs...) where {fmt} =
	throw(ArgumentError("Unsupported format: $(fmt)"))

@inline LineParameters(fmt::Symbol, args...; kwargs...) =
	LineParameters(Val(fmt), args...; kwargs...)
