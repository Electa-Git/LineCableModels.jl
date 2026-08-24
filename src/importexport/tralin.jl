const _TRALIN_COMP = ("CORE", "SHEATH", "ARMOUR")

function export_data(::Val{:tralin},
        cable_system::LineCableSystem,
        earth_props::EarthModel;
        freq = 50.0,
        file_name::Union{String, Nothing} = nothing
)::String

    # -- helpers ---------------------------------------------------------------
    _freqs(x) = x isa AbstractVector ? collect(x) : [x]
    _fmt(x) = string(round(Float64(nominal(x)); digits = 6))
    _maybe(x) = (x === nothing) ? "" : _fmt(x)

    # Prefix the output filename with "tr_", matching the XML exporter.
    if isnothing(file_name)
        file_name = joinpath(@__DIR__, "tr_$(cable_system.system_id).f05")
    else
        dir = dirname(file_name)
        fname = basename(file_name)
        # Add the prefix while preserving the supplied filename.
        prefixed_fname = startswith(fname, "tr_") ? fname : "tr_" * fname
        # Resolve relative paths beside this exporter.
        file_name = isabspath(file_name) ? joinpath(dir, prefixed_fname) :
                    joinpath(@__DIR__, dir, prefixed_fname)
    end

    num_phases = length(cable_system.cables)
    freqs = map(f -> nominal(f), _freqs(freq))

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
        soil = earth_props.layers[end]
        rho = _fmt(soil.rho)
        mu_r = _fmt(soil.mu_r)
        eps_r = _fmt(soil.eps_r)
        push!(lines, "SOIL-TYPE")
        push!(lines, "UNIFORM,$rho,$mu_r,$eps_r")
    else
        # [AIR, TOP, (CENTRAL...), BOTTOM] => HORIZONTAL
        push!(lines, "SOIL-TYPE")
        push!(lines, "HORIZONTAL")

        # AIR: no thickness -> explicit empty field `,,`
        push!(lines, "    LAYER,AIR,1e+18,,1,1")

        n_earth = nlayers - 1
        names = n_earth == 1 ? ["TOP"] :
                n_earth == 2 ? ["TOP", "BOTTOM"] :
                vcat("TOP", fill("CENTRAL", n_earth - 2), "BOTTOM")

        for (eidx, (lname, layer)) in enumerate(zip(names, earth_props.layers[2:end]))
            rho = _fmt(layer.rho)
            mu_r = _fmt(layer.mu_r)
            eps_r = _fmt(layer.eps_r)

            if eidx == n_earth
                # BOTTOM: no thickness -> explicit empty field `,,`
                push!(lines, "    LAYER,$lname,$rho,,$mu_r,$eps_r")
            else
                # TOP/CENTRAL: include thickness if available; otherwise leave it empty to keep the slot
                thk = (
                    isfinite(layer.thickness)
                ) ?
                      _fmt(layer.thickness) : ""
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
        outer_R = nominal(comps_vec[end].insulator_group.r_ex)
        push!(lines, "CABLE,CA-$(pidx),$(_fmt(outer_R))")

        # Strict connection vector
        conn = getproperty(cable, :conn)
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
            comp_id = String(getproperty(comp, :id))  # <-- component name from your datamodel

            conn_val = Int(conn[i])  # 0 or 1..N phases

            cond_group = comp.conductor_group
            ins_group = comp.insulator_group
            cond_props = comp.conductor_props
            ins_props = comp.insulator_props

            rin = _fmt(cond_group.r_in)
            rex = _fmt(cond_group.r_ex)
            rho = _fmt(cond_props.rho / 1.724e-8)
            muC = _fmt(cond_props.mu_r)
            epsI = _fmt(ins_props.eps_r)  # coating εr

            # COMPONENT,<id-string>,<conn-int>,<Rout>,<Rin>,<rho>,<mu_r>,0,<eps>
            push!(lines, "$label,$comp_id,$conn_val,$rex,$rin,$rho,$muC,0,$epsI")
        end
    end

    push!(lines, "ENDPROGRAM")

    open(file_name, "w") do fid
        for ln in lines
            write(fid, ln)
            write(fid, '\n')
        end
    end
    @info "TRALIN file saved to: $(_display_path(file_name))"
    return file_name
end

# --- internal utility: slice a block between an anchor and the next page header ---
# Finds the first line that contains `anchor` and returns the lines up to (but not including)
# the next "TRALIN package - PAGE" header. Throws if not found.
function _block_after_anchor(fileLines::Vector{String}, anchor::AbstractString)
    start_idx = findfirst(l -> occursin(anchor, l), fileLines)
    start_idx === nothing && throw(ArgumentError("Anchor not found: $anchor"))

    # Stop before the next page header.
    page_hdr = "TRALIN package - PAGE"
    stop_idx = findnext(l -> occursin(page_hdr, l), fileLines, start_idx + 1)
    stop_idx === nothing && (stop_idx = length(fileLines) + 1)

    # drop the anchor line itself and the terminating page header (if any)
    return fileLines[(start_idx + 1):(stop_idx - 1)]
end

function _infer_tralin_order(file_or_lines)::Int
    fileLines = file_or_lines isa AbstractString ? readlines(String(file_or_lines)) :
                file_or_lines

    block = _block_after_anchor(
        fileLines,
        "CHARACTERISTICS OF ALL CONDUCTORS"
    )

    # Table rows look like:
    #    1      1     1     1     1     core      0.00000   0.01885  ...
    # Columns (first 5 numbers): CONDUCTOR, GROUP, CABLE, COAX, PHASE
    # Read the fifth integer (PHASE) and retain distinct nonzero values.
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
$(TYPEDSIGNATURES)

Read operating frequencies from a TRALIN output section \\[Hz\\].

# Arguments

- `file_or_lines`: TRALIN output path or preloaded file lines.

# Returns

- Operating frequencies through the next page header \\[Hz\\].

# Errors

- Throws `ArgumentError` when the frequency section is absent or empty.
"""
function _extract_tralin_frequencies(file_or_lines)::Vector{Float64}
    fileLines = file_or_lines isa AbstractString ? readlines(String(file_or_lines)) :
                file_or_lines

    block = _block_after_anchor(
        fileLines,
        "FREQUENCY OF HARMONIC CURRENT:"
    )

    # Data lines look like:
    #       1       1.00
    #       6      0.215E+04
    # Read the second column as a floating-point value with optional E notation.
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
$(TYPEDSIGNATURES)

Parse frequency-indexed impedance, admittance, and potential-coefficient
matrices from a TRALIN output file.

# Arguments

- `filename`: TRALIN output path.

# Returns

- A tuple containing frequencies \\[Hz\\], series impedance \\[Ω/m\\], shunt
  admittance \\[S/m\\], and potential coefficients \\[Ω·m\\].
"""
function parse_tralin_file(filename)
    fileLines = readlines(filename)

    ord = _infer_tralin_order(fileLines)
    freqs = _extract_tralin_frequencies(fileLines)

    # Get all occurrences of "GROUND WIRES ELIMINATED"
    limited_str = "GROUND WIRES ELIMINATED"
    all_idx = findall(row -> occursin(limited_str, row), fileLines)

    # Initialise matrices for all frequency samples.
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
            "SHUNT ADMITTANCES (microsiemens/kilometer)"
        ),
        )
        Y_matrices[k] = Complex{Float64}.(
            extract_tralin_variable(
            block_lines,
            ord,
            "SHUNT ADMITTANCES (microsiemens/kilometer)",
            "SERIES ADMITTANCES (siemens.kilometer)"
        ),
        )
        P_matrices[k] = Complex{Float64}.(
            extract_tralin_variable(
            block_lines,
            ord,
            "POTENTIAL COEFFICIENTS (meghoms.kilometer)",
            "SERIES IMPEDANCES - (ohms/kilometer)"
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
$(TYPEDSIGNATURES)

Parse one upper-triangular complex matrix between two TRALIN section headers.

# Arguments

- `fileLines`: TRALIN output lines beginning before `str_init`.
- `order`: Matrix order.
- `str_init`: Header that begins the matrix section.
- `str_final`: Header that ends the matrix section.

# Returns

- A symmetric `order × order` complex matrix. Missing headers produce a zero
  matrix after writing a diagnostic to standard output.
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
    for line in fileLines[(variable_init + 15):(variable_final - 1)]
        numbers = take_complex_list(line)
        if !isempty(numbers)
            push!(variable_list_number, numbers)
        end
    end

    # Process, clean, and arrange data into matrix form
    variable_list_number = clean_variable_list(variable_list_number, order)

    # Initialise and fill the matrix, padding incomplete rows when necessary.
    matrix = zeros(ComplexF64, order, order)
    for (i, row) in enumerate(variable_list_number)
        matrix[i, 1:length(row)] = row
    end

    # Make symmetric by filling lower triangle
    matrix += tril(matrix, -1)'

    return matrix
end

"""
$(TYPEDSIGNATURES)

Parse the leading row index and `a + j b` values from one TRALIN matrix row.
"""
function take_complex_list(s)
    numbers = []

    # Match the first real number (decimal, integer, or scientific notation)
    first_real_pattern = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)"
    first_real_match = match(first_real_pattern, s)
    if !isnothing(first_real_match)
        real_part_str = strip(first_real_match.match)
        real_value = occursin(r"[Ee]", real_part_str) ? parse(Float64, real_part_str) :
                     parse(Float64, real_part_str) * 1
        push!(numbers, real_value)
    end

    # Match complex numbers (handles scientific notation or regular float, allowing extra whitespace before 'j')
    complex_pattern = r"([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)\s*\+\s*j\s*([-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?|\d+)"
    for m in eachmatch(complex_pattern, s)
        real_part_str, imag_part_str = m.captures
        real_value = occursin(r"[Ee]", real_part_str) ? parse(Float64, real_part_str) :
                     parse(Float64, real_part_str) * 1
        imag_value = occursin(r"[Ee]", imag_part_str) ? parse(Float64, imag_part_str) :
                     parse(Float64, imag_part_str) * 1
        push!(numbers, Complex(real_value, imag_value))
    end

    return numbers
end

"""
$(TYPEDSIGNATURES)

Remove row labels and pad parsed TRALIN rows to `order × order`.
"""
function clean_variable_list(data, order)
    # Remove entries that lack values, filter short lists
    filter!(lst -> length(lst) > 1, data)

    # Remove row labels.
    data = [lst[2:end] for lst in data]

    # Pad each row to the requested order.
    data_padded = [vcat(lst, fill(0.0 + 0.0im, order - length(lst))) for lst in data]

    # Add zero rows to reach the requested order.
    if length(data_padded) < order
        for _ in 1:(order - length(data_padded))
            push!(data_padded, fill(0.0 + 0.0im, order))
        end
    end

    return data_padded
end

# -- Direct TRALIN constructor
function LineParameters(::Val{:tralin}, file_name::AbstractString)
    f, Z_tralin, Y_tralin, _ = parse_tralin_file(file_name)

    # Convert parsed values to the requested numeric types.
    Z = ComplexF64.(Z_tralin)
    Y = ComplexF64.(Y_tralin)
    fv = Float64.(f)

    return LineParameters(SeriesImpedance(Z), ShuntAdmittance(Y), fv)
end

# Select the parser from the requested format.
function LineParameters(file_name::AbstractString; format::Symbol = :auto)
    fmt = format === :auto ? (endswith(lowercase(file_name), ".f09") ? :tralin : :unknown) :
          format
    if fmt === :tralin
        return LineParameters(Val(:tralin), file_name)
    else
        throw(
            ArgumentError("Unknown/unsupported format for '$file_name' (format=$format)."),
        )
    end
end

# Report unsupported format symbols as argument errors.
function LineParameters(::Val{fmt}, args...; kwargs...) where {fmt}
    throw(ArgumentError("Unsupported format: $(fmt)"))
end

@inline LineParameters(fmt::Symbol, args...; kwargs...) = LineParameters(Val(fmt), args...; kwargs...)
