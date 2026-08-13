
# ─────────────────────────── Mappings ───────────────────────────

const MAP = Dict{Symbol, Dict{String, String}}(
	:designation => Dict(
		"N"   => "DIN VDE standard",
		"(N)" => "similar to DIN VDE standard",
	),

	# Conductor material (omitted ⇒ copper)
	:conductor => Dict(
		"A" => "aluminium conductor",
		"-" => "copper conductor",
	),

	# Insulation (omitted ⇒ paper)
	:insulation => Dict(
		"2X" => "cross-linked PE (XLPE)",
		"Y"  => "PVC",
		"H"  => "LSOH compound",
		"-"  => "impregnated paper",
	),

	# Screen / concentric conductor
	:metallic_screen => Dict(
		"CW" => "concentric conductor of copper in waveconal formation",
		"CE" => "concentric conductor of copper over each individual core",
		"SE" => "screen of copper wires over each individual core",
		"C"  => "concentric conductor of copper",
		"S"  => "screen of copper wires",
		"H"  => "conductive layers",
	),

	# Water blocking right after screen, inside parentheses
	:waterblocking => Dict(
		"FL" => "longitudinally and radially water-proof protection",
		"F"  => "longitudinally water-proof protection",
		"L"  => "radially water-proof protection",
	),

	# Inner sheath (e.g., …XSH… : H before outer sheath)
	:inner_sheath => Dict(
		"H" => "LSOH compound inner sheath",
	),

	# Armouring
	:armouring => Dict(
		"B" => "steel tape armouring",
		"F" => "armour of galvanised flat steel wires",
		"G" => "counter helix of galvanised steel tape",
		"R" => "armour of galvanised round steel wires",
	),

	# Sheath
	:sheath => Dict(
		"KL" => "aluminium sheath",
		"K"  => "lead sheath",
	),

	# Outer sheath
	:outer_sheath => Dict(
		"A"  => "outer sheath made of fibrous material",
		"2Y" => "PE outer sheath",
		"Y"  => "PVC outer sheath",
	),

	# Grounding / protective conductor suffix
	:grounding => Dict(
		"I" => "with grounding (protective) conductor",
		"J" => "with grounding (protective) conductor",
		"O" => "without grounding (protective) conductor",
	),
)

# Canonical order of appearance within the stub
const ORDER = [
	:designation,
	:conductor,
	:insulation,
	:metallic_screen,
	:waterblocking,   # immediately after :metallic_screen, parenthesized
	:inner_sheath,    # H before metallic sheath (e.g., …XSH…)
	:sheath,
	:armouring,
	:outer_sheath,
	:grounding,
]

# ─────────────────────── Regex builders ────────────────────────

escape_for_rx(s) = replace(s, r"([.^$|?*+\[\]{}\\])" => "\\\\\1")

# Return the non-capturing alternation **as a string**, e.g. "(?:FL|F|L)"
function union_pat_str(tokens::Vector{String})
	ts = sort(tokens; by = length, rev = true) # longest first
	"(?:" * join(escape_for_rx.(ts), "|") * ")"
end

# If you really want a Regex anchored at start, wrap union_pat_str
union_pat(tokens::Vector{String}) = Regex("^" * union_pat_str(tokens))

const RXS = let rxs = Dict{Symbol, Regex}()
	for fld in ORDER
		fld_keys = collect(keys(MAP[fld]))
		if fld == :waterblocking
			# parenthesized immediately after :metallic_screen (e.g., "(FL)"), anchored
			core = union_pat_str(fld_keys)          # "(?:FL|F|L)"
			rxs[fld] = Regex("^\\(" * core * "\\)") # "^\((?:FL|F|L)\)"
		else
			rxs[fld] = union_pat(fld_keys)          # e.g. "^(?:CE|SE|CW|S|C|H)"
		end
	end
	rxs
end

# Trailing specs (anchored at start of tail)
const RX_CORES_X_CSA = r"^(\d+)\s*[x×]\s*(\d+(?:\.\d+)?)(?:\s*/\s*(\d+(?:\.\d+)?))?"
const RX_VOLTAGE     = r"^(\d+(?:\.\d+)?)\s*/\s*(\d+(?:\.\d+)?)\s*(?:kV|KV|kv)\b"
const RX_TYPE        = r"^([RSEMOH]{1,3})(?:\s*/\s*V)?\b"   # RM, SE, OH, … + optional /V

# KISS conductor type mapping
const TYPE_MAP = Dict(
	'R' => "round", 'S' => "sector", 'O' => "oval",
	'E' => "solid", 'M' => "stranded", 'H' => "hollow",
	'V' => "compact",
)

function decode_type(code::AbstractString; has_compact::Bool = false)
	words = String[]
	for c in code
		if haskey(TYPE_MAP, c)
			push!(words, TYPE_MAP[c])
		else
			@warn "Unknown conductor type letter ignored." letter=String(c)
		end
	end
	if has_compact
		push!(words, "compact")
	end
	# de-dup preserving order
	seen = Set{String}();
	uniq = String[]
	for w in words
		if !(w in seen)
			;
			push!(uniq, w);
			push!(seen, w);
		end
	end
	return join(uniq, ", ")
end

# ─────────────────────────── Parser ────────────────────────────

"""
	vdeparse(code::AbstractString) -> Dict{Symbol,String}

Parses VDE/DIN 0271/0276 cable codes:
- **stub** (first non-space token): designation → conductor_material (default copper) → insulation (default paper) → screen → waterblocking → inner_sheath → armouring → sheath → grounding
- **tail**: cores × cross-section (optional screen csa), voltage, conductor type (R/S/O + E/M/H, optional `/V` ⇒ compact)

Only parsed keys are returned; defaults are materialized when omitted.
"""
function vdeparse(code::AbstractString)::Dict{Symbol, String}
	# normalize spaces (NBSP -> space) and trim
	s = replace(code, '\u00A0' => ' ')
	s = strip(s)

	# split into stub token (first non-space chunk) + tail
	m = match(r"^\S+", s)
	if m === nothing
		return Dict{Symbol, String}()  # empty / whitespace line
	end
	stub = m.match                   # e.g., "2XS(F)2Y"
	tail = strip(s[(length(stub)+1):end])

	# parse stub left-to-right
	out = Dict{Symbol, String}()
	rest = stub
	for fld in ORDER
		rx = RXS[fld]
		mm = match(rx, rest)
		if mm !== nothing
			tok = mm.match
			key = (fld == :waterblocking) ? tok[2:(end-1)] : tok
			out[fld] = MAP[fld][key]
			rest = rest[(length(tok)+1):end]  # consume
		else
			# materialize normative omissions
			if fld == :conductor
				out[fld] = "copper conductor"
			elseif fld == :insulation
				out[fld] = "impregnated paper"
			end
		end
	end
	# Any non-empty rest means unknown extra within the stub itself
	if !isempty(rest)
		out[:unparsed_stub] = rest
		@warn "Unparsed stub residue." residue=rest
	end

	# parse tail in anchored passes
	if !isempty(tail)
		if (mm = match(RX_CORES_X_CSA, tail)) !== nothing
			out[:cores] = mm.captures[1]
			out[:conductor_cross_section] = mm.captures[2]
			if mm.captures[3] !== nothing
				out[:metallic_screen_cross_section] = mm.captures[3]
			end
			tail = strip(tail[(length(mm.match)+1):end])
		end
	end

	if !isempty(tail)
		if (mm = match(RX_VOLTAGE, tail)) !== nothing
			out[:voltage] = "$(mm.captures[1])/$(mm.captures[2]) kV"
			tail = strip(tail[(length(mm.match)+1):end])
		end
	end

	if !isempty(tail)
		if (mm = match(RX_TYPE, tail)) !== nothing
			raw = mm.captures[1]
			has_compact = occursin(r"/\s*V\b", mm.match)
			out[:conductor_type] = decode_type(raw; has_compact = has_compact)
			tail = strip(tail[(length(mm.match)+1):end])
		end
	end

	if !isempty(tail)
		out[:unparsed] = tail
		@warn "Unparsed trailing token(s)." residue=tail
	end

	return out
end

# # ───────────────────────── Smoke tests ─────────────────────────
# exs = [
# 	"2XS(F)2Y 3x240/25 12/20kV RM/V",
# 	"A2XS(FL)2Y 1x630 18/30kV RE",
# 	"N2XSY 1x150 6/10kV SE",
# 	"N2XSH2YI 3x95 0.6/1kV OM",
# 	"2XSY 1x240 RM",
# 	"A2XS(L)KL2Y 3x300/35 20/35kV OH/V",
# ]
# for ex in exs
# 	println(ex, " → ", vdeparse(ex))
# end


