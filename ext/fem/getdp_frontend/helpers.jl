# Helper functions for GetDP.jl

"""
	add_raw_code(s, raw_code, newline=true)

Add raw code to a string, optionally with a newline.
"""
function add_raw_code(s, raw_code, newline = true)
	nl = newline ? "\n" : ""
	return s * nl * raw_code
end

"""
	comment(s, style="short", newline=false)

Add a comment to a string, optionally with a newline.
"""
function comment(s; style = "short", newline = false)
	nl = newline ? "\n" : ""

	if style == "short"
		if length(s) > 80
			# Split long comments into multiple lines
			lines = split_long_string(s)
			return comment(join(lines, "\n"), style = "long", newline = newline)
		else
			return nl * "// " * s
		end
	elseif style == "long"
		return nl * "/*  " * s * " */"
	end
end

"""
	split_long_string(s, width=80)

Split a long string into multiple lines.
"""
function split_long_string(s, width = 80)
	lines = String[]
	words = split(s)
	current_line = ""

	for word in words
		if length(current_line) + length(word) + 1 <= width
			if isempty(current_line)
				current_line = word
			else
				current_line *= " " * word
			end
		else
			push!(lines, current_line)
			current_line = word
		end
	end

	if !isempty(current_line)
		push!(lines, current_line)
	end

	return lines
end

"""
	make_args(glist, sep=",", list_char=false)

Format arguments for GetDP code.
"""

function make_args(glist; sep::String = ", ", list_char::Bool = false)
	sep = sep * " "  # Add space after separator if needed

	if isa(glist, Array)
		if length(glist) == 1
			return "{" * string(glist[1]) * "}"
		else
			formatted = join([string(g) for g in glist], sep)
			return list_char ? "#{$formatted}" : "{$formatted}"
		end
	else
		return string(glist)
	end
end

# """
#     get_getdp_exe()

# Get the path to the GetDP executable.
# """
# function get_getdp_exe()
#     macos_getdp_location = "/Applications/Getdp.app/Contents/MacOS/getdp"
#     return isfile(macos_getdp_location) ? macos_getdp_location : "getdp"
# end

"""
	array2getdplist(l)

Convert a Julia array to a GetDP list.
"""
function array2getdplist(l)
	return "{" * join([string(item) for item in l], ",") * "}"
end


function get_getdp_major_version(getdp_exe = get_getdp_executable())
	out = read(`$getdp_exe --version`, String)
	ex = split(strip(out), ".")
	return parse(Int, ex[1])
end

# """
#     replace_formula(str_in, to_replace, replacement)

# Replace formulas in a string.
# """
# function replace_formula(str_in, to_replace, replacement)
#     str_in = replace(join(split(str_in)), "{" => "", "}" => "")

#     to_replace = [replace(join(split(r)), "{" => "", "}" => "") for r in to_replace]

#     # This is a simplified version of the Python tokenize approach
#     # In a real implementation, we might want to use a proper Julia tokenizer
#     for (rold, rnew) in zip(to_replace, replacement)
#         str_in = replace(str_in, rold => rnew)
#     end

#     return str_in
# end

"""
	get_getdp_executable()

Get the path to the GetDP executable.
"""
function get_getdp_executable()
	# Environment variable takes precedence
	env_path = get(ENV, "GETDP_EXECUTABLE", nothing)
	if env_path !== nothing
		isfile(env_path) || throw(
			ArgumentError(
				"GETDP_EXECUTABLE points to a missing file: $(repr(env_path)).",
			),
		)
		return abspath(env_path)
	end

	# Then system installation
	system_getdp = Sys.which("getdp")
	system_getdp !== nothing && return system_getdp

	throw(
		ArgumentError(
			"GetDP executable not found. Set GETDP_EXECUTABLE or add getdp to PATH.",
		),
	)
end
