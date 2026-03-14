# ==============================================================================
# THE AST JANITOR
# Drills through blocks, escapes, and un-evaluated macrocalls to find and replace 
# the struct definition, ensuring multiple macros stack flawlessly.
# ==============================================================================

# Recursively strips out :escape nodes.
function _strip_escapes(ex)
	if ex isa Expr && ex.head === :escape
		return _strip_escapes(ex.args[1])
	elseif ex isa Expr
		return Expr(ex.head, map(_strip_escapes, ex.args)...)
	else
		return ex
	end
end

# Drills into blocks and macrocalls to find the actual :struct node.
function _get_struct_node(ex)
	ex isa Expr || return nothing
	ex.head === :struct && return ex
	for arg in ex.args
		node = _get_struct_node(arg)
		node !== nothing && return node
	end
	return nothing
end

# Replaces the old struct node with the new one, leaving macrocalls intact.
function _replace_struct(ex, new_struct)
	if ex isa Expr && ex.head === :struct
		return new_struct
	elseif ex isa Expr
		return Expr(ex.head, map(arg -> _replace_struct(arg, new_struct), ex.args)...)
	else
		return ex
	end
end

# Universal field parser. Handles `field::T = default` for BOTH macros.
function _parse_fields(struct_body)
	fields = Any[]
	clean_body_args = Any[]

	for arg in struct_body.args
		if arg isa LineNumberNode || arg isa String
			push!(clean_body_args, arg)
			continue
		end

		has_default = false
		local field_name, default_val, clean_field

		# Detect `field::T = default`
		if arg isa Expr && arg.head === :(=)
			has_default = true
			clean_field = arg.args[1]
			default_val = arg.args[2]
		else
			clean_field = arg
		end

		if clean_field isa Expr && clean_field.head === :(::)
			field_name = clean_field.args[1]
			push!(clean_body_args, clean_field)
		elseif clean_field isa Symbol
			field_name = clean_field
			push!(clean_body_args, clean_field)
		else
			# Inner constructors / garbage passthrough
			push!(clean_body_args, arg)
			continue
		end

		push!(
			fields,
			(
				name = field_name,
				has_default = has_default,
				default_val = has_default ? default_val : nothing,
			),
		)
	end
	return fields, clean_body_args
end

# Rebuilds the AST block natively around whatever macrocalls are left.
function _rebuild_ast(ex, new_struct, new_funcs)
	replaced = _replace_struct(ex, new_struct)
	if replaced isa Expr && replaced.head === :block
		return Expr(:block, replaced.args..., new_funcs...)
	else
		return Expr(:block, replaced, new_funcs...)
	end
end

function _extract_struct_name(struct_node)
	sig = struct_node.args[2]
	if sig isa Symbol
		return sig
	elseif sig.head === :(<:)
		left = sig.args[1]
		return left isa Symbol ? left : left.args[1]
	elseif sig.head === :curly
		return sig.args[1]
	else
		error("Malformed struct signature.")
	end
end


# ==============================================================================
# MACRO 1: @gridspace (The Combinatorics Hook)
# ==============================================================================
macro gridspace(expr)
	raw_ast = _strip_escapes(expr)
	struct_node = _get_struct_node(raw_ast)
	struct_node === nothing && error("@gridspace must be applied to a struct.")

	struct_name = _extract_struct_name(struct_node)
	fields, clean_body_args = _parse_fields(struct_node.args[3])

	# We strictly enforce the removal of `= default` from the struct definition,
	# otherwise Julia will crash when it tries to compile the dumb struct.
	is_mutable = struct_node.args[1]
	struct_sig = struct_node.args[2]
	clean_struct = Expr(:struct, is_mutable, struct_sig, Expr(:block, clean_body_args...))

	kw_params = Any[]
	grid_calls = Any[]
	for f in fields
		if f.has_default
			push!(kw_params, Expr(:kw, f.name, f.default_val))
		else
			push!(kw_params, f.name)
		end
		push!(grid_calls, :(Grid($(f.name))))
	end

	params = Expr(:parameters, kw_params...)
	call_sig = Expr(:call, struct_name, params)
	grid_tuple = Expr(:tuple, grid_calls...)

	kw_func = Expr(:function, call_sig, quote
		return Gridspace{$struct_name}($grid_tuple)
	end)

	final_ast = _rebuild_ast(raw_ast, clean_struct, Any[kw_func])
	return esc(final_ast)
end


# ==============================================================================
# MACRO 2: @relax (The Type Promoter & Converter)
# ==============================================================================
macro relax(expr)
	raw_ast = _strip_escapes(expr)
	struct_node = _get_struct_node(raw_ast)
	struct_node === nothing && error("@relax must be applied to a struct.")

	struct_name = _extract_struct_name(struct_node)

	# We use your universal parser. If gridspace ran first, defaults are already gone.
	# If relax is running first (standalone or wrong order), we catch the illegal syntax here.
	fields, clean_body_args = _parse_fields(struct_node.args[3])

	if any(f -> f.has_default, fields)
		error(
			"""
	  Syntax Error: @relax encountered default values in struct '$struct_name'.
	  @relax only handles strict type promotion, not kwargs. 
	  If you want combinatorics and defaults, stack the macros as: `@gridspace @relax struct...` so @gridspace processes and strips the defaults before @relax sees them.
	  """,
		)
	end

	field_names = [f.name for f in fields]

	typeof_calls = [:(typeof($n)) for n in field_names]
	convert_args = [:(convert(T, $n)) for n in field_names]
	convert_m_args = [:(convert(T, m.$n)) for n in field_names]

	# If an idiot bypasses kwargs and feeds strings in here, promote_type yields Any.
	# The Target{Any} instantiation will throw a MethodError, which is exactly what they deserve.
	pos_func = quote
		@inline function $struct_name($(field_names...))
			T = promote_type($(typeof_calls...))
			return $struct_name{T}($(convert_args...))
		end
	end

	conv_func = quote
		@inline function Base.convert(
			::Type{$struct_name{T}},
			m::$struct_name,
		) where {T <: Real}
			return $struct_name{T}($(convert_m_args...))
		end
	end

	# Rebuild the AST exactly as it came in, appending the new functions.
	# Since we already validated no defaults exist, we can safely reuse struct_node.
	final_ast = _rebuild_ast(raw_ast, struct_node, Any[pos_func, conv_func])
	return esc(final_ast)
end
# macro relax(expr)
# 	raw_ast = _strip_escapes(expr)
# 	struct_node = _get_struct_node(raw_ast)
# 	struct_node === nothing && error("@relax must be applied to a struct.")

# 	struct_name = _extract_struct_name(struct_node)
# 	fields, _ = _parse_fields(struct_node.args[3])

# 	# Note: @relax doesn't care if defaults are present or not. It just reads the names.
# 	# It leaves the struct node unmodified. Let @gridspace handle the cleanup.
# 	field_names = [f.name for f in fields]

# 	typeof_calls = [:(typeof($n)) for n in field_names]
# 	convert_args = [:(convert(T, $n)) for n in field_names]
# 	convert_m_args = [:(convert(T, m.$n)) for n in field_names]

# 	# If a user bypasses kwargs and feeds strings in here, promote_type yields Any.
# 	# The Target{Any} instantiation will throw a MethodError, which is exactly what they deserve.
# 	pos_func = quote
# 		function $struct_name($(field_names...))
# 			T = promote_type($(typeof_calls...))
# 			return $struct_name{T}($(convert_args...))
# 		end
# 	end

# 	conv_func = quote
# 		function Base.convert(::Type{$struct_name{T}}, m::$struct_name) where {T <: Real}
# 			return $struct_name{T}($(convert_m_args...))
# 		end
# 	end

# 	final_ast = _rebuild_ast(raw_ast, struct_node, Any[pos_func, conv_func])
# 	return esc(final_ast)
# end
