# FunctionSpace: defining function spaces

using .GetDPFrontend: add_raw_code, comment, make_args
using .GetDPFrontend: AbstractBase_

# FunctionSpace struct
mutable struct FunctionSpace <: AbstractBase_
	name::String
	content::String
	items::Vector{Dict}
	comment::Union{String, Nothing}
	indent::String

	function FunctionSpace()
		new("FunctionSpace", "", [], nothing, " "^4)
	end
end

# Add a function space item
function add!(
	functionspace::FunctionSpace,
	name,
	operation,
	constraint;
	comment = nothing,
	kwargs...,
)
	item = Dict(
		:name => name,
		:basis_functions => [],
		:global_quantities => [],
		:constraints => [],
		:comment => comment,
		:kwargs => Dict(kwargs),
	)
	for (k, v) in kwargs
		if v !== nothing
			item[:kwargs][k] = v
		end
	end
	push!(functionspace.items, item)
	functionspace.content = code(functionspace)
	return name
end

# Add a basis function
function add_basis_function!(
	functionspace::FunctionSpace,
	name,
	nameOfCoef,
	func;
	Support,
	Entity,
	comment = nothing,
	kwargs...,
)
	if isempty(functionspace.items)
		Base.error("No FunctionSpace items defined. Call add! first.")
	end
	item = functionspace.items[end]  # Target the last added item
	bf = Dict(
		:name => name,
		:nameOfCoef => nameOfCoef,
		:function => func,
		:Support => Support,
		:Entity => Entity,
		:kwargs => Dict(kwargs),
		:comment => comment,
	)
	push!(item[:basis_functions], bf)
	functionspace.content = code(functionspace)
	return name
end

# Add a constraint
function add_constraint!(
	functionspace::FunctionSpace,
	nameOfCoef,
	entityType,
	nameOfConstraint;
	comment = nothing,
	kwargs...,
)
	if isempty(functionspace.items)
		Base.error("No FunctionSpace items defined. Call add! first.")
	end
	item = functionspace.items[end]  # Target the last added item
	c = Dict(
		:nameOfCoef => nameOfCoef,
		:entityType => entityType,
		:nameOfConstraint => nameOfConstraint,
		:kwargs => Dict(kwargs),
		:comment => comment,
	)
	push!(item[:constraints], c)
	functionspace.content = code(functionspace)
	return nameOfCoef
end

# Add a global quantity
function add_global_quantity!(
	functionspace::FunctionSpace,
	name,
	type;
	NameOfCoef,
	comment = nothing,
	kwargs...,
)
	if isempty(functionspace.items)
		Base.error("No FunctionSpace items defined. Call add! first.")
	end
	item = functionspace.items[end]  # Target the last added item
	gq = Dict(
		:name => name,
		:type => type,
		:NameOfCoef => NameOfCoef,
		:kwargs => Dict(kwargs),
		:comment => comment,
	)
	push!(item[:global_quantities], gq)
	functionspace.content = code(functionspace)
	return name
end

function add_raw_code!(functionspace::FunctionSpace, raw_code, newline = true)
	functionspace.content = add_raw_code(functionspace.content, raw_code, newline)
end

function add_comment!(functionspace::FunctionSpace, comment_text, newline = true)
	add_raw_code!(functionspace, comment(comment_text; newline = false), newline)
end

function code(functionspace::FunctionSpace)
	code_lines = String[]
	push!(code_lines, "\nFunctionSpace {")

	for item in functionspace.items
		c = "{ Name $(item[:name])"
		for (k, v) in item[:kwargs]
			if v !== nothing
				c *= "; $k $(make_args(v, sep=","))"
			end
		end
		c *= ";"
		if item[:comment] !== nothing
			c *= " " * comment(item[:comment], newline = false)
		end
		push!(code_lines, "  " * c)

		# Basis Functions
		if !isempty(item[:basis_functions])
			push!(code_lines, "    BasisFunction {")
			for bf in item[:basis_functions]
				bfc = "{ Name $(bf[:name]); NameOfCoef $(bf[:nameOfCoef]); Function $(bf[:function]); Support $(bf[:Support]); Entity $(bf[:Entity])"
				for (k, v) in bf[:kwargs]
					if k ∉ [:condition, :endCondition] && v !== nothing
						bfc *= "; $k $(make_args(v, sep=","))"
					end
				end
				bfc *= "; }"  # Ensure each basis function definition closes with }
				if bf[:comment] !== nothing
					bfc *= " " * comment(bf[:comment], newline = false)
				end
				if haskey(bf[:kwargs], :condition)
					push!(code_lines, "    $(bf[:kwargs][:condition])")
					push!(code_lines, "      " * bfc)
					if haskey(bf[:kwargs], :endCondition)
						push!(code_lines, "    $(bf[:kwargs][:endCondition])")
					end
				else
					push!(code_lines, "      " * bfc)
				end
			end
			push!(code_lines, "    }")
		end

		# Global Quantities
		if !isempty(item[:global_quantities])
			push!(code_lines, "    GlobalQuantity {")
			for gq in item[:global_quantities]
				gqc = "{ Name $(gq[:name]); Type $(gq[:type]); NameOfCoef $(gq[:NameOfCoef]);"
				for (k, v) in gq[:kwargs]
					if k ∉ [:condition, :endCondition] && v !== nothing
						gqc *= "; $k $(make_args(v, sep=","))"
					end
				end
				gqc *= "}"
				if gq[:comment] !== nothing
					gqc *= " " * comment(gq[:comment], newline = false)
				end
				push!(code_lines, "      " * gqc)
			end
			push!(code_lines, "    }")
		end

		# Constraints
		if !isempty(item[:constraints])
			push!(code_lines, "    Constraint {")
			for c in item[:constraints]
				cc = "{ NameOfCoef $(c[:nameOfCoef]); EntityType $(c[:entityType]); NameOfConstraint $(c[:nameOfConstraint]);"
				for (k, v) in c[:kwargs]
					if k ∉ [:condition, :endCondition] && v !== nothing
						cc *= "; $k $(make_args(v, sep=","))"
					end
				end
				cc *= "}"
				if c[:comment] !== nothing
					if haskey(c[:kwargs], :condition)
						push!(code_lines, "    $(c[:kwargs][:condition])")
						push!(code_lines, "      " * comment(c[:comment], newline = false))
						push!(code_lines, "      " * cc)
						if haskey(c[:kwargs], :endCondition)
							push!(code_lines, "    $(c[:kwargs][:endCondition])")
						end
					else
						push!(
							code_lines,
							"      " * cc * " " * comment(c[:comment], newline = false),
						)
					end
				else
					if haskey(c[:kwargs], :condition)
						push!(code_lines, "    $(c[:kwargs][:condition])")
						push!(code_lines, "      " * cc)
						if haskey(c[:kwargs], :endCondition)
							push!(code_lines, "    $(c[:kwargs][:endCondition])")
						end
					else
						push!(code_lines, "      " * cc)
					end
				end
			end
			push!(code_lines, "    }")
		end

		push!(code_lines, "  }")
	end

	push!(code_lines, "}")
	if functionspace.comment !== nothing
		return comment(functionspace.comment) * "\n" * join(code_lines, "\n") * "\n"
	else
		return join(code_lines, "\n") * "\n"
	end
end
