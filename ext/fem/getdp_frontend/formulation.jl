# Formulation: building equations
using .GetDPFrontend: add_raw_code, comment, make_args

"""
	Quantity

A quantity in a formulation.
"""
mutable struct Quantity
	name::String
	Type::String
	NameOfSpace::String
	comment::Union{String, Nothing}

	function Quantity(name::String; Type, NameOfSpace, comment = nothing)
		new(name, Type, NameOfSpace, comment)
	end
end

function code(q::Quantity)
	code_str = "{ Name $(q.name); Type $(q.Type); NameOfSpace $(q.NameOfSpace); }"
	if q.comment !== nothing
		code_str = comment(q.comment, newline = false) * "\n      " * code_str
	end
	code_str
end

"""
	EquationTerm

A term in an equation.
"""
mutable struct EquationTerm
	term_type::String
	term::String
	options::Dict
	comment::Union{String, Nothing}

	function EquationTerm(term_type, term; comment = nothing, kwargs...)
		options = Dict(kwargs)
		new(term_type, term, options, comment)
	end
end

function code(term::EquationTerm)
	code_lines = String[]
	if term.comment !== nothing
		push!(code_lines, comment(term.comment, newline = false))
	end
	if term.term_type == "Galerkin"
		push!(code_lines, "Galerkin { $(term.term);")
		# Order: In, Jacobian, Integration
		for key in [:In, :Jacobian, :Integration]
			if haskey(term.options, key)
				push!(code_lines, "    $key $(make_args(term.options[key], sep=","));")
			end
		end
		# Include any other options
		for (k, v) in term.options
			if k ∉ [:In, :Jacobian, :Integration]
				push!(code_lines, "    $k $(make_args(v, sep=","));")
			end
		end
		push!(code_lines, "}")
	elseif term.term_type == "GlobalTerm"
		push!(code_lines, "GlobalTerm {  $(term.term) ;")
		for (k, v) in term.options
			push!(code_lines, "    $k $(make_args(v, sep=","));")
		end
		push!(code_lines, "}")
	else
		Base.error("Unsupported term type: $(term.term_type)")
	end
	join(code_lines, "\n")
end

"""
	Equation

An equation in a formulation.
"""
mutable struct Equation
	items::Vector{EquationTerm}

	function Equation(; comment = nothing, kwargs...)
		new(EquationTerm[])
	end
end

function code(equation::Equation)
	code_lines = ["Equation {"]
	for term in equation.items
		term_code = code(term)
		for line in split(term_code, '\n')
			push!(code_lines, "  $line")
		end
	end
	push!(code_lines, "}")
	join(code_lines, "\n")
end

"""
	add!(equation::Equation, term_type::String, term::String; kwargs...)

Add an equation term to an equation.
"""
function add!(equation::Equation, term_type::String, term::String; kwargs...)
	item = EquationTerm(term_type, term; kwargs...)
	push!(equation.items, item)
	item
end

"""
	FormulationItem

An item in a formulation.
"""
mutable struct FormulationItem
	Name::String
	Type::String
	comment::Union{String, Nothing}
	items::Vector{Any}  # Will be Vector{Union{Quantity,Equation}}

	function FormulationItem(Name, Type; comment = nothing, kwargs...)
		new(Name, Type, comment, [])
	end
end

function code(item::FormulationItem)
	code_lines = ["{ Name $(item.Name); Type $(item.Type);"]
	if item.comment !== nothing
		push!(code_lines, "  $(comment(item.comment, newline=false))")
	end

	# Quantities
	quantities = filter(x -> isa(x, Quantity), item.items)
	if !isempty(quantities)
		push!(code_lines, "    Quantity {")
		for q in quantities
			push!(code_lines, "      $(code(q))")
		end
		push!(code_lines, "    }")
	end

	# Equation
	equations = filter(x -> isa(x, Equation), item.items)
	if !isempty(equations)
		equation = equations[1]  # Assuming only one equation
		equation_code = code(equation)
		for line in split(equation_code, '\n')
			push!(code_lines, "    $line")
		end
	end

	push!(code_lines, "}")
	join(code_lines, "\n")
end

"""
	add_quantity!(item::FormulationItem, name::String; kwargs...)

Add a quantity to a formulation item.
"""
function add_quantity!(
	item::FormulationItem,
	name::String;
	Type,
	NameOfSpace,
	comment = nothing,
)
	case = Quantity(name; Type = Type, NameOfSpace = NameOfSpace, comment = comment)
	push!(item.items, case)
	case
end

"""
	add_equation!(item::FormulationItem, args...; kwargs...)

Add an equation to a formulation item.
"""
function add_equation!(item::FormulationItem, args...; kwargs...)
	case = Equation(args...; kwargs...)
	push!(item.items, case)
	case
end

"""
	Formulation

Building equations.
"""
mutable struct Formulation <: AbstractGetDPObject
	name::String
	content::String
	comment::Union{String, Nothing}
	indent::String
	items::Vector{FormulationItem}

	function Formulation()
		new("Formulation", "", nothing, " "^4, FormulationItem[])
	end
end

"""
	content(formulation::Formulation)

Get the content of a formulation.
"""
function content(formulation::Formulation)
	code_lines = String[]
	for item in formulation.items
		push!(code_lines, code(item))
	end
	join(code_lines, "\n")
end

"""
	add!(formulation::Formulation, Name, Type; kwargs...)

Add a formulation item to a formulation.
"""
function add!(formulation::Formulation, Name, Type; kwargs...)
	o = FormulationItem(Name, Type; kwargs...)
	push!(formulation.items, o)
	o
end

"""
	add_raw_code!(formulation::Formulation, raw_code, newline=true)

Add raw code to the Formulation object.
"""
function add_raw_code!(formulation::Formulation, raw_code, newline = true)
	formulation.content = add_raw_code(formulation.content, raw_code, newline)
end

"""
	add_comment!(formulation::Formulation, comment_text, newline=true)

Add a comment to the Formulation object.
"""
function add_comment!(formulation::Formulation, comment_text, newline = true)
	add_raw_code!(formulation, comment(comment_text, newline = false), newline)
end

"""
	code(formulation::Formulation)

Generate GetDP code for a Formulation object.
"""
function code(formulation::Formulation)
	code_lines = ["\nFormulation {"]
	formulation_content = content(formulation)
	if !isempty(formulation_content)
		for line in split(formulation_content, '\n')
			if !isempty(line)
				push!(code_lines, "  $line")
			end
		end
	end
	if !isempty(formulation.content)
		for line in split(formulation.content, '\n')
			if !isempty(line)
				push!(code_lines, "  $line")
			end
		end
	end
	push!(code_lines, "}")
	join(code_lines, "\n") * "\n"
end
