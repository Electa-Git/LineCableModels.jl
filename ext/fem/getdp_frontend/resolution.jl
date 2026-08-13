# Resolution: defining how to solve the problem
using .GetDPFrontend: add_raw_code, comment, make_args

"""
	SystemItem

A system item in a resolution.
"""
mutable struct SystemItem
	name::String
	formulation::String
	Type::Union{String, Nothing}
	Frequency::Union{String, Nothing}
	comment::Union{String, Nothing}
	kwargs::Dict

	function SystemItem(
		name,
		formulation;
		Type = nothing,
		Frequency = nothing,
		comment = nothing,
		kwargs...,
	)
		new(name, formulation, Type, Frequency, comment, Dict(kwargs))
	end
end
# mutable struct SystemItem
#     name::String
#     formulation::String
#     Type::String
#     Frequency::String
#     comment::Union{String,Nothing}
#     kwargs::Dict

#     function SystemItem(name, formulation; Type, Frequency, comment=nothing, kwargs...)
#         new(name, formulation, Type, Frequency, comment, Dict(kwargs))
#     end
# end

function code(item::SystemItem)
	c = "{ Name $(item.name); NameOfFormulation $(item.formulation); Type $(item.Type); Frequency $(item.Frequency);"
	for (k, v) in item.kwargs
		c *= " $k $(make_args(v, sep=","));"
	end
	c *= " }"
	if item.comment !== nothing
		c = comment(item.comment, newline = false) * "\n" * c
	end
	c
end

"""
	Operation

An operation in a resolution.
"""
mutable struct Operation
	operations::Vector{String}
	comment::Union{String, Nothing}

	function Operation(; comment = nothing)
		new(String[], comment)
	end
end

function add_operation!(op::Operation, operation::String)
	push!(op.operations, operation)
end

function code(op::Operation)
	code_lines = String[]
	if op.comment !== nothing
		push!(code_lines, comment(op.comment, newline = false))
	end
	for operation in op.operations
		push!(code_lines, operation)
	end
	join(code_lines, "\n")
end

"""
	Resolution

Defining how to solve the problem.
"""
mutable struct Resolution <: AbstractGetDPObject
	name::String
	content::String
	comment::Union{String, Nothing}
	indent::String
	systems::Vector{SystemItem}
	operation::Operation

	function Resolution()
		new("Resolution", "", nothing, " "^4, SystemItem[], Operation())
	end
end

"""
	add!(resolution::Resolution, id, system_name; kwargs...)

Add a resolution with system and operation to the Resolution object.
"""
function add!(
	resolution::Resolution,
	id::String,
	system_name::String;
	NameOfFormulation = nothing,
	Type,
	Frequency,
	Operation,
	comment = nothing,
	kwargs...,
)
	@warn "In future versions this function will be deprecated, causing breaking changes. Use a vector of SystemItem instead"
	# Use id as resolution name, NameOfFormulation for system formulation
	resolution.name = id
	formulation = NameOfFormulation !== nothing ? NameOfFormulation : id
	system = SystemItem(
		system_name,
		formulation;
		Type = Type,
		Frequency = Frequency,
		comment = comment,
		kwargs...,
	)
	push!(resolution.systems, system)

	# Add operations to Operation
	for op in Operation
		add_operation!(resolution.operation, op)
	end

	resolution.content = code(resolution)
end

"""
	add!(resolution::Resolution, id::String, systems::Vector{SystemItem}; Operation::Vector{String})

Add a resolution with one or more systems and an operation to the Resolution object.
"""
function add!(
	resolution::Resolution,
	id::String,
	systems::Vector{SystemItem};
	Operation::Vector{String},
)
	# Use id as resolution name
	resolution.name = id

	# Add all provided system items
	for system in systems
		push!(resolution.systems, system)
	end

	# Add operations to Operation
	for op in Operation
		add_operation!(resolution.operation, op)
	end

	resolution.content = code(resolution)
end

# """
# 	add!(resolution::Resolution, id, system_name; kwargs...)

# Add a resolution with system and operation to the Resolution object.
# """
# function add!(
# 	resolution::Resolution,
# 	id,
# 	system_name;
# 	NameOfFormulation = nothing,
# 	Type,
# 	Frequency,
# 	Operation,
# 	comment = nothing,
# 	kwargs...,
# )
# 	# Use id as resolution name, NameOfFormulation for system formulation
# 	resolution.name = id
# 	formulation = NameOfFormulation !== nothing ? NameOfFormulation : id
# 	system = SystemItem(
# 		system_name,
# 		formulation;
# 		Type = Type,
# 		Frequency = Frequency,
# 		comment = comment,
# 		kwargs...,
# 	)
# 	push!(resolution.systems, system)

# 	# Add operations to Operation
# 	for op in Operation
# 		add_operation!(resolution.operation, op)
# 	end

# 	resolution.content = code(resolution)
# end

"""
	add_raw_code!(resolution::Resolution, raw_code, newline=true)

Add raw code to the Resolution object.
"""
function add_raw_code!(resolution::Resolution, raw_code, newline = true)
	resolution.content = add_raw_code(resolution.content, raw_code, newline)
end

"""
	add_comment!(resolution::Resolution, comment_text, newline=true)

Add a comment to the Resolution object.
"""
function add_comment!(resolution::Resolution, comment_text, newline = true)
	add_raw_code!(resolution, comment(comment_text, newline = false), newline)
end

"""
	code(resolution::Resolution)

Generate GetDP code for a Resolution object.
"""
function code(resolution::Resolution)
	code_lines = ["Resolution {"]
	push!(code_lines, "  { Name $(resolution.name);")
	push!(code_lines, "    System {")
	for system in resolution.systems
		system_code = code(system)
		for line in split(system_code, '\n')
			push!(code_lines, "      " * line)
		end
	end
	push!(code_lines, "    }")
	push!(code_lines, "    Operation {")
	operation_code = code(resolution.operation)
	for line in split(operation_code, '\n')
		push!(code_lines, "      " * line * ";")
	end
	push!(code_lines, "    }")
	push!(code_lines, "  }")
	push!(code_lines, "}")
	join(code_lines, "\n") * "\n"
end
# System struct
mutable struct System
	items::Vector{SystemItem}
	code::String
	comment::Union{String, Nothing}

	function System(; comment = nothing, kwargs...)
		c = "System { "
		if comment !== nothing
			c *= GetDP.comment(comment)
		end
		c *= "\n        }"
		new(SystemItem[], c, comment)
	end
end

function add!(system::System, name, formulation; kwargs...)
	item = SystemItem(name, formulation; kwargs...)
	s = system.code
	n = 10  # Length of "\n        }"
	system.code = s[1:(end-n)] * "\n       " * item.code * s[(end-n+1):end]
	push!(system.items, item)
	return item
end

# Operation methods
function generate(op::Operation, system_id)
	return "Generate[$system_id]"
end

function solve(op::Operation, system_id)
	return "Solve[$system_id]"
end

function solve_again(op::Operation, system_id)
	return "SolveAgain[$system_id]"
end

function set_global_solver_options(op::Operation, char_expression)
	return "SetGlobalSolverOptions['$char_expression']"
end

function generate_jac(op::Operation, system_id)
	return "GenerateJac[$system_id]"
end

function solve_jac(op::Operation, system_id)
	return "SolveJac[$system_id]"
end

function generate_separate(op::Operation, system_id)
	return "GenerateSeparate[$system_id]"
end

function generate_only(op::Operation, system_id, expression_cst_list)
	return "GenerateOnly[$system_id, $(to_getdp_list(expression_cst_list))]"
end

function generate_only_jac(op::Operation, system_id, expression_cst_list)
	return "GenerateOnlyJac[$system_id, $(to_getdp_list(expression_cst_list))]"
end

function generate_group(op::Operation)
	return "GenerateGroup[]"
end

function generate_right_hand_side_group(op::Operation)
	return "GenerateRightHandSdeGroup[]"
end

function update(op::Operation, system_id; expression = nothing)
	if expression !== nothing
		return "Update[$system_id, $expression]"
	else
		return "Update[$system_id]"
	end
end

function update_constraint(op::Operation, system_id, group_id, constraint_type)
	return "UpdateConstraint[$system_id, $group_id, $constraint_type]"
end

function get_residual(op::Operation, system_id, variable_id)
	return "GetResidual[$system_id, $variable_id]"
end

function get_norm_solution(op::Operation, system_id, variable_id)
	return "GetNormSolution[$system_id, $variable_id]"
end

function get_norm_right_hand_side(op::Operation, system_id, variable_id)
	return "GetNormRightHandSide[$system_id, $variable_id]"
end

function get_norm_residual(op::Operation, system_id, variable_id)
	return "GetNormResidual[$system_id, $variable_id]"
end

function get_norm_increment(op::Operation, system_id, variable_id)
	return "GetNormIncrement[$system_id, $variable_id]"
end

function swap_solution_and_residual(op::Operation, system_id)
	return "SwapSolutionAndResidual[$system_id]"
end

function swap_solution_and_right_hand_side(op::Operation, system_id)
	return "SwapSolutionAndRightHandSide[$system_id]"
end

function init_solution(op::Operation, system_id)
	return "InitSolution[$system_id]"
end

function init_solution1(op::Operation, system_id)
	return "InitSolution1[$system_id]"
end

function create_solution(op::Operation, system_id; expression_cst = nothing)
	if expression_cst !== nothing
		return "CreateSolution[$system_id, $expression_cst]"
	else
		return "CreateSolution[$system_id]"
	end
end

function apply(op::Operation, system_id)
	return "Apply[$system_id]"
end

function set_solution_as_right_hand_side(op::Operation, system_id)
	return "SetSolutionAsRightHandSide[$system_id]"
end

function set_right_hand_side_as_solution(op::Operation, system_id)
	return "SetRightHandSideAsSolution[$system_id]"
end

function residual(op::Operation, system_id)
	return "Residual[$system_id]"
end

function copy_solution(op::Operation, system_id, expression; reverse = false)
	expr = isa(expression, Array) ? "$(expression[1])()" : "'$expression'"
	if reverse
		return "CopySolution[$expr, $system_id]"
	else
		return "CopySolution[$system_id, $expr]"
	end
end

function copy_right_hand_side(op::Operation, system_id, expression; reverse = false)
	expr = isa(expression, Array) ? "$(expression[1])()" : "'$expression'"
	if reverse
		return "CopyRightHandSide[$expr, $system_id]"
	else
		return "CopyRightHandSide[$system_id, $expr]"
	end
end

function copy_residual(op::Operation, system_id, expression; reverse = false)
	expr = isa(expression, Array) ? "$(expression[1])()" : "'$expression'"
	if reverse
		return "CopyResidual[$expr, $system_id]"
	else
		return "CopyResidual[$system_id, $expr]"
	end
end

function save_solution(op::Operation, system_id)
	return "SaveSolution[$system_id]"
end

function save_solutions(op::Operation, system_id)
	return "SaveSolutions[$system_id]"
end

function remove_last_solution(op::Operation, system_id)
	return "RemoveLastSolution[$system_id]"
end

function transfer_solution(op::Operation, system_id)
	return "TransferSolution[$system_id]"
end

function transfer_init_solution(op::Operation, system_id)
	return "TransferInitSolution[$system_id]"
end

function evaluate(op::Operation, expression_list)
	expr =
		isa(expression_list, Array) ? join(string.(expression_list), ", ") :
		string(expression_list)
	return "Evaluate[$expr]"
end

function set_time(op::Operation, expression)
	return "SetTime[$expression]"
end

function set_time_step(op::Operation, expression)
	return "SetTimeStep[$expression]"
end

function set_dtime(op::Operation, expression)
	return "SetDTime[$expression]"
end

function set_frequency(op::Operation, system_id, expression)
	return "SetFrequency[$system_id, $expression]"
end

function system_command(op::Operation, expression_char)
	return "SystemCommand['$expression_char']"
end

function error(expression_char)
	return "Error['$expression_char']"
end

function test(op::Operation, expression, resolution_op_1; resolution_op_2 = nothing)
	if resolution_op_2 !== nothing
		return "Test[$expression {$resolution_op_1} {$resolution_op_2}]"
	else
		return "Test[$expression {$resolution_op_1}]"
	end
end

function while_loop(op::Operation, expression, resolution_op)
	return "While[$expression {$resolution_op}]"
end

function Break(op::Operation)
	return "Break[]"
end

function sleep(op::Operation, expression)
	return "Sleep[$expression]"
end

function set_extrapolation_order(op::Operation, expression_cst)
	return "SetExtrapolationOrder[$expression_cst]"
end

function print_expr(op::Operation, expression_list; file = nothing, format = nothing)
	expr = to_getdp_list(expression_list)
	if format !== nothing && file !== nothing
		return "Print[$expr, File $file, Format $format]"
	else
		return "Print[$expr]"
	end
end

function print_sys(
	op::Operation,
	system_id;
	file = nothing,
	dof_list = nothing,
	timestep = nothing,
)
	s = "Print[$system_id"
	for (k, v) in Dict("File" => file, "dof_list" => dof_list, "TimeStep" => timestep)
		if v !== nothing
			if k == "dof_list"
				s *= ", $(to_getdp_list(v))"
			else
				s *= ", $k $v"
			end
		end
	end
	s *= "]"
	return s
end

function eigen_solve(op::Operation, system_id, neig, shift_re, shift_im; filter = nothing)
	if filter !== nothing
		return "EigenSolve[$system_id, $neig, $shift_re, $shift_im, $filter]"
	else
		return "EigenSolve[$system_id, $neig, $shift_re, $shift_im]"
	end
end

function fourier_transform(op::Operation, system_id, system_id_dest, freq_list)
	return "FourierTransform[$system_id, $system_id_dest, $(to_getdp_list(freq_list))]"
end

function post_operation(op::Operation, post_operation_id)
	return "PostOperation[$post_operation_id]"
end

function gmsh_read(op::Operation, filename; tag = nothing)
	if tag !== nothing
		return "GmshRead['$filename', $tag]"
	else
		return "GmshRead['$filename']"
	end
end

function gmsh_write(op::Operation, filename, field)
	return "GmshWrite['$filename', $field]"
end

function gmsh_clear_all(op::Operation)
	return "GmshClearAll[]"
end

function delete_file(op::Operation, filename)
	return "DeleteFile['$filename']"
end

function rename_file(op::Operation, filename, field)
	return "RenameFile['$filename', '$field']"
end

function create_directory(op::Operation, dirname)
	return "CreateDirectory['$dirname']"
end

function mpi_set_comm_self(op::Operation)
	return "MPI_SetCommSelf[]"
end

function mpi_set_comm_world(op::Operation)
	return "MPI_SetCommWorld[]"
end

function mpi_barrier(op::Operation)
	return "MPI_Barrier[]"
end

function mpi_broadcast_fields(op::Operation; file_list = nothing)
	if file_list !== nothing
		return "MPI_BroadcastFields[$(to_getdp_list(file_list))]"
	else
		return "MPI_BroadcastFields[]"
	end
end

function mpi_broadcast_variables(op::Operation)
	return "MPI_BroadcastVariables[]"
end

# Placeholder for loop methods (to be expanded based on additional details)
function add_time_loop_theta(op::Operation; kwargs...)
	loop_code = "TimeLoopTheta { // Placeholder }"
	push!(op.loops, loop_code)
	return loop_code
end

function add_time_loop_newmark(op::Operation; kwargs...)
	loop_code = "TimeLoopNewmark { // Placeholder }"
	push!(op.loops, loop_code)
	return loop_code
end

function add_iterative_loop(op::Operation; kwargs...)
	loop_code = "IterativeLoop { // Placeholder }"
	push!(op.loops, loop_code)
	return loop_code
end

# function code(op::Operation) #CHANGE
#     code_lines = ["Operation {"]
#     for operation in op.operations
#         push!(code_lines, "    $operation;")
#     end
#     push!(code_lines, "}")
#     join(code_lines, "\n")
# end


# ResolutionItem struct
mutable struct ResolutionItem
	name::String
	code::String
	items::Vector{Union{Operation, System}}
	comment::Union{String, Nothing}
	kwargs::Dict

	function ResolutionItem(name; comment = nothing, kwargs...)
		c = "{ Name $name"
		for (k, v) in kwargs
			c *= "; $k $(make_args(v, sep=","))"
		end
		c *= "; "
		if comment !== nothing
			c *= GetDP.comment(comment)
		end
		c *= "\n}"
		new(name, c, Union{Operation, System}[], comment, Dict(kwargs))
	end
end

function add_operation!(item::ResolutionItem; kwargs...)
	op = Operation(; kwargs...)
	item.code = item.code[1:(end-2)] * "\n     " * op.code * item.code[(end-1):end]
	push!(item.items, op)
	return op
end

function add_system!(item::ResolutionItem; comment = nothing, kwargs...)
	sys = System(; comment, kwargs...)
	item.code = item.code[1:(end-2)] * "\n     " * sys.code * item.code[(end-1):end]
	push!(item.items, sys)
	return sys
end

function code(item::ResolutionItem)
	s = item.code
	for case in item.items
		s = s[1:(end-2)] * "\n     " * code(case) * s[(end-1):end]
	end
	return s
end
