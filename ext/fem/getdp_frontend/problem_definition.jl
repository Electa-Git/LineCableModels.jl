# Problem: defining the problem

using .GetDPFrontend: add_raw_code, comment, code
using .GetDPFrontend: Group, Function, Constraint, FunctionSpace, Jacobian, Integration
using .GetDPFrontend: Formulation, Resolution, PostProcessing, PostOperation
"""
    Problem

The main problem definition class that brings together all the components.
"""
mutable struct Problem
    _GETDP_CODE::Vector{String}
    filename::Union{String,Nothing}
    group::Group
    function_obj::Vector{Function}
    constraint::Constraint
    functionspace::FunctionSpace
    jacobian::Jacobian
    integration::Integration
    formulation::Formulation
    resolution::Resolution
    postprocessing::PostProcessing
    postoperation::PostOperation
    objects::Vector{String}

    function Problem(; gmsh_major_version=nothing)
        _GETDP_CODE = ["// File created with GetDP.jl: https://github.com/Electa-Git/GetDP.jl.\n"]

        group = Group()
        function_obj = Function[]  # Initialize as empty vector
        constraint = Constraint()
        functionspace = FunctionSpace()
        jacobian = Jacobian()
        integration = Integration()
        formulation = Formulation()
        resolution = Resolution()
        postprocessing = PostProcessing()
        postoperation = PostOperation()

        objects = [
            "group",
            "function_obj",
            "constraint",
            "functionspace",
            "jacobian",
            "integration",
            "formulation",
            "resolution",
            "postprocessing",
            "postoperation"
        ]

        new(
            _GETDP_CODE,
            nothing,
            group,
            function_obj,
            constraint,
            functionspace,
            jacobian,
            integration,
            formulation,
            resolution,
            postprocessing,
            postoperation,
            objects
        )
    end
end

"""
    get_code(problem::Problem)

Returns properly formatted GetDP code.
"""
function get_code(problem::Problem)
    return join(problem._GETDP_CODE, "")
end

"""
    add_raw_code!(problem::Problem, raw_code, newline=true)

Add raw code to the Problem object.
"""
function add_raw_code!(problem::Problem, raw_code, newline=true)
    problem._GETDP_CODE = [add_raw_code(get_code(problem), raw_code, newline)]
end

"""
    add_comment!(problem::Problem, comment_text, newline=true)

Add a comment to the Problem object.
"""
function add_comment!(problem::Problem, comment_text, newline=true)
    add_raw_code!(problem, comment(comment_text; newline=false), newline)
end

"""
    make_problem!(problem::Problem)

Generate the GetDP code for all objects in the Problem, including only non-empty components.
"""
function make_problem!(problem::Problem)
    for attr in problem.objects
        if attr == "function_obj"
            for func in problem.function_obj  # Iterate over all Function objects
                if !isempty(func.content)  # Check if functions are defined
                    push!(problem._GETDP_CODE, code(func))
                end
            end
        elseif attr == "group"
            p = getfield(problem, :group)
            if !isempty(p.content)  # Check if groups are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "constraint"
            p = getfield(problem, :constraint)
            if !isempty(p.constraints)  # Check if constraints are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "functionspace"
            p = getfield(problem, :functionspace)
            if !isempty(p.content)  # Check if items are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "jacobian"
            p = getfield(problem, :jacobian)
            if !isempty(p.content)  # Check if jacobians are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "integration"
            p = getfield(problem, :integration)
            if !isempty(p.content)  # Check if integrations are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "formulation"
            p = getfield(problem, :formulation)
            if !isempty(p.items)  # Check if formulations are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "resolution"
            p = getfield(problem, :resolution)
            if !isempty(p.content)  # Check if resolutions are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "postprocessing"
            p = getfield(problem, :postprocessing)
            if !isempty(p.content)  # Check if postprocessings are defined
                push!(problem._GETDP_CODE, code(p))
            end
        elseif attr == "postoperation"
            p = getfield(problem, :postoperation)
            if !isempty(p.content)  # Check if postoperations are defined
                push!(problem._GETDP_CODE, code(p))
            end
        end
    end
end

"""
    write_file(problem::Problem)

Write the GetDP code to a file.
"""
function write_file(problem::Problem)
    if problem.filename === nothing
        problem.filename = tempname()
    end

    # Remove file if it exists to ensure a new file is created
    isfile(problem.filename) && rm(problem.filename)

    open(problem.filename, "w") do f
        write(f, get_code(problem))
    end
end

"""
    include!(problem::Problem, incl_file)

Include another GetDP file.
"""
function include!(problem::Problem, incl_file)
    push!(problem._GETDP_CODE, "\nInclude \"$(incl_file)\";")
end


# """
#     write_multiple_problems(problems::Vector{Problem}, filename::String)

# Write the GetDP code for multiple Problem instances to a single file.
# The version comment is included only once at the top, and each problem's code
# is separated by a comment indicating its index.
# """
# function write_multiple_problems(problems::Vector{Problem}, filename::String)
#     if isempty(problems)
#         Base.error("No problems to write.")
#     end

#     # Define the version comment once, assuming GetDP.VERSION is accessible
#     version = GetDP.VERSION
#     version_comment = "// This code was created by GetDP.jl v$(version).\n"

#     open(filename, "w") do f
#         # Write the version comment at the top
#         write(f, version_comment)

#         # Process each problem
#         for (i, problem) in enumerate(problems)
#             # Generate the GetDP code for this problem
#             make_file!(problem)
#             # Skip the version comment (first element) and join the rest
#             component_code = join(problem._GETDP_CODE[2:end], "")
#             # Add a separator with problem index
#             write(f, "\n// Problem $i\n")
#             write(f, component_code)
#         end
#     end
# end
