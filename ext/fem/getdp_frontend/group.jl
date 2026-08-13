# Group: defining topological entities

using .GetDPFrontend: make_args

"""
	Group

Defining topological entities.
"""
mutable struct Group <: AbstractGetDPObject
	name::String
	content::String
	comment::Union{String, Nothing}
	indent::String
	idlist::Vector{String}

	function Group()
		new("Group", "", nothing, " "^4, String[])
	end
end

"""
	define!(group::Group, id="domain")

Define a group.
"""
function define!(group::Group, id = "domain")
	if isa(id, Array)
		id = join([string(g) for g in id], ", ")
	end

	c = "DefineGroup[$(id)];"
	group.content *= c * "\n"
	return nothing
end
"""
	add_raw_code!(group::group, raw_code, newline=true)

Add raw code to the group object.
"""
function add_raw_code!(group::Group, raw_code, newline = true)
	group.content = add_raw_code(group.content, raw_code, newline)
end
"""
	add_comment!(group::group, comment_text, newline=true)

Add a comment to the group object.
"""
function add_comment!(group::Group, comment_text; newline = true)
	add_raw_code!(group, comment(comment_text, newline = false), newline)
end
"""
	add!(group::Group, id="domain", glist=[1], gtype="Region"; operation="=", comment=nothing, kwargs...)

Add an expression to the Group object. The default group type is "Region". Use operation="=" to define a group or operation="+=" to extend an existing group.
"""
function add!(
	group::Group,
	id = "domain",
	glist = [1],
	gtype = "Region";
	operation = "=",
	comment = nothing,
	kwargs...,
)
	if operation == "=" && id in group.idlist
		Base.error("Identifier $(id) already in use.")
	elseif operation == "+=" && !(id in group.idlist)
		Base.error("Cannot extend non-existing group $(id).")
	end

	glist_str = make_args(glist)

	c = "$(id) $(operation) $(gtype)[ $(glist_str) "

	for (k, v) in kwargs
		if v !== nothing
			c *= ", " * string(k) * " " * make_args(v)
		end
	end

	c *= "];"
	group.content *= c

	if comment !== nothing
		add_comment!(group, comment; newline = false)
	end

	group.content *= "\n"
	if operation == "="
		push!(group.idlist, id)
	end

	return id
end

"""
	add_space!(group::Group, num_spaces::Int=1)

Add a specified number of empty lines to the Group object's content for spacing in the output.
"""
function add_space!(group::Group, num_spaces::Int = 1)
	group.content *= "\n" ^ num_spaces
	return nothing
end

"""
	code(group::Group)

Generate GetDP code for a Group object.
"""
function code(group::Group)
	code_lines = String[]
	push!(code_lines, "\nGroup{")

	if !isempty(group.content)
		for line in split(group.content, '\n')
			if !isempty(line)
				push!(code_lines, "  " * line)
			else
				push!(code_lines, "")  # Intended to preserve blank lines
			end
		end
	end

	push!(code_lines, "}")
	return join(code_lines, "\n")
end

"""
	Region(group::Group, args...; kwargs...)

Regions in R1.
"""
function Region(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "Region", kwargs...)
end

"""
	Global(group::Group, args...; kwargs...)

Regions in R1 (variant of Region used with global BasisFunctions BF_Global and BF_dGlobal).
"""
function Global(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "Global", kwargs...)
end

"""
	NodesOf(group::Group, args...; Not=nothing, kwargs...)

Nodes of elements of R1 (Not: but not those of R2).
"""
function NodesOf(group::Group, args...; Not = nothing, kwargs...)
	return add!(group, args...; gtype = "NodesOf", Not = Not, kwargs...)
end

"""
	EdgesOf(group::Group, args...; Not=nothing, kwargs...)

Edges of elements of R1 (Not: but not those of R2).
"""
function EdgesOf(group::Group, args...; Not = nothing, kwargs...)
	return add!(group, args...; gtype = "EdgesOf", Not = Not, kwargs...)
end

"""
	FacetsOf(group::Group, args...; Not=nothing, kwargs...)

Facets of elements of R1 (Not: but not those of R2).
"""
function FacetsOf(group::Group, args...; Not = nothing, kwargs...)
	return add!(group, args...; gtype = "FacetsOf", Not = Not, kwargs...)
end

"""
	VolumesOf(group::Group, args...; Not=nothing, kwargs...)

Volumes of elements of R1 (Not: but not those of R2).
"""
function VolumesOf(group::Group, args...; Not = nothing, kwargs...)
	return add!(group, args...; gtype = "VolumesOf", Not = Not, kwargs...)
end

"""
	ElementsOf(group::Group, args...; OnOneSideOf=nothing, OnPositiveSideOf=nothing, Not=nothing, kwargs...)

Elements of regions in R1.

- OnOneSideOf: only elements on one side of R2 (non-automatic, i.e., both sides if both in R1)
- OnPositiveSideOf: only elements on positive (normal) side of R2
- Not: but not those touching only its skin R3 (mandatory for free skins for correct separation of side layers)
"""
function ElementsOf(
	group::Group,
	args...;
	OnOneSideOf = nothing,
	OnPositiveSideOf = nothing,
	Not = nothing,
	kwargs...,
)
	return add!(
		group,
		args...;
		gtype = "ElementsOf",
		OnOneSideOf = OnOneSideOf,
		OnPositiveSideOf = OnPositiveSideOf,
		Not = Not,
		kwargs...,
	)
end

"""
	GroupsOfNodesOf(group::Group, args...; kwargs...)

Groups of nodes of elements of R1 (a group is associated with each region).
"""
function GroupsOfNodesOf(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "GroupsOfNodesOf", kwargs...)
end

"""
	GroupsOfEdgesOf(group::Group, args...; InSupport=nothing, kwargs...)

Groups of edges of elements of R1 (a group is associated with each region).
< InSupport: in a support R2 being a group of type ElementOf, i.e., containing elements >.
"""
function GroupsOfEdgesOf(group::Group, args...; InSupport = nothing, kwargs...)
	return add!(group, args...; gtype = "GroupsOfEdgesOf", InSupport = InSupport, kwargs...)
end

"""
	GroupsOfEdgesOnNodesOf(group::Group, args...; Not=nothing, kwargs...)

Groups of edges incident to nodes of elements of R1 (a group is associated with each node).
< Not: but not those of R2) >.
"""
function GroupsOfEdgesOnNodesOf(group::Group, args...; Not = nothing, kwargs...)
	return add!(group, args...; gtype = "GroupsOfEdgesOnNodesOf", Not = Not, kwargs...)
end

"""
	GroupOfRegionsOf(group::Group, args...; kwargs...)

Single group of elements of regions in R1
(with basis function BF_Region just one DOF is created for all elements of R1).
"""
function GroupOfRegionsOf(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "GroupOfRegionsOf", kwargs...)
end

"""
	EdgesOfTreeIn(group::Group, args...; StartingOn=nothing, kwargs...)

Edges of a tree of edges of R1
< StartingOn: a complete tree is first built on R2 >.
"""
function EdgesOfTreeIn(group::Group, args...; StartingOn = nothing, kwargs...)
	return add!(group, args...; gtype = "EdgesOfTreeIn", StartingOn = StartingOn, kwargs...)
end

"""
	FacetsOfTreeIn(group::Group, args...; StartingOn=nothing, kwargs...)

Facets of a tree of facets of R1
< StartingOn: a complete tree is first built on R2 >.
"""
function FacetsOfTreeIn(group::Group, args...; StartingOn = nothing, kwargs...)
	return add!(
		group,
		args...;
		gtype = "FacetsOfTreeIn",
		StartingOn = StartingOn,
		kwargs...,
	)
end

"""
	DualNodesOf(group::Group, args...; kwargs...)

Dual nodes of elements of R1.
"""
function DualNodesOf(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "DualNodesOf", kwargs...)
end

"""
	DualEdgesOf(group::Group, args...; kwargs...)

Dual edges of elements of R1.
"""
function DualEdgesOf(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "DualEdgesOf", kwargs...)
end

"""
	DualFacetsOf(group::Group, args...; kwargs...)

Dual facets of elements of R1.
"""
function DualFacetsOf(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "DualFacetsOf", kwargs...)
end

"""
	DualVolumesOf(group::Group, args...; kwargs...)

Dual volumes of elements of R1.
"""
function DualVolumesOf(group::Group, args...; kwargs...)
	return add!(group, args...; gtype = "DualVolumesOf", kwargs...)
end
