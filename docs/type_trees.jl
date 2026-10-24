module DocumentationTrees

using AbstractTrees
using InteractiveUtils: subtypes

struct TypeTree
    value::Type
end

function AbstractTrees.children(node::TypeTree)
    descendants = sort!(subtypes(node.value); by = string)
    return TypeTree.(descendants)
end

function AbstractTrees.printnode(io::IO, node::TypeTree; kwargs...)
    name = replace(string(node.value), "LineCableModels." => "")
    return print(io, name)
end

function type_tree(root::Type)
    io = IOBuffer()
    AbstractTrees.print_tree(io, TypeTree(root); maxdepth = typemax(Int))
    return String(take!(io))
end

end
