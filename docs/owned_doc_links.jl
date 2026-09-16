module OwnedDocLinks

import Documenter

const MarkdownAST = Documenter.MarkdownAST

export OwnedDocLinker

"""
    OwnedDocLinker(roots...)

Turn unlinked inline-code names owned by `roots` into ordinary Documenter `@ref`
links. Only bare or module-qualified names that Documenter has included in the
current build are linked. Existing links, code snippets, values, missing names,
and ambiguous names are left unchanged.
"""
mutable struct OwnedDocLinker <: Documenter.Plugin
    roots::Vector{Module}
    linked::Int

    function OwnedDocLinker(roots::Module...)
        isempty(roots) && throw(ArgumentError("at least one owning module is required"))
        return new(collect(roots), 0)
    end
end

abstract type LinkOwnedDocNames <: Documenter.Builder.DocumentPipeline end

# Template expansion registers the docstring inventory. CrossReferences then
# resolves the @ref nodes created here using Documenter's normal machinery.
Documenter.Selectors.order(::Type{LinkOwnedDocNames}) = 2.9

function Documenter.Selectors.runner(
        ::Type{LinkOwnedDocNames}, doc::Documenter.Document)
    linker = get(doc.plugins, OwnedDocLinker, nothing)
    linker === nothing && return
    Documenter.is_doctest_only(doc, "LinkOwnedDocNames") && return

    linker.linked = 0
    for page in values(doc.blueprint.pages)
        meta = copy(doc.user.meta)
        link_owned_names!(page.mdast, meta, doc, linker)
    end
    @info "LinkOwnedDocNames: linked $(linker.linked) owned inline-code references."
    return
end

function link_owned_names!(node, meta, doc, linker)
    element = node.element

    if element isa Documenter.MetaNode
        merge!(meta, element.dict)
    elseif element isa Documenter.DocsNode
        for (docstring, docmeta) in zip(element.mdasts, element.metas)
            docstring_meta = copy(meta)
            module_ = get(docmeta, :module, nothing)
            isnothing(module_) || (docstring_meta[:CurrentModule] = module_)
            link_owned_names!(docstring, docstring_meta, doc, linker)
        end
        return
    elseif element isa Union{MarkdownAST.Link, MarkdownAST.Image}
        # Existing links and image descriptions must never acquire nested links.
        return
    elseif element isa MarkdownAST.Code
        code = element.code
        target = documented_owned_target(code, meta, doc, linker.roots)
        if target !== nothing
            node.element = MarkdownAST.Link("@ref $target", "")
            push!(node.children, MarkdownAST.Node(element))
            linker.linked += 1
        end
        return
    end

    for child in node.children
        link_owned_names!(child, meta, doc, linker)
    end
    return
end

function is_name(code::AbstractString)
    expression = try
        Meta.parse(code)
    catch
        return false
    end
    return is_name(expression)
end

is_name(::Symbol) = true
function is_name(expression::Expr)
    Meta.isexpr(expression, :., 2) || return false
    field = expression.args[2]
    return is_name(expression.args[1]) &&
           field isa QuoteNode && field.value isa Symbol
end
is_name(_) = false

function documented_owned_target(code, meta, doc, roots)
    is_name(code) || return nothing
    expression = Meta.parse(code)
    current_module = get(meta, :CurrentModule, Main)
    modules = current_module === Main ? (Main,) : (current_module, Main)

    for module_ in (modules..., roots...)
        binding = try
            Documenter.DocSystem.binding(module_, expression)
        catch
            continue
        end
        object = Documenter.find_object(doc, binding, Union{})
        object === nothing && continue
        is_owned(object.binding, roots) || continue
        return Documenter.bindingstring(object.binding)
    end

    # Public-but-unexported names are common in this package. A bare name is
    # still safe to link when the generated inventory has exactly one owned
    # binding with that name; duplicate names remain deliberately unlinked.
    expression isa Symbol || return nothing
    targets = Set{String}()
    for binding in keys(doc.internal.bindings)
        binding.var === expression || continue
        object = Documenter.find_object(doc, binding, Union{})
        object === nothing && continue
        is_owned(object.binding, roots) || continue
        push!(targets, Documenter.bindingstring(object.binding))
    end
    return length(targets) == 1 ? only(targets) : nothing
end

function is_owned(binding, roots)
    any(root -> is_descendant(binding.mod, root), roots) && return true

    # A module's binding belongs to its parent, so handle an owning root module
    # (and its documented submodules) explicitly.
    isdefined(binding.mod, binding.var) || return false
    value = getfield(binding.mod, binding.var)
    return value isa Module && any(root -> is_descendant(value, root), roots)
end

function is_descendant(module_::Module, root::Module)
    while true
        module_ === root && return true
        parent = parentmodule(module_)
        parent === module_ && return false
        module_ = parent
    end
end

end
