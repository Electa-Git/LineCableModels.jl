using DocStringExtensions: DocStringExtensions, SIGNATURES, TYPEDSIGNATURES, TYPEDEF,
                           TYPEDFIELDS, FIELDS, FUNCTIONNAME, IMPORTS, EXPORTS

export SIGNATURES, TYPEDSIGNATURES, TYPEDEF, TYPEDFIELDS, FIELDS, FUNCTIONNAME,
       METHODLIST, IMPORTS, EXPORTS

#! explicit-imports: off
# METHODLIST adapts dependency-internal DocStringExtensions machinery. Keep this
# exception local to the adapter so the rest of the package remains strict.
struct _CleanMethodList <: DocStringExtensions.Abbreviation end

"`METHODLIST` abbreviation with package-local, CI-safe source paths."
const METHODLIST = _CleanMethodList()
const _PACKAGE_ROOT = normpath(dirname(@__DIR__))

function _method_path(method)
    file = string(method.file)
    isempty(file) && return "unknown"
    normalized = normpath(file)
    return startswith(normalized, _PACKAGE_ROOT) ?
           relpath(normalized, _PACKAGE_ROOT) : basename(normalized)
end

function DocStringExtensions.format(::_CleanMethodList, buffer, doc)
    binding = doc.data[:binding]
    typesig = doc.data[:typesig]
    module_name = doc.data[:module]
    function_value = Docs.resolve(binding)
    groups = DocStringExtensions.methodgroups(
        function_value,
        typesig,
        module_name;
        exact = false
    )
    isempty(groups) && return nothing

    println(buffer)
    for group in groups
        println(buffer, "```julia")
        for method in group
            DocStringExtensions.printmethod(buffer, binding, function_value, method)
            println(buffer)
        end
        println(buffer, "```\n")
        if !isempty(group)
            method = first(group)
            url = DocStringExtensions.url(method)
            if isempty(url) || startswith(url, "file:")
                println(
                    buffer,
                    "defined at `$(_method_path(method)):$(method.line)`."
                )
            else
                println(
                    buffer,
                    "defined at [`$(_method_path(method)):$(method.line)`]($url)."
                )
            end
        end
        println(buffer)
    end
    println(buffer)
    return nothing
end
#! explicit-imports: on
