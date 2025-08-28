using DocStringExtensions, Pkg
"""
Override `DocStringExtensions.format` for `_CLEANMETHODLIST`.
"""
struct _CleanMethodList <: DocStringExtensions.Abbreviation end
"Modified `_CLEANMETHODLIST` abbreviation with sanitized file paths."
const _CLEANMETHODLIST = _CleanMethodList()
function DocStringExtensions.format(::_CleanMethodList, buf, doc)
    local binding = doc.data[:binding]
    local typesig = doc.data[:typesig]
    local modname = doc.data[:module]
    local func = Docs.resolve(binding)
    local groups = DocStringExtensions.methodgroups(func, typesig, modname; exact=false)
    if !isempty(groups)
        println(buf)
        local pkg_root = Pkg.pkgdir(modname) # Use Pkg.pkgdir here
        if pkg_root === nothing
            @warn "Could not determine package root for module $modname using _CLEANMETHODLIST. Paths will be shown as basenames."
        end
        for group in groups
            println(buf, "```julia")
            for method in group
                DocStringExtensions.printmethod(buf, binding, func, method)
                println(buf)
            end
            println(buf, "```\n")
            if !isempty(group)
                local method = group[1]
                local file = string(method.file)
                local line = method.line
                local path =
                    if pkg_root !== nothing && !isempty(file) &&
                       startswith(file, pkg_root)
                        basename(file) # relpath(file, pkg_root)
                    # elseif !isempty(file) && isfile(file)
                    # 	basename(file)
                    else
                        string(method.file) # Fallback
                    end
                local URL = DocStringExtensions.url(method)
                isempty(URL) || println(buf, "defined at [`$path:$line`]($URL).")
            end
            println(buf)
        end
        println(buf)
    end
    return nothing
end