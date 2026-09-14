using Documenter
using LineCableModels

DocMeta.setdocmeta!(
    LineCableModels,
    :DocTestSetup,
    quote
        using LineCableModels
        using LineCableModels.DataModel.BaseParams
    end;
    recursive = true
)
doctest(LineCableModels)
