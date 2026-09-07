@testitem "Engine / literature assimilation / all matched source tables" begin
    using Base.Docs: meta
    using LineCableModels
    root=pkgdir(LineCableModels)
    survey=joinpath(root,"docs","theory")
    records=readlines(joinpath(survey,"assimilation.tsv"))
    @test first(records) == "Record\tScope\tImplementation\tModule\tIdentifier\tFormula source"
    @test length(records)==145
    seen=Set{String}()
    coverage=Dict{String,Set{Symbol}}()
    indexed=Set{String}()
    for page in ("internal_impedance.md","insulation_parameters.md",
                 "earth_return_impedance.md","earth_return_admittance.md")
        content=read(joinpath(survey,page),String)
        for target in eachmatch(r"\]\(([^)#]+\.md)#identification-and-source\)",content)
            @test target[1] ∉ indexed
            push!(indexed,target[1])
        end
    end
    for line in Iterators.drop(records,1)
        fields=split(line,'\t';keepempty=true)
        @test length(fields)==6
        path,scope,status,modules,ids,files=fields
        @test path ∉ seen
        push!(seen,path)
        @test isfile(joinpath(survey,path))
        @test scope in ("Candidate","Deferred")
        @test status in ("Pending","Deferred","Existing","Implemented","Verified","Equivalent")
        @test (scope=="Deferred") == (status=="Deferred")
        @test path in indexed
        status in ("Existing","Implemented","Verified","Equivalent") || continue
        source=read(joinpath(survey,path),String)
        frontmatter=match(r"(?s)(## Identification and source\n.*?)(?=\n\n\*\*Description\.\*\*)",source)[1]
        for (module_name,id,file) in zip(split(modules,';'),split(ids,';'),split(files,';'))
            owner=getproperty(LineCableModels.Engine,Symbol(module_name))
            identifier=Symbol(id)
            push!(get!(Set{Symbol},coverage,module_name),identifier)
            @test identifier in owner.formulas()
            descriptions=String[]
            for (binding,multidoc) in meta(owner)
                binding.var===:description || continue
                for (typesig,docstring) in multidoc.docs
                    first(Base.unwrap_unionall(only(typesig.parameters)).parameters)===identifier || continue
                    push!(descriptions,join(filter(x->x isa String,collect(docstring.text))))
                    @test normpath(String(docstring.data[:path]))==joinpath(root,file)
                end
            end
            @test length(descriptions)==1
            if status == "Equivalent"
                @test occursin(path, read(joinpath(survey,"deduplication.md"),String))
            else
                @test occursin(frontmatter,only(descriptions))
            end
        end
    end
    for (name,ids) in coverage
        owner=getproperty(LineCableModels.Engine,Symbol(name))
        unmatched=setdiff(Set(owner.formulas()),ids)
        expected=name=="EarthAdmittance" ? Set((:IdealGround,)) : Set{Symbol}()
        @test unmatched == expected
    end
    @test indexed == seen
end
