@testitem "Utils / legacy macros / parameterization and numeric promotion" tags = [:unit] begin
    macro_source = joinpath(pkgdir(LineCableModels), "src", "utils", "macros.jl")
    sandbox_name = gensym(:LegacyMacroSandbox)
    sandbox = Module(sandbox_name)
    Core.eval(sandbox, :(using MacroTools, DocStringExtensions))
    Base.include(sandbox, macro_source)
    Core.eval(
        sandbox,
        quote
            const SUPPORTED = Union{Float32, Float64}

            _coerce_args_to_T(values...) = promote_type(map(typeof, values)...)
            _coerce_scalar_to_T(value, ::Type{T}) where {T} = convert(T, value)
            _coerce_array_to_T(values, ::Type{T}) where {T} = T.(values)

            struct EarthModel{T}
                value::T
            end
            Base.convert(::Type{EarthModel{T}}, model::EarthModel) where {T} = EarthModel{T}(convert(
                T, model.value))
        end
    )

    parameterized = Core.eval(
        sandbox,
        :(@parameterize(AbstractArray{_, 2}, SUPPORTED))
    )
    @test parameterized == Union{AbstractMatrix{Float32}, AbstractMatrix{Float64}}
    @test_throws LoadError Core.eval(sandbox, :(@parameterize(Vector{_}, Float64)))
    @test_throws LoadError Core.eval(sandbox, :(@parameterize(Vector{_}, NOT_DEFINED)))

    Core.eval(
        sandbox,
        quote
            @measurify function promoted_sum(
                    left::T,
                    right::T,
                    values::Vector{T};
                    scale::T = one(T)
            ) where {T <: AbstractFloat}
                return (
                    typeof(left),
                    typeof(right),
                    eltype(values),
                    typeof(scale),
                    left + right + sum(values) * scale
                )
            end

            @measurify function anchored_sum(
                    model::EarthModel{T},
                    offset::T
            ) where {T <: AbstractFloat}
                return (typeof(model), typeof(offset), model.value + offset)
            end
        end
    )

    promoted = sandbox.promoted_sum(
        Float32(1),
        2.0,
        Float32[3, 4];
        scale = 2.0
    )
    @test promoted[1:4] == (Float64, Float64, Float64, Float64)
    @test promoted[5] == 17.0

    same_precision = sandbox.promoted_sum(
        Float32(1),
        Float32(2),
        Float32[3];
        scale = Float32(2)
    )
    @test same_precision[1:4] == (Float32, Float32, Float32, Float32)
    @test same_precision[5] == Float32(9)

    anchored = sandbox.anchored_sum(sandbox.EarthModel{Float32}(1), 2.0)
    @test anchored == (sandbox.EarthModel{Float64}, Float64, 3.0)
    @test_throws LoadError Core.eval(sandbox, :(@measurify 1 + 2))
end

@testitem "Utils / legacy macros / public-name export policy" tags = [:unit] begin
    macro_source = joinpath(pkgdir(LineCableModels), "src", "utils", "macros.jl")
    sandbox_name = gensym(:AutoExportSandbox)
    sandbox = Module(sandbox_name)
    Core.eval(sandbox, :(using MacroTools, DocStringExtensions))
    Base.include(sandbox, macro_source)
    Core.eval(
        sandbox,
        quote
            public_function() = :public
            _private_function() = :private
            struct PublicType end
            struct _PrivateType end
        end
    )
    Core.eval(sandbox, :(@autoexport))

    exported = names(sandbox)
    @test :public_function in exported
    @test :PublicType in exported
    @test !(:_private_function in exported)
    @test !(:_PrivateType in exported)
    @test sandbox.public_function() === :public
end
