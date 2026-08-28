using Test
using LineCableModels
using Bonito
using WGLMakie
using CairoMakie

CairoMakie.activate!()
include(joinpath(@__DIR__, "..", "app.jl"))

@testset "Bonito cable showcase" begin
    @testset "domain construction" begin
        cable = design(12.5, 8.0)
        values = geometry(cable)

        @test cable isa CableDesign
        @test cable.root isa Stack
        @test length(cable.root.items) == 2
        @test cable.root.items[1] isa Group
        @test cable.root.items[2] isa Region
        @test length(cable.geometry.regions) == 2
        @test cable.terminal_order == [:core]
        @test values.core_r_in_m == 0.0
        @test values.core_r_ex_m ≈ 12.5e-3
        @test values.insulation_r_in_m ≈ values.core_r_ex_m
        @test values.insulation_r_ex_m ≈ 20.5e-3
        @test values.core_area_m2 > 0
        @test values.core_area_m2 ≈ π * values.core_r_ex_m^2

        wider_core = design(15.0, 8.0)
        thicker_insulation = design(12.5, 10.0)
        @test wider_core !== cable
        @test thicker_insulation !== cable
        @test geometry(wider_core).core_r_ex_m ≈ 15.0e-3
        @test geometry(thicker_insulation).insulation_r_ex_m ≈ 22.5e-3

        for invalid in (0.0, -1.0, Inf, NaN)
            @test_throws DomainError design(invalid, 8.0)
            @test_throws DomainError design(12.5, invalid)
        end
    end

    @testset "plot construction" begin
        state = Observable(design(12.5, 8.0))
        plot = cable_figure(state)
        figure = plot.figure
        insulation_plot = plot.insulation_plot
        core_plot = plot.core_plot

        @test figure isa Figure
        @test length(plot.axis.scene.plots) == 2
        @test plot.insulation_circle[].r ≈ 20.5f0
        @test plot.core_circle[].r ≈ 12.5f0

        state[] = design(20.0, 3.0)
        @test plot.figure === figure
        @test plot.insulation_plot === insulation_plot
        @test plot.core_plot === core_plot
        @test plot.insulation_circle[].r ≈ 23.0f0
        @test plot.core_circle[].r ≈ 20.0f0

        mktempdir() do directory
            output = joinpath(directory, "cable.svg")
            Makie.save(output, figure)
            @test isfile(output)
            @test filesize(output) > 0
        end
    end

    @testset "application DOM" begin
        app = cable_app()
        @test app isa App

        theme = read(joinpath(@__DIR__, "..", "assets", "theme.css"), String)
        @test occursin("width: 100vw;", theme)
        @test occursin("height: 100vh;", theme)

        parent = Session(NoConnection(); asset_server = NoServer())
        io = IOBuffer()
        rendered = Bonito.show_html(io, app; parent)
        html = String(take!(io))
        try
            for marker in (
                "linecable-app",
                "linecable-slides",
                "core-radius-control",
                "insulation-thickness-control",
                "cable-plot",
                "reveal lc-deck"
            )
                @test occursin(marker, html)
            end
        finally
            close(rendered)
            close(parent)
        end
    end
end
