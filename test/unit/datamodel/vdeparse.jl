@testitem "DataModel / VDE parser / ordered tokens and residue" tags=[:unit, :datamodel] setup=[
    UseDataModelSupport
] begin
    parser=LineCableModels.DataModel.vdeparse

    cases=(
        (
            "2XS(F)2Y 3x240/25 12/20kV RM/V",
            Dict(
                :conductor=>"copper conductor",
                :insulation=>"cross-linked PE (XLPE)",
                :waterblocking=>"longitudinally water-proof protection",
                :outer_sheath=>"PE outer sheath",
                :cores=>"3",
                :conductor_cross_section=>"240",
                :metallic_screen_cross_section=>"25",
                :voltage=>"12/20 kV",
                :conductor_type=>"round, stranded, compact"
            )
        ),
        (
            "A2XS(FL)KL2Y 1x630 18/30kV RE",
            Dict(
                :conductor=>"aluminium conductor",
                :sheath=>"aluminium sheath",
                :cores=>"1",
                :conductor_type=>"round, solid"
            )
        ),
        ("N2XSY 1x150", Dict(:designation=>"DIN VDE standard", :cores=>"1"))
    )

    for (code, expected) in cases
        parsed=parser(code)
        @test all(key -> parsed[key] == expected[key], keys(expected))
    end

    @test isempty(parser(" \u00a0 "))
    @test LineCableModels.DataModel.decode_type("RRM") == "round, stranded"
    @test_logs (:warn, r"Unknown conductor type") begin
        @test LineCableModels.DataModel.decode_type("RZ") == "round"
    end
    @test_logs (:warn, r"Unparsed stub residue") begin
        parsed = parser("2XSQ 1x10")
        @test parsed[:unparsed_stub] == "Q"
    end
    @test_logs (:warn, r"Unparsed trailing") begin
        parsed = parser("2XSY 1x10 unsupported")
        @test parsed[:unparsed] == "unsupported"
    end
end
