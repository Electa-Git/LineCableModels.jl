@testitem "ModalAnalysis / armored study fixture declares nine independent terminals" tags=[:unit] setup=[ModalStudyFixtures] begin
    using LineCableModels
    for layout in (:trefoil,:horizontal)
        fixture=ModalStudyFixtures.study_problem(layout;frequencies=[50.0])
        @test fixture.design.terminal_order==[:core,:screen,:armor]
        @test length(fixture.system.terminal_order)==9
        @test sort(collect(v for mapping in fixture.connections for v in values(mapping)))==collect(1:9)
        @test line_length(fixture.problem.system)==600.0
        @test length(fixture.problem.frequencies)==1
    end
end
