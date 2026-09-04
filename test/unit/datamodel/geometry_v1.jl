@testitem "DataModel / v1 geometry / primitive and resolved contracts" tags=[:unit] begin
    const DM=LineCableModels.DataModel

    @test DM.Disk(2) isa DM.AbstractPrimitive{Float64}
    @test DM.Rectangle(4.0f0, 2) isa DM.Rectangle{Float32}
    @test DM.Ellipse(3, 2) isa DM.Ellipse{Float64}
    @test DM.Sector(span = π / 2, r_base = 0.2, r_back = 2, fillet = 0.1) isa DM.Sector
    @test DM.Annulus(1, 2) isa DM.Annulus
    @test DM.Shell(0.5) isa DM.Shell
    @test_throws DomainError DM.Disk(0)
    @test_throws DomainError DM.Rectangle(-1, 2)
    @test_throws DomainError DM.Ellipse(1, Inf)
    @test_throws DomainError DM.Sector(span = 0, r_base = 0.2, r_back = 2)
    @test_throws DomainError DM.Sector(span = 1, r_base = -1, r_back = 2)
    @test_throws DomainError DM.Annulus(2, 1)
    @test_throws DomainError DM.Shell(0)

    disk=@inferred DM.resolve(DM.EmptyBoundary(), DM.Disk(2.0))
    @test disk isa DM.Disk
    @test DM.boundary(disk) === disk
    @test DM.area(disk) ≈ 4π
    @test DM.centroid(disk) == (0.0, 0.0)
    @test DM.support(disk, 0.37) == 2.0
    @test DM.r_in(disk) == 0.0
    @test DM.r_ex(disk) == 2.0
    @test DM.thickness(disk) == 2.0

    annulus=@inferred DM.resolve(disk, DM.Shell(0.5))
    @test annulus == DM.Annulus(2.0, 2.5)
    @test DM.boundary(annulus) == DM.Disk(2.5)
    @test DM.area(annulus) ≈ π * (2.5^2 - 2.0^2)
    @test DM.r_in(annulus) == 2.0
    @test DM.r_ex(annulus) == 2.5
    @test DM.thickness(annulus) == 0.5
    @test_throws MethodError DM.resolve(
        DM.resolve(DM.EmptyBoundary(), DM.Rectangle(2, 1)),
        DM.Shell(0.1)
    )

    rectangle=DM.resolve(DM.EmptyBoundary(), DM.Rectangle(4, 2))
    @test DM.area(rectangle) == 8.0
    @test DM.centroid(rectangle) == (0.0, 0.0)
    @test DM.support(rectangle, 0) == 2.0
    @test DM.support(rectangle, π / 2) ≈ 1.0
    @test !hasmethod(DM.r_in, Tuple{typeof(rectangle)})
    @test_throws MethodError DM.r_in(rectangle)

    ellipse=DM.resolve(DM.EmptyBoundary(), DM.Ellipse(3, 2))
    @test DM.area(ellipse) ≈ 6π
    @test DM.perimeter(ellipse) ≈ 15.865439589290588
    @test DM.support(ellipse, 0) == 3.0
    @test DM.support(ellipse, π / 2) ≈ 2.0
    ellipse_shell=DM.resolve(ellipse, DM.Shell(0.5))
    @test ellipse_shell isa DM.ShellShape
    @test ellipse_shell.outer isa DM.EllipseOffset
    @test DM.support(ellipse_shell.outer, 0) ≈ 3.5
    @test DM.support(ellipse_shell.outer, π / 2) ≈ 2.5
    @test DM.area(ellipse_shell.outer) ≈
          DM.area(ellipse) + DM.perimeter(ellipse) * 0.5 + π * 0.5^2
    @test DM.perimeter(ellipse_shell.outer) ≈ DM.perimeter(ellipse) + π
    @test DM.centroid(ellipse_shell.outer) == DM.centroid(ellipse)
    @test length(DM.tessellate(ellipse_shell.outer; points_per_arc = 64)) == 64

    sector=DM.resolve(DM.EmptyBoundary(),
        DM.Sector(span = π / 2, r_base = 0.2, r_back = 2, fillet = 0.1))
    @test DM.area(sector) > 0
    @test DM.r_in(sector) == 0.2
    @test DM.r_ex(sector) == 2.0
    @test DM.support(sector, 0) ≈ 2.0
    @test DM.centroid(sector)[1] > 0
    @test abs(DM.centroid(sector)[2]) <= 16eps(Float64)

    polygon=DM.Polygon(((0, 0), (2, 0), (2, 1), (0, 1)))
    resolved_polygon=DM.resolve(DM.EmptyBoundary(), polygon)
    @test DM.area(resolved_polygon) == 2.0
    @test DM.centroid(resolved_polygon) == (1.0, 0.5)
    @test DM.support(resolved_polygon, 0) == 2.0
    @test_throws DomainError DM.Polygon(((0, 0), (1, 0), (2, 0)))

    converted=convert(
        DM.AbstractPrimitive{Float32},
        DM.Rectangle(4.0, 2.0)
    )
    @test converted isa DM.Rectangle{Float32}
    @test converted.w === 4.0f0
end

@testitem "DataModel / v1 geometry / pose and Gridspace materialization" tags=[:unit] begin
    using Random
    const DM=LineCableModels.DataModel
    const PB=LineCableModels.ParametricBuilder

    parent=DM.Pose2(1, 2, π / 2)
    child=DM.Pose2(3, 0, π / 2)
    composed=parent * child
    @test composed.x ≈ 1.0
    @test composed.y ≈ 5.0
    @test composed.φ ≈ π

    placed=DM.resolve(parent, DM.Rectangle(4, 2))
    @test DM.area(placed) == 8.0
    @test DM.centroid(placed) == (1.0, 2.0)
    @test DM.support(placed, 0) ≈ 2.0

    poses=PB.Pose2(x = PB.Grid((0.0, 1.0)), y = 2.0, φ = PB.Grid((0.0, π / 2)))
    @test poses isa PB.Gridspace{DM.Pose2}
    @test length(poses) == 4
    @test all(pose -> pose isa DM.Pose2, poses)

    disks=PB.Disk(PB.Grid((1.0f0, 2.0f0)))
    @test disks isa PB.Gridspace{DM.Disk}
    @test collect(disks) == collect(DM.Disk.((1.0f0, 2.0f0)))
    @test rand(MersenneTwister(7), disks) in collect(disks)

    rectangles=PB.Rectangle(PB.Grid((1.0, 2.0)), 3.0)
    @test rectangles isa PB.Gridspace{DM.Rectangle}
    @test collect(rectangles) ==
          [DM.Rectangle(1.0, 3.0), DM.Rectangle(2.0, 3.0)]
end
