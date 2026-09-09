@testset "SQL snapshots keep one concrete shape across cold query transitions" begin
    db = RT.SQLite.DB()
    saved = nothing
    try
        empty = RT.sql_rows(db, "SELECT 1 AS count WHERE 0")
        populated = RT.sql_rows(db, "SELECT 1 AS count")
        nullable = RT.sql_rows(db, "SELECT NULL AS value UNION ALL SELECT 'present' AS value")
        saved = only(RT.sql_rows(db, "SELECT 'retained' AS columns, x'0102' AS bytes, 2.5 AS quantity"))
        @test typeof(empty) === typeof(populated) === typeof(nullable) === Vector{RT.SQLRow}
        @test isempty(empty) && only(populated).count == 1
        @test ismissing(nullable[1].value) && nullable[2].value == "present"
        @test typeof(nullable[1]) === typeof(nullable[2]) === typeof(saved)
        @test Set(propertynames(saved)) == Set((:columns, :bytes, :quantity))
        @test hasproperty(saved, :quantity) && !hasproperty(saved, :absent)
        @test_throws KeyError saved.absent
        @test_throws Exception RT.sql_rows(db, "SELECT missing_column")
        @test only(RT.sql_rows(db, "SELECT 3 AS count")).count == 3
    finally
        RT.DBInterface.close!(db)
    end
    @test saved.columns == "retained" # A column may share the storage-field name.
    @test saved.bytes == UInt8[1, 2] && saved.quantity == 2.5
end
