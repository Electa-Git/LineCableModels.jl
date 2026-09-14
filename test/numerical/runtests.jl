include("references.jl")

isempty(ARGS) || throw(ArgumentError("numerical reference checks read the reviewed manifest; no live-run options are accepted"))
rows = NumericalReferences.check()
for row in rows
    println(row.id, '\t', row.quantity, '[', row.row, ',', row.column, ']',
        '\t', row.passed ? "pass" : "FAIL",
        "\tabs_rms=", row.absolute, "\trel_rms=", row.relative)
end
all(row -> row.passed, rows) || error("numerical-reference regression; stored arrays were not updated")
