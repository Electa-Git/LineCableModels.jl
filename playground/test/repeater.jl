using Bonito
using UUIDs

function test_repeat_entry(values=(
        name="Cable 01",
        enabled=true,
        model="Solid",
        taps=[1.0, 2.0],
    ))
    fields = (
        name=Observable(values.name),
        enabled=Observable(values.enabled),
        model=Observable(values.model),
        taps=Observable(copy(values.taps)),
    )
    return RepeatEntry(fields, DOM.div(values.name))
end

@testset "generic repeater" begin
    first_entry = test_repeat_entry()
    repeater = Repeater(
        test_repeat_entry;
        initial=[first_entry],
        duplicate=test_repeat_entry,
        min_entries=1,
        max_entries=3,
        label="Conductors"
    )

    @test length(repeater) == 1
    @test snapshot(repeater) == [(
        name="Cable 01",
        enabled=true,
        model="Solid",
        taps=[1.0, 2.0],
    )]

    first_entry.fields.name[] = "Cable A"
    first_entry.fields.enabled[] = false
    @test snapshot(first_entry).name == "Cable A"
    @test snapshot(first_entry).enabled == false

    second_entry = duplicate!(repeater, first_entry.id)
    @test length(repeater) == 2
    @test second_entry.id != first_entry.id
    @test snapshot(second_entry) == snapshot(first_entry)
    second_entry.fields.name[] = "Cable B"
    push!(second_entry.fields.taps[], 3.0)
    @test first_entry.fields.name[] == "Cable A"
    @test first_entry.fields.taps[] == [1.0, 2.0]

    third_entry = add!(repeater)
    @test length(repeater) == 3
    @test allunique(entry.id for entry in repeater)
    @test_throws ArgumentError add!(repeater)
    @test_throws ArgumentError duplicate!(repeater, first_entry.id)

    move!(repeater, third_entry.id, 1)
    @test repeater[1] === third_entry
    @test snapshot(repeater)[2].name == "Cable A"
    removed = remove!(repeater, second_entry.id)
    @test removed === second_entry
    @test length(repeater) == 2

    remove!(repeater, first_entry.id)
    @test length(repeater) == 1
    @test_throws ArgumentError remove!(repeater, third_entry.id)
    @test_throws KeyError move!(repeater, uuid4(), 1)
    @test_throws BoundsError move!(repeater, third_entry.id, 2)
end

@testset "repeater validation and extraction" begin
    @test_throws ArgumentError Repeater(test_repeat_entry; initial=-1)
    @test_throws ArgumentError Repeater(test_repeat_entry; initial=0, min_entries=1)
    @test_throws ArgumentError Repeater(
        test_repeat_entry;
        initial=1,
        min_entries=2,
        max_entries=1
    )
    entry = test_repeat_entry()
    @test_throws ArgumentError Repeater(
        test_repeat_entry;
        initial=[entry, entry]
    )

    custom = RepeatEntry(
        (real=Observable(3.0), ignored=Observable("metadata")),
        DOM.div();
        extract=fields -> fields.real[]^2
    )
    @test snapshot(custom) == 9.0

    without_duplication = Repeater(test_repeat_entry; initial=1)
    @test_throws ArgumentError duplicate!(without_duplication, without_duplication[1].id)
end
