@testitem "ObservedResult / retained re-expression preserves scientific records" tags=[:unit] begin
    using Measurements
    using LineCableModels.Grammar: observation_product,gridpoint_id
    using LineCableModels.Engine: retain_gridpoint,compare
    U=LineCableModels.Units
    shared=measurement(1.,0.1)
    line=LineParameters(fill(complex(shared,2shared),2,2,3),fill(3.0+4im,2,2,3),[10.,100.,1000.])
    candidate=retain_gridpoint(line,gridpoint_id())
    reference=retain_gridpoint(line,gridpoint_id())
    comparisons=compare(reference,candidate,[R];bands=(:all,))
    timings=(candidate_id=LineCableModels.Grammar.observation_gridpoint(candidate).id,seconds=0.5)
    observed=ObservedResult(candidate;comparisons,timings)
    repeated=observables(observed)
    @test isequal(repeated.quantities,observed.quantities)
    @test repeated.quantities[1].values!==observed.quantities[1].values
    @test LineCableModels.Grammar.observation_quantity(observed,R)==first(observed.quantities)
    @test repeated.errors==observed.errors
    @test repeated.timings==observed.timings
    @test only(observables([observed])).timings==timings
    @test only(observables([observed])).errors==comparisons
    selected=ObservedResult(observed,((R,[2,1],[2,1],[3,1]),);length_unit=:base,frequency_unit=:kilo)
    product=only(selected.quantities)
    @test nominal.(product.values)==ones(2,2,2)
    prefix_only=ObservedResult(selected;quantity_units=:milli)
    @test only(prefix_only.quantities).unit==U.units(:milli,:ohm;per=(:base,:meter))
    length_only=ObservedResult(prefix_only;length_unit=:kilo)
    @test only(length_only.quantities).unit==U.units(:milli,:ohm;per=(:kilo,:meter))
    @test product.coordinates.rows==[2,1]
    @test product.coordinates.columns==[2,1]
    @test product.coordinates.samples==[3,1]
    @test product.coordinates.frequencies==[1.,0.01]
    restored=ObservedResult(selected;length_unit=:kilo,frequency_unit=:base)
    original=observation_product(observed,(R,[2,1],[2,1],[3,1]))
    @test only(restored.quantities).values≈original.values
    @test only(restored.quantities).coordinates==original.coordinates
    @test only(restored.quantities).thresholds.values≈original.thresholds.values
    @test only(restored.quantities).thresholds.unit==original.thresholds.unit
    for field in (:available,:engineering_zero,:clipped,:missing_reason)
        @test isequal(getproperty(only(restored.quantities),field),getproperty(original,field))
    end
    @test restored.errors==comparisons && restored.timings==timings
    @test uncertainty(only(restored.quantities).values[1]-observe(observed,R)[1])==0
    @test uncertainty(only(restored.quantities).values[1])>0
    @test_throws ArgumentError ObservedResult(observed;atol=0)
    @test_throws ArgumentError ObservedResult(observed;clip=false)
    @test_throws ArgumentError ObservedResult(observed,(L,))
    @test_throws ArgumentError ObservedResult(observed;frequencies=[1.,2.,3.])
    @test_throws ArgumentError ObservedResult(observed;timings=(seconds=5.,))
    @test_throws ArgumentError ObservedResult(observed,(R,);units=(U.units(:base,:farad),))
    polar=ObservedResult(LineParameters(fill(3.0+4im,1,1,2),fill(1.0+2im,1,1,2),[1.,2.]),((Z,abs),(Z,angle)))
    angle_unit=only(filter(p -> p.quantity==U.quantity(Z,angle),polar.quantities)).unit
    radians=ObservedResult(polar,((Z,angle),);units=(U.units(:base,:radian),))
    degrees=ObservedResult(radians;units=(angle_unit,))
    @test observe(degrees,Z,angle)≈observe(polar,Z,angle)
    @test only(radians.quantities).thresholds==only(degrees.quantities).thresholds
    phase_product=only(radians.quantities)
    malformed=merge(phase_product,(thresholds=merge(phase_product.thresholds,(unit=phase_product.unit,)),))
    @test_throws ArgumentError ObservedResult(radians.gridpoint,[malformed],[],(;))
end

@testitem "ObservedResult / constructor rejects malformed scientific records" tags=[:unit] begin
    using LineCableModels.Grammar: observation_product
    using LineCableModels.ReportBuilder: tabulate
    source=ObservedResult(LineParameters(fill(1.0+2im,2,2,3),fill(3.0+4im,2,2,3),[1.,2.,3.]))
    q=first(source.quantities)
    make(product)=ObservedResult(source.gridpoint,[product],source.errors,source.timings)
    @test_throws ArgumentError make(Base.structdiff(q,(family=nothing,)))
    @test_throws ArgumentError make(merge(q,(family="Z",)))
    @test_throws ArgumentError make(merge(q,(unit=LineCableModels.Units.units(:base,:farad),)))
    @test_throws DimensionMismatch make(merge(q,(values=ones(2,2,2),)))
    @test_throws DimensionMismatch make(merge(q,(values=ones(3,2,2),)))
    @test_throws DimensionMismatch make(merge(q,(available=trues(1,1,3),)))
    @test_throws ArgumentError make(merge(q,(coordinates=merge(q.coordinates,(rows=[1,1],)),)))
    @test_throws ArgumentError make(merge(q,(coordinates=merge(q.coordinates,(rows=[1,3],)),)))
    @test_throws DimensionMismatch make(merge(q,(coordinates=merge(q.coordinates,(frequencies=[1.,2.],)),)))
    @test_throws ArgumentError make(merge(q,(coordinates=merge(q.coordinates,(frequency_unit=q.unit,)),)))
    @test_throws ArgumentError ObservedResult(source.gridpoint,[q,q],[],(;))
    @test_throws ArgumentError make(merge(q,(thresholds=merge(q.thresholds,(values=-1.,)),)))
    @test_throws DimensionMismatch make(merge(q,(thresholds=merge(q.thresholds,(values=[0.,0.],)),)))
    @test_throws ArgumentError make(merge(q,(coordinates=merge(q.coordinates,(indices=([2,1],:,:),)),)))
    encoded=LineCableModels.ImportExport.serialize_value(source)
    restored=LineCableModels.ImportExport.deserialize_value(encoded)
    @test tabulate(restored).Z.R==tabulate(source).Z.R
    # Restoration has the same validation boundary as native construction.
    fields=encoded["payload"]["fields"]["values"]
    record=fields[2]["values"][1]
    coordinate=record["values"][findfirst(==("coordinates"),record["names"])]
    coordinate["values"][findfirst(==("kind"),coordinate["names"])]=LineCableModels.ImportExport.serialize_value(:unknown)
    @test_throws ArgumentError LineCableModels.ImportExport.deserialize_value(encoded)
end

@testitem "ObservedResult / coordinate and unit interpretation is shared" tags=[:unit] begin
    using LineCableModels.Grammar: observation_product
    z=reshape(complex.(1.:12.,101.:112.),2,2,3)
    line=LineParameters(z,2z,[1.,10.,100.])
    a=ObservedResult(line;length_unit=:base)
    b=ObservedResult(line,((R,[2,1],[2,1],[1,3]),(X,[2,1],[2,1],[1,3]),
        (G,[2,1],[2,1],:),(B,[2,1],[2,1],:));length_unit=:kilo,frequency_unit=:kilo)
    products=observation_product((a,b),R)
    @test products[2].unit==products[1].unit
    @test products[2].coordinates.rows==products[2].coordinates.columns==[1,2]
    @test products[2].coordinates.samples==[1,3]
    @test products[2].coordinates.frequencies==[1.,100.]
    @test products[2].values≈products[1].values[:,:,[1,3]]
    @test b.quantities[1].coordinates.rows==[2,1]
    @test b.quantities[1].unit!=a.quantities[1].unit
    incompatible=ObservedResult(b,((R,[1],:,:),))
    @test_throws DimensionMismatch observation_product((a,incompatible),R)
end

@testitem "ObservedResult / table definitions retain one hierarchy" tags=[:unit] begin
    using DataFrames
    using LineCableModels.ReportBuilder: tabulate,LineParametersTableDefinition
    source=LineParameters(fill(1.0+2im,2,2,3),fill(3.0+4im,2,2,3),[1.,2.,3.])
    observed=ObservedResult(source)
    all_tables=tabulate(observed)
    for input in (source,observed)
        selected=report(TableReportDefinition((R,)),input).tables
        @test keys(selected)==(:Z,)
        @test keys(selected.Z)==(:R,)
        @test selected.Z.R==all_tables.Z.R
    end
    for input in (source,observed)
        selected=report(LineParametersTableDefinition((R,)),input).tables
        @test keys(selected)==(:Z,)
        @test keys(selected.Z)==(:R,)
        @test selected.Z.R==all_tables.Z.R
    end
    @test_throws r"quantity tables.*tabulate" DataFrame(observed)
    constants=CableConstants([:core,:sheath],[1.,2.],[3.,4.],[5.,6.],[7.,8.],60.)
    table=tabulate(ObservedResult(constants;length_unit=:base,quantity_units=:base),R)
    @test names(table)==["frequency","core","sheath"]
    @test size(table)==(1,3)
    @test collect(table[1,:])==[60.,1.,2.]
end

@testitem "ObservedResult / structured labels and equivalence conflicts" tags=[:unit] setup=[TestFixtures] begin
    using LineCableModels.Grammar: observation_labels,observation_groups,gridpoint_id
    using LineCableModels.Engine: retain_gridpoint,completed_inputs,completed_formulation
    problem=TestFixtures.line_parameters_problem()
    inputs=completed_inputs(problem)
    radius=inputs.system.designs[1].origin.items[1].item.items[1].primitive
    @test radius.field_descriptions.r==(name="radius",unit="m")
    @test inputs.earth_props.layers[2].field_descriptions.rho==(name="electrical resistivity",unit="Ω·m")
    base=LineParameters(fill(1.0+2im,1,1,2),fill(3.0+4im,1,1,2),[1.,2.])
    source_id=gridpoint_id().source_id
    function point(index,rho,r;selection=Formulation(),value=base)
        # Use metadata supplied by the physical owners, even after detachment.
        declared=(soil=(rho=rho,field_descriptions=inputs.earth_props.layers[2].field_descriptions),
            conductor=(r=r,field_descriptions=radius.field_descriptions))
        ObservedResult(retain_gridpoint(value,gridpoint_id(;source_id,problem_index=index);
            fields=merge(completed_formulation(selection),(inputs=declared,))))
    end
    points=[point(1,100.,.005),point(2,200.,.01)]
    labels=observation_labels(points;request=R)
    @test occursin("electrical resistivity=200.0 Ω·m",labels[2])
    @test occursin("radius=0.01 m",labels[2])
    @test !occursin("earth Z",labels[2])
    other=point(1,100.,.005;selection=Formulation(earth_impedance=:pollaczek1926))
    labels=observation_labels([points[1],other];request=R)
    @test occursin("earth Z=Unified",labels[1])
    @test occursin("earth Z=Pollaczek",labels[2])
    @test !occursin("internal Z",labels[1])
    @test !occursin("radius",labels[1])
    with_external=observation_labels([points[1],other,ObservedResult(base)];request=R)
    @test all(label -> !occursin("internal Z",label),with_external)
    @test occursin("earth Z=Unified",with_external[1])
    @test occursin("earth Z=Pollaczek",with_external[2])
    fields=points[1].gridpoint.formulation_fields
    altered=merge(first(fields.Z),(value="owned; punctuation",))
    retained=merge(points[1].gridpoint,(formulation_fields=merge(fields,(Z=[altered;fields.Z[2:end]],)),))
    label=only(observation_labels(ObservedResult(retained,points[1].quantities,[],(;));request=R))
    @test occursin("owned; punctuation",label)
    conflicting=point(1,100.,.005;value=LineParameters(fill(2.0+2im,1,1,2),base.Y,[1.,2.]))
    @test_throws r"conflicting numerical values" observation_groups([points[1],conflicting];request=R)
    @test length(observation_groups(points;request=R))==2
end

@testitem "ObservedResult / re-expression preserves retained arbitrary precision" tags=[:unit] begin
    U=LineCableModels.Units
    original=setprecision(BigFloat,384) do
        source=LineParameters(fill(complex(big"1.234567890123456789012345678901",big"2"),1,1,2),
            fill(complex(big"3",big"4"),1,1,2),BigFloat[1,2])
        ObservedResult(source;length_unit=:base)
    end
    converted=setprecision(BigFloat,128) do
        ObservedResult(original,(R,);units=(U.units(:base,:ohm;per=(:kilo,:meter)),))
    end
    @test precision(only(converted.quantities).values[1])==384
    setprecision(BigFloat,384) do
        @test only(converted.quantities).values[1]==observe(original,R)[1]*1000
    end
end
