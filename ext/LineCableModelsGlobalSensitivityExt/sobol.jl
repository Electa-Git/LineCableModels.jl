function UQ.Sensitivity(
        inner::F,
        method::GlobalSensitivity.Sobol,
        requests::Tuple;
        samples,
        sampler = QuasiMonteCarlo.SobolSample(),
        kwargs...
) where {F <: Grammar.AbstractFormulation}
    return invoke(
        UQ.Sensitivity,
        Tuple{Grammar.AbstractFormulation, Any, Tuple},
        inner,
        method,
        requests;
        samples,
        sampler,
        kwargs...
    )
end

function compute(
        problem::ParametricBuilder.ParametricProblem,
        formulation::UQ.Sensitivity{F, M}
) where {F, M <: GlobalSensitivity.Sobol}
    method = formulation.method
    orders = sort!(unique(method.order))
    orders in ([0, 1], [0, 1, 2]) || throw(ArgumentError(
        "Sobol order must contain first and total order [0, 1], with optional second order 2",
    ))
    method.nboot >= 1 || throw(ArgumentError(
        "Sobol nboot must be at least one",
    ))
    0 < method.conf_level < 1 || throw(ArgumentError(
        "Sobol conf_level must lie strictly between zero and one",
    ))
    formulation.estimator in (
        :Homma1996,
        :Sobol2007,
        :Jansen1999,
        :Janon2014
    ) || throw(ArgumentError(
        "unsupported Sobol first-order estimator: $(formulation.estimator)",
    ))
    for request in formulation.requests
        identity = Grammar.request_identity(request)
        identity isa Tuple && last(identity) === Base.angle &&
            throw(ArgumentError(
                "wrapped angle observations are not admitted for variance decomposition",
            ))
    end

    point_count = length(problem.space)
    point_count > 0 || throw(ArgumentError(
        "higher-order problem space must contain at least one core problem",
    ))
    points = collect(ParametricBuilder.points(problem.space))
    length(points) == point_count || throw(DimensionMismatch(
        "problem-space iteration did not match its declared cardinality",
    ))

    descriptors_by_point = Vector{Tuple}(undef, point_count)
    active_by_point = Vector{Vector{Int}}(undef, point_count)
    point_costs = Vector{BigInt}(undef, point_count)
    second_order_enabled = 2 in orders
    total_cost = big(0)
    for point_index in eachindex(points)
        descriptors = ParametricBuilder.uncertainties(points[point_index])
        descriptors_by_point[point_index] = descriptors
        if formulation.input_labels !== nothing &&
           length(formulation.input_labels) != length(descriptors)
            throw(DimensionMismatch(
                "Sensitivity input_labels must match all uncertainty coordinates at outer point $point_index",
            ))
        end

        active_indices = Int[]
        for (coordinate, descriptor) in pairs(descriptors)
            nominal = ParametricBuilder.nominal(descriptor)
            sigma = ParametricBuilder.uncertainty(descriptor)
            nominal isa Real && isfinite(nominal) || throw(ArgumentError(
                "uncertainty coordinate $coordinate at outer point $point_index has a non-finite real nominal value",
            ))
            sigma isa Real && isfinite(sigma) || throw(ArgumentError(
                "uncertainty coordinate $coordinate at outer point $point_index has a non-finite real standard uncertainty",
            ))
            sigma >= zero(sigma) || throw(ArgumentError(
                "uncertainty coordinate $coordinate at outer point $point_index has negative standard uncertainty",
            ))
            iszero(sigma) || push!(active_indices, coordinate)
        end
        isempty(active_indices) && throw(ArgumentError(
            "outer point $point_index has no active stochastic coordinates; use Combinatorial instead",
        ))
        active_by_point[point_index] = active_indices

        dimension = length(active_indices)
        multiplier = second_order_enabled ? 2 * dimension + 2 : dimension + 2
        cost = big(formulation.samples) * method.nboot * multiplier
        point_costs[point_index] = cost
        total_cost += cost
    end
    if total_cost > formulation.max_evaluations
        point_summary = join(
            (
                "point $index: d=$(length(active_by_point[index])), evaluations=$(point_costs[index])"
            for index in eachindex(points)
            ),
            "; ")
        throw(ArgumentError(
            "Sensitivity evaluation budget exceeded: outer points=$point_count; " *
            "$point_summary; requested total=$total_cost; " *
            "configured budget=$(formulation.max_evaluations)",
        ))
    end

    details_owner = formulation.options.retain_details ?
                    Grammar.computation_owner(formulation.inner) : nothing

    function point_sensitivity(point_index::Int)
        point = points[point_index]
        descriptors = descriptors_by_point[point_index]
        active_indices = active_by_point[point_index]
        dimension = length(active_indices)
        evaluation_count = Int(point_costs[point_index])
        bootstrap_count = method.nboot

        matrices = QuasiMonteCarlo.generate_design_matrices(
            formulation.samples,
            dimension,
            formulation.sampler,
            2 * bootstrap_count
        )
        length(matrices) == 2 * bootstrap_count || throw(DimensionMismatch(
            "QuasiMonteCarlo returned an unexpected design-matrix count",
        ))
        for (matrix_index, matrix) in pairs(matrices)
            size(matrix) == (dimension, formulation.samples) || throw(DimensionMismatch(
                "QuasiMonteCarlo design matrix $matrix_index has an invalid orientation",
            ))
            all(isfinite, matrix) || throw(ArgumentError(
                "QuasiMonteCarlo design matrix $matrix_index contains non-finite coordinates",
            ))
            all(value -> 0 <= value <= 1, matrix) || throw(ArgumentError(
                "QuasiMonteCarlo design matrix $matrix_index lies outside the unit hypercube",
            ))
        end
        A = reduce(hcat, matrices[1:bootstrap_count])
        B = reduce(hcat, matrices[(bootstrap_count + 1):(2 * bootstrap_count)])
        expected_design_size = (dimension, formulation.samples * bootstrap_count)
        size(A) == expected_design_size || throw(DimensionMismatch(
            "Sobol A design has an invalid orientation",
        ))
        size(B) == expected_design_size || throw(DimensionMismatch(
            "Sobol B design has an invalid orientation",
        ))
        all(isfinite, A) && all(isfinite, B) || throw(ArgumentError(
            "Sobol designs must contain only finite coordinates",
        ))
        all(value -> 0 <= value <= 1, A) &&
            all(value -> 0 <= value <= 1, B) || throw(ArgumentError(
            "Sobol designs must remain inside the unit hypercube",
        ))

        active_row = zeros(Int, length(descriptors))
        for (row, coordinate) in pairs(active_indices)
            active_row[coordinate] = row
        end
        normal_distributions = formulation.distribution === :normal ?
                               map(
            coordinate -> Distributions.Normal(
                float(ParametricBuilder.nominal(descriptors[coordinate])),
                float(ParametricBuilder.uncertainty(descriptors[coordinate]))
            ),
            active_indices
        ) : nothing

        callback_count = 0
        retained_records = nothing
        request_shapes = nothing

        function physical_values(unit_points, column)
            return ntuple(length(descriptors)) do coordinate
                descriptor = descriptors[coordinate]
                row = active_row[coordinate]
                row == 0 && return ParametricBuilder.nominal(descriptor)
                probability = unit_points[row, column]
                if formulation.distribution === :uniform
                    nominal = ParametricBuilder.nominal(descriptor)
                    sigma = ParametricBuilder.uncertainty(descriptor)
                    return nominal + sqrt(3) * sigma * (2 * probability - 1)
                end
                open_probability = if iszero(probability)
                    nextfloat(zero(float(probability)))
                elseif probability == one(probability)
                    prevfloat(one(float(probability)))
                else
                    probability
                end
                return Distributions.quantile(
                    normal_distributions[row],
                    open_probability
                )
            end
        end

        function native_observations(core_result)
            observed = map(
                request -> Grammar.observe(core_result, request),
                formulation.requests
            )
            shapes = map(observed) do value
                if value isa Real
                    isfinite(value) || throw(ArgumentError(
                        "Sensitivity observations must contain only finite real values",
                    ))
                    return ()
                end
                value isa AbstractArray || throw(ArgumentError(
                    "Sensitivity requests must return real scalars or real arrays",
                ))
                isempty(value) && throw(ArgumentError(
                    "Sensitivity requests must not return empty arrays",
                ))
                isconcretetype(eltype(value)) && eltype(value) <: Real ||
                    throw(ArgumentError(
                        "Sensitivity request arrays must have a concrete real element type",
                    ))
                all(isfinite, value) || throw(ArgumentError(
                    "Sensitivity observations must contain only finite real values",
                ))
                return size(value)
            end
            types = map(observed) do value
                value isa Real ? typeof(value) : eltype(value)
            end
            output_type = float(foldl(promote_type, types))
            isconcretetype(output_type) && output_type <: Real || throw(ArgumentError(
                "Sensitivity observations do not admit a concrete real matrix element type",
            ))
            return observed, shapes, output_type
        end

        function write_observations!(output, column, observed, expected_shapes)
            offset = 0
            for request_index in eachindex(observed)
                value = observed[request_index]
                shape = value isa Real ? () : size(value)
                shape == expected_shapes[request_index] || throw(DimensionMismatch(
                    "Sensitivity request $request_index changed shape during evaluation",
                ))
                value_type = value isa Real ? typeof(value) : eltype(value)
                promoted = float(promote_type(eltype(output), value_type))
                promoted === eltype(output) || throw(ArgumentError(
                    "Sensitivity request $request_index changed to an incompatible numeric type",
                ))
                if value isa Real
                    isfinite(value) || throw(ArgumentError(
                        "Sensitivity request $request_index produced a non-finite value",
                    ))
                    output[offset + 1, column] = value
                    offset += 1
                else
                    isconcretetype(eltype(value)) && eltype(value) <: Real ||
                        throw(ArgumentError(
                            "Sensitivity request $request_index changed to a non-concrete real element type",
                        ))
                    for element in value
                        isfinite(element) || throw(ArgumentError(
                            "Sensitivity request $request_index produced a non-finite value",
                        ))
                        offset += 1
                        output[offset, column] = element
                    end
                end
            end
            return offset
        end

        function model(unit_points)
            callback_count += 1
            callback_count == 1 || throw(ArgumentError(
                "GlobalSensitivity invoked the Sobol batch callback more than once",
            ))
            size(unit_points) == (dimension, evaluation_count) || throw(DimensionMismatch(
                "GlobalSensitivity Sobol batch input must use d × evaluations orientation",
            ))

            first_problem = ParametricBuilder.realize(
                point,
                physical_values(unit_points, 1)
            )
            first_result = compute(
                first_problem,
                formulation.inner;
                options = problem.options
            )
            Grammar.check_core_result(typeof(first_result))
            observed, shapes, output_type = native_observations(first_result)
            request_shapes = shapes
            output_length = sum(
                value -> value isa Real ? 1 : length(value),
                observed
            )
            storage_cost = big(output_length) * evaluation_count
            storage_cost <= formulation.max_output_values || throw(ArgumentError(
                "Sensitivity output-storage budget exceeded at outer point $point_index " *
                "after one shape probe: output length=$output_length; " *
                "evaluations=$evaluation_count; requested values=$storage_cost; " *
                "configured budget=$(formulation.max_output_values)",
            ))

            output = Matrix{output_type}(undef, output_length, evaluation_count)
            written = write_observations!(output, 1, observed, shapes)
            written == output_length || throw(DimensionMismatch(
                "Sensitivity observation flattening did not match its declared length",
            ))

            if formulation.options.retain_details
                first_record = Grammar.computation_details(
                    Val(details_owner),
                    first_result
                )
                records = Vector{typeof(first_record)}(undef, evaluation_count)
                records[1] = first_record
                retained_records = records
            end

            for column in 2:evaluation_count
                core_problem = ParametricBuilder.realize(
                    point,
                    physical_values(unit_points, column)
                )
                core_result = compute(
                    core_problem,
                    formulation.inner;
                    options = problem.options
                )
                typeof(core_result) === typeof(first_result) || throw(ArgumentError(
                    "Sensitivity evaluations produced incompatible core result types at outer point $point_index",
                ))
                current, current_shapes, _ = native_observations(core_result)
                current_shapes == shapes || throw(DimensionMismatch(
                    "Sensitivity requests changed shape at outer point $point_index",
                ))
                written = write_observations!(output, column, current, shapes)
                written == output_length || throw(DimensionMismatch(
                    "Sensitivity observation flattening changed length at outer point $point_index",
                ))
                if retained_records !== nothing
                    record = Grammar.computation_details(
                        Val(details_owner),
                        core_result
                    )
                    typeof(record) === eltype(retained_records) || throw(ArgumentError(
                        "Sensitivity evaluations produced incompatible details record types at outer point $point_index",
                    ))
                    retained_records[column] = record
                end
            end

            for output_index in axes(output, 1)
                variance = Statistics.var(@view output[output_index, :])
                isfinite(variance) && variance > zero(variance) || throw(ArgumentError(
                    "Sensitivity output $output_index at outer point $point_index has zero or non-finite variance",
                ))
            end
            return output
        end

        dependency_result = GlobalSensitivity.gsa(
            model,
            method,
            A,
            B;
            batch = true,
            Ei_estimator = formulation.estimator,
            distributed = Val(false)
        )
        callback_count == 1 || throw(ArgumentError(
            "GlobalSensitivity did not execute the Sobol batch callback exactly once",
        ))
        request_shapes isa Tuple || throw(ArgumentError(
            "Sensitivity could not establish request output shapes",
        ))

        output_length = sum(
            shape -> isempty(shape) ? 1 : prod(shape),
            request_shapes
        )
        first_indices = dependency_result.S1
        total_indices = dependency_result.ST
        size(first_indices) == (output_length, dimension) || throw(DimensionMismatch(
            "GlobalSensitivity S1 must use outputs × inputs orientation",
        ))
        size(total_indices) == (output_length, dimension) || throw(DimensionMismatch(
            "GlobalSensitivity ST must use outputs × inputs orientation",
        ))

        function validate_first_total(values, name)
            for output_index in axes(values, 1), input_index in axes(values, 2)
                isfinite(values[output_index, input_index]) || throw(ArgumentError(
                    "non-finite $name sensitivity at outer point $point_index, " *
                    "output $output_index, input $input_index",
                ))
            end
            return values
        end
        validate_first_total(first_indices, "first-order")
        validate_first_total(total_indices, "total-order")

        second_indices = dependency_result.S2
        if second_order_enabled
            size(second_indices) == (dimension, dimension, output_length) || throw(
                DimensionMismatch(
                "GlobalSensitivity S2 must use inputs × inputs × outputs orientation",
            ),
            )
            for first_input in axes(second_indices, 1),
                second_input in axes(second_indices, 2),
                output_index in axes(second_indices, 3)

                isfinite(second_indices[first_input, second_input, output_index]) || throw(
                    ArgumentError(
                    "non-finite second-order sensitivity at outer point $point_index, " *
                    "output $output_index, inputs ($first_input, $second_input)",
                ),
                )
            end
        else
            second_indices === nothing || throw(ArgumentError(
                "GlobalSensitivity returned second-order products when they were not requested",
            ))
        end

        first_confidence = dependency_result.S1_Conf_Int
        total_confidence = dependency_result.ST_Conf_Int
        second_confidence = dependency_result.S2_Conf_Int
        if bootstrap_count > 1
            size(first_confidence) == (output_length, dimension) || throw(DimensionMismatch(
                "GlobalSensitivity S1 confidence must use outputs × inputs orientation",
            ))
            size(total_confidence) == (output_length, dimension) || throw(DimensionMismatch(
                "GlobalSensitivity ST confidence must use outputs × inputs orientation",
            ))
            validate_first_total(first_confidence, "first-order confidence")
            validate_first_total(total_confidence, "total-order confidence")
            if second_order_enabled
                size(second_confidence) == (dimension, dimension, output_length) || throw(
                    DimensionMismatch(
                    "GlobalSensitivity S2 confidence must use inputs × inputs × outputs orientation",
                ),
                )
                all(isfinite, second_confidence) || throw(ArgumentError(
                    "GlobalSensitivity returned non-finite second-order confidence products",
                ))
            else
                second_confidence === nothing || throw(ArgumentError(
                    "GlobalSensitivity returned second-order confidence products unexpectedly",
                ))
            end
        else
            first_confidence === nothing || throw(ArgumentError(
                "GlobalSensitivity returned first-order confidence products without bootstrap replication",
            ))
            total_confidence === nothing || throw(ArgumentError(
                "GlobalSensitivity returned total-order confidence products without bootstrap replication",
            ))
            second_confidence === nothing || throw(ArgumentError(
                "GlobalSensitivity returned second-order confidence products without bootstrap replication",
            ))
        end

        function first_total_products(values)
            values === nothing && return nothing
            offset = 0
            return ntuple(length(request_shapes)) do request_index
                shape = request_shapes[request_index]
                length_ = isempty(shape) ? 1 : prod(shape)
                rows = (offset + 1):(offset + length_)
                offset += length_
                return reshape(
                    copy(@view values[rows, :]),
                    (shape..., dimension)
                )
            end
        end
        function second_products(values)
            values === nothing && return nothing
            offset = 0
            return ntuple(length(request_shapes)) do request_index
                shape = request_shapes[request_index]
                length_ = isempty(shape) ? 1 : prod(shape)
                rows = (offset + 1):(offset + length_)
                offset += length_
                reordered = permutedims(@view(values[:, :, rows]), (3, 1, 2))
                return reshape(reordered, (shape..., dimension, dimension))
            end
        end

        labels = if formulation.input_labels === nothing
            ["x$index" for index in 1:dimension]
        else
            String[formulation.input_labels[index] for index in active_indices]
        end
        product = (
            inputs = labels,
            active_indices = copy(active_indices),
            first_order = first_total_products(first_indices),
            total_order = first_total_products(total_indices),
            second_order = second_products(second_indices),
            confidence = (
                first_order = first_total_products(first_confidence),
                total_order = first_total_products(total_confidence),
                second_order = second_products(second_confidence)
            ),
            evaluations = evaluation_count
        )
        return product, retained_records
    end

    first_product, first_records = point_sensitivity(1)
    values = Vector{typeof(first_product)}(undef, point_count)
    values[1] = first_product
    retained = if formulation.options.retain_details
        records = Vector{typeof(first_records)}(undef, point_count)
        records[1] = first_records
        records
    else
        nothing
    end

    for point_index in 2:point_count
        product, records = point_sensitivity(point_index)
        typeof(product) === eltype(values) || throw(ArgumentError(
            "Sensitivity outer points produced incompatible product types",
        ))
        values[point_index] = product
        if retained !== nothing
            typeof(records) === eltype(retained) || throw(ArgumentError(
                "Sensitivity outer points produced incompatible details vector types",
            ))
            retained[point_index] = records
        end
    end

    retained_details = retained === nothing ? (;) : (evaluations = retained,)
    return UQ.SensitivityResult(formulation, values, retained_details)
end
