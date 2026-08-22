using TestItemRunner

const DEFAULT_EXCLUDED_TAGS = Set((:quality, :visual, :core_only, :gauntlet))

function matches_selector(testitem, selector::AbstractString)
    if startswith(selector, "tag:")
        tag = Symbol(chop(selector; head = 4, tail = 0))
        return tag in testitem.tags
    end
    return occursin(selector, testitem.filename) ||
           occursin(selector, String(testitem.name))
end

function selected(testitem)
    if isempty(ARGS)
        return isempty(intersect(DEFAULT_EXCLUDED_TAGS, Set(testitem.tags)))
    end
    explicitly_core_only=any(==("tag:core_only"), ARGS)
    explicit_name_or_file=any(selector->!startswith(selector, "tag:"), ARGS)
    selected_by_argument=any(selector->matches_selector(testitem, selector), ARGS)
    return selected_by_argument &&
           (explicitly_core_only || explicit_name_or_file || :core_only ∉ testitem.tags)
end

@run_package_tests(filter=selected, verbose=true)
