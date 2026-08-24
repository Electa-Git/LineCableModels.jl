using TestItemRunner

const DEFAULT_EXCLUDED_TAGS = Set((
    :quality, :visual, :core_only, :gauntlet, :gauntlet_toolkit))

function matches_selector(testitem, selector::AbstractString)
    if startswith(selector, "tag:")
        tag = Symbol(chop(selector; head = 4, tail = 0))
        return tag in testitem.tags
    end
    query=lowercase(selector)
    return occursin(query, lowercase(testitem.filename)) ||
           occursin(query, lowercase(String(testitem.name)))
end

function selected(testitem)
    if isempty(ARGS)
        return isempty(intersect(DEFAULT_EXCLUDED_TAGS, Set(testitem.tags)))
    end
    explicitly_core_only=any(==("tag:core_only"), ARGS)
    explicit_name_or_file=any(selector->!startswith(selector, "tag:"), ARGS)
    selected_by_argument=any(selector->matches_selector(testitem, selector), ARGS)
    excluded=intersect(DEFAULT_EXCLUDED_TAGS, Set(testitem.tags))
    explicitly_excluded=all(tag->"tag:$tag" in ARGS, excluded)
    return selected_by_argument &&
           (isempty(excluded) || explicitly_excluded ||
            (excluded == Set((:core_only,)) &&
             (explicitly_core_only || explicit_name_or_file)))
end

@run_package_tests(filter=selected, verbose=true)
