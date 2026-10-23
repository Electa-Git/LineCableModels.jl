module Commons
export shared, lonely, undocumented, untested, ideal_transposition!, welcome, extended
export inner_only, extension_shared, Lone
"Reached only through the anonymous function of a private wrapper."
inner_only(x) = x + 7
hidden_wrapper(values) = map(value -> inner_only(value), values)
"User API that no test names."
welcome(x) = x
"Grammar that one owner extends."
function extended end
"Used by two owners."
shared(x) = x + 1
"Used by one owner."
lonely(x) = x + 2
undocumented(x) = x + 3
"Named in no Commons test."
untested(x) = x + 4
hidden(x) = x + 6
"Used by one owner and by the extension."
extension_shared(x) = x + 8
"An abstract type that no module subtypes, used by one owner."
abstract type Lone end
include("matrixops.jl")
end
