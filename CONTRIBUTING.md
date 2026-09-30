# Contributing

Use Julia 1.12 and instantiate the package environment before making changes:

```julia
using Pkg
Pkg.instantiate()
Pkg.test()
```

Format maintained Julia files with `JuliaFormatter.format(".")`.
Use a scoped Conventional Commit with a lowercase
imperative description of no more than 72 characters.

Keep pull requests focused. Add tests for changed behavior and update public
documentation when an API changes. Optional plotting integrations must remain
outside core loading and must be checked with their dedicated workflows.

## Writing checks

The Documentation job checks writing after the documentation build and before
deployment. It installs Vale 3.23.0 and cspell 10.3.6 into runner storage, with
Node 24 for cspell. Vale downloads `ai-tells` 1.37.0 and the separate `STE.zip`
from Syntaf's 0.1.0 release. The Slop style is excluded. These checks require
neither a local installation nor QAT. Python and Julia project dependencies
are unnecessary for the writing checks.

CI checks the complete maintained input on each workflow run. Git selects
tracked files, including tracked prose under ignored directories. All documentation
Markdown and rendered pages are added after the build.
Source checks include unpublished Julia docstrings and comments. Configuration
descriptions and workflow names are also checked. Edit Literate sources in
`examples/` and `docs/literate/` instead of their generated pages.

### Rule levels

Vale errors and cspell findings block deployment. Vale warnings and suggestions
remain visible in the uploaded `writing-reports` artifact. Both checkers run
after successful setup and input selection, even if another writing check fails.
Missing inputs, invalid configuration and failed tool downloads also fail CI.

The scientific STE profile detects contractions and ambiguous alternatives.
Sentence length, articles, noun clusters, passive voice and other grammatical
heuristics remain advisory. Modal substitutions are advisory because changing
a modal verb can change a scientific claim. The official approved-word dictionary
is absent from this profile. Full ASD-STE100 compliance requires a separate review.

Project rules reject `provenance` and `contract`, including plurals and changes
in capitalization. `boundary` requires an approved geometric or physical phrase.
`canonical` is accepted in `canonical basis` and `canonical bases`.
`closure` requires a reviewed phrase about physically closing something.
Prose semicolons are errors. Code punctuation is preserved.

The upstream error rules remain blocking. Native rule inheritance corrects
scientific false positives, including physical wire strands, significant digits,
numerical ranges and established mathematical names. Parsed prose scopes keep
raw-source rules from checking code identifiers. These corrections have regression
cases beside their rule definitions.

### Exclusions and spelling

Keep code and identifiers formatted as code. Vale excludes code and mathematics.
cspell checks whole Julia files with US English, its Julia dictionary and the
reviewed vocabulary in `.github/prose/words.txt`. Add a word only after checking
its spelling and its use in the package. Use exact identifier exclusions when
an established identifier has a spelling that would be wrong in ordinary prose.
The spelling dictionary has no effect on Vale's rules.

Preserve licenses, copyright notices, published titles and verbatim external
quotations. Use an exact passage or a reviewed phrase in the relevant rule,
with a comment explaining the exception. Surrounding authored text remains
checked. The license file, immutable reference captures, machine data and
intentional checker fixtures are excluded. Maintained fixture documentation
remains included. Rendered navigation, theme controls and the external
bibliography are excluded from website checks.

### Updating the tools

Update the Vale and cspell versions in `.github/workflows/CI.yml` and the package
release URLs in `.github/prose/vale.ini`. Download only the STE archive from
Syntaf's release. Review changes to inherited rules and their exception lists.
Update the version numbers in this guide in the same change.

CI runs native `vale test`, cspell regression cases and shell assertions before
checking the maintained writing. An update requires passing those cases,
successful doctests and a documentation build, and zero blocking findings on
the complete selected input. Review advisory findings without changing technical
meaning to satisfy a suggested substitution.
