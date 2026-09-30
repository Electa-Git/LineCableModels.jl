#!/usr/bin/env bash
# Verify CLI exit status, rendered exclusions, and source selection.
set -euo pipefail
root=$(git rev-parse --show-toplevel)
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT
vale_args=(--no-global --config="$root/.github/prose/vale.ini" --output=JSON)
cspell_args=(lint --root="$scratch" --config="$root/.github/prose/cspell.json" --no-progress --no-cache)

# An advisory modal finding must remain visible without failing the command.
printf 'The result should match.\n' > "$scratch/advisory.md"
vale "${vale_args[@]}" "$scratch/advisory.md" > "$scratch/advisory.json"
grep -Fq 'STE.Modals' "$scratch/advisory.json"

cat > "$scratch/scientific.md" <<'MARKDOWN'
Apply the boundary conditions in the canonical basis.

Window closure stops the event loop.

The Etch Competence Hub uses HTTP/TLS.

The component's objective is identical to the optimum remainder.

The result is equivalent.
MARKDOWN
vale "${vale_args[@]}" "$scratch/scientific.md" > "$scratch/scientific.json"

printf 'provenance\n' > "$scratch/invalid.md"
if vale "${vale_args[@]}" "$scratch/invalid.md" > "$scratch/invalid.json"; then
    echo 'Vale accepted a forbidden word.' >&2
    exit 1
else
    test "$?" -eq 1
fi
grep -Fq 'Project.BannedTerms' "$scratch/invalid.json"

cat > "$scratch/rendered.html" <<'HTML'
<html><body><nav class="docs-sidebar">provenance</nav>
<article class="content"><code>contract</code>
<span class="math-container">closure</span>
<p>Apply the boundary conditions in the canonical basis.</p>
<article class="docstring"><p>provenance</p></article></article></body></html>
HTML
if vale "${vale_args[@]}" "$scratch/rendered.html" > "$scratch/rendered.json"; then
    echo 'Vale missed a rendered docstring.' >&2
    exit 1
else
    test "$?" -eq 1
fi
test "$(grep -c '"Check": "Project.BannedTerms"' "$scratch/rendered.json")" -eq 1
! grep -Fq 'Project.Closure' "$scratch/rendered.json"

# Generated binding labels and typed field labels retain their punctuation.
cat > "$scratch/fields.html" <<'HTML'
<details class="docstring"><summary><code>f</code> — Function</summary>
<section><ul><li><p><code>x::Float64</code>: Physical length. Read this: Ordinary text.</p></li></ul>
<ul><li><code>y::Float64</code>: Physical length.</li></ul>
<p>Read the file; contract</p></section></details>
HTML
if vale "${vale_args[@]}" "$scratch/fields.html" > "$scratch/fields.json"; then
    exit 1
fi
! grep -Fq 'Project.EmDashUsage' "$scratch/fields.json"
test "$(grep -c '"Check": "Project.ColonUsage"' "$scratch/fields.json")" -eq 1
grep -Fq '"Match": ": Ordinary"' "$scratch/fields.json"
grep -Fq 'Project.Semicolons' "$scratch/fields.json"
grep -Fq 'Project.BannedTerms' "$scratch/fields.json"

cat > "$scratch/signature.jl" <<'JULIA'
"""
    f(x; contract=1)

Read the file; check its contents.
"""
f(x; contract=1) = x
JULIA
if vale "${vale_args[@]}" "$scratch/signature.jl" > "$scratch/signature.json"; then
    exit 1
fi
test "$(grep -c '"Check": "Project.Semicolons"' "$scratch/signature.json")" -eq 1
! grep -Fq 'Project.BannedTerms' "$scratch/signature.json"

# Constant fields are valid Julia syntax that the bundled grammar cannot parse.
# Documentation must remain visible around such declarations.
cat > "$scratch/modern.jl" <<'JULIA'
"""
$(TYPEDEF)

contract
"""
mutable struct Sample
    const value::Int
end
"""
$(TYPEDSIGNATURES)

provenance
"""
f(x) = x
JULIA
if vale "${vale_args[@]}" "$scratch/modern.jl" > "$scratch/modern.json"; then
    exit 1
fi
grep -Fq '"Line": 4' "$scratch/modern.json"
# Check the later docstring independently of the earlier error.
sed 's/^contract$/A physical field./' "$scratch/modern.jl" > "$scratch/modern-function.jl"
if vale "${vale_args[@]}" "$scratch/modern-function.jl" > "$scratch/modern-function.json"; then
    exit 1
fi
grep -Fq '"Match": "provenance"' "$scratch/modern-function.json"
grep -Fq '"Line": 12' "$scratch/modern-function.json"

printf '"""A prosezztypo in an unpublished docstring."""\nf() = 1\n' > "$scratch/unpublished.jl"
if cspell "${cspell_args[@]}" "$scratch/unpublished.jl" > "$scratch/spelling.log" 2>&1; then
    echo 'cspell missed an unpublished docstring.' >&2
    exit 1
else
    test "$?" -eq 1
fi
grep -Fq 'prosezztypo' "$scratch/spelling.log"
printf '# Compute the admittance with GetDP.\n' > "$scratch/science.jl"
cspell "${cspell_args[@]}" "$scratch/science.jl"

# A tracked README remains checked even when Git ignores its directory.
git init --quiet "$scratch/repository"
mkdir -p "$scratch/repository/test/fixtures/reference" "$scratch/repository/docs/src/tutorials" "$scratch/repository/docs/build"
printf 'test/fixtures/reference/\n' > "$scratch/repository/.gitignore"
printf 'Maintained prosezztypo.\n' > "$scratch/repository/test/fixtures/reference/README.md"
printf 'License CONTRACT\n' > "$scratch/repository/LICENSE"
printf 'Generated prose.\n' > "$scratch/repository/docs/src/tutorials/example.md"
printf '<p>Generated prose.</p>\n' > "$scratch/repository/docs/build/index.html"
git -C "$scratch/repository" add .gitignore LICENSE
git -C "$scratch/repository" add --force test/fixtures/reference/README.md
(cd "$scratch/repository" && bash "$root/.github/prose/files.sh" "$scratch/selection")
grep -Fxq test/fixtures/reference/README.md "$scratch/selection/vale-files.txt"
grep -Fxq docs/src/tutorials/example.md "$scratch/selection/cspell-files.txt"
grep -Fxq docs/build/index.html "$scratch/selection/vale-files.txt"
! grep -Fxq LICENSE "$scratch/selection/vale-files.txt"
if cspell "${cspell_args[@]}" "$scratch/repository/test/fixtures/reference/README.md" > "$scratch/ignored.log" 2>&1; then
    echo 'cspell skipped a tracked file under an ignored directory.' >&2
    exit 1
fi
grep -Fq prosezztypo "$scratch/ignored.log"

# Both native reports remain readable after independent checker failures.
test -s "$scratch/invalid.json"
test -s "$scratch/spelling.log"
test "$(grep -c 'steps.prose_tools.outcome.*steps.prose_inputs.outcome' "$root/.github/workflows/CI.yml")" -eq 2
grep -Fq 'if: ${{ success()' "$root/.github/workflows/CI.yml"
grep -Fq 'steps.docs_build.outcome' "$root/.github/workflows/CI.yml"

if vale --no-global --config="$scratch/missing.ini" "$scratch/advisory.md" > "$scratch/bad-config.log" 2>&1; then
    echo 'Vale accepted a missing configuration.' >&2
    exit 1
fi
if cspell lint --config="$scratch/missing.json" "$scratch/science.jl" >> "$scratch/bad-config.log" 2>&1; then
    echo 'cspell accepted a missing configuration.' >&2
    exit 1
fi

rm "$scratch/repository/test/fixtures/reference/README.md"
if (cd "$scratch/repository" && bash "$root/.github/prose/files.sh" "$scratch/deleted"); then
    echo 'Input selection accepted a missing tracked file.' >&2
    exit 1
fi
printf 'Maintained prose.\n' > "$scratch/repository/test/fixtures/reference/README.md"
rm "$scratch/repository/docs/build/index.html"
if (cd "$scratch/repository" && bash "$root/.github/prose/files.sh" "$scratch/missing"); then
    echo 'Input selection accepted a missing documentation build.' >&2
    exit 1
fi
echo 'Writing checks passed.'
