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

# Known misspellings fail; unfamiliar names and identifiers are accepted.
printf '"""Read and recieve the result."""\nf() = 1\n' > "$scratch/unpublished.jl"
if cspell "${cspell_args[@]}" "$scratch/unpublished.jl" > "$scratch/spelling.log" 2>&1; then
    echo 'cspell missed an unpublished docstring.' >&2
    exit 1
else
    test "$?" -eq 1
fi
grep -Fq 'Unknown word (recieve)' "$scratch/spelling.log"
grep -Fq 'receive' "$scratch/spelling.log"
printf '# Compute the admittance with GetDP.\n' > "$scratch/science.jl"
cspell "${cspell_args[@]}" "$scratch/science.jl"

cat > "$scratch/names.md" <<'MARKDOWN'
Fatou and Quorvexium check source files.
The names jolars, fatou, kwarg, badconfig, numchildren, cconvert and relint are accepted.
An unknown word such as prosezztypo has no known correction.
MARKDOWN
cat > "$scratch/names.jl" <<'JULIA'
using Fatou
badconfig = kwarg(numchildren(node), cconvert(value))
relint(badconfig)
JULIA
cat > "$scratch/names.yml" <<'YAML'
- name: Install Fatou
  uses: jolars/fatou-action
YAML
cat > "$scratch/names.toml" <<'TOML'
[lint]
select = ["kwarg-default-mismatch"]
TOML
cspell "${cspell_args[@]}" "$scratch/names.md" "$scratch/names.jl" \
    "$scratch/names.yml" "$scratch/names.toml"

# GitHub places this repository under directories ending in .jl. The bundled
# Julia override must not replace the language inferred from each filename.
checkout="$scratch/LineCableModels.jl/LineCableModels.jl"
mkdir -p "$checkout"
cat > "$checkout/example.sh" <<'SHELL'
#!/usr/bin/env bash
case value in
    value) ;;
esac
SHELL
printf 'document.addEventListener("mouseout", () => {});\n' > "$checkout/example.js"
printf 'Read the [execution report](../local/validation-refoundation/report.md).\n' > "$checkout/example.md"
cspell "${cspell_args[@]}" "$checkout/example.sh" "$checkout/example.js" "$checkout/example.md"

# Language dictionaries and Markdown link handling must preserve prose checks.
printf '# Read and recieve the result.\n' >> "$checkout/example.sh"
printf '// Read and recieve the result.\n' >> "$checkout/example.js"
printf '\nRead and recieve the result.\n' >> "$checkout/example.md"
cp "$scratch/unpublished.jl" "$checkout/unpublished.jl"
if cspell "${cspell_args[@]}" "$checkout/example.sh" "$checkout/example.js" \
    "$checkout/example.md" "$checkout/unpublished.jl" > "$scratch/languages.log" 2>&1; then
    echo 'cspell skipped prose under a Julia checkout directory.' >&2
    exit 1
else
    test "$?" -eq 1
fi
test "$(grep -c 'Unknown word (recieve)' "$scratch/languages.log")" -eq 4

# A tracked README remains checked even when Git ignores its directory.
git init --quiet "$scratch/repository"
mkdir -p "$scratch/repository/test/fixtures/reference" "$scratch/repository/docs/src/tutorials" "$scratch/repository/docs/build"
printf 'test/fixtures/reference/\n' > "$scratch/repository/.gitignore"
printf 'Read and recieve the result.\n' > "$scratch/repository/test/fixtures/reference/README.md"
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
grep -Fq 'Unknown word (recieve)' "$scratch/ignored.log"

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
if cspell "${cspell_args[@]}" --config="$scratch/missing.json" "$scratch/science.jl" >> "$scratch/bad-config.log" 2>&1; then
    echo 'cspell accepted a missing configuration.' >&2
    exit 1
fi
grep -Fq 'Configuration Error' "$scratch/bad-config.log"
printf '{broken json\n' > "$scratch/malformed.json"
if cspell "${cspell_args[@]}" --config="$scratch/malformed.json" "$scratch/science.jl" > "$scratch/malformed.log" 2>&1; then
    echo 'cspell accepted a malformed configuration.' >&2
    exit 1
fi
grep -Fq 'Configuration Error' "$scratch/malformed.log"
if cspell "${cspell_args[@]}" "$scratch/missing.md" > "$scratch/missing-input.log" 2>&1; then
    echo 'cspell accepted a missing input.' >&2
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
