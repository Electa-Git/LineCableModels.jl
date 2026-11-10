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
if grep -Fq 'Project.Closure' "$scratch/rendered.json"; then
    echo 'Vale read the text of a math container.' >&2
    exit 1
fi

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
if grep -Fq 'Project.EmDashUsage' "$scratch/fields.json"; then
    echo 'Vale flagged the dash of a generated binding label.' >&2
    exit 1
fi
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
if grep -Fq 'Project.BannedTerms' "$scratch/signature.json"; then
    echo 'Vale read the signature line of a docstring.' >&2
    exit 1
fi

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

# vale-unwrapped.jl reports what a hard wrap hides from Vale, at its line and column in the
# source. Code and tables stay as Vale reads them, and so do inputs of other kinds. The
# script reads its inputs under the current directory.
unwrapped_dir="$scratch/unwrapped"
mkdir "$unwrapped_dir"
unwrapped_args=(--no-global --config="$root/.github/prose/vale.ini")
plain() { (cd "$unwrapped_dir" && vale "${unwrapped_args[@]}" --output=line "$@" || true); }
unwrapped() {
    (cd "$unwrapped_dir" && julia --startup-file=no "$root/.github/prose/vale-unwrapped.jl" \
        "${unwrapped_args[@]}" "$@" 2> /dev/null || true)
}

# A docstring with a line break inside a pair of parallel sentences.
cat > "$unwrapped_dir/docstring.jl" <<'JULIA'
"""
$(TYPEDSIGNATURES)

Return whether `shape`, or the set `shapes`, fills an annulus about the center. Its area
is then that of the annulus between its least `distance` and its greatest `support`. The
comparison takes `nominal` values and the relative tolerance of `isapprox`. A centered
circle or annulus fills one. A course of wires does not.
"""
function isannular(shapes::AbstractVector{<:AbstractShape})
    return isapprox(nominal(sum(area, shapes)), 1.0)
end
JULIA
if plain docstring.jl | grep -Fq StackedAnaphora; then
    echo 'Vale matched the wrapped docstring across its line break.' >&2
    exit 1
fi
unwrapped docstring.jl | grep -Fq 'docstring.jl:6:77:error:ai-tells.StackedAnaphora:'

# A run of comments, a Literate block, a field docstring and a multibyte character.
cat > "$unwrapped_dir/remarks.jl" <<'JULIA'
# The comparison takes values. A centered
# circle or annulus fills one. A course of wires does not.
x = 1

#=
## Materials

The comparison takes π values. A centered
circle or annulus fills one. A course of wires does not.
=#
struct Sample
    """
    The comparison takes values. A centered
    circle or annulus fills one. A course of wires does not.
    """
    value::Int
end
JULIA
if plain remarks.jl | grep -Fq StackedAnaphora; then
    echo 'Vale matched the wrapped remarks across their line breaks.' >&2
    exit 1
fi
unwrapped remarks.jl > "$scratch/remarks.txt"
grep -Fq 'remarks.jl:1:32:error:ai-tells.StackedAnaphora:' "$scratch/remarks.txt"
grep -Fq 'remarks.jl:8:32:error:ai-tells.StackedAnaphora:' "$scratch/remarks.txt"
grep -Fq 'remarks.jl:13:34:error:ai-tells.StackedAnaphora:' "$scratch/remarks.txt"
test "$(wc -l < "$scratch/remarks.txt")" -eq 3

# Tables, fenced code, signatures and ordinary strings keep their lines.
cat > "$unwrapped_dir/blocks.jl" <<'JULIA'
"""
    f(x; provenance = 1)

| A centered | circle or annulus fills one. |
|---|---|
| A course of wires | does not. |

```julia
A centered
circle or annulus fills one. A course of wires does not.
```
"""
f(x; provenance = 1) = x
template = """
A centered
circle or annulus fills one. A course of wires does not.
"""
JULIA
test -z "$(unwrapped blocks.jl)"

# A list item in Markdown, as GitHub renders it.
cat > "$unwrapped_dir/list.md" <<'MARKDOWN'
- **Cable modeling:** Detailed representation of solid and stranded conductors,
  with round wires or rectangular strands laid on helices.
MARKDOWN
if plain list.md | grep -Fq ColonUsage; then
    echo 'Vale matched the wrapped list item across its line break.' >&2
    exit 1
fi
unwrapped list.md | grep -Fq 'list.md:1:19:error:Project.ColonUsage:'

# Vale reads an input of another kind as it is, with the alerts of a plain run.
printf '<html><body><p>The value has provenance.</p></body></html>\n' > "$unwrapped_dir/page.html"
test "$(unwrapped page.html | cut -d: -f1-3,5)" = "$(plain page.html | cut -d: -f1-4)"

# The exit status is that of Vale.
printf 'The result matches.\n' > "$unwrapped_dir/clean.md"
(cd "$unwrapped_dir" && julia --startup-file=no "$root/.github/prose/vale-unwrapped.jl" \
    "${unwrapped_args[@]}" clean.md 2> /dev/null)
if (cd "$unwrapped_dir" && julia --startup-file=no "$root/.github/prose/vale-unwrapped.jl" \
    "${unwrapped_args[@]}" page.html > /dev/null 2>&1); then
    echo 'An error alert left the exit status at 0.' >&2
    exit 1
else
    test "$?" -eq 1
fi
# In a stream that joins both outputs, as the log of CI does, every alert comes before the
# summary.
for i in $(seq 1 3000); do printf 'It has provenance.\n\n'; done > "$unwrapped_dir/many.md"
(cd "$unwrapped_dir" && julia --startup-file=no "$root/.github/prose/vale-unwrapped.jl" \
    "${unwrapped_args[@]}" many.md > "$scratch/many.txt" 2>&1 || true)
test "$(tail -n 1 "$scratch/many.txt")" = '3000 errors, 0 warnings and 0 suggestions in 1 files.'
test "$(grep -cxE 'many.md:[0-9]+:[0-9]+:error:Project.BannedTerms:.*' "$scratch/many.txt")" -eq 3000
# The path of the temporary directory can hold a parent step. The alerts keep the path of
# each input and its columns.
(cd "$unwrapped_dir" && TMPDIR="$unwrapped_dir/.." julia --startup-file=no \
    "$root/.github/prose/vale-unwrapped.jl" "${unwrapped_args[@]}" docstring.jl 2> /dev/null || true) |
    grep -q '^docstring.jl:6:77:error:ai-tells.StackedAnaphora:'
grep -Fq 'julia --startup-file=no .github/prose/vale-unwrapped.jl --no-global --config=.github/prose/vale.ini --no-color "${files[@]}"' \
    "$root/.github/workflows/CI.yml"

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
mkdir -p "$scratch/repository/test/fixtures/reference" "$scratch/repository/docs/src/examples" "$scratch/repository/docs/build"
printf 'test/fixtures/reference/\n' > "$scratch/repository/.gitignore"
printf 'Read and recieve the result.\n' > "$scratch/repository/test/fixtures/reference/README.md"
printf 'License CONTRACT\n' > "$scratch/repository/LICENSE"
printf 'Generated prose.\n' > "$scratch/repository/docs/src/examples/example.md"
printf '<p>Generated prose.</p>\n' > "$scratch/repository/docs/build/index.html"
git -C "$scratch/repository" add .gitignore LICENSE
git -C "$scratch/repository" add --force test/fixtures/reference/README.md
(cd "$scratch/repository" && bash "$root/.github/prose/files.sh" "$scratch/selection")
grep -Fxq test/fixtures/reference/README.md "$scratch/selection/vale-files.txt"
grep -Fxq docs/src/examples/example.md "$scratch/selection/cspell-files.txt"
grep -Fxq docs/build/index.html "$scratch/selection/vale-files.txt"
if grep -Fxq LICENSE "$scratch/selection/vale-files.txt"; then
    echo 'Input selection listed LICENSE.' >&2
    exit 1
fi
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
