# Development of the package

This page records the way LineCableModels 0.1.0 was built and the way further work
continues. It says who decides what and how a change moves from a specification to a commit.
It also lists the checks every change passes and where the records are kept. A contributor
who follows it works the way the project works.

## Roles

The release had four roles. The author decides. The others propose and execute.

| Role | Decides | Produces | Never |
|---|---|---|---|
| the author | scope, physics, names, rules, every commit, every decision | decisions as dated quotes | nothing |
| the architect | nothing, and recommends | the specifications, the audits at milestones, the amendments, the questions for the author | edits the tree or writes an executor prompt, launches a process or states a decision as the author's without a dated quote |
| the coordinator | how a step is cut and run, what passes review | the step prompts, the review, the report, the commit on confirmation, the state file | designs a type, a slot or a formula the specification does not give, or accepts an executor deviation, or commits without the author's word, or reopens a decision |
| the executor | nothing | the files a prompt names, a report of at most 25 lines | proposes or refines anything. It never adds a test, a name, a check or a file the prompt does not name |

In the 0.1.0 release the architect and the coordinator were two separate language-model
threads, and the executor was a third tool, used under strict bounds. The executor is an
available tool, paid for, and it executes: it never designs or decides. Prototypes are the
executor's too. Prototype code
is never written or run by the coordinator. The architect specifies and audits, and measures
before specifying with scratch runs that never enter the tree.

A decision is taken once, by the role that owns it, and stated once. Only the decisions that
belong to the author go to the author. Each fact is written to a file in the turn it
appears. A conversation is never the only record of a fact.

## The protocol of a work package

1. **Specification.** The architect writes the specification of the work package. It lists
   every type, slot, method signature, file and test, the sources to port with file and
   line, and the success criteria. A specification cites the author's decision for
   everything it rests on.
2. **Decisions.** The coordinator takes the open decisions to the author, at most five at a
   time. The author decides. The coordinator records the dated quote. The architect adds the
   decision to the specification. No step starts before its decisions are recorded.
3. **Steps.** For each step, the coordinator writes the prompt and launches the executor. It
   reviews the result against the specification and the checklist below. Its report to the
   author gives the executor's report verbatim, the verdict with evidence and the commit
   block. The author confirms, and the coordinator commits. The state file is rewritten.
4. **Audit.** At each milestone the architect reads the committed tree against the
   specification and the physics sources. The audit says accepted, or lists numbered
   amendments. Each amendment gives the file, the line, what is wrong and what the
   specification says. Amendments become a step applied by a second executor. The work
   package closes when the audit accepts.
5. **Records.** Every fact, decision, number and in-flight process is written in the file
   that owns it, in the turn it appears.

Without waiting for a milestone, the architect receives two things. One is a step that needs
a type, slot, field or method the specification does not give. The other is a target missed
beyond its tolerance. Everything that is the author's goes to the author, such as a name, a
physics choice or a scope question.

## The prompt modes

A step runs in one of two modes.

- **Patch mode**, for structure and wiring. The prototyping executor produces a patch on an
  export of the tree. The coordinator reviews it, and a second executor applies it in the
  worktree, formats, runs the checks and reports.
- **Implementation mode**, for self-contained mathematics with targets. The prompt gives the
  signatures, the registered names, the file to create and the source lines to port. It
  gives a test file whose targets the prototyping executor computed from the physics source.
  A second executor, which never sees the prototype's implementation, writes the code until
  the test passes. The executor that implements never generates a target.

A prompt never contains a request for a design, a name, a study or a proposal. When
something does not fit, the executor stops and reports.

## Questions for the author

Before a question goes to the author, it passes five checks.

1. Is it already decided? Then act on the recorded decision and ask nothing.
2. Does any computation read the thing? If only presentation reads it, its fate is not a
   decision.
3. Is the premise in this layer? A convenience keeps the scope the author gave it.
4. No option keeps a superficial result by adding machinery.
5. The brief states what the thing is, with file and line, who made it and why, and who
   reads it. It gives each option as example calls, with what happens inside and its effect
   on the public and private interfaces. It ends with a recommendation derived from a quoted
   decision of the author.

## The review checklist

Tests first. Each test calls feature code and checks a feature result. The review rejects an
assertion of what cannot fail and a check that a constructor or `validate` already blocks.
It rejects a defensive branch and a test input that does not change a checked result. An
executor that reports thousands of passing assertions has written tests that do not check
anything.

Then the diff. Only the files named. No private function, sentinel, forwarder, `Val`-first
method, decision by absence, Symbol switch or `isa` outside the admitted forms. Every name
in the register with the author's confirmation.

Then the numbers. The checks rerun independently in a capped process. The achieved values
compared with the targets and recorded.

Then the records. `git status --porcelain` equals the paths of the step. The commit paths
dry run in a scratch index. The message passes gitlint. The formatter is idempotent. Vale
and cspell pass as the CI job runs them.

Rejected on sight: assertion-of-assertion tests, imaginary scenarios, refinement loops,
anything outside the files of the prompt, an executor deviation accepted as known.

## The checks of every change

Every change passes these before it is reported.

- `julia --project=. -e 'using Pkg; Pkg.test()'`, with every guard at zero. The
  finite-element tests run when `LINECABLEMODELS_FEM` is set.
- The formatter: `using JuliaFormatter; format(".")` in the SciML style, idempotent.
- Vale 3.23.0 and cspell 10.3.6 on every maintained file, as `.github/prose/files.sh`
  lists them: zero errors, zero unknown words. A finding is fixed by rewriting the sentence.
- gitlint on the commit message: `type(scope): description`, an imperative subject of at
  most 72 characters, without an attribution line.

A commit is one step, with exactly the paths of the step. The author confirms each commit by
name. An agent never pushes, merges or publishes.

## Continuous integration

- `CI` runs the whole suite on Linux, macOS and Windows with Julia 1.12. It uploads coverage
  to Codecov and runs the writing checks on the sources and the rendered pages.
- `FEM` runs the suite with the finite-element cases on pull requests that are not drafts.
- Codecov requires 95 % of `src/` and `ext/` on pull requests, over the merged report of the
  four jobs. The number comes only from tests that check feature results, and the coverage
  report has an empty ignore list. An audit of the whole test tree for gamed tests preceded
  the release, and every later test is read against the same rule.
- `CompatHelper` proposes updates of the compat bounds. `TagBot` tags registered versions.

## Records

The working records live under the untracked folder `local/` of the author's checkout. The
state file lists the goal and the work packages with the state of each step. It says what is
next, what is in flight and which follow-ups are open. It records the way a question goes to
the author and every decision as a dated quote. The architect's folder keeps the
specifications, one audit per milestone, and the answers to the coordinator's briefs. A
restart prompt for each role is rewritten at every handoff, so that a new thread resumes
from the files.

The repository keeps what the records produced, namely the code, the tests, the register,
the frozen references with their README and this documentation.

## For a contributor

1. Use Julia 1.12. Clone the repository and run `Pkg.test()` from the package directory.
   Set `LINECABLEMODELS_FEM=true` to run the finite-element cases, which download GetDP on
   first use.
2. Read [Architecture](architecture.md) and [Rules and patterns](rules.md) before changing
   anything. A change that adds a name needs the author's confirmation of the name and its
   meaning, and the matching entries in `test/architecture/names.toml`.
3. Format with JuliaFormatter in the SciML style. Install gitlint and its hook once per
   clone: `uv tool install gitlint-core` and `gitlint install-hook`.
4. Run the writing checks with Vale 3.23.0 and cspell 10.3.6, configured in
   `.github/prose/`, on every file the change touches. Never add a disable directive.
5. Keep a pull request to one step. Add tests that check feature results. Update the
   documentation when the public API changes.
