# Repository working rules

## Protected showcase branch

Before changing files, dependencies, generated artifacts, or Git state in this
working tree, run:

```bash
git branch --show-current
```

Proceed with repository changes only when the output is exactly:

```text
feat/codespaces-showcase
```

If another branch is checked out, inspect the worktree for uncommitted changes,
preserve them, and switch to `feat/codespaces-showcase` before changing any
other file. Do not rely on the user to switch branches manually. If Git cannot
carry the existing worktree changes across safely, stop before making any other
change and report the conflict. Read-only inspection needed to establish the
branch and worktree state is allowed.
