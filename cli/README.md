# Package command line

`cli/lcm` is the package entrypoint. It selects an application and executes its
launcher; each application owns its arguments, environment and process lifecycle.
The router preserves the caller's working directory, argument boundaries, exit
status and signals. Global help, version and path inspection start no Julia process.

Install the command from the checkout that should own the global entrypoint:

```sh
./cli/install.sh
lcm --help
lcm --version
lcm --paths
lcm gauntlet --help
lcm playground --help
```

Installation creates `~/.local/bin/lcm -> CHECKOUT/cli/lcm`. It updates an existing
symlink, including the old playground link, and refuses to replace a regular file
or directory. It neither installs application dependencies nor starts applications.
Keep `~/.local/bin` on `PATH`. `--bin-dir DIRECTORY` or `LCM_INSTALL_DIR` selects another
installation directory.

## Applications in separate worktrees

Applications in the launcher's own checkout take precedence. An application that
is absent there can be explicitly located in another worktree:

```sh
./cli/install.sh \
  --application playground /path/to/LineCableModels-playground/playground

lcm --paths
lcm gauntlet run --definition /path/to/benchmark.jl --directory /path/to/drafts
lcm playground start
lcm runtime --help
```

The installer records ordinary directory symlinks under
`${XDG_CONFIG_HOME:-~/.config}/lcm/applications`. `LCM_APPLICATIONS_DIR` selects a
different directory for both installation and invocation. No configuration file is
executed, no Git branches are searched, and the current working directory does not
change application selection. Relative input/output arguments retain their meaning
at the caller's location. A missing application fails with exit status 2.

`lcm --paths` shows the actual package and application directories. After the
applications are merged into the selected checkout, its local directories take
precedence and the corresponding external symlinks can be removed. Run the installer
from another checkout to deliberately change the global entrypoint.

Playground also owns the existing `runtime`, `worker`, `presentation`, `nats`,
`container` and `demo` commands. The router passes those commands to its launcher
unchanged. Its bootstrap installs application dependencies; global command
installation belongs to `cli/install.sh`.

## Application launcher contract

An application supplies an executable `APPLICATION/lcm`. Its name uses lowercase
letters, digits, `_` or `-`, beginning with a letter. `cli`, `help` and the existing
playground command names are reserved. A new application directory is sufficient
for ordinary routing and help discovery; it needs no Julia registry or root switch.

The launcher receives the complete command vector, including the application name:
`lcm gauntlet run ...` executes `gauntlet/lcm gauntlet run ...`. Gauntlet consumes its
namespace and starts its Julia CLI with its own project. Playground retains its
existing command parser and server supervision. `LCM_JULIA` selects the Julia
executable used by either application; it is one executable path, not shell code.

The package router owns only command selection and process handoff. The installer
owns the global symlink and explicit application locations, using ordinary shell
filesystem operations. Neither duplicates application parsing or dependency setup.

Test these boundaries without numerical packages or running services:

```sh
python3 -m unittest discover -s test/cli -v
```
