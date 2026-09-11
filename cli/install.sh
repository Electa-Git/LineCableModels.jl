#!/usr/bin/env bash
set -euo pipefail

lcm_installer=$(readlink -f -- "${BASH_SOURCE[0]}")
lcm_root=$(cd -- "$(dirname -- "$lcm_installer")/.." && pwd -P)
lcm_bin=${LCM_INSTALL_DIR:-$HOME/.local/bin}
lcm_applications=${LCM_APPLICATIONS_DIR:-${XDG_CONFIG_HOME:-$HOME/.config}/lcm/applications}
lcm_names=()
lcm_directories=()

while (( $# )); do
    case "$1" in
        --help|-h)
            cat <<'EOF'
Usage: cli/install.sh [--bin-dir DIRECTORY] [--application NAME DIRECTORY]...

Install this checkout's cli/lcm as the global lcm command.
Optionally locate an application in another checkout using a directory symlink.
The directory must contain its executable lcm launcher. No dependencies are installed.

Existing symlinks are updated; regular files and directories are never replaced.
Defaults: ~/.local/bin and ${XDG_CONFIG_HOME:-~/.config}/lcm/applications.
LCM_APPLICATIONS_DIR and LCM_INSTALL_DIR may override those locations.
EOF
            exit 0
            ;;
        --bin-dir)
            if (( $# < 2 )) || [[ -z "$2" ]]; then
                printf 'install: --bin-dir requires a directory\n' >&2
                exit 2
            fi
            lcm_bin=$2
            shift 2
            ;;
        --application)
            if (( $# < 3 )) || [[ ! "$2" =~ ^[a-z][a-z0-9_-]*$ ]]; then
                printf 'install: --application requires NAME DIRECTORY\n' >&2
                exit 2
            fi
            case "$2" in
                cli|help|runtime|worker|presentation|nats|container|demo)
                    printf 'install: %s is not an application directory name; playground owns its related commands\n' "$2" >&2
                    exit 2
                    ;;
            esac
            if [[ ! -x "$3/lcm" || -d "$3/lcm" ]]; then
                printf 'install: application launcher is missing or not executable: %s/lcm\n' "$3" >&2
                exit 2
            fi
            if [[ " ${lcm_names[*]} " == *" $2 "* ]]; then
                printf 'install: repeated application: %s\n' "$2" >&2
                exit 2
            fi
            lcm_names+=("$2")
            lcm_directories+=("$(cd -- "$3" && pwd -P)")
            shift 3
            ;;
        *)
            printf 'install: unknown argument: %s\n' "$1" >&2
            exit 2
            ;;
    esac
done

lcm_public="$lcm_bin/lcm"
lcm_destinations=("$lcm_public")
for lcm_name in "${lcm_names[@]}"; do
    lcm_destinations+=("$lcm_applications/$lcm_name")
done
for lcm_destination in "${lcm_destinations[@]}"; do
    if [[ -e "$lcm_destination" && ! -L "$lcm_destination" ]]; then
        printf 'install: refusing to replace non-symlink path: %s\n' "$lcm_destination" >&2
        exit 1
    fi
done

mkdir -p -- "$lcm_bin"
if (( ${#lcm_names[@]} )); then
    mkdir -p -- "$lcm_applications"
fi
for lcm_index in "${!lcm_names[@]}"; do
    lcm_destination="$lcm_applications/${lcm_names[$lcm_index]}"
    ln -sfn -- "${lcm_directories[$lcm_index]}" "$lcm_destination"
    printf '%s -> %s\n' "$lcm_destination" "${lcm_directories[$lcm_index]}"
done
ln -sfn -- "$lcm_root/cli/lcm" "$lcm_public"
printf '%s -> %s\n' "$lcm_public" "$lcm_root/cli/lcm"
case ":$PATH:" in
    *":$lcm_bin:"*) ;;
    *) printf 'Add %s to PATH to invoke lcm from any directory.\n' "$lcm_bin" ;;
esac
