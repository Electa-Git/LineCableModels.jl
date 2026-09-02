#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
DESTINATION=${LCM_INSTALL_DIR:-"$HOME/.local/bin"}
mkdir -p "$DESTINATION"
ln -sfn "$ROOT/cli/lcm" "$DESTINATION/lcm"
printf 'Installed %s\n' "$DESTINATION/lcm"
