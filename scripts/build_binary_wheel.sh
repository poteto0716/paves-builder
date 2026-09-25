#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SOURCE=${POLYPAVES_SOURCE:-$ROOT/../paves}
VERSION=${POLYPAVES_VERSION:-0.3.0}
PYTHON=${PYTHON:-python3}
STAGE=$(mktemp -d)
trap 'rm -rf -- "$STAGE"' EXIT

PACKAGE=$SOURCE/python/polypaves
API=$PACKAGE/api.py
# A native extension built with -ffile-prefix-map, so no build path survives in it.
NATIVE_DIR=${POLYPAVES_NATIVE_DIR:-$SOURCE/build/python-wheel313}
EXT_SUFFIX=$($PYTHON -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')

[[ -f "$API" ]] || {
    echo "polypaves API source not found: $API" >&2
    exit 2
}

NATIVE=$(find "$NATIVE_DIR" -type f -name "_native$EXT_SUFFIX" -print -quit)
[[ -n "$NATIVE" ]] || {
    echo "built polypaves native extension not found below: $NATIVE_DIR" >&2
    exit 2
}

mkdir -p "$STAGE/polypaves"
cp "$ROOT/packaging/pyproject.toml" "$ROOT/packaging/setup.py" "$STAGE/"
cp "$ROOT/packaging/polypaves_init.py" "$STAGE/polypaves/__init__.py"
for module in "$PACKAGE"/*.py; do
    [[ $(basename "$module") == __init__.py ]] && continue
    cp "$module" "$STAGE/polypaves/"
done
cp "$NATIVE" "$STAGE/polypaves/"
strip --strip-unneeded "$STAGE/polypaves/$(basename "$NATIVE")"
cp "$ROOT/LICENSE" "$STAGE/LICENSE"

mkdir -p "$ROOT/dist"
POLYPAVES_VERSION=$VERSION "$PYTHON" -m pip wheel \
    --no-deps --no-cache-dir --wheel-dir "$ROOT/dist" "$STAGE"

WHEEL=$(find "$ROOT/dist" -maxdepth 1 -type f -name "polypaves-${VERSION}-*.whl" -print -quit)
[[ -n "$WHEEL" ]] || {
    echo "wheel was not created" >&2
    exit 2
}

"$ROOT/scripts/audit_binary_wheel.sh" "$WHEEL"
echo "$WHEEL"
