#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SOURCE=${POLYPAVES_SOURCE:-$ROOT/../paves}
VERSION=${POLYPAVES_VERSION:-0.4.0}
PYTHON=${PYTHON:-python3}
ABI=$($PYTHON -c 'import sys; print(f"cp{sys.version_info.major}{sys.version_info.minor}")')
STAGE=$(mktemp -d)
trap 'rm -rf -- "$STAGE"' EXIT

PACKAGE=$SOURCE/python/polypaves
API=$PACKAGE/api.py
# A native extension built with -ffile-prefix-map, so no build path survives in it.
NATIVE_DIR=${POLYPAVES_NATIVE_DIR:-$SOURCE/build/public-wheel${ABI#cp}}
EXT_SUFFIX=$($PYTHON -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')

[[ -f "$API" ]] || {
    echo "polypaves API source not found: $API" >&2
    exit 2
}

if ! grep -q '^POLYPAVES_AMOEBA09:BOOL=OFF$' "$NATIVE_DIR/CMakeCache.txt"; then
    echo 'Public native builds must use -DPOLYPAVES_AMOEBA09=OFF' >&2
    exit 2
fi

NATIVE=$(find "$NATIVE_DIR" -type f -name "_native$EXT_SUFFIX" -print -quit)
[[ -n "$NATIVE" ]] || {
    echo "built polypaves native extension not found below: $NATIVE_DIR" >&2
    exit 2
}

mkdir -p "$STAGE/polypaves"
cp "$ROOT/packaging/pyproject.toml" "$ROOT/packaging/setup.py" "$STAGE/"
cp "$ROOT/packaging/polypaves_init.py" "$STAGE/polypaves/__init__.py"
# Copy modules recursively so reaction's implementation is also compiled. Only
# package initializers remain as public facades; the wheel audit checks this.
while IFS= read -r -d '' module; do
    relative=${module#"$PACKAGE/"}
    [[ $relative == __init__.py ]] && continue
    mkdir -p "$STAGE/polypaves/$(dirname "$relative")"
    cp "$module" "$STAGE/polypaves/$relative"
done < <(find "$PACKAGE" -type f -name '*.py' -print0)
cp "$NATIVE" "$STAGE/polypaves/"
strip --strip-unneeded "$STAGE/polypaves/$(basename "$NATIVE")"
cp "$ROOT/LICENSE" "$STAGE/LICENSE"

mkdir -p "$ROOT/dist"
POLYPAVES_VERSION=$VERSION "$PYTHON" -m pip wheel \
    --no-deps --no-cache-dir --wheel-dir "$ROOT/dist" "$STAGE"

WHEEL=$(find "$ROOT/dist" -maxdepth 1 -type f -name "polypaves-${VERSION}-${ABI}-${ABI}-*.whl" -print -quit)
[[ -n "$WHEEL" ]] || {
    echo "wheel was not created" >&2
    exit 2
}

"$ROOT/scripts/audit_binary_wheel.sh" "$WHEEL"
echo "$WHEEL"
