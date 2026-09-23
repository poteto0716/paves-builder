#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
WHEEL=${1:-}

if [[ -z "$WHEEL" ]]; then
    WHEEL=$(find "$ROOT/dist" -maxdepth 1 -type f -name 'polyse-*.whl' -print -quit)
fi
[[ -n "$WHEEL" && -f "$WHEEL" ]] || {
    echo "polyse wheel not found: ${WHEEL:-<none>}" >&2
    exit 2
}

for command in unzip readelf nm strings file; do
    command -v "$command" >/dev/null || {
        echo "required audit command not found: $command" >&2
        exit 2
    }
done

STAGE=$(mktemp -d)
trap 'rm -rf -- "$STAGE"' EXIT
unzip -q "$WHEEL" -d "$STAGE"

mapfile -t MEMBERS < <(unzip -Z1 "$WHEEL")
for member in "${MEMBERS[@]}"; do
    case "$member" in
        polyse/__init__.py|*.so|*.dist-info/*) ;;
        *.py|*.pyi|*.pyx|*.pxd|*.c|*.cc|*.cpp|*.h|*.hh|*.hpp|*.o|*.a|*.debug)
            echo "source or development artifact found in wheel: $member" >&2
            exit 1
            ;;
    esac
done

mapfile -t EXTENSIONS < <(find "$STAGE/polyse" -maxdepth 1 -type f -name '*.so' -print | sort)
[[ ${#EXTENSIONS[@]} -ge 2 ]] || {
    echo "expected the native engine and compiled modules, found ${#EXTENSIONS[@]} extensions" >&2
    exit 1
}
[[ -n $(find "$STAGE/polyse" -maxdepth 1 -name '_native*.so' -print -quit) ]] || {
    echo "native engine extension missing" >&2
    exit 1
}

for extension in "${EXTENSIONS[@]}"; do
    file "$extension" | grep -q ', stripped' || {
        echo "native extension is not stripped: $(basename "$extension")" >&2
        exit 1
    }
    if readelf -W -S "$extension" | grep -Eq '\.(debug|symtab)([._[:space:]]|$)'; then
        echo "debug or static symbol section found: $(basename "$extension")" >&2
        exit 1
    fi
    if readelf -W -d "$extension" | grep -Eq 'RPATH|RUNPATH'; then
        echo "RPATH or RUNPATH found: $(basename "$extension")" >&2
        exit 1
    fi
    if strings -a "$extension" | grep -Eq '/home/|/Users/|/private/var/|/tmp/[^[:space:]]*(polyse|pip-)'; then
        echo "private build path found: $(basename "$extension")" >&2
        exit 1
    fi

    module=$(basename "$extension")
    module=${module%%.*}
    mapfile -t EXPORTS < <(nm -D --defined-only "$extension" | awk '{print $3}')
    [[ ${#EXPORTS[@]} -eq 1 && ${EXPORTS[0]} == "PyInit_${module}" ]] || {
        echo "unexpected exported symbols in $(basename "$extension")" >&2
        printf '  %s\n' "${EXPORTS[@]}" >&2
        exit 1
    }
done

echo "PASS: $(basename "$WHEEL") contains only the public Python facade and hidden, stripped extensions"
