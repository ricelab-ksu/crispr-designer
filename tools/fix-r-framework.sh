#!/bin/bash
set -euo pipefail
FW_RES="${1}"
LIB="$FW_RES/lib"

if [ ! -d "$LIB" ]; then
  echo "ERROR: R.framework Resources/lib directory not found at: $LIB" >&2
  exit 1
fi

RESOURCES_DIR="$FW_RES/.."
FW="$RESOURCES_DIR"

if [ ! -d "$FW" ]; then
  echo "ERROR: R.framework root not found at: $FW" >&2
  exit 1
fi

rel_to() {
  python3 - "$1" "$2" <<'PY'
import os, sys
print(os.path.relpath(sys.argv[2], sys.argv[1]))
PY
}

echo "Relocating R.framework: $FW"

OLD_PREFIX="/Library/Frameworks/R.framework/Versions"
patch_count=0
sign_count=0

# rewrite dylib ids
while IFS= read -r -d '' f; do
  case "$f" in
    *.dSYM*|*.TBD*) continue ;;
  esac

  idline=""
  if [[ "$f" == *.dylib ]]; then
    idline=$(otool -D "$f" 2>/dev/null | sed -n 2p || true)
  fi

  args=()

  if [[ -n "$idline" && "$idline" == "$OLD_PREFIX"* ]]; then
    rel=$(rel_to "$(dirname "$f")" "$LIB")
    args+=( -id "@loader_path/${rel}/$(basename "$idline")" )
    patch_count=$((patch_count + 1))
  fi

  while read -r dep; do
    [[ -z "$dep" ]] && continue
    if [[ "$dep" == "$OLD_PREFIX"* ]]; then
      base=$(basename "$dep")
      if [ -e "$LIB/$base" ] || [ -L "$LIB/$base" ]; then
        rel=$(rel_to "$(dirname "$f")" "$LIB")
        args+=( -change "$dep" "@loader_path/${rel}/${base}" )
        patch_count=$((patch_count + 1))
      else
        echo "  WARN: skipping dependency not in bundle lib/: $dep" >&2
      fi
    fi
  done < <(otool -L "$f" 2>/dev/null | sed -n '2,$p' | awk '{print $1}')

  if [ ${#args[@]} -gt 0 ]; then
    install_name_tool "${args[@]}" "$f" >/dev/null 2>&1 || echo "  FAIL: $f" >&2
  fi
done < <(find "$FW/Resources" -type f \( -name '*.dylib' -o -name '*.so' \) ! -path '*/.dSYM/*' ! -name '*.TBD' -print0 2>/dev/null)

# exec/R rewrite
for exe in "$FW/Resources/bin/exec/R" "$FW/Resources/bin/exec/Rscript"; do
  if [ -f "$exe" ]; then
    rel=$(rel_to "$(dirname "$exe")" "$LIB")
    while read -r dep; do
      [[ -z "$dep" ]] && continue
      if [[ "$dep" == "$OLD_PREFIX"* ]]; then
        base=$(basename "$dep")
        if [ -e "$LIB/$base" ] || [ -L "$LIB/$base" ]; then
          install_name_tool -change "$dep" "@loader_path/${rel}/${base}" "$exe" >/dev/null 2>&1 || echo "  FAIL: $exe" >&2
          patch_count=$((patch_count + 1))
        fi
      fi
    done < <(otool -L "$exe" 2>/dev/null | sed -n '2,$p' | awk '{print $1}')
  fi
done

echo "Patched $patch_count load-command paths."

# re-sign
while IFS= read -r -d '' f; do
  case "$f" in
    *.dSYM*|*.TBD*) continue ;;
  esac
  if [[ "$f" == *.so || "$f" == *.dylib || -x "$f" ]]; then
    codesign --force --sign - "$f" >/dev/null 2>&1 || true
    sign_count=$((sign_count + 1))
  fi
done < <(find "$FW" -type f \( -name '*.so' -o -name '*.dylib' -o -perm +111 \) ! -path '*/.dSYM/*' -print0 2>/dev/null)

echo "Re-signed $sign_count binaries (ad hoc)."

# Sign the R.framework itself last, so its seal is computed after its contents are finalized
if [ -d "$FW/Versions/A" ]; then
  codesign --force --sign - "$FW/Versions/A" >/dev/null 2>&1 || true
  sign_count=$((sign_count + 1))
fi
if [ -d "$FW/Versions/Current" ]; then
  codesign --force --sign - "$FW/Versions/Current" >/dev/null 2>&1 || true
  sign_count=$((sign_count + 1))
fi

# sanity check
leftover=$(find "$FW/Resources" -type f \( -name '*.dylib' -o -name '*.so' \) ! -path '*/.dSYM/*' -print0 2>/dev/null |
  xargs -0 otool -L 2>/dev/null | grep -c "/Library/Frameworks/R.framework" || true)
echo "Remaining /Library/Frameworks R references: $leftover"

if [ "$leftover" -gt 0 ]; then
  echo "WARNING: some binaries still reference the system R.framework" >&2
else
  echo "OK: R.framework is fully relocatable."
fi
