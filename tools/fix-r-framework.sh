#!/bin/bash
set -euo pipefail
FW_RES="${1}"
LIB="$FW_RES/lib"

if [ ! -d "$LIB" ]; then
  echo "ERROR: R.framework Resources/lib directory not found at: $LIB" >&2
  exit 1
fi

# The separate package-library tree that sits beside the framework:
#   <r>/R.framework/Versions/<v>/Resources   (= $FW_RES)
#   <r>/library                              (= $PKGLIB, when present)
# Its package .so files still reference the system R paths. If they are left
# unpatched, dyld cannot match their libR dependency against the bundled
# libR (whose install id was rewritten), so it loads a SECOND copy of libR
# and R segfaults at dyn.load time. Patch them with the same rewrites.
PKGLIB="$FW_RES/../../../../library"

rel_to() {
  python3 - "$1" "$2" <<'PY'
import os, sys
print(os.path.relpath(sys.argv[2], sys.argv[1]))
PY
}

OLD_PREFIX="/Library/Frameworks/R.framework/Versions"

# Patch every .dylib/.so under $1: rewrite install ids and dependency load
# commands that point into the system R (/Library/Frameworks/...) to
# @loader_path-relative locations inside this bundle's lib/ directory.
# Prints the number of patched load-command paths.
patch_tree() {
  local root="$1"
  local count=0
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
      count=$((count + 1))
    fi

    while read -r dep; do
      [[ -z "$dep" ]] && continue
      if [[ "$dep" == "$OLD_PREFIX"* ]]; then
        base=$(basename "$dep")
        if [ -e "$LIB/$base" ] || [ -L "$LIB/$base" ]; then
          rel=$(rel_to "$(dirname "$f")" "$LIB")
          args+=( -change "$dep" "@loader_path/${rel}/${base}" )
          count=$((count + 1))
        else
          echo "  WARN: skipping dependency not in bundle lib/: $dep" >&2
        fi
      fi
    done < <(otool -L "$f" 2>/dev/null | sed -n '2,$p' | awk '{print $1}')

    if [ ${#args[@]} -gt 0 ]; then
      install_name_tool "${args[@]}" "$f" >/dev/null 2>&1 || echo "  FAIL: $f" >&2
    fi
  done < <(find "$root" -type f \( -name '*.dylib' -o -name '*.so' \) ! -path '*/.dSYM/*' ! -name '*.TBD' -print0 2>/dev/null)
  echo "$count"
}

# Ad-hoc re-sign every binary under $1 (install_name_tool invalidates the
# existing signature of every file it edits). Prints the signed file count.
sign_tree() {
  local root="$1"
  local count=0
  while IFS= read -r -d '' f; do
    case "$f" in
      *.dSYM*|*.TBD*) continue ;;
    esac
    if [[ "$f" == *.so || "$f" == *.dylib || -x "$f" ]]; then
      codesign --force --sign - "$f" >/dev/null 2>&1 || true
      count=$((count + 1))
    fi
  done < <(find "$root" -type f \( -name '*.so' -o -name '*.dylib' -o -perm +111 \) ! -path '*/.dSYM/*' -print0 2>/dev/null)
  echo "$count"
}

echo "Relocating R.framework: $FW_RES"

n=$(patch_tree "$FW_RES")
echo "Patched $n load-command paths in the framework."

if [ -d "$PKGLIB" ]; then
  n=$(patch_tree "$PKGLIB")
  echo "Patched $n load-command paths in the package library."
fi

# exec/R + exec/Rscript are plain executables — patch their deps separately.
for exe in "$FW_RES/bin/exec/R" "$FW_RES/bin/exec/Rscript"; do
  if [ -f "$exe" ]; then
    rel=$(rel_to "$(dirname "$exe")" "$LIB")
    while read -r dep; do
      [[ -z "$dep" ]] && continue
      if [[ "$dep" == "$OLD_PREFIX"* ]]; then
        base=$(basename "$dep")
        if [ -e "$LIB/$base" ] || [ -L "$LIB/$base" ]; then
          install_name_tool -change "$dep" "@loader_path/${rel}/${base}" "$exe" >/dev/null 2>&1 || echo "  FAIL: $exe" >&2
        fi
      fi
    done < <(otool -L "$exe" 2>/dev/null | sed -n '2,$p' | awk '{print $1}')
  fi
done

# re-sign
n=$(sign_tree "$FW_RES/..")
echo "Re-signed $n framework binaries (ad hoc)."
if [ -d "$PKGLIB" ]; then
  n=$(sign_tree "$PKGLIB")
  echo "Re-signed $n package-library binaries (ad hoc)."
fi

# sanity check: no /Library/Frameworks load commands should remain anywhere.
leftover=0
for root in "$FW_RES" "$PKGLIB"; do
  [ -d "$root" ] || continue
  c=$(find "$root" -type f \( -name '*.dylib' -o -name '*.so' \) ! -path '*/.dSYM/*' -print0 2>/dev/null |
    xargs -0 otool -L 2>/dev/null | grep -c "/Library/Frameworks/R.framework" || true)
  leftover=$((leftover + c))
done
echo "Remaining /Library/Frameworks R references: $leftover"

if [ "$leftover" -gt 0 ]; then
  echo "WARNING: some binaries still reference the system R.framework" >&2
else
  echo "OK: R runtime is fully relocatable."
fi
