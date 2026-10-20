#!/usr/bin/env bash
# ============================================================================
# build-r-runtime.sh — assemble a self-contained R runtime for packaging.
#
# Produces shared/r/:
#   R.framework/Versions/<v>/Resources/…   ← base R + recommended packages
#   library/<pkg>                          ← the packages the apps need
#   janelia_atalanta_plasmids/             ← the bundled pJAT collection
#
# Only the base R and the apps' dependency closure are copied, which keeps
# the runtime to roughly 200 MB instead of the ~2 GB of a full R install.
#
# Usage:  tools/build-r-runtime.sh [--force]
# ============================================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="$REPO_ROOT/shared/r"
FRAMEWORK_SRC="/Library/Frameworks/R.framework"
FORCE=0
[ "${1:-}" = "--force" ] && FORCE=1

if [ ! -d "$FRAMEWORK_SRC" ]; then
  echo "ERROR: no system R framework at $FRAMEWORK_SRC" >&2
  exit 1
fi

if [ -d "$OUT_DIR/R.framework" ] && [ "$FORCE" -eq 0 ]; then
  echo "R runtime already present at $OUT_DIR (use --force to rebuild)"
  exit 0
fi

echo "==> Determining the R version and package closure"
# The framework's version directory is named like "4.6-x86_64", which is not
# the same as R.version.string ("4.6.0"). Use the directory name.
VERSION="$(ls -1 "$FRAMEWORK_SRC/Versions" | grep -E '^[0-9]' | sort -V | tail -1)"
if [ -z "$VERSION" ]; then
  echo "ERROR: could not determine the R framework version" >&2
  exit 1
fi
echo "    R framework version directory: $VERSION"

# The packages the two apps need, plus their recursive dependency closure.
PKGS="$(Rscript --vanilla "$REPO_ROOT/tools/r-closure.R" 2>/dev/null)"

PKG_COUNT="$(printf '%s\n' "$PKGS" | grep -c .)"
echo "    $PKG_COUNT packages in the closure"

echo "==> Copying the R framework (base + recommended only)"
rm -rf "$OUT_DIR/R.framework" "$OUT_DIR/R.framework.tmp_versions"
mkdir -p "$OUT_DIR"

# Rebuild the framework skeleton by hand so we only carry what we need.
FRAMEWORK_DST="$OUT_DIR/R.framework"
mkdir -p "$FRAMEWORK_DST/Versions"

VER_SRC="$FRAMEWORK_SRC/Versions/$VERSION"
VER_DST="$FRAMEWORK_DST/Versions/$VERSION"
mkdir -p "$VER_DST/Resources"

# "Current" must point at the version directory, exactly as the system
# framework does. (Pointing it at itself creates a symlink loop.)
ln -sfn "$VERSION" "$FRAMEWORK_DST/Versions/Current"

# The executable, its libraries and the headers are needed at runtime.
# `doc` and `tests` are documentation only and are skipped to save space.
for item in bin lib include share etc fontconfig modules; do
  if [ -e "$VER_SRC/Resources/$item" ]; then
    cp -R "$VER_SRC/Resources/$item" "$VER_DST/Resources/"
  fi
done
# R itself and the version marker live one level up.
for item in R Rscript COPYING Info.plist SVN-REVISION; do
  [ -e "$VER_SRC/$item" ] && cp -R "$VER_SRC/$item" "$VER_DST/" || true
done

echo "==> Copying base + recommended packages"
LIB_DST="$VER_DST/Resources/library"
mkdir -p "$LIB_DST"
for pkg in $(Rscript --vanilla -e '
  cat(rownames(installed.packages(priority = "base")), sep = "\n")
  cat(rownames(installed.packages(priority = "recommended")), sep = "\n")
'); do
  if [ -d "$VER_SRC/Resources/library/$pkg" ]; then
    cp -R "$VER_SRC/Resources/library/$pkg" "$LIB_DST/"
  fi
done

echo "==> Copying the app package closure"
LIB_APP="$OUT_DIR/library"
mkdir -p "$LIB_APP"
for pkg in $PKGS; do
  if [ -d "$VER_SRC/Resources/library/$pkg" ]; then
    cp -R "$VER_SRC/Resources/library/$pkg" "$LIB_APP/"
  else
    echo "    WARNING: package '$pkg' is not installed in the source library" >&2
  fi
done

# The app package tree must also be visible as the framework's own library,
# because the R runtime resolves R_LIBS against Resources/library. The
# packages are copied rather than symlinked: a symlink inside a signed
# application bundle breaks code-signing, and a relative symlink is easy to
# break when the bundle is moved.
for pkg in $PKGS; do
  if [ -d "$LIB_APP/$pkg" ]; then
    rm -rf "$LIB_DST/$pkg"
    cp -R "$LIB_APP/$pkg" "$LIB_DST/$pkg"
  fi
done

rm -rf "$OUT_DIR/R.framework.tmp_versions"

echo "==> Bundling the pJAT plasmid collection"
rm -rf "$OUT_DIR/janelia_atalanta_plasmids"
mkdir -p "$OUT_DIR/janelia_atalanta_plasmids"
if [ -d "$REPO_ROOT/shared/reference/janelia_atalanta_plasmids" ]; then
  cp "$REPO_ROOT/shared/reference/janelia_atalanta_plasmids"/*.gb* \
     "$OUT_DIR/janelia_atalanta_plasmids/" 2>/dev/null || true
  cp "$REPO_ROOT/shared/reference/janelia_atalanta_plasmids"/*.gbff \
     "$OUT_DIR/janelia_atalanta_plasmids/" 2>/dev/null || true
fi

echo "==> Verifying the runtime"
VERIFY_R="$VER_DST/Resources/bin/exec/R"
VERIFY_LIBS="$LIB_APP:$LIB_DST"
R_HOME="$VER_DST/Resources" \
R_LIBS="$VERIFY_LIBS" \
R_LIBS_USER="$LIB_APP" \
R_ENVIRON_USER='' R_PROFILE_USER='' \
  "$VERIFY_R" --no-save --no-restore --slave -e '
  pkgs <- c("shiny", "Biostrings", "data.table", "DT", "plasmidplot",
            "httpuv", "httr", "readr", "png", "yaml", "jsonlite", "testthat")
  bad <- character(0)
  for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-12s %s\n", p, if (ok) "ok" else "MISSING"))
    if (!ok) bad <- c(bad, p)
  }
  if (length(bad)) {
    stop("Verification failed for: ", paste(bad, collapse = ", "))
  }
  cat("All required packages load.\n")
'

SIZE="$(du -sh "$OUT_DIR" | cut -f1)"
echo "==> R runtime ready at $OUT_DIR ($SIZE)"
