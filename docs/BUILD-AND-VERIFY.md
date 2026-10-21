# Build & Verify Runbook

Build both apps and prove they work — with observed evidence at every gate.

## One-time setup

```bash
npm install
tools/build-r-runtime.sh     # assembles shared/r/ from the system R (~cached)
```

`build-r-runtime.sh` also runs `fix-r-framework.sh`, which must finish with
`Remaining /Library/Frameworks R references: 0`. Never skip this check.

## Build an app

```bash
node tools/build-app.js plasmid-designer
node tools/build-app.js crispr-designer
```

Output: `dist/<id>/<App Name>-darwinwin-x64/<App Name>.app` (~1.0 GB).
There is no signing step — the bundle is intentionally unsigned (HCR model).

## Gate 1 — bundle content checks (seconds)

```bash
A='dist/crispr-designer/CRISPR Designer-darwin-x64/CRISPR Designer.app'
diff -q shared/electron/bootstrap.js "$A/Contents/Resources/app/shared/electron/bootstrap.js"
otool -L "$A/Contents/Resources/app/r/library/httpuv/libs/httpuv.so" | grep libR
# expect: @loader_path/...  (NOT /Library/Frameworks)
codesign -dv "$A" 2>&1 | head -1
# expect: "code object is not signed at all"
ls "$A/Contents/Resources/app/payload" "$A/Contents/Resources/app/preload.js"
```

## Gate 2 — R simulation (the exact runner spawn, no GUI)

Replicates what bootstrap.spawnR() does. Save as /tmp/rsim.sh:

```bash
R="<repo>"; RROOT="$R/shared/r"; RES="$RROOT/R.framework/Versions/4.6-x86_64/Resources"
env HOME="$HOME" PATH=/usr/bin:/bin:/usr/sbin:/sbin \
  R_HOME="$RES" R_SHARE_DIR="$RES/share" R_INCLUDE_DIR="$RES/include" R_DOC_DIR="$RES/doc" \
  R_LIBS="$RROOT/library:$RES/library" R_LIBS_USER="$RROOT/library:$RES/library" \
  R_ENVIRON_USER= R_PROFILE_USER= \
  R_DEFAULT_PACKAGES='methods,stats,graphics,grDevices,utils,datasets' \
  APP_SCRIPT='Cas9_sgRNA_design.R' \
  "$RES/bin/exec/R" --no-save --no-restore --no-environ --no-site-file --no-init-file --slave \
  -f "$R/shared/electron/runner.R" --args "$R/apps/crispr-designer/payload"
```

Run it in the background. **Pass criteria:** stdout contains
`APP_PORT=<n>` and `Listening on http://127.0.0.1:<n>`, then
`curl -s -o /dev/null -w '%{http_code}' http://127.0.0.1:<n>/` returns **200**.

If it segfaults at `dyn.load`, the package library tree is unpatched —
re-run `tools/fix-r-framework.sh` (see ARCHITECTURE.md for the mechanism).

## Gate 3 — real launch

`open "dist/<id>/…/<App Name>.app"`, then within ~30 s:

- the app process is alive (`ps ax | grep '<App Name>.app'`)
- an R child is alive (`ps ax | grep runner.R`, path inside the app bundle)
- `<userData>/main.log` shows `spawning R:` … then no fatal error
- the window shows the Shiny UI

Any failure text appears in the window AND in `<userData>/main.log`.

## Troubleshooting (all observed 2026-10-20)

| Symptom | Root cause | Fix |
|---|---|---|
| dyld crash "different Team IDs" / "Check with the developer" | ad-hoc signing with `--options runtime` enables library validation | remove all Electron-side signing (unsigned, HCR model) |
| R not spawned, "No bundled R runtime found" | resolver expected `Resources/r`, bundle has `Resources/app/r` | dual-layout `resolveRRoot()` |
| "App script not found" / payload missing | payload at `Resources/app/payload/` not resolved; `APP_SCRIPT` unset | `resolveAppDir()` payload branch + set `APP_SCRIPT` |
| R segfault at dyn.load loading shiny | second libR loaded from unpatched `r/library` `.so` deps | `fix-r-framework.sh` patches BOTH trees |
| "TypeError: startApp is not a function" | mangled bootstrap.js (module.exports in dead scope) | clean bootstrap.js rewrite |
| "CLI mode requires --gene" | payload's standalone dispatch fired when sourced by runner | `if (sys.nframe() == 0L)` guard |

## Testing

```bash
npm run test:r                    # payload testthat suites
node tools/verify-references.js   # reference-store parity check
```
