# Architecture

How the two desktop apps (Plasmid Designer, CRISPR Designer) work end to end.

Reference model: [HCR Probe Designer](https://github.com/ricelab-ksu/HCR-probe-design)
— our runtime contract and signing posture deliberately match it.

---

## Launch chain

```
App.app
 └─ Contents/MacOS/<name>            Electron main
     └─ Resources/app/main.js        thin per-app config
         └─ shared/electron/bootstrap.js   startApp(config)
             1. splash window (data: URL, shows fatal errors in-page)
             2. resolveAppDir()  → Resources/app/payload
             3. resolveRRoot()   → Resources/app/r   (both layouts supported)
             4. spawnR():
                  exec/R --no-save --no-restore --no-environ --no-site-file
                         --no-init-file --slave
                         -f shared/electron/runner.R --args <payload dir>
                  env: R_HOME, R_LIBS (both library trees), APP_SCRIPT
             5. runner.R:
                  - prepends the bundled library to .libPaths()
                  - sources <payload>/<APP_SCRIPT>   ← payload only DEFINES ui/server
                  - httpuv::randomPort(30000, 60000)
                  - prints "APP_PORT=<n>" to stdout  ← handshake
                  - shiny::runApp(host 127.0.0.1, port, launch.browser=FALSE)
             6. bootstrap watches R stdout for APP_PORT=, then loads
                http://127.0.0.1:<port> in the window (with connection retries)
```

The payload R script must NOT call `runApp`/`shinyApp` at top level when
sourced. Both payloads guard their standalone dispatch with
`if (sys.nframe() == 0L)` — true only for `Rscript` / `R -f` execution,
false when `source()`d by the runner. (Plasmid_Designer.R learned this the
hard way; see docs/DEBUG-LOG-2026-10-20.md.)

## Bundle layouts

Two layouts are supported by every resolver (`r-runtime.js`):

| Layout | R location | Payload | Produced by |
|---|---|---|---|
| packager (current) | `Resources/app/r/` | `Resources/app/payload/` | `@electron/packager` via tools/build-app.js |
| HCR / electron-builder | `Resources/r/` | `Resources/app/` (flat) | planned Phase 3 restructure |

## R runtime (`shared/r/`)

Two trees, both required:

- `R.framework/` — R 4.6 base + recommended packages
- `library/` — the app package closure (shiny stack, Bioconductor, …)

### Relocation: why BOTH trees must be patched

The CRAN R.framework and every CRAN binary package `.so` reference
`/Library/Frameworks/R.framework/Versions/<v>/Resources/lib/libR.dylib`
by absolute path. `tools/fix-r-framework.sh` rewrites those load commands
to `@loader_path`-relative paths and re-signs each edited file ad-hoc.

If a tree is left unpatched, dyld cannot match its libR dependency against
the bundled (rewritten) libR, so it loads a **second copy of libR** into the
process — R then **segfaults at dyn.load()** the moment shiny's dependencies
load (observed: `requireNamespace("shiny")` → SIGSEGV). The script therefore
patches and re-signs both the framework tree and the sibling `library/`
tree, and its sanity check must report **0 remaining
`/Library/Frameworks` references** across both.

macOS dyld dedups dependencies by install name / resolved path — this is
also why the reference HCR app works with unpatched package binaries: its
bundled libR *kept* the system install name, so dyld matches instead of
double-loading.

## Signing model (matches HCR exactly)

- **Electron side: completely unsigned.** No ad-hoc signing, no hardened
  runtime. `codesign -dv` reports "code object is not signed at all" —
  identical to the HCR Probe Designer.app. Do NOT add `--options runtime`
  ad-hoc signing: the hardened-runtime flag enables dyld library validation,
  which rejects the unsigned Electron frameworks ("different Team IDs") and
  makes the app unlaunchable.
- **R side: per-file ad-hoc only**, applied by `fix-r-framework.sh` because
  `install_name_tool` invalidates CRAN's original signatures. No runtime flag.

Distributed copies need right-click → Open on first launch (Gatekeeper),
same as HCR's released DMGs.

## Diagnostics

The bootstrap mirrors all main-process console output (including R stderr,
which it prefixes) to `<userData>/main.log` — e.g.
`~/Library/Application Support/CRISPR Designer/main.log`. Packaged-app
stdout is otherwise detached (LaunchServices), so this file is the first
place to look when a launch fails.

## Reference data

See `shared/electron/resolve-references.js` and the README. Genome store:
shared with HCR (`~/Library/Application Support/hcr-probe-designer/
reference_files`), resolved read-only. pJAT plasmid collection: bundled at
`<r root>/janelia_atalanta_plasmids`.

## Payload copy model

`apps/<id>/payload/` files are **adapted copies** of the originals in
`/Volumes/Backup_Plus/projects/bioinformatics/projects/<project>/`, with:
portable-path preamble (loads `R/app_paths.R`), runner-compatible dispatch
(`sys.nframe()` guard), and no launch scripts. Changes to the originals
must be re-copied AND re-adapted; changes made in the electron copies
should be ported back when they fix generally-applicable bugs
(example: the Plasmid `sys.nframe()` guard).
