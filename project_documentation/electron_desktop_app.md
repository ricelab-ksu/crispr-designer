# Electron Desktop App — CRISPR Designer (Cas9_sgRNA_design)

Status record for the standalone macOS desktop build of this tool.

**Status (2026-10-20): WORKING.** The packaged app launches to its Shiny UI,
verified end to end (R sidecar, port handshake, HTTP 200 on the Shiny server,
and a user-confirmed UI launch). This was the first of the two apps to work
after the launch-fix session.

## Where things live

| What | Where |
|---|---|
| Electron repo (build system + adapted payload) | `/Users/work/Desktop/projects/bioinformatics-electron` |
| Packaged app | `dist/crispr-designer/CRISPR Designer-darwin-x64/CRISPR Designer.app` (~1.0 GB, x64, unsigned) |
| Reference app (architecture model) | HCR Probe Designer (`projects/HCR-probe-design-main` + the shipped .app) |
| Docs | `docs/ARCHITECTURE.md`, `docs/BUILD-AND-VERIFY.md`, `docs/DEBUG-LOG-2026-10-20.md` in the electron repo |

## How the desktop app works (summary)

Electron spawns the bundled R runtime (R 4.6 + full package closure, ~1 GB,
relocated to be fully self-contained) with
`exec/R -f runner.R --args <payload dir>`. runner.R sources
`Cas9_sgRNA_design.R`, which only *defines* `ui`/`server` in this mode —
this script already guards its standalone dispatch with
`if (sys.nframe() == 0L)`, which is exactly why it worked first; the same
idiom was subsequently applied to Plasmid_Designer.R. The runner picks a
random port, prints `APP_PORT=<n>`, and hosts the Shiny app on 127.0.0.1;
Electron loads that URL in the app window. Full detail:
`docs/ARCHITECTURE.md` in the electron repo.

## Payload copy model

The payload in the electron repo (`apps/crispr-designer/payload/`) is an
**adapted copy** of this project's `Cas9_sgRNA_design.R`:

- a portable-path preamble was added (loads `R/app_paths.R` first),
- the standalone dispatch retains the `sys.nframe()` guard.

### Sync rules

- Changes here (features, bug fixes) → re-copy into the electron payload and
  re-adapt (preamble), then rebuild: `node tools/build-app.js crispr-designer`.
- Changes made in the electron payload that fix general bugs → port back
  here.

## Building / verifying the desktop app

See `docs/BUILD-AND-VERIFY.md` in the electron repo — three verification
gates (bundle content checks, an exact-replica R simulation that must print
`APP_PORT=` and serve HTTP 200, and a real launch). The packaged app writes
diagnostics to `~/Library/Application Support/CRISPR Designer/main.log`.

## Packaging posture

Unsigned, exactly like the shipped HCR Probe Designer.app (the reference
this build intentionally matches). Distributed copies need
right-click → Open on first launch. The app is x64-only for now; arm64
requires the planned CRAN-download runtime bundler (Phase 2 in the electron
repo's roadmap).
