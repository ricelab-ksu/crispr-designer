# Bioinformatics Electron Apps

Standalone macOS desktop applications for the **Plasmid Designer** and
**CRISPR Designer** R/Shiny tools, built on a shared Electron + R framework.

Both apps are the same application with a different R payload. They share:

- a single Electron main process (`shared/electron/`),
- a single bundled R runtime (`shared/r/`),
- a single reference-data resolver,
- a single reference store on disk.

---

## Architecture

```
┌─────────────────────┐        ┌─────────────────────┐
│  Plasmid Designer   │        │  CRISPR Designer    │
│  (Electron + R)     │        │  (Electron + R)     │
└──────────┬──────────┘        └──────────┬──────────┘
           │                              │
           └──────────┬───────────────────┘
                      ▼
        shared/electron/          shared/r/
        ├── bootstrap.js           ├── R.framework/   (base + recommended)
        ├── r-runtime.js           ├── library/       (85 packages)
        ├── resolve-references.js  └── janelia_atalanta_plasmids/
        ├── preload.template.js
        └── runner.R
                      │
                      ▼
        reference store (shared, read-only)
        └── <userData>/reference_files/d_<species>/
            ├── genomic.fna.gz
            ├── genomic.gff.gz
            └── rna.fna.gz
```

The R runtime is spawned as a child process. It picks a free TCP port, prints
`APP_PORT=<n>` on stdout, and serves the Shiny app on `127.0.0.1`. Electron
waits for that handshake and then points its window at the server.

---

## Reference data

Two independent stores are involved.

| Store | Contents | Where it lives |
|---|---|---|
| **Genome store** | 7 Drosophila species assemblies (`.fna.gz` / `.gff.gz`) | Shared with HCR Probe Designer; downloaded on demand |
| **pJAT store** | 24 Janelia pJAT plasmid GenBank files (616 KB) | **Bundled inside each app** |

### Genome store resolution order

1. `CRISPR_REFERENCE_DIR` / `PLASMID_REFERENCE_DIR` (explicit override)
2. `HCR_REFERENCE_DIR`
3. `~/Library/Application Support/hcr-probe-designer/reference_files`
4. `~/Downloads/reference_files`
5. `<app userData>/reference_files`
6. a packaged payload, if a build ever ships genomes

The layout is `reference_files/d_<species>/{genomic.fna.gz,genomic.gff.gz,rna.fna.gz}`
— byte-for-byte what HCR Probe Designer writes, so **one download serves all
three applications**.

Resolution is strictly read-only. Nothing is symlinked into an application
bundle, and no file is written inside a bundle: both would break macOS code
signing.

---

## Repository layout

```
apps/
  plasmid-designer/
    app.config.json      app metadata (id, title, colours, env prefix)
    payload/             the R application (unchanged from the original repo)
  crispr-designer/
    app.config.json
    payload/
shared/
  electron/              shared Electron main process + R sidecar
  r/app_paths.R          portable R path resolution (shared by both payloads)
  reference/             the bundled pJAT plasmid collection
  icons/
tools/
  stage-app.js           assemble apps/<id>/build/
  build-app.js           package + ad-hoc sign a .app
  build-r-runtime.sh     assemble shared/r/
  r-closure.R            compute the R dependency closure
```

---

## Building

```bash
npm install                 # Electron + @electron/packager
npm run stage:all           # stage both apps
npm run build:mac:all       # package + sign both apps
```

The R runtime is built separately and is cached:

```bash
tools/build-r-runtime.sh    # → shared/r/  (~330 MB)
```

Results land in `dist/<app-id>/<App Name>.app`.

### Gatekeeper

The apps are **ad-hoc signed**, not notarized. On first launch macOS will
refuse to open them. To run one:

1. Right-click the app in Finder → **Open**.
2. Confirm **Open** in the dialog.

(Or run `xattr -dr com.apple.quarantine "/path/to/App.app"`.)

---

## Development

Run against the system R without packaging:

```bash
cd apps/crispr-designer
electron --dev .
```

The Electron main process sets the reference directories in the environment
before spawning R:

| Variable | Meaning |
|---|---|
| `<PREFIX>_REFERENCE_DIR` | genome store |
| `<PREFIX>_PJAT_DIR` | pJAT plasmid collection |
| `APP_SCRIPT` | R file defining `ui` / `server` |
| `APP_ROOT` | application root (optional override) |

`PLASMID` is the prefix for the Plasmid Designer, `CRISPR` for the CRISPR
Designer.

---

## Testing

```bash
npm run test:r              # the payloads' own testthat suites
node tools/verify-references.js   # reference-store parity check
```
