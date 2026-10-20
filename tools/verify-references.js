#!/usr/bin/env node
/**
 * Verify that the reference stores resolve correctly and that the R side
 * agrees with the Electron side.
 *
 * Checks:
 *   1. the Electron resolver finds the genome store and all 7 species,
 *   2. the R resolver (shared/r/app_paths.R) finds the same directory,
 *   3. every species has the three expected files,
 *   4. the pJAT collection is present and complete,
 *   5. gene lookup works against the compressed GFF3 for a known gene.
 *
 * Usage: node tools/verify-references.js [--json]
 */

const path = require('path');
const fs = require('fs');
const { execFileSync } = require('child_process');

const { resolveReferences } = require('../shared/electron/resolve-references');
const REPO_ROOT = path.join(__dirname, '..');
const AS_JSON = process.argv.includes('--json');

const EXPECTED_SPECIES = [
  'ananassae', 'erecta', 'eugracilis', 'melanogaster',
  'simulans', 'teissieri', 'virilis',
];
const EXPECTED_FILES = ['genomic.fna.gz', 'genomic.gff.gz', 'rna.fna.gz'];

const results = [];
let failed = 0;

/** Record a check. */
function check(name, ok, detail) {
  results.push({ name, ok, detail });
  if (!ok) failed += 1;
  if (!AS_JSON) {
    console.log(`${ok ? '  ok  ' : ' FAIL '} ${name}${detail ? ` — ${detail}` : ''}`);
  }
}

// ── 1. Electron-side resolution ───────────────────────────────────────────

const resolved = resolveReferences({
  userData: path.join(REPO_ROOT, '.verify-userdata'),
  resourcesPath: path.join(REPO_ROOT, '.verify-resources'),
  env: process.env,
  dev: true,
});

check('genome store found', resolved.genome.status === 'ok',
  resolved.genome.dir || 'not found');
check('genome store is the shared HCR store',
  (resolved.genome.dir || '').includes('hcr-probe-designer')
  || (resolved.genome.dir || '').includes('Downloads'),
  resolved.genome.source);

const missingSpecies = EXPECTED_SPECIES.filter((s) => !resolved.genome.species.includes(s));
check('all 7 species present', missingSpecies.length === 0,
  missingSpecies.length ? `missing: ${missingSpecies.join(', ')}`
    : resolved.genome.species.join(', '));

// ── 2. Per-species files ───────────────────────────────────────────────────

if (resolved.genome.dir) {
  for (const sp of EXPECTED_SPECIES) {
    const dir = path.join(resolved.genome.dir, `d_${sp}`);
    const absent = EXPECTED_FILES.filter(
      (f) => !fs.existsSync(path.join(dir, f)));
    check(`d_${sp} has all three files`, absent.length === 0,
      absent.length ? `missing: ${absent.join(', ')}` : 'ok');
  }
}

// ── 3. pJAT collection ─────────────────────────────────────────────────────

const pjatDir = resolved.pjat.dir
  || path.join(REPO_ROOT, 'shared', 'reference', 'janelia_atalanta_plasmids');
const pjatFiles = fs.existsSync(pjatDir)
  ? fs.readdirSync(pjatDir).filter((f) => /^pJAT\d+\./i.test(f))
  : [];
check('pJAT collection present', pjatFiles.length > 0,
  `${pjatFiles.length} files in ${pjatDir}`);

// ── 4. R-side agreement ────────────────────────────────────────────────────

function runR(code) {
  return execFileSync('Rscript', ['--vanilla', '-e', code], {
    cwd: path.join(REPO_ROOT, 'apps', 'crispr-designer', 'payload'),
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', 'pipe'],
  });
}

try {
  const out = runR(`
    source("R/app_paths.R")
    cat(get_app_reference_dir("CRISPR"))
  `).trim();
  check('R resolver agrees with Electron',
    out === resolved.genome.dir,
    out === resolved.genome.dir ? 'same directory' : `R=${out} vs JS=${resolved.genome.dir}`);
} catch (err) {
  check('R resolver agrees with Electron', false, err.message.split('\n')[0]);
}

// ── 5. Gene lookup against the compressed GFF3 ─────────────────────────────

try {
  const out = runR(`
    suppressWarnings(suppressMessages({
      source("R/genome_registry.R")
      source("R/gene_lookup.R")
    }))
    res <- lookup_gene("white", species = "melanogaster")
    if (is.null(res) || !length(res)) {
      cat("NO_RESULT")
    } else {
      ids <- if (!is.null(res$fbgn)) res$fbgn else
             if (!is.null(res$gene_id)) res$gene_id else "unknown"
      cat(as.character(ids)[1])
    }
  `).trim();
  check('gene lookup works (white / melanogaster)',
    out && out !== 'NO_RESULT',
    out === 'NO_RESULT' ? 'no result' : `resolved ${out}`);
} catch (err) {
  check('gene lookup works (white / melanogaster)', false,
    err.message.split('\n')[0]);
}

// ── 6. Genomic sequence extraction ─────────────────────────────────────────

try {
  const out = runR(`
    suppressWarnings(suppressMessages({
      source("R/genome_registry.R")
    }))
    seq <- extract_genomic_subseq("melanogaster", "X", 1, 60)
    cat(if (is.null(seq) || !length(seq)) "EMPTY" else nchar(as.character(seq)[1]))
  `).trim();
  check('genomic sequence extraction works',
    out && out !== 'EMPTY' && Number(out) > 0,
    `${out} bp`);
} catch (err) {
  check('genomic sequence extraction works', false, err.message.split('\n')[0]);
}

// ── Summary ────────────────────────────────────────────────────────────────

if (AS_JSON) {
  console.log(JSON.stringify({ results, failed }, null, 2));
} else {
  console.log(`\n${results.length - failed}/${results.length} checks passed.`);
  if (failed) {
    console.log('FAILED: ' + results.filter((r) => !r.ok).map((r) => r.name).join('; '));
  }
}
process.exit(failed ? 1 : 0);
