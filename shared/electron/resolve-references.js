/**
 * Reference-data resolution — shared by every app in this repo.
 *
 * Two distinct stores are involved:
 *
 *   1. GENOME store — per-species Drosophila assemblies
 *                      (genomic.fna.gz / genomic.gff.gz / rna.fna.gz)
 *   2. PJAT store   — the Janelia pJAT plasmid GenBank collection
 *
 * Resolution is environment-first, then well-known locations, then a
 * user-writable fallback. Nothing here ever writes into an application
 * bundle, and no symlink is placed inside a bundle: bundles are signed and
 * symlinked payloads break code-signing on macOS.
 */

const path = require('path');
const fs = require('fs');
const os = require('os');

/** Species folders we know about (mirrors R/genome_registry.R). */
const KNOWN_SPECIES = [
  'melanogaster', 'simulans', 'eugracilis',
  'erecta', 'teissieri', 'ananassae', 'virilis',
];

/** Name of the per-species directory inside the genome store. */
function speciesDirName(species) {
  return `d_${String(species).toLowerCase()}`;
}

/**
 * True when `dir` looks like a usable genome store: at least one species
 * folder containing a genome file.
 */
function isGenomeStore(dir) {
  if (!dir || !fs.existsSync(dir)) return false;
  let entries;
  try {
    entries = fs.readdirSync(dir, { withFileTypes: true });
  } catch {
    return false;
  }
  return entries
    .filter((e) => e.isDirectory() && /^d_[a-z0-9_]+$/i.test(e.name))
    .some((e) => {
      try {
        return fs.readdirSync(path.join(dir, e.name))
          .some((f) => /^genomic\.fna(\.gz)?$/i.test(f));
      } catch {
        return false;
      }
    });
}

/** True when `dir` looks like the pJAT plasmid collection. */
function isPjatsStore(dir) {
  if (!dir || !fs.existsSync(dir)) return false;
  let files;
  try {
    files = fs.readdirSync(dir);
  } catch {
    return false;
  }
  return files.some((f) => /^pJAT\d+\.(gbk|gb|gbff|genbank)$/i.test(f));
}

/**
 * Candidate locations for the genome store, most specific first.
 */
function genomeStoreCandidates({ home, userData, resourcesPath, env, dev }) {
  const c = [];

  // 1. Explicit override — always wins.
  const override = env.CRISPR_REFERENCE_DIR || env.HCR_REFERENCE_DIR;
  if (override) c.push(override);

  // 2. HCR Probe Designer's own user-data store. HCR downloads the
  //    per-species sets into <userData>/reference_files/d_<species>/, which is
  //    exactly the layout these apps consume. Sharing it means one download
  //    serves all three applications.
  c.push(path.join(home, 'Library', 'Application Support',
    'hcr-probe-designer', 'reference_files'));

  // 3. HCR's documented download fallback.
  c.push(path.join(home, 'Downloads', 'reference_files'));

  // 4. Our own user-data store (per-app, so the two new apps do not fight).
  if (userData) c.push(path.join(userData, 'reference_files'));

  // 5. Packaged payload, if a build ever ships genomes.
  if (resourcesPath) {
    c.push(path.join(resourcesPath, 'reference_files'));
    c.push(path.join(resourcesPath, 'reference_files', 'drosophila_species_genomes'));
  }

  // 6. Development: a store checked out beside the repo.
  if (dev) {
    c.push(path.join(process.cwd(), 'reference_files'));
    c.push(path.join(process.cwd(), 'reference_files', 'drosophila_species_genomes'));
  }

  // 7. A store relative to the repo root, when running from a worktree.
  c.push(path.join(__dirname, '..', '..', '..', 'reference_files'));
  c.push(path.join(__dirname, '..', '..', '..', 'reference_files',
    'drosophila_species_genomes'));

  return [...new Set(c)];
}

/**
 * Candidate locations for the pJAT plasmid store.
 */
function pjatStoreCandidates({ home, userData, resourcesPath, env }) {
  const c = [];

  const override = env.PJAT_REFERENCE_DIR;
  if (override) c.push(override);

  if (userData) c.push(path.join(userData, 'janelia_atalanta_plasmids'));

  if (resourcesPath) {
    c.push(path.join(resourcesPath, 'janelia_atalanta_plasmids'));
    c.push(path.join(resourcesPath, 'reference_files', 'janelia_atalanta_plasmids'));
  }

  c.push(path.join(home, 'Downloads', 'janelia_atalanta_plasmids'));
  c.push(path.join(home, 'Library', 'Application Support',
    'hcr-probe-designer', 'janelia_atalanta_plasmids'));

  c.push(path.join(__dirname, '..', '..', '..', 'reference_files',
    'janelia_atalanta_plasmids'));

  return [...new Set(c)];
}

/** First candidate that satisfies `predicate`. */
function firstMatch(candidates, predicate) {
  for (const c of candidates) {
    try {
      if (predicate(c)) return c;
    } catch {
      /* ignore unreadable candidates */
    }
  }
  return null;
}

/**
 * Human-readable label for where a resolved path came from.
 *
 * Classified by inspecting the path itself rather than by position, because
 * the candidate list is conditional (the env override is only present when
 * the user actually set it).
 */
function describeSource(candidates, found) {
  if (!found) return 'none';
  const p = found;
  if (p.includes('hcr-probe-designer')) return 'HCR Probe Designer user data';
  if (p.includes('Downloads')) return 'downloads folder';
  if (p.includes('Application Support')) return 'application user data';
  if (candidates.indexOf(p) === 0) return 'environment override';
  if (p.includes(process.cwd())) return 'development checkout';
  return 'repository checkout';
}

/**
 * Resolve both reference stores.
 *
 * @returns {{
 *   genome: {dir: string|null, source: string, species: string[], status: string},
 *   pjat:   {dir: string|null, source: string, status: string}
 * }}
 */
function resolveReferences({ userData, resourcesPath, env = {}, dev = false } = {}) {
  const home = env.HOME || os.homedir();

  const genomeCandidates = genomeStoreCandidates(
    { home, userData, resourcesPath, env, dev });
  const genomeDir = firstMatch(genomeCandidates, isGenomeStore);

  let species = [];
  if (genomeDir) {
    try {
      species = fs.readdirSync(genomeDir, { withFileTypes: true })
        .filter((e) => e.isDirectory() && /^d_[a-z0-9_]+$/i.test(e.name))
        .map((e) => e.name.replace(/^d_/, '').toLowerCase())
        .filter((s) => KNOWN_SPECIES.includes(s))
        .sort();
    } catch {
      species = [];
    }
  }

  const pjatCandidates = pjatStoreCandidates({ home, userData, resourcesPath, env });
  const pjatDir = firstMatch(pjatCandidates, isPjatsStore);

  return {
    genome: {
      dir: genomeDir,
      source: genomeDir ? describeSource(genomeCandidates, genomeDir) : 'none',
      species,
      status: genomeDir ? 'ok' : 'missing',
      candidates: genomeCandidates,
    },
    pjat: {
      dir: pjatDir,
      source: pjatDir ? describeSource(pjatCandidates, pjatDir) : 'none',
      status: pjatDir ? 'ok' : 'missing',
      candidates: pjatCandidates,
    },
  };
}

module.exports = {
  resolveReferences,
  isGenomeStore,
  isPjatsStore,
  speciesDirName,
  KNOWN_SPECIES,
};
