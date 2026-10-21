#!/usr/bin/env node
/**
 * Package an application into a macOS .app bundle.
 */

const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');

const { stageApp } = require('./stage-app');

const REPO_ROOT = path.join(__dirname, '..');
const DIST_DIR = path.join(REPO_ROOT, 'dist');
const R_ROOT = path.join(REPO_ROOT, 'shared', 'r');

/** Recursively copy a directory tree (preserves relative symlink targets). */
function copyTree(src, dest) {
  fs.mkdirSync(dest, { recursive: true });
  for (const entry of fs.readdirSync(src, { withFileTypes: true })) {
    if (entry.name.startsWith('._')) continue;
    const s = path.join(src, entry.name);
    const d = path.join(dest, entry.name);
    if (entry.isDirectory()) {
      copyTree(s, d);
    } else if (entry.isSymbolicLink()) {
      const rel = fs.readlinkSync(s);          // keep the relative target as-is
      try {
        fs.symlinkSync(rel, d);                // re-create with same relative target
      } catch (e) {
        fs.copyFileSync(fs.realpathSync(s), d); // fallback: copy the resolved file
      }
    } else {
      fs.copyFileSync(s, d);
    }
  }
}

async function buildApp(id, { sign = true } = {}) {
  const cfg = JSON.parse(
    fs.readFileSync(path.join(REPO_ROOT, 'apps', id, 'app.config.json'), 'utf8'));

  // 1. Stage.
  const buildDir = stageApp(id);

  // 2. Attach the R runtime.
  if (!fs.existsSync(path.join(R_ROOT, 'R.framework'))) {
    throw new Error(
      'No R runtime found. Run tools/build-r-runtime.sh before packaging.');
  }
  console.log('[build] copying the R runtime…');
  copyTree(R_ROOT, path.join(buildDir, 'r'));

  // 3. Package.
  const { packager } = require('@electron/packager');
  const outDir = path.join(DIST_DIR, id);
  fs.mkdirSync(outDir, { recursive: true });

  const iconPath = path.join(REPO_ROOT, 'shared', 'icons', `${id}.icns`);
  const icon = fs.existsSync(iconPath)
    ? iconPath
    : path.join(REPO_ROOT, 'shared', 'icons', 'icon.icns');

  console.log('[build] running @electron/packager…');
  const paths = await packager({
    dir: buildDir,
    out: outDir,
    appVersion: cfg.version || '1.0.0',
    appBundleId: `edu.ricelab.${cfg.id}`,
    appCategoryType: 'public.app.category.education',
    icon: fs.existsSync(icon) ? icon : undefined,
    overwrite: true,
    asar: false,
    name: cfg.title,
    quiet: true,
  });

  const appPath = paths[0];
  console.log(`[build] packaged: ${appPath}`);

  // packager can emit either a bare `.app` bundle or a platform directory
  // containing one. Resolve to the `.app` so signing and size reporting are
  // correct regardless.
  let resolvedApp = appPath;
  if (fs.existsSync(appPath) && !fs.statSync(appPath).isDirectory()) {
    // shouldn't happen, but be defensive
  } else if (fs.existsSync(appPath)) {
    const cands = fs.readdirSync(appPath)
      .map((n) => path.join(appPath, n))
      .filter((p) => fs.existsSync(p) && p.endsWith('.app') && fs.statSync(p).isDirectory());
    if (cands.length) resolvedApp = cands[0];
  }

  // 4. Skip Electron-side signing (matching HCR Probe Designer: unsigned
  // Electron side; only per-file ad-hoc on relocated R binaries via
  // fix-r-framework.sh, which runs inside copyTree → build-r-runtime.sh).
  // No --deep, no --options runtime, no Team-ID: the app launches clean.
  // if (sign) { … } block removed entirely.

  // Report the final size.
  const size = execFileSync('du', ['-sh', resolvedApp]).toString().split('	')[0];
  console.log(`[build] done: ${resolvedApp} (${size})`);

  return resolvedApp;
}

if (require.main === module) {
  const args = process.argv.slice(2);
  const id = args.find((a) => !a.startsWith('--')) || 'plasmid-designer';
  const sign = !args.includes('--no-sign');
  buildApp(id, { sign })
    .then(() => process.exit(0))
    .catch((err) => {
      console.error(`[build] FAILED: ${err.message}`);
      process.exit(1);
    });
}

module.exports = { buildApp, copyTree };
