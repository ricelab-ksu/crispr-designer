#!/usr/bin/env node
/**
 * Stage an application for packaging.
 *
 * Reads apps/<id>/app.config.json plus apps/<id>/payload/ and writes a
 * self-contained build directory at apps/<id>/build/ that @electron/packager
 * can turn into a .app:
 *
 *   build/
 *     package.json      app metadata (name, id, version, entry point)
 *     main.js           thin wrapper around shared/electron/bootstrap.js
 *     preload.js        shared preload, namespaced to this app
 *     shared/electron/  the shared Electron + R bootstrap code
 *     payload/          the R app: <AppScript>.R, R/, www/, donor_templates/ …
 *
 * The R payload is copied (not symlinked) so the resulting bundle is
 * self-contained and survives code-signing.
 */

const fs = require('fs');
const path = require('path');

const REPO_ROOT = path.join(__dirname, '..');
const SHARED_ELECTRON = path.join(REPO_ROOT, 'shared', 'electron');

/** Files and folders in the payload that are never packaged. */
const PAYLOAD_EXCLUDE = [
  'Plasmid Designer.app',
  'CRISPR Designer.app',
  '.git',
  '.gitignore',
  '.cline',
  '.DS_Store',
  'launch_app.sh',
  'subagent_logs',
  'logs',
];

/** Recursively copy `src` into `dest`, applying the exclusion list. */
function copyTree(src, dest, exclude) {
  fs.mkdirSync(dest, { recursive: true });
  for (const entry of fs.readdirSync(src, { withFileTypes: true })) {
    if (exclude.includes(entry.name)) continue;
    if (entry.name.startsWith('._')) continue; // AppleDouble
    const s = path.join(src, entry.name);
    const d = path.join(dest, entry.name);
    if (entry.isDirectory()) {
      copyTree(s, d, exclude);
    } else if (entry.isSymbolicLink()) {
      // Recreate the symlink rather than resolving it through copyFileSync,
      // which fails on special files such as sockets.
      const target = fs.realpathSync(s);
      try {
        fs.symlinkSync(target, d);
      } catch (e) {
        fs.copyFileSync(target, d);
      }
    } else {
      fs.copyFileSync(s, d);
    }
  }
}

/** Stage one app. */
function stageApp(id) {
  const appDir = path.join(REPO_ROOT, 'apps', id);
  const cfgPath = path.join(appDir, 'app.config.json');
  if (!fs.existsSync(cfgPath)) {
    throw new Error(`Missing app config: ${cfgPath}`);
  }
  const config = JSON.parse(fs.readFileSync(cfgPath, 'utf8'));
  const payloadDir = path.join(appDir, 'payload');
  if (!fs.existsSync(payloadDir)) {
    throw new Error(`Missing payload: ${payloadDir}`);
  }
  const buildDir = path.join(appDir, 'build');

  // Clean slate — a stale file in build/ would silently ship.
  fs.rmSync(buildDir, { recursive: true, force: true });
  fs.mkdirSync(buildDir, { recursive: true });

  // 1. The R payload.
  copyTree(payloadDir, path.join(buildDir, 'payload'), PAYLOAD_EXCLUDE);

  // 2. The shared Electron code.
  const sharedOut = path.join(buildDir, 'shared', 'electron');
  fs.mkdirSync(sharedOut, { recursive: true });
  for (const f of fs.readdirSync(SHARED_ELECTRON)) {
    fs.copyFileSync(path.join(SHARED_ELECTRON, f), path.join(sharedOut, f));
  }

  // 3. The pJAT plasmid collection (616 KB) is bundled with the app so the
  //    plasmid designer is fully self-contained. It is also placed inside the
  //    payload so a source checkout behaves like a packaged build.
  const pjatSrc = path.join(REPO_ROOT, 'shared', 'r', 'janelia_atalanta_plasmids');
  if (fs.existsSync(pjatSrc)) {
    copyTree(pjatSrc, path.join(buildDir, 'payload', 'janelia_atalanta_plasmids'), []);
    copyTree(pjatSrc, path.join(buildDir, 'r', 'janelia_atalanta_plasmids'), []);
  }

  // 3. main.js — a thin wrapper that supplies the per-app configuration.
  const mainJs = `/**
 * ${config.title} — Electron main process.
 *
 * All of the real work lives in shared/electron/bootstrap.js; this file only
 * supplies the per-app configuration.
 */

const { startApp } = require('./shared/electron/bootstrap');

startApp({
  id: ${JSON.stringify(config.id)},
  title: ${JSON.stringify(config.title)},
  productName: ${JSON.stringify(config.productName)},
  appScript: ${JSON.stringify(config.appScript)},
  accent: ${JSON.stringify(config.accent)},
  splashSub: ${JSON.stringify(config.splashSub)},
  envPrefix: ${JSON.stringify(config.envPrefix)},
  outputDirName: ${JSON.stringify(config.outputDirName || null)},
});
`;
  fs.writeFileSync(path.join(buildDir, 'main.js'), mainJs);

  // 4. preload.js — the shared template, namespaced to this app.
  const template = fs.readFileSync(path.join(SHARED_ELECTRON, 'preload.template.js'), 'utf8');
  fs.writeFileSync(
    path.join(buildDir, 'preload.js'),
    template.replace(/__APP_ID__/g, config.id));

  // 5. The Electron-facing package.json.
  const pkg = {
    name: config.id,
    productName: config.title,
    version: config.version || '1.0.0',
    description: config.description || config.title,
    main: 'main.js',
    author: config.author || 'Rice Lab',
    license: 'MIT',
  };
  fs.writeFileSync(path.join(buildDir, 'package.json'), `${JSON.stringify(pkg, null, 2)}\n`);

  console.log(`[stage] ${id} -> ${path.relative(REPO_ROOT, buildDir)}`);
  return buildDir;
}

if (require.main === module) {
  const id = process.argv[2] || 'plasmid-designer';
  try {
    stageApp(id);
  } catch (err) {
    console.error(`[stage] FAILED: ${err.message}`);
    process.exit(1);
  }
}

module.exports = { stageApp, copyTree, PAYLOAD_EXCLUDE };
