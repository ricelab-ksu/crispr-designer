#!/usr/bin/env node
// sign-app.js - bottom-up ad-hoc signer for Electron + R bundles
// Signs every framework / helper / dylib from the inside out, then the outer .app last.
// Uses codesign with NO --deep and NO --options runtime.
// Hardened Runtime on ad-hoc causes dyld to reject different Team IDs; avoiding --options runtime sidesteps that.

const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');

const APP = process.argv[2];
if (!APP) { console.error('Usage: node tools/sign-app.js /path/to/MyApp.app'); process.exit(1); }
if (!fs.existsSync(APP)) { console.error('Not found: ' + APP); process.exit(1); }

// Sign inner dylibs / .so / executables inside frameworks
function signDylibs(rootDir) {
  const results = [];

  function walk(dir) {
    const entries = fs.readdirSync(dir, { withFileTypes: true });
    for (const entry of entries) {
      if (entry.name.startsWith('.')) continue;
      const full = path.join(dir, entry.name);
      if (entry.isDirectory()) {
        walk(full);
      } else if (entry.isFile()) {
        const basename = entry.name;
        if (basename.endsWith('.dylib') || basename.endsWith('.so') || basename.endsWith('.node')) {
          results.push({ kind: 'dylib', path: full });
        }
      }
    }
  }

  const fwDir = path.join(rootDir, 'Contents', 'Frameworks');
  if (fs.existsSync(fwDir)) {
    const fwEntries = fs.readdirSync(fwDir, { withFileTypes: true });
    for (const entry of fwEntries) {
      if (entry.name.endsWith('.framework')) {
        walk(path.join(fwDir, entry.name));
      }
    }
  }

  const resourcesBin = path.join(rootDir, 'Contents', 'Resources', 'bin');
  if (fs.existsSync(resourcesBin)) {
    walk(resourcesBin);
  }

  for (const item of results) {
    try {
      console.log('signing ' + item.kind + ' ' + item.path);
      execFileSync('codesign', ['--force', '--sign', '-', '--timestamp=none', item.path], { encoding: 'utf8' });
    } catch (e) {
      console.error('abort: could not sign ' + item.path);
      process.exit(1);
    }
  }
}

// Sign the Electron Framework itself (after its dylibs are signed)
function signElectronFramework(rootDir) {
  const fwPath = path.join(rootDir, 'Contents', 'Frameworks', 'Electron Framework.framework');
  if (!fs.existsSync(fwPath)) return;

  const fwBinary = path.join(fwPath, 'Versions', 'A', 'Electron Framework');
  if (fs.existsSync(fwBinary)) {
    try {
      console.log('signing framework ' + fwBinary);
      execFileSync('codesign', ['--force', '--sign', '-', '--timestamp=none', fwBinary], { encoding: 'utf8' });
    } catch (e) { console.error('abort: could not sign ' + fwBinary); }
  }

  const fwBundle = path.join(fwPath, 'Versions', 'A');
  if (fs.existsSync(fwBundle)) {
    try {
      console.log('signing framework bundle ' + fwBundle);
      execFileSync('codesign', ['--force', '--sign', '-', '--timestamp=none', fwBundle], { encoding: 'utf8' });
    } catch (e) { console.error('abort: could not sign ' + fwBundle); }
  }
}

// Sign nested Electron helper .app bundles inside Frameworks
function signNestedHelpers(appPath) {
  const fwDir = path.join(appPath, 'Contents', 'Frameworks');
  if (!fs.existsSync(fwDir)) return;

  function walkDir(d) {
    const entries = fs.readdirSync(d, { withFileTypes: true });
    for (const e of entries) {
      if (e.name.startsWith('.')) continue;
      const full = path.join(d, e.name);
      if (e.isDirectory()) {
        if (e.name.endsWith('.app')) {
          try { console.log('signing nested helper .app ' + full); execFileSync('codesign', ['--force', '--sign', '-', '--timestamp=none', full], { encoding: 'utf8' }); }
          catch (e) { process.exit(1); }
        } else {
          walkDir(full);
        }
      }
    }
  }
  walkDir(fwDir);
}

// Sign the main executable
function signMainExecutable(appPath) {
  const mePath = path.join(appPath, 'Contents', 'MacOS');
  if (!fs.existsSync(mePath)) return;
  const entries = fs.readdirSync(mePath, { withFileTypes: true });
  for (const e of entries) {
    if (e.name.startsWith('.')) continue;
    const full = path.join(mePath, e.name);
    if (e.isFile() && !e.name.endsWith('.app')) {
      try { console.log('signing executable ' + full); execFileSync('codesign', ['--force', '--sign', '-', '--timestamp=none', full], { encoding: 'utf8' }); }
      catch (e) { process.exit(1); }
    }
  }
}

// Sign the outer .app bundle (final step)
function signOuterApp(appPath) {
  try { console.log('signing outer app ' + appPath); execFileSync('codesign', ['--force', '--sign', '-', '--timestamp=none', appPath], { encoding: 'utf8' }); }
  catch (e) { process.exit(1); }
}

// Execution order (bottom-up, then the enclosing app last)
signDylibs(APP);
signElectronFramework(APP);
signNestedHelpers(APP);
signMainExecutable(APP);
signOuterApp(APP);

console.log('SIGNING COMPLETED for ' + APP);
console.log('done — all signed');
