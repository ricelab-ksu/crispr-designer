/**
 * R runtime resolution — shared by every app in this repo.
 *
 * Layout inside a packaged .app:
 *   Contents/Resources/r/R.framework/Versions/<v>/Resources/…   ← R itself
 *   Contents/Resources/r/library/<pkg>                          ← our packages
 *
 * The framework is a plain copy of the system R (base + recommended only),
 * so the bundled runtime is self-contained and does not depend on whatever
 * happens to be installed on the user's machine.
 */

const path = require('path');
const fs = require('fs');

/** Highest framework version directory under a given root. */
function newestVersion(versionsDir) {
  if (!fs.existsSync(versionsDir)) return null;
  const vers = fs.readdirSync(versionsDir).filter((v) => /^\d/.test(v)).sort().reverse();
  return vers.length ? vers[0] : null;
}

/**
 * Root of the bundled R runtime inside a packaged .app, or null.
 *
 * Two layouts are supported:
 *   Contents/Resources/r       ← electron-builder extraResources (HCR style)
 *   Contents/Resources/app/r   ← @electron/packager (everything under app/)
 */
function resolveRRoot(resourcesPath) {
  if (!resourcesPath) return null;
  const candidates = [
    path.join(resourcesPath, 'r'),
    path.join(resourcesPath, 'app', 'r'),
  ];
  for (const c of candidates) {
    if (fs.existsSync(c)) return c;
  }
  return null;
}

/**
 * Locate an R runtime to spawn.
 *
 * @param {object} opts
 * @param {boolean} opts.dev         true in unpackaged development mode
 * @param {string}  opts.resourcesPath  Electron's process.resourcesPath
 * @param {string}  opts.env         the process environment
 * @returns {{kind:string, executable:string, args:string[], env:object}|null}
 */
function resolveRRuntime({ dev, resourcesPath, env }) {
  // ── Development: prefer an explicit override, then system Rscript ────────
  if (dev) {
    const override = env.APP_RSCRIPT || env.RSCRIPT;
    if (override && fs.existsSync(override)) {
      return { kind: 'system', executable: override, args: [], env: {} };
    }
    return { kind: 'system', executable: 'Rscript', args: [], env: {} };
  }

  // ── Packaged: bundled R.framework ───────────────────────────────────────
  const rRoot = resolveRRoot(resourcesPath);
  const versionsDir = rRoot ? path.join(rRoot, 'R.framework', 'Versions') : null;
  const version = newestVersion(versionsDir);
  if (rRoot && version) {
    const res = path.join(versionsDir, version, 'Resources');
    const execR = path.join(res, 'bin', 'exec', 'R');
    if (fs.existsSync(execR)) {
      const lib = path.join(rRoot, 'library');
      const env = {
        R_HOME: res,
        R_SHARE_DIR: path.join(res, 'share'),
        R_INCLUDE_DIR: path.join(res, 'include'),
        R_DOC_DIR: path.join(res, 'doc'),
        // Our package tree first, so the bundled versions always win over
        // anything the user happens to have in their personal library.
        R_LIBS: lib,
        R_LIBS_USER: lib,
        R_ENVIRON_USER: '',
        R_PROFILE_USER: '',
      };
      if (fs.existsSync(lib)) {
        env.R_LIBS = `${lib}${path.delimiter}${path.join(res, 'library')}`;
      }
      return { kind: 'framework', executable: execR, args: [], env };
    }
  }

  // ── Packaged: Windows-style R-x.y layout (future platform support) ──────
  if (rRoot && fs.existsSync(rRoot)) {
    const sub = fs.readdirSync(rRoot).find((d) => /^R-\d/.test(d));
    if (sub) {
      const rHome = path.join(rRoot, sub);
      const script = path.join(rHome, 'bin', 'x64', 'Rscript.exe');
      if (fs.existsSync(script)) {
        return {
          kind: 'windows',
          executable: script,
          args: [],
          env: {
            R_HOME: rHome,
            R_LIBS: path.join(rRoot, 'library'),
            R_LIBS_USER: path.join(rRoot, 'library'),
            R_ENVIRON_USER: '',
            R_PROFILE_USER: '',
          },
        };
      }
    }
  }

  return null;
}

/** Path to the bundled package library, or null when there is none. */
function bundledLibrary(resourcesPath) {
  const rRoot = resolveRRoot(resourcesPath);
  if (!rRoot) return null;
  const lib = path.join(rRoot, 'library');
  return fs.existsSync(lib) ? lib : null;
}

module.exports = { resolveRRuntime, resolveRRoot, bundledLibrary, newestVersion };
