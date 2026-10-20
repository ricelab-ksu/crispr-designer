/**
 * Shared Electron main-process bootstrap.
 *
 * Every app in this repo has the same shape:
 *   1. show a splash window,
 *   2. spawn the R runtime, which picks a free port and announces it on stdout,
 *   3. point the window at the Shiny server once it answers,
 *   4. surface R's stderr in the splash if anything goes wrong.
 *
 * apps/<name>/main.js is a thin wrapper around startApp(config).
 */

const { app, BrowserWindow, dialog, Menu, Notification, session, shell, ipcMain } =
  require('electron');
const { spawn } = require('child_process');
const path = require('path');
const fs = require('fs');

const { resolveRRuntime } = require('./r-runtime');
const { resolveReferences } = require('./resolve-references');

const REPO_ROOT = path.join(__dirname, '..', '..');

/**
 * @typedef {object} AppConfig
 * @property {string}   id             app id, e.g. 'plasmid-designer'
 * @property {string}   title          window + notification title
 * @property {string}   productName    human-readable product name
 * @property {string}   appScript      R file defining ui/server
 * @property {string}   accent         splash accent colour
 * @property {string}   splashSub      splash subtitle line
 * @property {string}   envPrefix      prefix for reference env vars, e.g. 'PLASMID'
 * @property {string=}  outputDirName  folder under Documents for outputs
 */

/**
 * Boot an application.
 *
 * @param {AppConfig} config
 */
function startApp(config) {
  const DEV_MODE = process.argv.includes('--dev') || !app.isPackaged;
  const APP_VERSION = app.isPackaged
    ? app.getVersion()
    : require(path.join(REPO_ROOT, 'package.json')).version;

  const IPC = `${config.id}:`;

  // ── Path resolution ──────────────────────────────────────────────────────

  /** Directory holding the R payload (the .R app script and its R/ folder). */
  function resolveAppDir() {
    if (DEV_MODE) {
      const dev = path.join(REPO_ROOT, 'apps', config.id, 'payload');
      if (fs.existsSync(dev)) return dev;
    }
    if (process.resourcesPath) {
      const p = path.join(process.resourcesPath, 'app');
      if (fs.existsSync(p)) return p;
    }
    return __dirname;
  }

  /** Directory that holds the shared R runtime and the pJAT plasmid store. */
  function resolveRRoot() {
    if (process.resourcesPath) return path.join(process.resourcesPath, 'r');
    return null;
  }

  // ── R child process + port handshake ─────────────────────────────────────

  let rProc = null;
  let rBufferedError = '';
  let rPortResolve = null;

  function spawnR(appDir, env) {
    const rt = resolveRRuntime({
      dev: DEV_MODE,
      resourcesPath: process.resourcesPath,
      env: process.env,
    });
    if (!rt) {
      throw new Error(
        'No bundled R runtime found and not running in --dev mode. '
        + 'Run `npm run build:mac` first to stage the R runtime.');
    }

    const runner = path.join(__dirname, 'runner.R');
    // exec/R framework takes -f <file> --args ...; Rscript takes the script positionally.
    const scriptArgs = rt.kind === 'framework'
      ? ['-f', runner, '--args', appDir]
      : [runner, appDir];
    const args = ['--no-save', '--no-restore', '--no-environ', '--no-site-file',
      '--no-init-file', '--slave', ...scriptArgs];

    const spawnEnv = {
      ...process.env,
      ...env,
      R_DEFAULT_PACKAGES: 'methods,stats,graphics,grDevices,utils,datasets',
      R_PROFILE_USER: '', R_ENVIRON_USER: '',
    };
    if (rt.env) Object.assign(spawnEnv, rt.env);

    console.log(`[${config.id}] spawning R: ${rt.executable} ${args.join(' ')}`);
    rProc = spawn(rt.executable, args, {
      env: spawnEnv, stdio: ['ignore', 'pipe', 'pipe'],
    });

    let stdoutBuf = '';
    rProc.stdout.on('data', (chunk) => {
      stdoutBuf += chunk.toString();
      const m = stdoutBuf.match(/(?:APP|HCR|CAS9)_PORT=(\d+)/);
      if (m && rPortResolve) {
        rPortResolve(parseInt(m[1], 10));
        rPortResolve = null;
      }
    });
    rProc.stderr.on('data', (chunk) => {
      const s = chunk.toString();
      rBufferedError += s;
      if (rBufferedError.length > 200000) rBufferedError = rBufferedError.slice(-200000);
      console.error(`[${config.id}] R >>`,
        s.trimEnd().split('\n').join(`\n[${config.id}] R >> `));
    });
  // ── Window management ────────────────────────────────────────────────────

  let mainWindow = null;

  function createLoadingWindow() {
    mainWindow = new BrowserWindow({
      width: 1280, height: 830, minWidth: 1000, minHeight: 700,
      title: config.title,
      backgroundColor: '#1e1e2e',
      show: false,
      webPreferences: {
        preload: path.join(__dirname, 'preload.js'),
        contextIsolation: true, nodeIntegration: false, sandbox: true,
      },
    });
    mainWindow.on('closed', () => { mainWindow = null; });
    mainWindow.webContents.setWindowOpenHandler(({ url }) => {
      if (url.startsWith('http')) shell.openExternal(url);
      return { action: 'deny' };
    });

    const bootstrapHtml = `<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<style>
  body{margin:0;background:#1e1e2e;color:#cdd6f4;font-family:-apple-system,Segoe UI,sans-serif;
    display:flex;align-items:center;justify-content:center;height:100vh;flex-direction:column;gap:18px;}
  .logo{font-size:30px;letter-spacing:1px;font-weight:600;color:${config.accent};}
  .sub{color:#9399b2;font-size:14px;}
  .bar{width:260px;height:4px;background:#313244;border-radius:2px;overflow:hidden;}
  .bar i{display:block;width:40%;height:100%;background:#89b4fa;border-radius:2px;
    animation:move 1.1s ease-in-out infinite;}
  @keyframes move{0%{transform:translateX(-100%)}100%{transform:translateX(650%)}}
  #err{display:none;color:#f38ba8;font-size:13px;max-width:620px;text-align:center;
    white-space:pre-wrap;font-family:ui-monospace,Menlo,monospace;}
  #hint{display:none;color:#9399b2;font-size:12px;margin-top:8px;}
</style></head><body>
  <div class="logo">${config.productName}</div>
  <div class="sub">${config.splashSub}</div>
  <div class="bar"><i></i></div>
  <div id="err"></div>
  <div id="hint"></div>
  <script>
    const err = document.getElementById('err');
    const hint = document.getElementById('hint');
    window.__APP__.onError((msg) => {
      err.style.display = 'block';
      err.textContent = msg;
    });
    window.__APP__.onReferences((info) => {
      if (info && info.genome && info.genome.status === 'missing') {
        hint.style.display = 'block';
        hint.textContent = 'No genome store found — you can download one from the app.';
      }
    });
  </script></body></html>`;

    mainWindow.loadURL('data:text/html;charset=utf-8,' + encodeURIComponent(bootstrapHtml));
    mainWindow.show();
    return mainWindow;
  }

  function onShinyReady(port) {
    if (!mainWindow || mainWindow.isDestroyed()) return;
    const url = `http://127.0.0.1:${port}`;
    let attempts = 0;

    // runner.R prints the port just before runApp() binds the socket, so the
    // first load may hit ERR_CONNECTION_REFUSED. Retry until the server is up.
    const load = () => {
      attempts += 1;
      mainWindow.loadURL(url).catch((err) => {
        if (attempts < 120 && err && err.errno === 'ERR_CONNECTION_REFUSED') {
          setTimeout(load, 400);
        }
      });
    };
    mainWindow.webContents.on('did-fail-load', (_e, code, desc, _url, isMainFrame) => {
      if (isMainFrame && code === -102 && attempts < 120) {
        setTimeout(load, 400);
      }
    });
    load();

    mainWindow.once('page-title-updated', (event, title) => {
      event.preventDefault();
      mainWindow.setTitle(`${config.title} - ${title}`);
    });
  }

  function showRFatalError(tail) {
    if (!mainWindow || mainWindow.isDestroyed()) return;
    const summary = tail.length > 8000 ? tail.slice(-8000) : tail;
    mainWindow.webContents.send(`${IPC}error`,
      `The R design engine could not start.\n\n${summary}\n\n`
      + 'Make sure the app resources are intact, then relaunch.');
  }

    rProc.on('error', (err) => {
      if (rPortResolve) { rPortResolve(null); rPortResolve = null; }
      console.error(`[${config.id}] failed to launch R:`, err.message);
    });
    rProc.on('exit', (code, signal) => {
      console.log(`[${config.id}] R exited code=${code} signal=${signal}`);
      if (rPortResolve) { rPortResolve(null); rPortResolve = null; }
    });
    return rProc;
  // ── Reference resolution ─────────────────────────────────────────────────

  /**
   * Resolve the reference stores and expose them to R via the environment.
   *
   * The pJAT plasmid collection is bundled with the app, so it is always
   * available; the genome store is large and shared, so it is located at
   * runtime and may legitimately be absent.
   */
  function resolveAndExportReferences() {
    const rRoot = resolveRRoot();
    const bundledPjats = rRoot
      ? path.join(rRoot, 'janelia_atalanta_plasmids')
      : path.join(REPO_ROOT, 'shared', 'reference', 'janelia_atalanta_plasmids');

    const env = { ...process.env };

    // Let an explicit override win, then fall back to the bundled collection.
    if (!env[`${config.envPrefix}_PJAT_DIR`] && fs.existsSync(bundledPjats)) {
      env[`${config.envPrefix}_PJAT_DIR`] = bundledPjats;
    }

    const resolved = resolveReferences({
      userData: app.getPath('userData'),
      resourcesPath: process.resourcesPath,
      env,
      dev: DEV_MODE,
    });

    if (resolved.genome.dir) {
      env[`${config.envPrefix}_REFERENCE_DIR`] = resolved.genome.dir;
    }
    if (resolved.pjat.dir) {
      env[`${config.envPrefix}_PJAT_DIR`] = resolved.pjat.dir;
    }

    // The bundled pJAT collection is the app's own payload, so it must win
    // over anything found in the wider environment.
    if (fs.existsSync(bundledPjats)) {
      env[`${config.envPrefix}_PJAT_DIR`] = bundledPjats;
    }

    return { resolved, env };
  }

  // ── App lifecycle ────────────────────────────────────────────────────────

  async function boot() {
    Menu.setApplicationMenu(null);
    createLoadingWindow();

    // Route downloads to the user's Downloads folder (needs app-ready session).
    session.defaultSession.on('will-download', (event, item) => {
      const downloads = app.getPath('downloads');
      item.setSavePath(path.join(downloads, item.getFilename()));
      item.once('done', (_e, state) => {
        if (state === 'completed' && Notification.isSupported()) {
          new Notification({ title: 'Download complete', body: item.getFilename() }).show();
        }
      });
    });

    const appDir = resolveAppDir();
    const { resolved, env } = resolveAndExportReferences();

    console.log(`[${config.id}] appDir=${appDir} dev=${DEV_MODE}`);
    console.log(`[${config.id}] genomes: ${resolved.genome.dir || '(none)'} `
      + `[${resolved.genome.source}] species=${resolved.genome.species.length}`);
    console.log(`[${config.id}] pJAT:    ${resolved.pjat.dir || '(none)'} `
      + `[${resolved.pjat.source}]`);

    if (mainWindow && !mainWindow.isDestroyed()) {
      mainWindow.webContents.send(`${IPC}references`, {
        genome: {
          status: resolved.genome.status,
          dir: resolved.genome.dir,
          source: resolved.genome.source,
          species: resolved.genome.species,
        },
        pjat: {
          status: resolved.pjat.status,
          dir: resolved.pjat.dir,
          source: resolved.pjat.source,
        },
      });
    }

    if (!fs.existsSync(path.join(appDir, config.appScript))) {
      dialog.showErrorBox(config.title,
        `App payload missing: ${config.appScript}\n\nExpected it in: ${appDir}`);
      app.quit();
      return;
    }

    const port = await new Promise((resolve) => {
      rPortResolve = resolve;
      try { spawnR(appDir, env); } catch (err) {
        console.error(`[${config.id}] spawn failed:`, err.message);
        resolve(null);
      }
    });

    if (port) {
      onShinyReady(port);
    } else {
      const tail = rBufferedError || '(no output)';
      console.error(`[${config.id}] R failed to report a port.\n`, tail);
      showRFatalError(tail);
    }
  }

  app.whenReady().then(boot);

  app.on('window-all-closed', () => app.quit());

  app.on('before-quit', () => {
    if (rProc && !rProc.killed) {
      try { rProc.kill('SIGTERM'); } catch (e) { /* ignore */ }
      setTimeout(() => {
        if (rProc && !rProc.killed) {
          try { rProc.kill('SIGKILL'); } catch (e2) { /* ignore */ }
        }
      }, 1500);
    }
  });

  // ── IPC ──────────────────────────────────────────────────────────────────

  ipcMain.on(`${IPC}quit`, () => app.quit());
  ipcMain.handle(`${IPC}version`, () => ({
    app: APP_VERSION,
    electron: process.versions.electron,
    node: process.versions.node,
  }));
  ipcMain.on(`${IPC}open-reference-folder`, () => {
    const { resolved } = resolveAndExportReferences();
    if (resolved.genome.dir) shell.openPath(resolved.genome.dir);
  });
  ipcMain.on(`${IPC}open-outputs-folder`, () => {
    const out = path.join(app.getPath('documents'), config.outputDirName || `${config.title} Outputs`);
    fs.mkdirSync(out, { recursive: true });
    shell.openPath(out);
  });
  ipcMain.handle(`${IPC}references`, () => {
    const { resolved } = resolveAndExportReferences();
    return {
      genome: {
        status: resolved.genome.status,
        dir: resolved.genome.dir,
        source: resolved.genome.source,
        species: resolved.genome.species,
      },
      pjat: {
        status: resolved.pjat.status,
        dir: resolved.pjat.dir,
        source: resolved.pjat.source,
      },
    };
  });
}

module.exports = { startApp };

  }
