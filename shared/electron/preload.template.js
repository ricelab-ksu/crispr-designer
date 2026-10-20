/**
 * Shared preload script.
 *
 * Exposes a minimal, safe IPC surface to the Shiny renderer. The channel
 * names are namespaced per app, so the preload reads the app id from the
 * package metadata that Electron injects.
 */

const { contextBridge, ipcRenderer } = require('electron');

// The app id is baked in at build time by tools/stage-app.js, which copies
// this file and rewrites the placeholder.
const APP_ID = '__APP_ID__';
const IPC = `${APP_ID}:`;

contextBridge.exposeInMainWorld('__APP__', {
  /** Open the writable reference-data folder (species genome store). */
  openReferenceFolder: () => ipcRenderer.send(`${IPC}open-reference-folder`),
  /** Open the folder where design outputs are written. */
  openOutputsFolder: () => ipcRenderer.send(`${IPC}open-outputs-folder`),
  /** Close the app. */
  quit: () => ipcRenderer.send(`${IPC}quit`),
  /** Version info for the splash screen. */
  version: () => ipcRenderer.invoke(`${IPC}version`),
  /** Current reference-store locations. */
  references: () => ipcRenderer.invoke(`${IPC}references`),
  /** Subscribe to fatal bootstrap errors. */
  onError: (cb) => {
    ipcRenderer.on(`${IPC}error`, (_e, msg) => cb(msg));
  },
  /** Subscribe to reference-store resolution updates. */
  onReferences: (cb) => {
    ipcRenderer.on(`${IPC}references`, (_e, info) => cb(info));
  },
});
