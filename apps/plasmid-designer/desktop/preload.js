const {contextBridge, ipcRenderer} = require('electron');

// Expose protected functions to the renderer process
contextBridge.exposeInMainWorld('electronAPI, {
  lookupReference: (query) => {
    ipcRenderer.send('reference-lookup', query);
    return true;
  },
  onRReady: (callback) => {
    ipcRenderer.on('r-status', (e, data) => callback(data));
  }
};
