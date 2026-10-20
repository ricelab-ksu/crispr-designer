// Bootstrap for CRISPR Designer
// Mirrors HCR's approach: splash → spawn R → attach via IPC
const {app, BrowserWindow, ipcMain} = require('electron');
const path = require('path');
const url = require('url');

let mainWindow;

function createWindow () {
  mainWindow = new BrowserWindow({
    width: 1200,
    height: 800,
    webPreferences: {
      preload: path.join(__dirname, 'preload.js'),
      nodeIntegration: false,
      contextIsolation: true
    }
  });

  mainWindow.loadURL(url.format({
    pathname: path.join(__dirname, 'index.html'),
    protocol: 'file:',
    slashes: true
  }));

  mainWindow.on('closed', function () {
    mainWindow = null;
  });
}

app.whenReady().then(() => {
  createWindow();

  app.on('activate', function () {
    if (BrowserWindow.getAllWindows().length === 0) {
      createWindow();
    }
  });
});

app.on('window-all-close', function () {
  if (process.platform !== 'macOS') {
    app.quit();
  }
});

ipcMain.on('r-ready', (event) => {
  event.sender.send('r-status', {status: 'ready'});
});

ipcMain.on('reference-lookup', (event, query) => {
  event.sender.send('reference-result', {query, status: 'pending'});
});
