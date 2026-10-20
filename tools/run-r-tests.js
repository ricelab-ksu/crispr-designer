#!/usr/bin/env node
/**
 * Run the payloads' own testthat suites (both apps).
 *
 * Usage: node tools/run-r-tests.js [app-id]   (omit app-id to run both)
 */

const { execFileSync } = require('child_process');
const path = require('path');

const REPO_ROOT = path.join(__dirname, '..');

function runTests(id) {
  const cwd = path.join(REPO_ROOT, 'apps', id, 'payload');
  console.log(`=== Running R tests for ${id} (cwd=${cwd}) ===`);
  execFileSync(
    'Rscript', ['--vanilla', '-e',
      'library(testthat); testthat::set_max_fails(Inf); testthat::test_dir("tests/testthat", package=NULL, stop_on_failure=FALSE)'],
    { cwd, stdio: 'inherit' });
}

const ids = process.argv.slice(2).filter((a) => !a.startsWith('--'));
if (ids.length) {
  ids.forEach(runTests);
} else {
  runTests('plasmid-designer');
  runTests('crispr-designer');
}
