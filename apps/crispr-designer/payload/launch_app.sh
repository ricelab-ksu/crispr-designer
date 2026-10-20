#!/bin/bash
# ============================================================================
# launch_app.sh — One-click launcher for the CRISPR Cas9 gRNA + Homology-Arm
#                 Designer Shiny app (Drosophila, v1.0.0)
# ============================================================================
# Double-click (or run from Terminal) to:
#   1. Start the Shiny app on a local server (http://127.0.0.1:<PORT>)
#   2. Auto-open your default web browser to the app
#      (Cas9_sgRNA_design.R already passes launch.browser = TRUE to shinyApp)
#   3. Keep running in this Terminal window (Ctrl+C to stop)
#
# This works because Cas9_sgRNA_design.R, when run via Rscript with NO --gene
# argument, falls through to interactive mode:
#     shiny::shinyApp(ui, server, options = list(port = PORT, launch.browser = TRUE))
# ============================================================================

set -e
cd "$(dirname "$0")"

echo "==============================================="
echo " CRISPR Cas9 gRNA + Homology-Arm Designer"
echo " for Drosophila (v1.0.0)"
echo "==============================================="
echo ""
echo "Starting Shiny app on a local server..."
echo "(This Terminal window must stay open while you use the app."
echo " Press Ctrl+C here to stop the app when you are done.)"
echo ""

exec Rscript Cas9_sgRNA_design.R

