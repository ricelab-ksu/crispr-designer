#!/bin/bash
# ===========================================================================
# launch_app.sh — One-click launcher for the Plasmid Designer Shiny app
#                 (vector design / Gateway donor-fragment assembly)
# ===========================================================================
# Double-click (or run from Terminal) to:
#   1. Start the Plasmid Designer Shiny app on a local server
#      (http://127.0.0.1:7790)
#   2. Auto-open your default web browser to the app
#   3. Keep running in this Terminal window (Ctrl+C to stop)
# ===========================================================================

# Preflight + launch via the bundled launcher
cd "$(dirname "$0")"
exec '/Volumes/Backup_Plus/projects/bioinformatics/projects/Plasmid_Designer/Plasmid Designer.app/Contents/MacOS/plasmid-designer.sh' "$@"
