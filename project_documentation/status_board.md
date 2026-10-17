# Status Board — Live Ticket Tracker

> Orchestrator updates this after every review. Workflow states: Ready / To do / In progress /
> In review / Done / Re-ticketed (defects listed).

| Ticket | Owner | Status | Notes |
|--------|-------|--------|-------|
| T1 — Tool survey & ranking | Research subagent | **Done** | 10 tools ranked in `tool_survey.md` from local `awesome-CRISPR/README.md` + GitHub |
| T2 — Genome inventory & loader | `research-drosophila` / `dev-backend` | **In progress** | Dispatched to Agent Team (`cas9-sgrna`) & Kanban repo |
| T3 — Gene-lookup index & API | `research-bioinformatics` / `dev-backend` | **In progress** | Dispatched to Agent Team (`cas9-sgrna`) & Kanban repo |
| T4 — gRNA generation & ranking | `research-bioinformatics` / `dev-backend` | **Queued (Wave 2)** | Depends on T2, T3 |
| T5 — Homology-arm designer | `research-drosophila` / `dev-backend` | **Queued (Wave 2)** | Depends on T4 |
| T6 — Shiny GUI & Electron shell | `dev-frontend` | **Queued (Wave 3)** | Depends on T3, T4, T5; reference app at `HCR Probe Designer.app` |
| T7 — End-to-end validation | `dev-reviewer` / `research-reviewer` | **Queued (Wave 3)** | Depends on T3–T6 |

## Review Log
- **2026-10-16:** Completed T1 (`tool_survey.md`) and initial T2 audit (`genome_inventory.md`); updated all paths to `/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design` and `/Volumes/Backup_Plus/projects/bioinformatics/reference_files/drosophila_species_genomes`.
- **2026-10-16:** Initialized git repository at `/Volumes/Backup_Plus/projects/bioinformatics/projects/Cas9_sgRNA_design`, launched Kanban board (`launch-kanban.sh`) and specialist Agent Team `cas9-sgrna` (`launch-team.sh` with roles `research-bioinformatics,research-drosophila,dev-backend,dev-frontend,dev-reviewer`) to execute tickets T2–T7.


