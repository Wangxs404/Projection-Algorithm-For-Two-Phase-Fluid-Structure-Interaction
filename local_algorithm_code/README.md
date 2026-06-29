# Local Algorithm Code Archive

Generated: 2026-06-29T10:20:54

This branch keeps the published repository root unchanged and adds local algorithm-code snapshots under `local_algorithm_code/`.

## Source sets included

- `old_phase_field_ibm_local_repo/`: local source from the older nested repository `Two-Phase-Fluid-Structure-Interaction-Using-Phase-Field-IBM`.
- `bubble_case_code/Unity VC-Public-IBM-Implicit-c/`: local bubble-case algorithm source.
- `bubble_case_code/Unity VC-Public-IBM-Implicit-DiracDensity/`: local Dirac-density bubble-case algorithm source.
- `code_snapshot_unique_sources/Unity VC-Public-IBM-Implicit-c/`: source files from the local snapshot that survived duplicate slimming.

## Source sets not duplicated here

The main projection-algorithm code at this branch root already matches the local nested repository:

`/Volumes/Elements/windows-docu/Research/Publications/[1] TwoPhase-FSI/05_Code/git_worktree_dirty/Unity VC-Public-IBM-Implicit-c - git/Projection-Algorithm-For-Two-Phase-Fluid-Structure-Interaction`

That local repository's committed tree matched `origin/main`; its working-tree modifications were CRLF/whitespace only.

## Exclusions

Simulation data and binary state files (`.dat`, `.mat`) were intentionally not added to this Git branch. This branch is for source code, not raw simulation outputs.

See `MANIFEST.tsv` for original paths and SHA256 checksums.
