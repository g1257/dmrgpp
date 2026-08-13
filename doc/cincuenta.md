# cincuenta documentation has moved

Design docs, diagrams, and input-format guides for `cincuenta` (the DMFT
driver) now live in [`cincuenta/doc/`](../cincuenta/doc/), next to its
source, rather than in this monorepo-wide `doc/` directory.

Start with [`cincuenta/doc/neq_dmft_ed_input.md`](../cincuenta/doc/neq_dmft_ed_input.md)
for the non-equilibrium DMFT ED input format, or
[`cincuenta/doc/cincuenta_design.md`](../cincuenta/doc/cincuenta_design.md)
for the overall (equilibrium) design.

**Backport note:** this branch is a synthetic/slice PR spun off `master`.
Its one piece of content that actually originated on `neq_dmft_merge_master_dev`
is `cincuenta_design.md` + `classes.mmd`/`component_arch.mmd`/`sequence_dmft.mmd`,
which on that dev branch still live at the old top-level `doc/` path (added
2026-05-15). When backporting this PR to `neq_dmft_merge_master_dev`, move
dev's copy of those four files to `cincuenta/doc/` (matching this PR) rather
than landing a second, duplicate copy at the old path. Everything else here
(`neq_dmft_ed_input.md`, this pointer file) is new — no backport conflict
for those. Expect further additions to `cincuenta/doc/` from later slice PRs
(the neq-DMFT/GBEK ED input design doc is planned to land alongside this
input-format guide), so re-check this note as those merge in.
