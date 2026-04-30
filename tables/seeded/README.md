# Seeded Run Tables

This directory contains the compact, tracked outputs needed to reproduce the
seeded-run figure inputs without committing the full Bonito run workspace.

- `seed_manifest.csv`: the ten RNG seeds reused across seeded model families.
- `run_manifest.csv`: seeded model run metadata, including model family,
  variant, seed, training input, and training hyperparameters.
- `threshold_long.csv`: cached threshold-level concordance metrics for the
  seeded model runs used by the final concordance ranking, heatmap, and ECDF
  figures.

The raw seeded-run workspace (`seeded_bonito_v1`) contains model checkpoints,
FASTQ/BAM files, and per-run test work products. It is intentionally excluded
from git because it is too large for the repository; the checked-in cache and
manifests are the repository-sized provenance for the final figures.
