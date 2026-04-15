# Fine-tuned Basecaller Models

Nine Bonito/Dorado-compatible basecaller checkpoints fine-tuned from
`dna_r10.4.1_e8.2_400bps_sup@v5.2.0` on the Fc amplicon datasets described
in Chapter 2 of the thesis. These are the models benchmarked in Chapter 3.

## Naming convention

Each model is identified by `<family><variant>`:

| Family | Training POD5 source                                 |
|--------|-------------------------------------------------------|
| Alpha  | Run A (Fc amplicon, dataset A)                       |
| Beta   | Run B (Fc amplicon, dataset B)                       |
| Gamma  | Combined Run A + Run B (dataset AB)                  |

| Variant | Training corpus                         | `min_accuracy_save_ctc` |
|---------|-----------------------------------------|--------------------------|
| (base)  | Full POD5 (length-only filtering)       | 0                        |
| QCL     | Quality-control low stringency pool     | 0                        |
| QCH     | Quality-control high stringency pool    | 0.99                     |

The nine resulting checkpoints are: `Alpha`, `AlphaQCL`, `AlphaQCH`,
`Beta`, `BetaQCL`, `BetaQCH`, `Gamma`, `GammaQCL`, `GammaQCH`.

Per-model training statistics (chunk counts, learning rate, loss, validation
mean/median) are tabulated in `tables/model.csv`.

## Excluded artefacts

Two classes of Bonito artefact are **not** checked in:

- `weights_1.tar` (~315 MB per model) exceeds GitHub's 100 MB per-file
  push limit and has been removed from each `<name>_bonito/` folder. The
  Dorado-ready package under `<name>_dorado/` is a complete, self-contained
  runtime equivalent, so inference reproducibility is preserved.
- `data/` training-array directories (`basecalls.bam`, `chunks.npy`,
  `references.npy`, `reference_lengths.npy`, `basecalls_summary.tsv`) are
  excluded for size. They can be regenerated from the raw POD5s via
  `src/bash/DMSPolishing.sh` plus the helpers in `src/python/`
  (`Pod5Splitter.py`, `Pod5Merger.py`, `MergeNumpy.py`).

Access to the original `weights_1.tar` checkpoints and training arrays —
needed only to resume Bonito fine-tuning or to re-export the Dorado
package — follows the same path as the raw sequencing data (see the
top-level `README.md` § Data availability).

## Folder layout

```
models/
├── README.md                       (this file)
└── <name>/                         (one per model, e.g. Alpha/, BetaQCH/)
    ├── <name>_bonito/              Bonito training record (weights excluded)
    │   ├── config.toml             Bonito training config incl. provenance
    │   ├── losses_1.csv            Per-batch training log
    │   └── training.csv            Per-epoch train/validation summary
    └── <name>_dorado/              Bonito-exported, Dorado-ready package
        ├── config.toml             Runtime config (no training section)
        └── *.tensor                157 weight tensors (conv / transformer / CRF / upsample)
```

## Reproducing a basecall run

### Dorado (inference)

Pass the `<name>_dorado/` directory directly as the model argument:

```bash
dorado basecaller models/Alpha/Alpha_dorado input.pod5 > basecalls.bam
```

### Bonito (requires `weights_1.tar`)

Bonito inference and resumed fine-tuning both read `weights_1.tar`, which
is not shipped in this repository. Once obtained, place it back into the
appropriate `<name>_bonito/` folder and run the usual Bonito commands,
e.g.:

```bash
bonito basecaller models/Alpha/Alpha_bonito input.pod5 > basecalls.fastq
bonito train --pretrained models/Alpha/Alpha_bonito new_training_dir
```

### Re-exporting the Dorado package

With `weights_1.tar` restored into `<name>_bonito/`, run `bonito export`
against that directory to regenerate `<name>_dorado/`; see
`bonito export --help` for the exact flag layout of your Bonito version.

Exact flags used to produce the checkpoints in this repository (including
`--min-accuracy-save-ctc`, chunk counts, and learning rates) are recorded
in each model's `<name>_bonito/config.toml` under the `[training]` section
and summarised in `tables/model.csv`.

## Provenance

- Pretrained source: `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` (ONT)
- Architecture: Transformer CRF, 18 layers, 8 heads, 512 d_model,
  2048 FFN, state_len=5, attn_window [127, 128]
- Sample rate: 5000 Hz (R10.4.1, 400 bps)
- Fine-tuning hardware: NVIDIA RTX 5090
- Bonito version: 1.0.1 (see `environment/README.md`)
