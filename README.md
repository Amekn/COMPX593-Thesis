# COMPX593 Thesis Toolkit

Source code, summary tables, figures, and reference inputs used for the
COMPX593 Master's thesis *"Fine-Tuning Oxford Nanopore Basecalling Models to Enable High-Fidelity Single-Read Variant Calling in Antibody Libraries"* (University
of Waikato, 2026).

The repository bundles:

- Custom C++ command-line tools for UMI extraction, DMS-aware polishing,
  and variant-concordance scoring against FASTQ/BAM inputs.
- Bash workflows that drive the basecalling → filtering → alignment →
  polishing → re-alignment → metrics pipeline described in Chapter 2.
- Python utilities for POD5 preprocessing and Bonito training-array
  construction.
- A single Jupyter notebook (`src/ipynb/plot.ipynb`) that renders every
  figure, panel, and summary table used in the thesis body and appendices.
- Summary tables, the Fc reference, and all rendered figure SVGs required
  to re-verify the published results.

## Repository Layout

```
COMPX593-Thesis/
├── CMakeLists.txt        # Build definition for the C++ tools
├── LICENSE               # GPL-3.0
├── README.md             # (this file)
├── environment/          # Python + external-tool dependency pins
├── figures/              # Rendered thesis figures (P1–P23, F1–F3, ELISA)
├── reference/            # Fc amplicon reference + DMS-window definition
├── src/
│   ├── bash/             # Pipeline drivers and table generators
│   ├── c++/              # 16 CLI tools + shared headers
│   │   └── include/ont_tools/
│   ├── python/           # POD5 + NumPy utilities
│   └── ipynb/
│       └── plot.ipynb    # Master notebook for all figures and tables
└── tables/               # Summary CSVs consumed by the notebook
```

Sequencing inputs (POD5, BAM, FASTQ), per-model test outputs, the nine
fine-tuned Bonito model weights, and derived caches are **not** tracked in
git (see [Data availability](#data-availability) below).

## Building the C++ tools

The C++ executables require a C++17 compiler, CMake ≥ 3.20, `pkg-config`,
and `htslib`. Build out-of-source from the project root.

### Linux

```bash
cmake -S . -B build -G Ninja
cmake --build build
```

### macOS

```bash
cmake -S . -B build -G Xcode
cmake --build build --config Release
```

### Windows

```powershell
cmake -S . -B build -G "Visual Studio 17 2022"
cmake --build build --config Release
```

Executable names match the source filenames under `src/c++`. Binaries are
written to `build/` and are not committed.

## Python environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
```

External bioinformatics tools expected on `PATH`: `samtools`, `minimap2`,
`dorado`, `fastplong`, and `bonito`. Exact versions are recorded in
`environment/README.md`.

## Reproducing the thesis figures

1. Build the C++ tools as described above.
2. Install the Python environment.
3. Obtain the sequencing inputs and per-model outputs (see
   [Data availability](#data-availability)) and place them under
   `src/ipynb/plot/{data,test,model}` (paths resolved relative to the repo
   root by the notebook's `ROOT` variable).
4. Open `src/ipynb/plot.ipynb` and run all cells. Rendered SVGs are written
   into `figures/` with filenames matching the panel titles printed in the
   thesis (e.g. `P1 - Read retention across filtering stages by dataset.svg`).
5. Summary tables consumed by the notebook live in `tables/`:
   `model.csv`, `data.csv`, `test.csv`, `correlation.csv`, `mutation.csv`,
   `umi.csv`.

## Thesis chapter ↔ code map

| Section                                | Code / artefact                                                                                                        |
|----------------------------------------|------------------------------------------------------------------------------------------------------------------------|
| §2.2.2 Basecalling, filtering, align.  | `src/bash/DMSPolishing.sh`, `src/c++/PrimerTrimmer.cpp`                                                                 |
| §2.2.3 UMI processing                  | `src/c++/UmiExtractor.cpp`, `UmiCloner.cpp`, `UmiSpliter.cpp`, `UmiFilter.cpp`, `UmiOverlap.cpp`, `SingleUmiOverlap.cpp`|
| §2.2.4 Model fine-tuning               | `src/python/Pod5Splitter.py`, `Pod5Merger.py`, `MergeNumpy.py`; `environment/README.md` (Bonito command template)       |
| §2.2.5 DMS-aware correction            | `src/c++/DualSiteDMSFilter.cpp`                                                                                        |
| §2.2.6 Benchmarking & concordance      | `src/c++/Benchmarker.cpp`, `VariantConcordance.cpp`, `VariantConcordanceCli.cpp`, `VariantKeyOverlap.cpp`, `MQAAR.cpp`, `ParseStats.cpp` |
| §3.1 Dataset QC tables                 | `src/bash/GenerateDataCsv.sh` → `tables/data.csv`                                                                      |
| §3.2–3.3 Per-model metrics             | `src/bash/CalcStats.sh`, `GenerateTestCsv.sh` → `tables/test.csv`                                                       |
| §3.4 Variant concordance               | `src/bash/GenerateCorrelationCsv.sh` → `tables/correlation.csv`                                                        |
| All figures & panels                   | `src/ipynb/plot.ipynb` → `figures/*.svg`                                                                               |
| Reference & DMS windows                | `reference/fc_reference.fa(.fai/.mmi)`, `reference/DMSZones.txt`                                                        |

## Data availability

The 687 bp Fc amplicon reference and the four DMS window coordinates are
checked into `reference/`. The raw nanopore POD5s, basecalled FASTQs, BAMs,
Bonito training arrays, per-model test artefacts, and fine-tuned model
weights are **excluded from the repository** (see `.gitignore`) because of
their size. Access can be arranged by contacting the author via the
thesis submission record held by the University of Waikato School of
Computing and Mathematical Sciences.

## License

Released under the terms of the GNU General Public License v3.0 (see
`LICENSE`).
