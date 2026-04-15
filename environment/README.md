# Environment

Runtime dependencies for reproducing the COMPX593 thesis pipeline.

## Python (`requirements.txt`)

Used by `src/ipynb/plot.ipynb` and by the helper scripts in `src/python/`.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
```

## External bioinformatics tools

These must be on `PATH` for the Bash workflows in `src/bash/` and for the
per-position mismatch-density cells in `plot.ipynb`:

| Tool                  | Version used in the thesis |
|-----------------------|----------------------------|
| `samtools`            | 1.19.2                     |
| `minimap2`            | 2.28                       |
| `dorado`              | 1.3.0                      |
| `fastplong`           | 0.2.2                      |
| `bonito`              | 1.0.1                      |

## C++ toolchain

The executables in `src/c++/` are built with CMake (≥ 3.20) against a C++17
compiler and `htslib`. See the top-level `README.md` for per-platform build
instructions.

## GPU / model weights

Fine-tuning and basecalling used the ONT release
`dna_r10.4.1_e8.2_400bps_sup@v5.2.0` on an Nvidia RTX 5090. Representative
Bonito commands are given in Appendix D of the thesis.
