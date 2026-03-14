# COMPX593 Thesis Toolkit

This repository contains the source code used for the COMPX593 thesis project.
It combines small C++ command-line tools for FASTQ and BAM processing with
Python and Bash utilities for POD5 handling, signal visualisation, and
end-to-end DMS polishing workflows.

## Repository Layout

- `src/c++`: standalone C++ executables built with CMake
- `src/bash`: shell pipelines and reporting helpers
- `src/python`: POD5, NumPy, and plotting utilities
- `build`: optional out-of-source build directory

## Build Requirements

The C++ tools require:

- a C++17 compiler
- CMake 3.20 or newer
- `pkg-config`
- `htslib`

Several scripts also expect command-line bioinformatics tools such as
`samtools`, `minimap2`, `dorado`, and `fastplong`. Python utilities depend on
packages such as `numpy`, `matplotlib`, `pod5`, and optionally `tqdm`.

## Building the C++ Tools

Perform an out-of-source build from the project root.

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

## Notes

- Executable names match their source filenames under `src/c++`.
- The Bash workflows assume the C++ tools are either on `PATH` or available in
  the repository `build` directory.
- The Python utilities are written as standalone entry-point scripts and can be
  executed directly with `python`.
