# COMPX593-Thesis

This repository will include all the source code for Kang's master project.

## Instructions

This code is distributed without an pre-build binary. Therefore it is necessary for the user to perform an out-of-source build on machines they want to execute the program.

See following command to perform out-of-source build (assume you have already cloned this directory and is located within the project root):

### For MacOS & Linux Builds

```bash
mkdir build
cmake -S . -B build -G Ninja
cmake --build build
```

### For Windows Builds

```powershell
New-Item -ItemType Directory -Name "build"
cmake -S . -B build -G "Visual Studio 17 2022"
cmake --build build --config Release
```



