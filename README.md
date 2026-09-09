# Wave Equation Project

Finite Element solver for the wave equation using deal.II and Newmark time integration.

## Source organization

All source files are contained in the `src` folder.

Main files:

- `src/main.cpp`
- `src/TimeADR.cpp`
- `src/TimeADR.hpp`
- `src/TestCases.hpp`
- `src/convergence.cpp`

## Build instructions

Load the required modules:

```bash
module load gcc-glibc dealii
```

Then compile with:

```bash
mkdir build
cd build
cmake ..
make
```

The executable is generated inside the `build` directory.

## Running

From inside the `build` directory, run:

```bash
./<executable-name>
```

Replace `<executable-name>` with the executable name defined in `CMakeLists.txt`.
