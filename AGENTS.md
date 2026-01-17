# Repository Guidelines

## Project Structure & Module Organization
- `testDiskGalaxyGrid.cpp` is the primary AMReX driver and grid-tagging logic.
- `GNUmakefile` and `Make.package` define the AMReX build configuration.
- `amrex/` is a git submodule with the AMReX source.
- `input.buf0` and `input.buf3` are runtime ParmParse inputs for grid setup.
- `plot_grids.py` generates plots; `outputs/` holds sample `*.grids` and `*.png`.
- `tmp_build_dir/`, `main3d.gnu.ex`, and `buf*.grids`/`buf*.png` are build/run artifacts.

## Build, Test, and Development Commands
- `git submodule update --init --recursive` fetch the AMReX dependency.
- `make -j` builds the executable `main3d.gnu.ex`.
- `./main3d.gnu.ex input.buf3` runs the grid generator using the specified inputs.
- `./main3d.gnu.ex input.buf0` reproduces the smaller-refined-volume behavior.
- `python3 plot_grids.py buf3.grids --inputs input.buf3 --output buf3.png --units kpc`
  visualizes the resulting grid hierarchy.

## Configuration & Inputs
- Refinement region is controlled by `agora_galaxy.refine_Rmax_kpc` and
  `agora_galaxy.refine_zmax_kpc` in the input files.
- `amr.n_error_buf` differs between `input.buf0` and `input.buf3` and affects tagging.
- Output filenames are set via `disk_galaxy_grid.gridfile`.

## Coding Style & Naming Conventions
- C++ only; keep standard headers before AMReX headers.
- Use CamelCase for classes/methods (e.g., `DiskGalaxyGrid`, `PrintGrids`) and
  snake_case for local variables (e.g., `refine_Rmax_kpc`).
- Match indentation in the file you edit; current code uses tabs in class methods
  and 2-space indentation in `main`.
- Preserve AMReX idioms (`amrex::ParallelFor`, `AMREX_GPU_DEVICE`).

## Testing Guidelines
- No automated tests; validate by running the binary and plotting grids.
- Use `input.buf3` to check expected refinement (R < 16 kpc, |z| < 4 kpc).
- Update reference plots in `outputs/` only when behavior changes intentionally.

## Commit & Pull Request Guidelines
- Follow existing history: short, imperative summaries (e.g., "Update README.md").
- In PRs, include the purpose, commands run, and before/after plots if grids change.
