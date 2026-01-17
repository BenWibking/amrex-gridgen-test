# AMReX grid generation iteration bug (DiskGalaxy reproducer)

## Summary
This repo reproduces a grid generation issue in AMReX where the regridding loop
stops before the grids converge. With zero error-buffer cells
(`amr.n_error_buf=0`), the refined region comes out too small even though the
tagging logic is correct. Increasing the iteration count in
`AmrMesh::MakeNewGrids` gives the grid generator enough passes to converge, and
the refined region matches the target cylinder (R < 16 kpc, |z| < 4 kpc).

## Reproducer
- Build and run:
  - `make -j`
  - `./main3d.gnu.ex input.buf0`
  - `./main3d.gnu.ex input.buf3`
  - `python3 plot_grids.py buf0.grids --inputs input.buf0 --output buf0.png --units kpc`
  - `python3 plot_grids.py buf3.grids --inputs input.buf3 --output buf3.png --units kpc`
- Expected: refined grids cover the cylindrical region defined in
  `testDiskGalaxyGrid.cpp` and the inputs.
- Observed:
  - `input.buf3` (`amr.n_error_buf=3`) produces the expected refinement volume.
  - `input.buf0` (`amr.n_error_buf=0`) yields a noticeably smaller refined
    region.

## Why this happens
- The tagger in `testDiskGalaxyGrid.cpp` marks cells whose corners intersect a
  cylinder. The geometry is fixed in physical space, but the exact set of tagged
  cells changes with resolution because the cell corners change.
- AMReX builds grids by:
  1) tagging cells on each level,
  2) buffering tags (`amr.n_error_buf`),
  3) clustering and enforcing proper nesting,
  4) iterating this process until the grid layout stops changing.
- With `amr.n_error_buf=0`, there is no padding to smooth out the boundary, so
  the first pass tends to under-cover the cylinder. Each subsequent iteration
  pushes the coarse-level grids outward a bit as fine grids are projected down
  and re-tagged, but the default limit (4 passes) is not enough for this
  `amr.max_level=8` case. The loop exits while the grids are still evolving, so
  the final grids are too small.

## Why increasing the iteration count fixes it
- The iteration cap lives in
  `amrex/Src/AmrCore/AMReX_AmrMesh.cpp` inside `AmrMesh::MakeNewGrids` (the
  `iterate_on_new_grids` loop).
- Increasing the cap from 4 to 10 gives the algorithm more passes to:
  - propagate fine-level tags back to coarser levels,
  - re-cluster with the updated tags and proper nesting,
  - converge to a stable BoxArray layout.
- With the higher limit, the grid layout stops changing before the cap and the
  refined volume matches the target geometry even with `amr.n_error_buf=0`.
  The `amr.n_error_buf=3` case already converges quickly because the tags are
  pre-expanded, which is why it does not show the bug.

## Notes
- The local `amrex/` checkout already contains the iteration-count change.
