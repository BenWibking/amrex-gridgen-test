# amrex-gridgen-test

* `inputs.buf3` gives a set of AMR grids that respect the desired tagging:
  R < 16 kpc and |z| < 4 kpc.
* `inputs.buf0` has too small of a refined volume.

```
make -j
./main3d.gnu.ex input.buf0
./main3d.gnu.ex input.buf3
python3 plot_grids.py buf0.grids --inputs input.buf0 --output buf0.png --units kpc
python3 plot_grids.py buf3.grids --inputs input.buf3 --output buf3.png --units kpc
```

Runtime debugging/tuning:
- `amr.dynamic_grid_iterations = 1` enables the dynamic convergence loop
  (cycle detection + 50-iteration fail-safe). Default is the original fixed
  4-iteration path for backward compatibility.
- `amr.debug_grid_iterations = 1` prints the iteration count and the reason
  the iteration loop stopped.

The problem is caused by this line:
https://github.com/AMReX-Codes/amrex/blob/8ae62c47e1d6660cb1a3926b9c2de4795b3d3c13/Src/AmrCore/AMReX_AmrMesh.cpp#L1013
