# Q4 cantilever convergence

Reference is the 32x8 mesh tip displacement.

| nx | ny | tip uy | rel error | observed order |
|---:|---:|---:|---:|---:|
| 2 | 1 | -1.0000000000e-01 | 1.650270e-01 | - |
| 4 | 2 | -1.8661600000e-01 | 7.841100e-02 | 1.0736 |
| 8 | 2 | -2.3717000000e-01 | 2.785700e-02 | 1.4930 |
| 16 | 4 | -2.5867900000e-01 | 6.348000e-03 | 2.1337 |

## ParaView post-processing

The 32x8 reference mesh is used for the report displacement and von Mises contour figures. Generate the deformed-coordinate VTK file with a visual scale factor of 3:

```bash
python3 make/q4_to_vtk.py data/q4_convergence/q4_cantilever_32x8.dat data/q4_convergence/q4_cantilever_32x8.out data/q4_convergence/q4_cantilever_32x8.vtk --scale 3
```
