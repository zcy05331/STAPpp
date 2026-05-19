# Q4 validation example

The validation case reuses a 4x2 uniaxial membrane with an exact independent continuum solution. Expected stress is sigma_x=1, sigma_y=0, tau_xy=0 and displacement field is u=x/E, v=-nu*y/E.

Max stress error: 3.682050e-16.
Max displacement error: 1.730770e-19.

## ParaView post-processing

`q4_validation_uniaxial.vtk` is a legacy VTK unstructured-grid file generated from the STAPpp `.dat/.out` pair and can be opened directly in ParaView. It contains deformed coordinates (scale 100), nodal displacement vectors, cell-average `sigma_x`, and cell-average plane-stress von Mises values.

Regenerate with:

```bash
python3 tools/q4_to_vtk.py data/q4_validation/q4_validation_uniaxial.dat data/q4_validation/q4_validation_uniaxial.out data/q4_validation/q4_validation_uniaxial.vtk --scale 100
```
