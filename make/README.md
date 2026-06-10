# Build, validation, and post-processing helpers

The source keeps CMake configuration in `src/CMakeLists.txt`, matching the original STAPpp layout. Use a fresh temporary build directory to avoid stale local cache files:

```bash
tmpbuild=$(mktemp -d /tmp/stappp-build-XXXXXX)
cmake -S src -B "$tmpbuild" -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build "$tmpbuild" -- -j2
```

Run all regression and Q4 course-project checks:

```bash
python3 make/validate_q4_cases.py --exe "$tmpbuild/stap++"
```

Regenerate the ParaView legacy VTK outputs:

```bash
python3 make/q4_to_vtk.py \
  data/q4_validation/q4_validation_uniaxial.dat \
  data/q4_validation/q4_validation_uniaxial.out \
  data/q4_validation/q4_validation_uniaxial.vtk \
  --scale 100

python3 make/q4_to_vtk.py \
  data/q4_patch_shear/q4_patch_shear.dat \
  data/q4_patch_shear/q4_patch_shear.out \
  data/q4_patch_shear/q4_patch_shear.vtk \
  --scale 80

python3 make/q4_to_vtk.py \
  data/q4_plate_hole/q4_plate_hole_tension.dat \
  data/q4_plate_hole/q4_plate_hole_tension.out \
  data/q4_plate_hole/q4_plate_hole_tension.vtk \
  --scale 20

python3 make/q4_to_vtk.py \
  data/q4_convergence/q4_cantilever_32x8.dat \
  data/q4_convergence/q4_cantilever_32x8.out \
  data/q4_convergence/q4_cantilever_32x8.vtk \
  --scale 3
```

Render the report figures with ParaView:

```bash
pvpython make/render_q4_paraview.py
```

On macOS, if ParaView is installed as an app bundle and `pvpython` is not on `PATH`, use the bundled executable path, for example:

```bash
/Applications/ParaView-6.1.1.app/Contents/bin/pvpython make/render_q4_paraview.py
```
