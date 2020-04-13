# Coord Mover

Visually place coordinates on a scout orientated t1 image and generate spectra from `siarray.1.1` MRS Images.

The simplest invocation is to clone and cd this directory
```matlab
coord_mover FF
```

![](../imgs/ml_cordmover.png?raw=True)
## Features

  * Can use FS GM nifti to report total GM voxels
     * The `GM1` button will move the current roi to hight GM in a 5 pixel radius
     * The `GM` button will do the same for all voxels
  * ROIs placed too close will get a red radius
  * `Spectrum` button will generate spectrum for all placed coordinates
  * arrows and mouse buttons have convenience bindings
    * right click to select closest roi
    * `u` to undo
    * arrow keys to nudge roi box
    * scroll wheel - 

## Prerequisites

  * `rorig.nii` is the T1 image registered to the scout image used to align the MRSI acquisition
  * to use GM counts, `rgm.nii` is FS gray matter mask also registered to the scout
  * to make a spectrum, the `siarray.1.1` file is needed
