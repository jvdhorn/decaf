A Python package for simulating diffuse scattering of X-rays in protein crystals

# Requirements
* Phenix/cctbx: runs in phenix.python (2.7.x) and depends on cctbx modules
* Numpy
* Scipy

# Usage
Make sure `phenix.python` is available.

Update `PATH` with
```
source source_env.sh
```

Command-line parameters can be provided using Phil syntax. Run any of the modules without arguments to get an overview of the avaiable parameters for that module.

# Overview of modules
* schimpy - SuperCell Hierarchical Model
* stimpy - Statistical Image-processing
* slicemtz - Viewing utility for mtz-files
* submtz - Subtract mtz-files
* plotmtz - Plot intensity distributions in mtz-files
* ccmtz - Calculate correlation coefficient and R1-values between mtz-files
* rmbragg - Remove supercell voxels around Bragg positions
* rmerge - Calculate R-int
* mtzstats - Show various statistics for mtz-files
* stimpy3d - Apply the stimpy-procedure to resolution bins in mtz-files
* pattsize - Estimate the size of the Patterson origin-peak
* mtz2txt - Extract raw intensities from mtz-files
* filter_mtz - Apply a kernel-filter to an mtz-file
* pdbrad - Estimate the size of a pdb-object

## schimpy
## stimpy
## slicemtz
## submtz
## plotmtz
## ccmtz
## rmbragg
## rmerge
## mtzstats
## stimpy3d
## pattsize
## mtz2txt
## filter_mtz
## pdbrad

# References
