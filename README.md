A Python package for simulating diffuse scattering of X-rays in protein crystals

# Requirements
* Phenix/cctbx: runs in phenix.python (2.7.x) and depends on cctbx modules
* Numpy
* Scipy

# Usage
Make sure `phenix.python` is available.

Add package to `PATH` with
```
source source_env.sh
```

Command-line parameters can be provided using Phil syntax. Run any of the modules without arguments to get an overview of the available parameters for that module.

# Overview of modules
* `schimpy` - SuperCell Hierarchical Model
* `stimpy` - Statistical Image-processing
* `slicemtz` - Viewing utility for mtz-files
* `submtz` - Subtract mtz-files
* `plotmtz` - Plot intensity distributions in mtz-files
* `ccmtz` - Calculate correlation coefficient and R1-values between mtz-files
* `rmbragg` - Remove supercell voxels around Bragg positions
* `rmerge` - Calculate R-int
* `mtzstats` - Show various statistics for mtz-files
* `stimpy3d` - Apply the stimpy-procedure to resolution bins in mtz-files
* `pattsize` - Estimate the size of the Patterson origin-peak
* `mtz2txt` - Extract raw intensities from mtz-files
* `filter_mtz` - Apply a kernel-filter to an mtz-file
* `pdbrad` - Estimate the size of a pdb-object

## Most important command-line parameters
### schimpy
* `pdb_in` - refined structure file (pdb)
* `tls_in` - result of the ECHT B-factor distribution (json)
* `sc_size` - size of the supercell (e.g. "5 5 10")
* `correlate` - enable or disable correlation of TLS-groups (default True)
* `stretch` - stretch parameter for the anchor points (default 0.25)
* `cutoff` - distance cutoff for intermolecular interactions (default 3)
* `weights` - power for the number of interactions (default 1.5)
* `max_level` - highest level of the TLS-hierarchy to include (default 1)
* `high_resolution` - high resolution cutoff (default 2)
* `k_sol` and `b_sol` - bulk solvent parameters (default 0.35 and 50)
* `processes` - number of parallel processes (default 1, watch memory usage!)
* `models` - number of models to simulate (default 128)

### stimpy
* `image` - raw image file
* `polar` - polarization-corrected intermediate background
* `bin_photons` - bin regions within this range of counts (default 1)
* `N` - expected number of independent rotations (default 1)
* `bragg_mask_median_filter` - kernel size of the median filter (default 9)
* `bragg_mask_dilation` - kernel size of the mask dilation (default 7)

### slicemtz
* `mtz_1` - input mtz-file
* `lbl_1` - array of interest (default IDFF)
* `slice` - desired slice to plot (default hk0)
* `depth` - additional depth of the slice on both sides (default 0)
* `sc_size` - size of the supercell for drawing Bragg positions (e.g. "5 5 10")
* `log` - plot log10 of the intensities (default False)
* `min` and `max` - minimum and maximum values to plot

### submtz
### plotmtz
### ccmtz
### rmbragg
### rmerge
### mtzstats
### stimpy3d
### pattsize
### mtz2txt
### filter_mtz
### pdbrad

# References
