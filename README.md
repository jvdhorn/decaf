A Python package for simulating diffuse scattering of X-rays in protein crystals

# Requirements
* cctbx/phenix: runs in cctbx.python (2.7.x / 3.x) and depends on cctbx modules
* numpy
* scipy
* matplotlib

# Usage
Make sure `cctbx.python` is available, which can be obtained [here](https://github.com/cctbx/cctbx_project?tab=readme-ov-file#installation).

Clone the repository
```
git clone https://github.com/jvdhorn/diffuse
```

Add package to `PATH` with
```
source source_env.sh
```

Optional: install package with
```
cctbx.python -m pip install .
```

Command-line parameters can be provided using Phil syntax. For example:
```
schimpy pdb_in=structure.pdb tls_in=optimised_model.json
```

Run any of the modules without arguments to get an overview of the available parameters for that module.

# Overview of modules
* `schimpy` - SuperCell Hierarchical Model
* `stimpy` - Statistical Image-processing
* `slicemtz` - Viewing utility for mtz-file
* `submtz` - Subtract two mtz-files
* `plotmtz` - Plot intensity distributions in mtz-files
* `ccmtz` - Calculate correlation coefficient and R1-values between mtz-files
* `rmbragg` - Remove supercell voxels around Bragg positions in mtz-file
* `mtzstats` - Show various statistics including R-int for mtz-file
* `stimpy3d` - Apply the stimpy-procedure to resolution bins in mtz-file
* `pattsize` - Estimate the size of the Patterson origin-peak
* `mtz2txt` - Extract raw intensities from mtz-file
* `map2mtz` - Convert map-file with experimental intensities to mtz
* `filter_mtz` - Apply a kernel-filter to an mtz-file
* `pdbrad` - Estimate the size of a pdb-object

## Most important command-line parameters
### schimpy
* `pdb_in` - refined structure file (pdb)
* `tls_in` - result of the ECHT B-factor distribution (json)
* `sc_size` - size of the supercell (e.g. `"5 5 10"`)
* `correlate` - enable or disable correlation of TLS-groups (default `True`)
* `stretch` - stretch parameter for the anchor points (default `0.25`)
* `cutoff` - distance cutoff for intermolecular interactions (default `3.0`)
* `weights` - power for the number of interactions (default `1.5`)
* `max_level` - highest level of the TLS-hierarchy to include (default `1`)
* `resolution` - resolution limit (default `2.0`)
* `k_sol` and `b_sol` - bulk solvent parameters (default `0.35` and `50.0`)
* `processes` - number of parallel simulations (default `1`, watch memory usage!)
* `interval` - number of seconds between consecutive simulations (default `1.0`)
* `n_models` - number of supercells to simulate (default `128`)

### stimpy
* `image` - raw image file
* `radial` - polarization-corrected radial average
* `bin_counts` - bin regions within this range of counts (default `1.0`)
* `N` - expected number of independent rotations (default `1.0`)
* `median_size` - kernel size of the median filter (default `9`)
* `dilation_size` - kernel size of the mask dilation (default `7`)

### slicemtz
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `slice` - desired slice to plot (default `hk0`)
* `depth` - additional depth of the slice on both sides (default `0`)
* `sc_size` - size of the supercell for drawing Bragg positions (e.g. `"5 5 10"`)
* `log` - plot log10 of the intensities (default `False`)
* `min` and `max` - minimum and maximum values to plot

### submtz
* `mtz_1` - first input mtz-file
* `mtz_2` - second input mtz-file
* `lbl_1` - first array of interest (default `IDFF`)
* `lbl_2` - second array of interest (default `IDFF`)
* `add` - add instead of subtract (default `False`)
* `scale` - scale factor for second mtz-file (default `1.0`)

### plotmtz
* `mtz` - input mtz-file (can be multiple)
* `lbl` - array of interest (default `IDFF`)
* `log` - plot logarithmic distributions (default `True`)
* `resolution` - set low and high resolution (eg `"3.6 3.4"`)

### ccmtz
* `mtz_1` - first input mtz-file
* `mtz_2` - second input mtz-file
* `lbl_1` - first array of interest (default `IDFF`)
* `lbl_2` - second array of interest (default `IDFF`)
* `bins` - number of resolution shells (default `10`)
* `hlim` and `klim` and `llim` - limit h, k, and l (e.g. `"-20 20"`)

### rmbragg
* `mtz` - input mtz-file
* `sc_size` - supercell size (e.g. `"5 5 10"`)
* `box` - number of voxels to remove around every Bragg position (e.g. `"1 1 1"`)
* `keep` - invert selection (default `False`)

### mtzstats
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `sg` - Space group for R-int (symbol or number, e.g. `P43212` or `96`)

### stimpy3d
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `bins` - number of resolution shells (default `1`)
* `N` - expected number of independent rotations (default `1.0`)

### pattsize
* `map` - input patterson map
* `binsize` - size of radial bins in Angstrom (default `1.0`)
* `sigma` - sigma cutoff to determine size (default `2.0`)

### mtz2txt
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `resolution` - set low and high resolution (eg `"3.6 3.4"`)

### map2mtz
* `map` - input map-file
* `resolution` - set resolution limit (eg `2.0`)

### filter_mtz
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `size` - filter size (default `1`)
* `filter` - filter type (`gaussian` or `uniform`, default `gaussian`)

### pdbrad
* `pdb` - input pdb-file

# References
