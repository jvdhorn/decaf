# DECaF - Diffuse: Extract & Calculate F<sup>2</sup>

A Python package for extracting and simulating diffuse scattering of X-rays in protein crystals

## How to cite
If you use this software, please cite the following paper:
TBC

# Requirements
* [phenix/cctbx](https://github.com/cctbx/cctbx_project?tab=readme-ov-file#computational-crystallography-toolbox): runs in cctbx.python (2.7.x / 3.x) and depends on cctbx modules
* [dials/dxtbx](https://github.com/dials/dxtbx?tab=readme-ov-file#diffraction-experiment-toolbox): used for reading crystallographic image files
* [numpy](https://github.com/numpy/numpy?tab=readme-ov-file)
* [scipy](https://github.com/scipy/scipy?tab=readme-ov-file)
* [matplotlib](https://github.com/matplotlib/matplotlib?tab=readme-ov-file)

# Usage
Make sure `cctbx.python` is available and the `dxtbx` package is installed.

Clone the repository
```
git clone https://github.com/jvdhorn/decaf
```

Update `PATH` and `PYTHONPATH`
```
source source_env.sh
```

Optional: install package
```
cctbx.python -m pip install .
```

Command-line parameters can be provided using Phil syntax. For example:
```
schimpy pdb_in=structure.pdb tls_in=optimised_model.json
```

Run any of the modules without arguments to get an overview of the available parameters for that module.

# Overview of modules
* [`schimpy`](#schimpy) - SuperCell Hierarchical Model
* [`stimpy`](#stimpy) - Statistical Image-processing
* [`slicemtz`](#slicemtz) - Viewing utility for mtz-file
* [`submtz`](#submtz) - Subtract two mtz-files
* [`plotmtz`](#plotmtz) - Plot intensity distributions in mtz-files
* [`ccmtz`](#ccmtz) - Calculate correlation coefficient and R1-values between mtz-files
* [`rmbragg`](#rmbragg) - Remove supercell voxels around Bragg positions in mtz-file
* [`mtzstats`](#mtzstats) - Show various statistics including R-int for mtz-file
* [`stimpy3d`](#stimpy3d) - Apply the stimpy-procedure to resolution bins in mtz-file
* [`patterson`](#patterson) - Calculate Patterson map from mtz-file
* [`pattsize`](#pattsize) - Estimate the size of the Patterson origin-peak
* [`mtz2txt`](#mtz2txt) - Extract raw intensities from mtz-file
* [`mtz2map`](#mtz2map) - Convert mtz to ccp4 map-file
* [`map2mtz`](#map2mtz) - Convert ccp4 map-file to mtz
* [`filter_mtz`](#filter_mtz) - Apply a kernel-filter to an mtz-file
* [`qdep`](#qdep) - Find power law for intensity decay around Bragg positions in mtz-file
* [`pdbrad`](#pdbrad) - Estimate the size of a pdb-object
* [`pdbdist`](#pdbdist) - Plot C-alpha covariance matrix for multistate pdb
* [`covbydist`](#covbydist) - Plot C-alpha covariances by distance and fit decay function
* [`ensemble2adp`](#ensemble2adp) - Convert multistate pdb to (anisotropic) ADPs
* [`subadp`](#subadp) - Subtract ADPs from two pdb-files
* [`btrace`](#btrace) - Plot C-alpha B-factor trace

## Most important command-line parameters
### schimpy
* `pdb_in` - refined structure file (pdb)
* `tls_in` - result of the [ECHT](https://pandda.bitbucket.io/pandemic/echt.html) B-factor distribution (json)
  - if not provided, TLS-matrices are extracted from `pdb_in` (if available)
* `sc_size` - size of the supercell (e.g. `"5 5 10"`)
* `correlate` - enable or disable correlation of TLS-groups (default `True`)
* `use_pbc` - enable or disable periodic boundary conditions (default `True`)
* `stretch` - stretch parameter for the anchor points (default `0.25`)
* `cutoff` - distance cutoff for intermolecular interactions (default `3.0`)
* `weights` - power for the number of interactions (default `1.5`)
* `max_level` - highest level of the TLS-hierarchy to include (default `1`)
* `resolution` - resolution limit (default `2.0`)
* `k_sol` and `b_sol` - bulk solvent parameters (default `0.35` and `50.0`)
* `tls_multipliers` - multipliers for all input tls matrices (e.g. `"1.0 0.0 0.0"`)
* `skip` - skip correlation (but not displacements!) for these levels (e.g. `"2 3 4"`)
* `randomize_after` - randomly swap positions of all molecules after this level (e.g. `2`)
* `reverse` - start procedure at highest level of the hierarchy (default `True`)
* `remove_waters` - remove all water molecules from the input (default `True`)
* `swap_frac` - end correlation prematurely after this fraction of swaps (default `1.0`)
* `energy_percentile` - only allow swaps for which local energy exceeds this percentile (default `0.0`)
* `single_mtz` - write MTZ file with phases after first supercell (default `False`)
* `processes` - number of parallel simulations (default `1`, watch memory usage!)
* `interval` - number of seconds between consecutive simulations (default `1.0`)
* `n_models` - number of supercells to simulate (default `128`)

### stimpy
* `image` - raw image file
* `radial` - polarization-corrected radial average
* `bin_counts` - bin regions within this range of counts (default `1.0`)
* `N` - expected number of independent rotations (default `1.0`)
* `mode` - Noisy Wilson distribution treatment (`continuous`, `discrete` or `dual`, default `continuous`)
* `median_size` - kernel size of the median filter (default `9`)
* `dilation_size` - kernel size of the mask dilation (default `5`)

### slicemtz
* `mtz` - input mtz or map-file
* `lbl` - array of interest (default `IDFF`)
* `slice` - desired slice to plot (default `hk0`)
* `save` - save image as high-res PNG instead of plotting (default `False`)
* `depth` - additional depth of the slice on both sides (default `0`)
* `sc_size` - size of the supercell for drawing Bragg positions (e.g. `"5 5 10"`)
* `overlay` - colour of axes and Bragg position indicators (default `black`)
* `log` - plot log10 of the intensities (default `False`)
* `min` and `max` - minimum and maximum values to plot
* `autoscale` - decrease upper bound for more detail in skewed slices (default `True`)
* `zoom` - zoom in around a given coordinate (`"[x] [y] [pad]"`, e.g. `"65 65 27"`)
* `inset` - inset zoom around a given coordinate (`"[x] [y] [pad]"`, e.g. `"65 65 27"`)
* `contours` - plot this many contours instead of raw values (e.g. `127`)
* `projection` - construct a cylindrical projection (`gp`, `eq` or `mercator`)
* `center` - shift grid to put the corner in the center (default `False`)

### submtz
* `mtz_1` - first input mtz-file
* `mtz_2` - second input mtz-file (can be multiple)
* `lbl_1` - first array of interest (default `IDFF`)
* `lbl_2` - second array of interest (default `IDFF`)
* `mode` - use a different operator (`sub`, `add`, `mul`, `div` or `merge`, default `sub`)
* `scale` - scale factor for second mtz-file (`0` for autoscale, default `1.0`)

### plotmtz
* `mtz` - input mtz-file (can be multiple)
* `lbl` - array of interest (default `IDFF`)
* `log` - plot logarithmic distributions (default `True`)
* `byres` - plot average intensity in resolution bins (default `False`)
* `resolution` - low and high resolution (e.g. `"3.6 3.4"`)

### ccmtz
* `mtz_1` - first input mtz-file
* `mtz_2` - second input mtz-file
* `lbl_1` - first array of interest (default `IDFF`)
* `lbl_2` - second array of interest (default `IDFF`)
* `bins` - number of resolution shells (default `10`)
* `hlim` and `klim` and `llim` - limit h, k, and l (e.g. `"-20 20"`)
* `resolution` - low and high resolution (e.g. `"3.6 3.4"`)

### rmbragg
* `mtz` - input mtz-file
* `sc_size` - supercell size (e.g. `"5 5 10"`)
* `box` - number of voxels to remove around every Bragg position (default `"1 1 1"`)
* `fraction` - remove this fraction of highest intensities in every box (e.g. `0.05`)
* `subtract` - subtract common intensities in resolution shells (`mean` or `min`)
* `keep` - invert selection (default `False`)

### mtzstats
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `sg` - space group for R-int (symbol or number, e.g. `P43212` or `96`)
* `resolution` - low and high resolution (e.g. `"3.6 3.4"`)

### stimpy3d
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `bins` - number of resolution shells (default `1`)
* `N` - expected number of independent rotations (default `1.0`)

### patterson
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `bins` - divide into this many shells and subtract the mean from each one (e.g. `50`)
* `sample` - sampling of the grid in 1/Angstrom, higher is finer (default `3.0`)
* `use_intensities` - calculate Patterson using intensities instead of F (default `False`)
* `center` - place the origin in the center of the map (default `True`)
* `limit` - real space limit in Angstrom of the output map (e.g. `10.0`)
* `resolution` - low and high resolution (e.g. `"3.6 3.4"`)

### pattsize
* `map` - input patterson map
* `binsize` - size of radial bins in Angstrom (default `1.0`)
* `sigma` - sigma cutoff to determine size (default `2.0`)

### mtz2txt
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `resolution` - low and high resolution (e.g. `"3.6 3.4"`)

### mtz2map
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `fill` - fill value for missing reflections (default `-1000`)
* `resolution` - low and high resolution (e.g. `"3.6 3.4"`)

### map2mtz
* `map` - input map-file
* `resolution` - resolution limit (e.g. `2.0`)

### filter_mtz
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `size` - filter size (default `1`)
* `filter` - filter type (`gaussian` or `uniform`, default `gaussian`)
* `interpolate` - interpolate missing intensities using a normalized convolution (default `False`)

### qdep
* `mtz` - input mtz-file
* `lbl` - array of interest (default `IDFF`)
* `sc_size` - supercell size (e.g. `"5 5 10"`)
* `strong` - consider only this number of strongest reflections (default `100`)

### pdbrad
* `pdb` - input multistate pdb-file

### pdbdist
* `pdb` - input multistate pdb-file (can be multiple)
* `mode` - ensemble treatment (`cov`, `cc`, `std`, `var` or `mean`, default `cov`)
* `combine` - multiple input treatment (`both`, `sub`, `div`, `add` or `mul`, default `both`)
* `lines` - plot lines at these x and y-positions (e.g. `"25.5 75.5"`)

### covbydist
* `pdb` - input multistate pdb-file
* `include_neighbours` - include neighbouring molecules in the analysis (default `True`)

### ensemble2adp
* `pdb` - input multistate pdb-file
* `models` - limit number of models from the input pdb (e.g. `100`)

### subadp
* `pdb_1` - first input pdb-file
* `pdb_2` - second input pdb-file (can be multiple)
* `mode` - use a different operator (`sub` or `add`, default `sub`)

### btrace
* `input` - input pdb or json (can be multiple)
* `lines` - plot vertical lines at these x-positions (e.g. `"25.5 75.5"`)


# References
This software relies on methods described in the following papers:
1. Chapman, Henry N., et al. "[Continuous diffraction of molecules and disordered molecular crystals.](https://doi.org/10.1107/S160057671700749X)" Journal of applied crystallography 50.4 (2017): 1084-1103.
2. Pearce, Nicholas M., and Piet Gros. "[A method for intuitively extracting macromolecular dynamics from structural disorder.](https://doi.org/10.1038/s41467-021-25814-x)" Nature communications 12.1 (2021): 5493.
3. Urzhumtsev, Alexandre, et al. "[From deep TLS validation to ensembles of atomic models built from elemental motions. Addenda and corrigendum.](https://doi.org/10.1107/S2059798316013048)" Acta Crystallographica Section D: Structural Biology 72.9 (2016): 1073-1075.
