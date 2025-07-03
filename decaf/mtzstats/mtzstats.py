from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
from cctbx import crystal
from . import phil
from scipy import stats, ndimage
import numpy as np

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().mtzstats
  print('Reading', p.input.mtz)
  obj    = mtz.object(p.input.mtz)
  lores  = max(p.input.resolution)
  hires  = min(p.input.resolution)
  if lores == hires: lores = float('inf')
  first  = obj.crystals()[0].miller_set(False).array(obj.get_column(p.input.lbl
           ).extract_values().as_double()).resolution_filter(lores, hires)
  npdata = first.data().select(first.data() > p.input.cutoff).as_numpy_array()

  # Determine number of local maxima
  ind    = first.indices().as_vec3_double().as_numpy_array().astype(int)
  offset = abs(ind).max(axis=0) + 1
  grid   = np.full(offset * 2 + 1, npdata.min()-1)
  grid[tuple((offset+ind).T)] = npdata
  grid[tuple((offset-ind).T)] = npdata
  ksize  = p.params.size if len(p.params.size) == 3 else p.params.size[0]
  fgrid  = ndimage.maximum_filter(grid, ksize, mode='constant', cval=grid.min())
  maxima = (grid == fgrid)[tuple((offset+ind).T)].sum()

  # Determine other stats 
  norm   = (npdata - npdata.min()) / (npdata.max() - npdata.min())
  rmsc   = norm.std()
  spcont = npdata.var() / npdata.mean() ** 2
  distr  = np.histogram(npdata, bins=p.params.bins)[0] / npdata.size
  distr  = distr[distr > 0]
  entropy= -sum(distr * np.log(distr))

  # Output
  strf = '%f' if p.params.floats else '%e'
  print('Space group     :', first.space_group_info().symbol_and_number())
  print('Resolution range:', *map('{:.5f}'.format,first.resolution_range()))
  print('Min max indices :', *first.min_max_indices())
  print('Unit cell param.:', *map('{:.5f}'.format,first.unit_cell().parameters()))
  print('Reciprocal cell :', *map('{:.5f}'.format,first.unit_cell().reciprocal_parameters()))
  print('Number of refl. :', npdata.size)
  print('Fraction neg.   :', (npdata < 0).sum() / npdata.size)
  print('No. local maxima:', maxima)
  print('Minimum         :', strf%npdata.min())
  print('Maximum         :', strf%npdata.max())
  print('Mean            :', strf%npdata.mean())
  print('Median          :', strf%np.median(npdata))
  print('Std. deviation  :', strf%npdata.std())
  print('Variance        :', strf%npdata.var())
  print('Skewness        :', float(stats.skew(npdata)))
  print('RMS contrast    :', rmsc)
  print('Var[I]/Mean[I]^2:', spcont)
  print('Shannon entropy :', entropy)

  if p.params.sg is not None:
    sg     = p.params.sg or str(first.space_group_info())
    symm   = crystal.symmetry(unit_cell=first.unit_cell(), space_group_symbol=sg)
    merged = first.customized_copy(crystal_symmetry = symm).merge_equivalents()
    symbol = ''.join(str(symm.space_group_info()).split())
    print('R_int {:10s}:'.format(symbol), merged.r_merge())
