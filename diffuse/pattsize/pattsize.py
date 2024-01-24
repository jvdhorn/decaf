from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx.ccp4_map import map_reader
from . import phil
import numpy as np
import boost_adaptbx.boost.python as bp
import matplotlib.pyplot as plt

cctbx_maptbx_ext = bp.import_ext('cctbx_maptbx_ext')

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().pattsize
  print('Reading', p.input.map)
  inp    = map_reader(p.input.map)
  origin = flex.vec3_double([[0,0,0]])
  cell   = inp.unit_cell()
  data   = inp.map_data().as_numpy_array()
  lim    = data.std() * p.input.sigma
  print('RMS cutoff:', lim)
  radius = min(cell.parameters()[:3]) / 2
  xy     = []
  r      = p.input.binsize
  value  = float('inf')
  pmask  = None
  print('Calculating radial averages')
  while r <= radius and (p.input.full or value > lim):
    mask = cctbx_maptbx_ext.mask(sites_frac = origin,
                                 unit_cell  = cell,
                                 n_real     = data.shape,
                                 mask_value_inside_molecule  = 1,
                                 mask_value_outside_molecule = 0,
                                 radii      = flex.double(1, r),
                                 wrapping   = True,
    ).as_numpy_array().astype(bool)
    shell = data[mask if pmask is None else mask ^ pmask]
    pmask = mask
    value = abs(shell).mean()
    xy.append((r, value))
    r    += p.input.binsize

  x, y = map(np.array, zip(*xy))
  i, j = next((i-1, i) for i, j in enumerate(y) if j <= lim)
  size = (y[i] - lim) / (y[i] - y[j]) * p.input.binsize + x[i]

  print('Estimated object size:', size)

  if p.input.full:
    plt.title('Patterson map intensity distribution')
    plt.plot(x, y, label='Mean absolute intensity in shell')
    plt.axhline(lim, color='C6', label='Intensity cutoff ({:.1f}$\sigma$)'.format(p.input.sigma))
    plt.vlines(size, 0, lim, linestyle='dashed')
    plt.yscale(p.input.scale)
    if p.input.scale == 'linear':
      plt.ylim(0, lim * 3)
    plt.xlim(0, max(x))
    plt.xlabel('Distance from origin ($\AA$)')
    plt.ylabel('Intensity (a.u.)')
    plt.legend()
    plt.savefig('pattsize.png', dpi=300)
  
    if p.input.show:
      plt.show()

