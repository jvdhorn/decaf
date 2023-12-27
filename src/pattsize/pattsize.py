from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx.map_manager import map_manager
from scipy.optimize import curve_fit, minimize_scalar
import phil
import numpy as np
import boost_adaptbx.boost.python as bp
import matplotlib.pyplot as plt

cctbx_maptbx_ext = bp.import_ext('cctbx_maptbx_ext')

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  p      = phil.phil_parse(args = args).pattsize
  print('Reading', p.input.map_1)
  mmm    = map_manager(p.input.map_1).as_map_model_manager()
  origin = flex.vec3_double([[0,0,0]])
  cell   = mmm.crystal_symmetry().unit_cell()
  data   = mmm.map_data().as_numpy_array()
  radius = min(cell.parameters()[:3]) / 2
  xy     = []
  pmask  = None
  print('Calculating radial distributions')
  for r in np.arange(0,radius,p.input.binsize) + p.input.binsize:
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

  x, y = map(np.array, zip(*xy))
  lim  = data.std() * p.input.sigma
  print('RMS cutoff:', lim)
  i, j = next((i-1, i) for i, j in enumerate(y) if j <= lim)
  size = (y[i] - lim) / (y[i] - y[j]) * p.input.binsize + x[i]

  func = None
  '''
  func = lambda x, a, b, c, d, e: a * np.exp(-x/b) + c * np.exp(-x/d) + e
  popt, pcov = curve_fit(func, x, np.log(y))
  size = minimize_scalar(lambda x: (func(x, *popt) - np.log(lim)) ** 2).x
  print('Fit:', *('{}={:.2f};'.format(n,v)
                  for n, v in zip(func.__code__.co_varnames[1:], popt)))
  '''

  print('Estimated object size:', size)

  plt.title('Patterson map intensity distribution')
  plt.plot(x, y, label='Mean absolute intensity in shell')
  if func: plt.plot(x, np.exp(func(x, *popt)), label='$f(x)=ae^{-x/b}+be^{-x/d}$ fit')
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

