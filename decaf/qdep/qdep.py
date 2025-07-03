from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
from matplotlib import pyplot as plt
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().qdep
  if not p.params.sc_size: print('Please provide sc_size'); return
  print('Reading', p.input.mtz)
  sc_size       = p.params.sc_size
  obj           = mtz.object(p.input.mtz)
  labels        = obj.column_labels()
  label         = p.input.lbl if p.input.lbl in labels else labels[3]
  hi, lo        = (sorted(p.params.resolution) + [float('inf')])[:2]
  arr           = obj.crystals()[0].miller_set(False).array(obj.get_column(label
                  ).extract_values()).expand_to_p1().resolution_filter(lo, hi)
  data          = arr.data().as_numpy_array()
  ind           = arr.indices().as_vec3_double().as_numpy_array().astype(int)
  origin        = abs(ind).max(axis=0)
  shape         = 2 * origin + 1
  grid          = np.zeros(shape=shape, dtype=float)
  grid[tuple((-ind + origin).T)] = data
  grid[tuple(( ind + origin).T)] = data
  unit_cell     = arr.crystal_symmetry().unit_cell().parameters()[:3]
  indices       = np.stack(np.meshgrid(*map(np.arange, grid.shape))).reshape(3,-1).T
  bragg         = indices[((indices - origin) % sc_size == 0).all(axis=1)]
  ranges        = [dim if dim>0 else (sc-1)//2 for sc, dim in zip(sc_size,p.params.range)]
  premask       = ~((bragg<ranges) ^ (-bragg+origin+origin<ranges)).any(axis=1)
  lower         = tuple(map(int, p.params.min or (bragg - origin).min(axis=0)))
  upper         = tuple(map(int, p.params.max or (bragg - origin).max(axis=0)))
  lower, upper  = zip(*map(sorted, zip(lower, upper)))
  premask      &= (lower<=bragg-origin).all(axis=1) & (bragg-origin<=upper).all(axis=1)
  realbragg     = bragg[premask]
  L, D          = realbragg.shape
  data_sub      = []
  for n, dim in enumerate(ranges):
    extent            = list(range(-1,~dim,-1)) + list(range(1,dim+1))
    selection         = realbragg[:,None].repeat(len(extent),axis=1)
    selection[:,:,n] += extent
    data_sub.append(grid[tuple(selection.reshape(-1,D).T)].reshape(L, -1))
  i_lim         = p.params.strong
  mesh          = np.mgrid[tuple(slice(-r,r+1) for r in ranges)].reshape(3,-1).T
  select        = grid[
                    tuple((realbragg[:,None].repeat(len(mesh),axis=1) + mesh).T)
                  ].sum(axis=0).argsort()
  select        = select[np.concatenate(data_sub, axis=1)[select].min(axis=1)>0.]
  select        = select[-i_lim:] if i_lim >=0 else select[:-i_lim]
  for n,arr in enumerate(data_sub):
    fold        = sum(np.split(arr[select], 2, axis=1)) / 2.
    q_log       = np.log(np.arange(1, arr.shape[1]//2 + 1) * unit_cell[n] * 2 * np.pi)
    fit         = np.polyfit(q_log, np.log(fold.T), 1)[0]
    bins        = int(np.ceil(np.log2(fit.shape[0]))+1)
    if p.params.plot:
      plt.hist(fit, bins=bins, histtype='stepfilled', alpha=0.3, label='abc'[n]+'*')
      plt.hist(fit, bins=bins, histtype='step', color='black', lw=.5)
    print('{}*: {} pts, range {}, mean exponent {:.3f} +/- {:.3f}'.format(
          'abc'[n], fit.shape[0], arr.shape[1]//2, fit.mean(), fit.std()))
  if p.params.plot:
    plt.legend();
    plt.title('Intensity decay around Bragg positions')
    plt.xlabel('Best fit exponent')
    plt.ylabel('Occurence')
    plt.show(block=True)
