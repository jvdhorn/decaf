from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
import numpy as np
from . import phil

def run(args):

  scope  = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p      = scope.extract().rmbragg
  print('Reading', p.input.mtz)
  obj     = mtz.object(p.input.mtz)
  arr     = obj.crystals()[0].miller_set(False).array(obj.get_column(p.input.lbl
            ).extract_values()).expand_to_p1()
  indices = arr.indices().as_vec3_double().as_numpy_array().astype(int)
  mask    = np.ones(indices.shape[0]).astype(bool)
  frac    = p.params.fraction and (1-p.params.fraction) * 100
  data    = arr.data().as_numpy_array()

  print('Removing voxels')
  # Remove intensities around Bragg positions
  if p.params.sc_size is not None:
    half   = np.array(p.params.box) // 2
    box    = half * 2 + 1

    offset = abs(indices).max(axis=0) + half + p.params.sc_size
    shape  = 2 * offset + 1
    grid   = np.full(shape, np.nan)
    grid[tuple((-indices + offset).T)] = data
    grid[tuple(( indices + offset).T)] = data
    shadow = grid.copy()

    # Find all indices
    bragg  = arr.miller_set(arr.indices(), False).complete_set().indices(
             ).as_vec3_double().as_numpy_array().astype(int)
    # Filter for Bragg
    bragg  = bragg[(bragg % p.params.sc_size == 0).all(axis=1)] + offset - half
    # Expand by one supercell size on all sides
    exp    = np.mgrid[-1:2,-1:2,-1:2].reshape(3,-1).T * p.params.sc_size
    bragg  = (np.repeat(bragg,27,axis=0).reshape(-1,27,3) + exp).reshape(-1, 3)
    bragg  = np.unique(bragg, axis=0)

    with np.errstate(all='ignore'):
      for ind in bragg:
        sel = tuple(map(slice,ind,ind+box))
        if frac is None:
          shadow[sel] = np.nan
        else:
          val = grid[sel]
          if not np.isnan(val).all():
            val[val >= np.nanpercentile(val,frac)] = np.nan
            shadow[sel] = val

    mask &= ~(np.isnan(shadow[tuple(( indices + offset).T)])
             |np.isnan(shadow[tuple((-indices + offset).T)]))

  # Remove fraction of whole map
  elif frac is not None:
    with np.errstate(all='ignore'):
      mask &= (data < np.nanpercentile(data,frac))


  # Select within limits
  h,k,l = indices.T
  hlim  = list(map(float, p.params.hlim.split()))
  klim  = list(map(float, p.params.klim.split()))
  llim  = list(map(float, p.params.llim.split()))
  hpos  = ( h >= min(hlim)) & ( h <= max(hlim))
  kpos  = ( k >= min(klim)) & ( k <= max(klim))
  lpos  = ( l >= min(llim)) & ( l <= max(llim))
  hneg  = (-h >= min(hlim)) & (-h <= max(hlim))
  kneg  = (-k >= min(klim)) & (-k <= max(klim))
  lneg  = (-l >= min(llim)) & (-l <= max(llim))
  mask &= (hpos & kpos & lpos) | (hneg & kneg & lneg)

  # Limit resolution
  hi,lo = (sorted(p.params.resolution) + [float('inf')])[:2]
  mask &= arr.resolution_filter_selection(lo, hi).as_numpy_array()

  # Invert selection
  if p.params.keep: mask = ~mask

  # Radial correction
  if p.params.subtract:
    print('Applying radial corrections')
    bins = arr.setup_binner(n_bins=p.params.bins)
    func = {'min': np.min, 'mean': np.mean}[p.params.subtract]
    for n in bins.range_all():
      selection = bins.array_indices(n).as_numpy_array()
      if selection.size: data[selection] -= func(data[selection][mask[selection]])
    obj.get_column(p.input.lbl).set_values(flex.float(data))

  # Apply
  remove = flex.size_t(np.where(~mask)[0])
  obj.delete_reflections(remove)

  if p.output.mtz_out:
    label = p.output.mtz_out.replace('.mtz','')
  else:
    label = '{}_reduced'.format(p.input.mtz.replace('.mtz',''))

  print('Writing new MTZ')
  obj.write(label+'.mtz')
