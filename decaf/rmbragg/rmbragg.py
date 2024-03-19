from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().rmbragg
  print('Reading', p.input.mtz)
  obj         = mtz.object(p.input.mtz)

  indices     = obj.extract_miller_indices().as_vec3_double().as_numpy_array()
  mask        = np.ones(indices.shape[0]).astype(bool)

  # Remove Bragg
  if p.params.sc_size is not None:
    ind_reduced = indices % p.params.sc_size
    lower       = [i//2 for i in p.params.box]
    upper       = [i-j for i,j in zip(p.params.sc_size, lower)]
    mask[((ind_reduced >= upper)|(ind_reduced <= lower)).all(axis=1)] = False

  # Select within limits
  h,k,l  = indices.T
  hlim   = list(map(float, p.params.hlim.split()))
  klim   = list(map(float, p.params.klim.split()))
  llim   = list(map(float, p.params.llim.split()))
  hpos   = ( h >= min(hlim)) & ( h <= max(hlim))
  kpos   = ( k >= min(klim)) & ( k <= max(klim))
  lpos   = ( l >= min(llim)) & ( l <= max(llim))
  hneg   = (-h >= min(hlim)) & (-h <= max(hlim))
  kneg   = (-k >= min(klim)) & (-k <= max(klim))
  lneg   = (-l >= min(llim)) & (-l <= max(llim))
  mask  &= (hpos & kpos & lpos) | (hneg & kneg & lneg)

  # Select intensity fraction
  data   = obj.get_column(p.input.lbl).extract_values().as_numpy_array()
  perc   = abs(p.params.fraction) * 100
  if p.params.fraction > 0:
    mask&= data < np.percentile(data[mask], perc)
  else:
    mask&= data > np.percentile(data[mask], 100 - perc)

  print('Removing voxels')
  if p.params.keep: mask = ~mask
  remove = flex.size_t(np.where(~mask)[0])
  obj.delete_reflections(remove)

  if p.output.mtz_out:
    label = p.output.mtz_out.replace('.mtz','')
  else:
    label = '{}_reduced'.format(p.input.mtz.replace('.mtz',''))

  print('Writing new MTZ')
  obj.write(label+'.mtz')
