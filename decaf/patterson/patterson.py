from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().patterson
  print('Reading', p.input.mtz)
  obj    = mtz.object(p.input.mtz)
  lores  = max(p.input.resolution)
  hires  = min(p.input.resolution)
  if lores == hires: lores = float('inf')
  first  = obj.crystals()[0].miller_set(False).array(obj.get_column(p.input.lbl
           ).extract_values().as_double()).resolution_filter(lores, hires)
  if p.input.fill is not None and first.completeness() < 1:
    first = first.complete_array(new_data_value=p.input.fill)

  if p.params.bins > 0:
    print('Subtracting resolution shell means')
    data   = first.data().as_numpy_array()
    binner = first.setup_binner(n_bins=p.params.bins)
    ind    = binner.bin_indices().as_numpy_array()
    for i in binner.range_used():
      sel   = data[ind == i]
      data[ind == i] -= sel.mean()
    first._data = flex.double(data)

  print('Calculating Patterson')
  patt   = first.f_sq_as_f().patterson_map(resolution_factor = 1./p.params.sample)

  if p.params.center:
    patt._real_map = pdata = patt._real_map.as_numpy_array()
    if not pdata[-1,:,:].any(): pdata = pdata[:-1,:,:]
    if not pdata[:,-1,:].any(): pdata = pdata[:,:-1,:]
    if not pdata[:,:,-1].any(): pdata = pdata[:,:,:-1]
    offset   = np.array([i//2 for i in pdata.shape])
    pdata    = np.roll(pdata, offset, (0,1,2))
    if pdata.shape[0] % 2 == 0:
      pdata = np.concatenate((pdata,pdata[:1,:,:]), axis=0)
    if pdata.shape[1] % 2 == 0:
      pdata = np.concatenate((pdata,pdata[:,:1,:]), axis=1)
    if pdata.shape[2] % 2 == 0:
      pdata = np.concatenate((pdata,pdata[:,:,:1]), axis=2)
    if p.params.limit is not None:
      lim    = np.array([p.params.limit * p.params.sample]*3).astype(int)
      pdata  = pdata[tuple(map(slice,offset-lim, offset+lim+1))]
      offset = lim
    patt._real_map = flex.double(pdata.ravel())
    new_grid = flex.grid(offset-pdata.shape+1, offset+1)
    patt._real_map.reshape(new_grid)

  print('Writing map')
  if p.output.map_out:
    label = p.output.map_out.replace('.map', '')
  else:
    label = p.input.mtz.replace('.mtz', '_patterson')
  patt.as_map_manager().write_map(label+'.map')
