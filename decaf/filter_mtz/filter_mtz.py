from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
from scipy.ndimage import gaussian_filter, uniform_filter
from scipy.interpolate import griddata
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().filter_mtz
  print('Reading', p.input.mtz)
  obj     = mtz.object(p.input.mtz)
  arr     = obj.crystals()[0].miller_set(False).array(obj.get_column(p.input.lbl
            ).extract_values()).expand_to_p1().complete_array(new_data_value=np.nan)
  size    = p.params.size if p.params.size[1:] else p.params.size[0]
  print('Applying {} filter with size {}'.format(p.params.filter, size))
  data    = arr.data().as_numpy_array()
  valid   = ~np.isnan(data)
  vdata   = data[valid]
  ind     = arr.indices().as_vec3_double().as_numpy_array().astype(int)
  vind    = ind[valid]
  offset  = abs(ind).max(axis=0)
  shape   = 2 * offset + 1
  # stackoverflow.com/a/36307291 for filtering approach
  filt    = {'gaussian':gaussian_filter, 'uniform':uniform_filter}[p.params.filter]
  agrid   = np.zeros(shape=shape)
  agrid[tuple(( vind + offset).T)] = vdata
  agrid[tuple((-vind + offset).T)] = vdata
  bgrid   = np.zeros(shape=shape)
  bgrid[tuple(( vind + offset).T)] = 1
  bgrid[tuple((-vind + offset).T)] = 1
  with np.errstate(all='ignore'):
    data = (filt(agrid, size, mode='constant')
           /filt(bgrid, size, mode='constant')
           )[tuple((ind + offset).T)]
  arr._data = flex.double(data)

  if p.params.interpolate:
    arr   = arr.select(flex.bool(~np.isnan(data)))
  else:
    arr   = arr.select(flex.bool(~np.isnan(data) & valid))

  print('Writing new MTZ')
  mtzobj  = arr.as_mtz_dataset(column_root_label=p.input.lbl, column_types='J')
  if p.output.mtz_out is None:
    name = p.input.mtz.replace('.mtz','_filt.mtz')
  else:
    name = p.output.mtz_out
  mtzobj.mtz_object().write(name)
