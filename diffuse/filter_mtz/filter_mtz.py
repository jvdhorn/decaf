from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
from scipy.ndimage import gaussian_filter, uniform_filter
import numpy as np
from . import phil

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().filter_mtz
  print('Reading', p.input.mtz_1)
  obj     = mtz.object(p.input.mtz_1)
  arr     = obj.crystals()[0].miller_set(False).array(obj.get_column(
            p.input.lbl_1).extract_values()).expand_to_p1()
  size    = p.params.size if p.params.size[1:] else p.params.size[0]
  filt    = gaussian_filter if p.params.filter == 'gaussian' else uniform_filter
  print('Applying {} filter with size {}'.format(p.params.filter, size))
  data    = arr.data().as_numpy_array()
  ind     = arr.indices().as_vec3_double().as_numpy_array().astype(int)
  offset  = abs(ind).max(axis=0)
  shape   = 2 * offset + 1
  # stackoverflow.com/a/36307291 for filtering approach
  agrid   = np.zeros(shape=shape)
  agrid[tuple(( ind + offset).T)] = data
  agrid[tuple((-ind + offset).T)] = data
  bgrid   = np.zeros(shape=shape)
  bgrid[tuple(( ind + offset).T)] = 1
  bgrid[tuple((-ind + offset).T)] = 1
  afilt   = filt(agrid, size)[tuple((ind + offset).T)]
  bfilt   = filt(bgrid, size)[tuple((ind + offset).T)]
  arr._data = flex.double(afilt/bfilt)
  print('Writing new MTZ')
  mtzobj  = arr.as_mtz_dataset(column_root_label=p.input.lbl_1, column_types='J')
  if p.output.mtz_out is None:
    name = p.input.mtz_1.replace('.mtz','_filt.mtz')
  else:
    name = p.output.mtz_out
  mtzobj.mtz_object().write(name)
