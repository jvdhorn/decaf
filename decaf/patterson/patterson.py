from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().patterson
  print('Reading', p.input.mtz)
  obj    = mtz.object(p.input.mtz)
  first  = obj.crystals()[0].miller_set(False).array(obj.get_column(
           p.input.lbl).extract_values().as_double())

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
  print('Writing map')
  if p.output.map_out:
    label = p.output.map_out.replace('.map', '')
  else:
    label = p.input.mtz.replace('.mtz', '_patterson')
  patt.as_ccp4_map(label+'.map')
