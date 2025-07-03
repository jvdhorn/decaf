from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().submtz
  print('Reading', p.input.mtz_1)
  obj     = mtz.object(p.input.mtz_1)
  first   = obj.crystals()[0].miller_set().array(obj.get_column(
            p.input.lbl_1).extract_values().as_double())
  first  = first.select(first.data() > p.params.cutoff)
  join   = '_' + p.params.mode
  name   = p.input.mtz_1.replace('.mtz', join)
  fscale = flex.sum(first.data()**2)**.5
  for file in p.input.mtz_2:
    name+= '_' + file.replace('.mtz', '')
    print('Reading', file)
    obj     = mtz.object(file)
    second  = obj.crystals()[0].miller_set().array(obj.get_column(
              p.input.lbl_2).extract_values().as_double())
    second = second.select(second.data() > p.params.cutoff)
    first, second = first.common_sets(second)
    if p.params.scale == 0:
      scale = fscale / flex.sum(second.data() ** 2)**.5
      print('Using scale factor', scale)
    else:
      scale = p.params.scale
    second._data *= scale
    if p.params.mode == 'sub':
      first._data -= second.data()
    elif p.params.mode == 'add':
      first._data += second.data()
    elif p.params.mode == 'mul':
      first._data *= second.data()
    elif p.params.mode == 'div':
      first._data /= second.data()

  print('Writing new MTZ')
  if p.output.mtz_out is not None:
    name = p.output.mtz_out
  else:
    name = name.replace('/','_')+'.mtz'
  mtzobj = first.as_mtz_dataset(column_root_label=p.input.lbl_1, column_types='J')
  mtzobj.mtz_object().write(name)
