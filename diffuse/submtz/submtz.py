from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
import numpy as np
from . import phil

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().submtz
  print('Reading', p.input.mtz_1)
  obj     = mtz.object(p.input.mtz_1)
  first   = obj.crystals()[0].miller_set().array(obj.get_column(
            p.input.lbl_1).extract_values().as_double())
  first  = first.select(first.data() > p.params.cutoff)
  print('Reading', p.input.mtz_2)
  obj     = mtz.object(p.input.mtz_2)
  second  = obj.crystals()[0].miller_set().array(obj.get_column(
            p.input.lbl_2).extract_values().as_double())
  second = second.select(second.data() > p.params.cutoff)
  first, second = first.common_sets(second)
  if p.params.scale == 0:
    print('Scaling arrays')
    scale = flex.sum(first.data()**2)**.5 / flex.sum(second.data() ** 2)**.5
  else:
    scale = p.params.scale
  print('Using scale factor', scale)
  second._data *= scale
  if p.params.add:
    print('Adding arrays')
    first._data += second.data()
    join = '_add_'
  else:
    print('Subtracting arrays')
    first._data -= second.data()
    join = '_sub_'
  print('Writing new MTZ')
  mtzobj = first.as_mtz_dataset(column_root_label=p.input.lbl_1, column_types='J')
  name   = (p.input.mtz_1.replace('.mtz', join) + p.input.mtz_2).replace('/','_')
  mtzobj.mtz_object().write(name)
