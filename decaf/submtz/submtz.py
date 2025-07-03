from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
import numpy as np
from . import phil

def scale_arrays(a, b):

  rms_a = flex.sum(a.data() ** 2.) ** 0.5
  rms_b = flex.sum(b.data() ** 2.) ** 0.5
  return rms_a / rms_b

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().submtz
  print('Reading', p.input.mtz_1)
  obj     = mtz.object(p.input.mtz_1)
  first   = obj.crystals()[0].miller_set().array(obj.get_column(
            p.input.lbl_1).extract_values().as_double())
  first  = first.select(first.data() > p.params.cutoff)
  copy   = first.copy()
  join   = '_' + p.params.mode
  name   = p.input.mtz_1.replace('.mtz', join)
  for file in p.input.mtz_2:
    name+= '_' + file.replace('.mtz', '')
    print('Reading', file)
    obj     = mtz.object(file)
    second  = obj.crystals()[0].miller_set().array(obj.get_column(
              p.input.lbl_2).extract_values().as_double())
    second  = second.select(second.data() > p.params.cutoff)

    if p.params.scale:
      scale = p.params.scale
    else:
      if p.params.scale_common:
        scale = scale_arrays(*copy.common_sets(second,
                                               assert_is_similar_symmetry=False))
      else:
        scale = scale_arrays(copy, second)
      print('Using scale factor', scale)
    second._data *= scale

    if p.params.mode == 'merge':
      first = first.concatenate(second, assert_is_similar_symmetry=False)
    else:
      first, second = first.common_sets(second, assert_is_similar_symmetry=False)
      if p.params.mode == 'sub':
        first._data -= second.data()
      elif p.params.mode == 'add':
        first._data += second.data()
      elif p.params.mode == 'mul':
        first._data *= second.data()
      elif p.params.mode == 'div':
        first._data /= second.data()

  first = first.merge_equivalents().array()

  print('Writing new MTZ')
  if p.output.mtz_out is not None:
    name = p.output.mtz_out
  else:
    name = name.replace('/','_')+'.mtz'
  mtzobj = first.as_mtz_dataset(column_root_label=p.input.lbl_1, column_types='J')
  mtzobj.mtz_object().write(name)
