from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
import numpy as np
from . import phil

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().mtz2txt
  print('Reading', p.input.mtz)
  obj     = mtz.object(p.input.mtz)
  lores   = max(p.input.resolution)
  hires   = min(p.input.resolution)
  if lores == hires: lores = 9e99
  arr     = obj.crystals()[0].miller_set(False).array(obj.get_column(
            p.input.lbl).extract_values()).resolution_filter(lores, hires)
  data    = arr.data().as_numpy_array()

  if p.output.txt_out:
    label = p.output.txt_out.replace('.txt', '')
  else:
    label = p.input.mtz.replace('.mtz', '')

  print('Writing TXT')
  np.savetxt(label+'.txt', data.ravel())
