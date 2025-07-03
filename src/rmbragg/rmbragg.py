from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
import numpy as np
import phil

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  p             = phil.phil_parse(args = args).rmbragg
  print('Reading', p.input.mtz_1)
  keep          = p.params.keep
  sc_size       = list(map(int, p.params.sc_size))
  lower         = [int(i)//2 for i in p.params.box]
  upper         = [i-j for i,j in zip(sc_size, lower)]
  mtz_object    = mtz.object(p.input.mtz_1)

  print('Removing voxels')
  indices       = mtz_object.extract_miller_indices()
  ind_reduced   = np.array(indices) % sc_size
  mask = ((ind_reduced >= upper)|(ind_reduced <= lower)).all(axis=1)
  if keep: mask = ~mask
  remove = flex.size_t(np.where(mask)[0])
  mtz_object.delete_reflections(remove)

  if p.output.mtz_out:
    label = p.output.mtz_out.replace('.mtz','')
  else:
    label = '{}_reduced'.format(p.input.mtz_1.replace('.mtz',''))

  print('Writing new MTZ')
  mtz_object.write(label+'.mtz')
