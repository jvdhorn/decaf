from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
import phil
import numpy as np

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  p      = phil.phil_parse(args = args).mtzstats
  print('Reading', p.input.mtz_1)
  obj    = mtz.object(p.input.mtz_1)
  first  = obj.crystals()[0].miller_set(False).array(obj.get_column(
           p.input.lbl_1).extract_values()).expand_to_p1()
  npdata = first.data().select(first.data() > p.input.cutoff).as_numpy_array()
  rms    = (npdata - npdata.mean()).std()
  norm   = (npdata - npdata.min()) / (npdata.max() - npdata.min())
  nrms   = (norm - norm.mean()).std()
  spcont = npdata.var() / npdata.mean() ** 2
  histo  = np.histogram(npdata, bins=1000)[0] / npdata.size
  histo  = histo[histo > 0]
  entropy= -sum(histo * np.log(histo))
  print('Min             :', npdata.min())
  print('Max             :', npdata.max())
  print('Mean            :', npdata.mean())
  print('Var             :', npdata.var())
  print('RMS contrast    :', rms)
  print('   >  normalized:', nrms)
  print('Var[I]/Mean[I]^2:', spcont)
  print('Shannon entropy :', entropy)
