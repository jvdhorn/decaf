from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
from cctbx import crystal
from . import phil
import numpy as np

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().mtzstats
  print('Reading', p.input.mtz)
  obj    = mtz.object(p.input.mtz)
  first  = obj.crystals()[0].miller_set(False).array(obj.get_column(
           p.input.lbl).extract_values().as_double())
  npdata = first.data().select(first.data() > p.input.cutoff).as_numpy_array()
  rms    = (npdata - npdata.mean()).std()
  norm   = (npdata - npdata.min()) / (npdata.max() - npdata.min())
  nrms   = (norm - norm.mean()).std()
  spcont = npdata.var() / npdata.mean() ** 2
  distr  = np.histogram(npdata, bins=p.params.bins)[0] / npdata.size
  distr  = distr[distr > 0]
  entropy= -sum(distr * np.log(distr))
  print('Number of refl. :', npdata.size)
  print('Min             :', npdata.min())
  print('Max             :', npdata.max())
  print('Mean            :', npdata.mean())
  print('Var             :', npdata.var())
  print('RMS contrast    :', rms)
  print('   >  normalized:', nrms)
  print('Var[I]/Mean[I]^2:', spcont)
  print('Shannon entropy :', entropy, '({} bins)'.format(p.params.bins))

  sg     = p.params.sg or str(first.space_group_info())
  symm   = crystal.symmetry(unit_cell=first.unit_cell(), space_group_symbol=sg)
  merged = first.customized_copy(crystal_symmetry = symm).merge_equivalents()
  symbol = symm.space_group_info().symbol_and_number()
  print('R_merge         :', merged.r_merge(), '(for {})'.format(symbol))
