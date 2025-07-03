from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
from matplotlib import pyplot as plt, ticker
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().plotmtz

  # Plot distribution
  plt.close()
  plt.figure(figsize=(16/3,4))
  scf = ticker.ScalarFormatter()
  scf.set_powerlimits((-2,2))
  scf.set_useMathText(True)
  plt.gca().xaxis.set_major_formatter(scf)
  low = high = None
  lores = max(p.params.resolution)
  hires = min(p.params.resolution)
  if lores == hires: lores = 9e99
  files = list(p.input.mtz)
  for file in files:
    print('Reading', file)
    obj    = mtz.object(file)
    first  = obj.crystals()[0].miller_set().array(obj.get_column(
             p.input.lbl).extract_values().as_double())
    first  = first.select(first.data() > p.params.cutoff)
    first  = first.resolution_filter(lores, hires)
    data   = first.data().as_numpy_array()
    if low == high == None:
      low  = data.min()
      high = data.max()
    label = file.replace('/','_').replace('.mtz','')
    plt.hist(data, bins=p.input.bins, range=(low, high), histtype='stepfilled',
             alpha=0.5, log=p.params.log, label=label, ec='black')
  plt.xlabel('Intensity')
  plt.ylabel('Counts')
  if p.params.legend:
    plt.legend()
    plt.title('Resolution range {:.2f} - {:.2f}'.format(*first.resolution_range()))
  plt.tight_layout()
  file = files[0].replace('/','_')
  name = 'plotmtz_{}_distribution.png'.format(file.replace('.mtz',''))
  if not p.params.log: name.replace('.png','_lin.png')
  name = p.output.png_out or name
  plt.savefig(name, dpi=300)

  if p.input.show:
    plt.show()
