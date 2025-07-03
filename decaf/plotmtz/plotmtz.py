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
  files  = list(p.input.mtz)
  print('Reading', files[0])
  obj    = mtz.object(files[0])
  first  = obj.crystals()[0].miller_set().array(obj.get_column(
           p.input.lbl).extract_values().as_double())
  first  = first.select(first.data() > p.params.cutoff)
  binner = first.setup_binner(n_bins=p.input.bins)
  ind    = binner.bin_indices()
  means  = [first.data().select(ind==i).min_max_mean().mean for i in binner.range_used()]
  mean   = first.data().min_max_mean().mean

  x = ['$\infty$'] + ['%.2f'%binner.bin_d_range(i)[1] for i in binner.range_used()]
  y = means
  xticks = [i-0.5 for i in range(len(x))]

  if len(y) >= 10:
    xsel   = [int(i//(len(y)/10)) for i in range(len(y))]+[10]
    xind   = [xsel.index(i) for i in range(11)]
    x      = [x[i] for i in xind]
    xticks = [xticks[i] for i in xind]
   
  plt.bar(range(len(y)), y, label='Binned')
  plt.plot(xticks, [mean] * len(xticks))
  plt.plot(xticks, [mean] * len(xticks), label='All')
  plt.xticks(xticks, x)
  plt.xlabel('Resolution shell limits ($\AA$)')
  plt.ylabel('Mean intensity')
  if p.params.legend: plt.legend()
  plt.tight_layout()
  name ='plotmtz_{}_byres.png'.format(files[0].replace('/','_').replace('.mtz',''))
  plt.savefig(name, dpi=300)

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
  for file in files:
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
  plt.savefig(name, dpi=300)

  if p.input.show:
    plt.show()
