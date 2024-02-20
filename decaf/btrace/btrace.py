from __future__ import print_function, division
from iotbx import pdb
from . import phil
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import json
import numpy as np

class ColorCycler():

  def __init__(self, cmap='tab10', start=0):

    self.cmap = plt.get_cmap(cmap)
    self.n    = start % len(self.cmap.colors)

  def __call__(self):

    color  = self.cmap(self.n)
    self.n = (self.n + 1) % len(self.cmap.colors)
    return color


def extract_b(file, chain=0):

  if file.endswith('.json'):
    with open(file, 'r') as inp:
      data  = json.load(inp)
      x     = [float(i) for i in data['x']]
      y     = [float(i) for i in data['y']]
      split = [0]+[n+1 for n,(i,j) in enumerate(zip(x,x[1:])) if j<i]+[len(x)]
      x     = x[split[chain]:split[chain+1]]
      y     = y[split[chain]:split[chain+1]]
      label = data['label']
    return x, y, label

  elif file.endswith('.pdb'):
    hier = pdb.input(file).construct_hierarchy()
    sel  = hier.atom_selection_cache().selection('name CA')
    calp = hier.select(sel).models()[0].chains()[chain].atoms()
    x, y = zip(*[(a.parent().parent().resseq_as_int(), a.b) for a in calp])
    x, y = zip(*[(n, y[i]) for i, n in enumerate(x) if x.index(n)==i])
    return x, y, 'PDB'


def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().btrace

  counter = 0
  colors  = ColorCycler(p.params.cmap, 0)
  plt.figure(figsize=(6,3))

  nlines  = max(map(len, (p.input.input, p.params.chain)))
  chains  = (p.params.chain * len(p.input.input))[:nlines]
  files   = (p.input.input * len(p.params.chain))[:nlines]

  # Plot all traces
  for chain, arg in zip(chains, files):
    if arg == 'skip': colors(); continue
    if '*' in arg:
      arg, mul = arg.split('*')
      mul      = float(mul)
    else:
      mul      = 1.

    x, y, label = extract_b(arg, chain)
    if 'ensemble' in label.lower():
      label = 'Simulated'
      ls    = (0,(6,2))
      color = colors()
      zord  = counter
    else:
      ls    = 'solid'
      color = colors()
      zord  = -counter
    plt.plot(x, np.array(y)*mul, label=label, ls=ls, color=color, zorder=zord)
    counter += 1

  # Plot vertical lines
  lo, hi = plt.ylim()
  plt.vlines(p.params.lines, lo, hi, linewidth=0.5, linestyles='dotted')

  plt.xlabel('Residue number')
  plt.ylabel('C$_{\\alpha}$ B-factor ($\mathrm{\AA}^{2}$)')
  plt.ylim((p.params.min, p.params.max))
  plt.tight_layout()
  if p.params.legend:
    plt.legend()
  plt.savefig(p.output.png_out, dpi=300)
  if p.params.show:
    plt.show()

