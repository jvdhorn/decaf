from __future__ import print_function, division
from iotbx import pdb
from cctbx.array_family import flex
from . import phil
import matplotlib.pyplot as plt
import numpy as np

funcs = {'std':np.std, 'var':np.var, 'mean':np.mean, 'add':np.add,
         'sub':np.subtract, 'div':np.divide, 'mul':np.multiply}
modes = {'std':'standard deviation', 'var':'variance', 'mean':'mean',
         'add':'sum', 'sub':'difference', 'div':'ratio', 'mul':'product'}

def extract_grid(file, mode, chain, resrng, modrng):

  print('Reading', file)
  hier = pdb.input(file).construct_hierarchy()
  sel  = hier.atom_selection_cache().selection('name CA')
  calp = hier.select(sel)
  resi = [atom.parent().parent().resseq_as_int() for atom in
          calp.models()[chain].chains()[chain].atoms()]
  rrng = resrng or [min(resi), max(resi)]
  mrng = modrng or [None]
  lo,hi= min(rrng), max(rrng)
  sel  = flex.size_t([resi.index(n) for n in sorted(set(resi)) if lo <= n <= hi])
  resi = [resi[i] for i in sel]
  grid = []

  for model in calp.models()[min(mrng):max(mrng)+1]:
    xyz  = model.chains()[chain].atoms().extract_xyz().select(sel)
    dmat = [(xyz - coord).dot()**0.5 for coord in xyz]
    grid.append(dmat)

  mode = funcs[mode]
  grid = mode(np.array(grid), axis=0)
  grid[np.triu_indices(grid.shape[0],1)] = np.nan

  for n in range(lo, hi):
    if n not in resi:
      grid = np.insert(grid, n-lo, np.nan, axis=0)
      grid = np.insert(grid, n-lo, np.nan, axis=1)

  return grid, lo, hi


def align_grids(grids):

  grids, lows, highs = map(list,zip(*grids))
  lo, hi = max(lows), min(highs)

  for n, (grid, low, high) in enumerate(zip(grids, lows, highs)):
    start    = lo - low or None
    end      = hi - high or None
    grids[n] = grid[start:end, start:end]

  return grids, (lo, hi)


def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().pdbdist

  # Extract and align grids
  grids = [extract_grid(file, p.params.ensemble, p.params.chain,
                        p.params.residue_range, p.params.model_range)
           for file in p.input.pdb[:4]]
  grids, (lo, hi) = align_grids(grids)

  # Combine grids
  if len(grids) == 1:
    grid = grids[0]
  elif p.params.mode == 'combine':
    grid, second = grids[:2]
    second = second.T
    valid  = ~np.isnan(second)
    grid[valid] = second[valid]
  else:
    pairs = zip(*(2*[iter(grids)]))
    for n, (a, b) in enumerate(pairs):
      with np.errstate(divide='ignore', invalid='ignore'):
        merge = funcs[p.params.mode](a, b)
      if n == 0:
        grid = merge
      else:
        merge = merge.T
        valid = ~np.isnan(merge)
        grid[valid] = merge[valid]
 
  # Construct label 
  base = 'C$_{\\alpha}$ distance'
  unit = '($\mathrm{\AA}'+'^{2}'*(p.params.ensemble in {'var'})+'$)'
  elbl = modes[p.params.ensemble]
  if p.params.mode == 'combine':
    label = [base, elbl, unit]
  else:
    mlbl = modes[p.params.mode]
    label = [base, elbl, mlbl, unit]
    if p.params.mode in {'div','mul'}: label.remove(unit)

  # Make plot
  fig, ax = plt.subplots(figsize=(5,4))
  im = ax.imshow(grid, origin='lower', vmin=p.params.min, vmax=p.params.max,
                 extent=(lo-.5,hi+.5,lo-.5,hi+.5), cmap=p.params.cmap)
  ax.xaxis.set_ticks_position('top')
  ax.xaxis.set_label_position('top')
  ax.set_frame_on(False)
  ax.set_xlabel('Residue number')
  ax.set_ylabel('Residue number')
  cb = fig.colorbar(im, label=' '.join(label))
  fig.tight_layout()
  ax.set_yticks([t for t in ax.get_xticks() if lo <= t <= hi])
  plt.savefig(p.output.png_out, dpi=300)

  if p.params.show:
    plt.show()
