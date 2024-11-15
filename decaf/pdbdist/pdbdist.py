from __future__ import print_function, division
from iotbx import pdb
from cctbx.array_family import flex
from . import phil
from matplotlib import pyplot as plt, colors
import numpy as np

funcs = {'std':np.std, 'var':np.var, 'mean':np.mean, 'add':np.add,
         'sub':np.subtract, 'div':np.divide, 'mul':np.multiply}

def extract_grid(file, mode, chain, resrng, states, ref=None):

  print('Reading', file)
  hier = pdb.input(file).construct_hierarchy(False,False,False)
  sel  = hier.atom_selection_cache().selection('name CA')
  calp = hier.select(sel)
  resi = [atom.parent().parent().resseq_as_int() for atom in
          calp.models()[0].chains()[chain].atoms()]
  rrng = resrng or [min(resi), max(resi)]
  mrng = states or [0, len(calp.models())-1]
  lo,hi= min(rrng), max(rrng)
  sel  = flex.size_t([resi.index(n) for n in sorted(set(resi)) if lo <= n <= hi])
  resi = [resi[i] for i in sel]
  grid = []

  if mode in {'cc', 'cov'}:
    for model in calp.models()[min(mrng):max(mrng)+1]:
      grid.append(model.chains()[chain].atoms().extract_xyz(
                  ).select(sel).as_numpy_array())

    mode = {'cc':np.corrcoef, 'cov':np.cov}[mode]
    grid = mode((grid - np.mean(grid, axis=0)).reshape(len(grid),-1).T)
    grid = np.array((
             grid[0::3,0::3], grid[1::3,1::3], grid[2::3,2::3],
#             grid[0::3,1::3], grid[0::3,2::3], grid[1::3,2::3],
           )).sum(axis=0)

  else:
    for model in calp.models()[min(mrng):max(mrng)+1]:
      xyz  = model.chains()[chain].atoms().extract_xyz().select(sel)
      dmat = [((xyz - coord).dot()**0.5).as_numpy_array() for coord in xyz]
      grid.append(dmat)
 
    grid = np.array(grid)
 
    if ref is not None and ref[1:] == (lo, hi):
      grid -= ref[0]
 
    mode = funcs[mode]
    grid = mode(grid, axis=0)

  grid[np.triu_indices(grid.shape[0],1)] = np.nan

  for n in range(lo, hi):
    if n not in resi:
      grid = np.insert(grid, n-lo, np.nan, axis=0)
      grid = np.insert(grid, n-lo, np.nan, axis=1)

  print('Mean matrix value:', np.nanmean(grid))

  return grid, lo, hi


def align_grids(grids):

  grids, lows, highs = map(list,zip(*grids))
  lo, hi = min(lows), max(highs)

  for n, (grid, low, high) in enumerate(zip(grids, lows, highs)):
    grid     = np.concatenate((
                  np.full((low-lo, grid.shape[1]),np.nan), grid,
                  np.full((hi-high, grid.shape[1]),np.nan)
               ), axis=0)
    grids[n] = np.concatenate((
                  np.full((grid.shape[0],low-lo),np.nan), grid,
                  np.full((grid.shape[0],hi-high),np.nan)
               ), axis=1)

  return grids, (lo, hi)


def write_pymol_bonds(file, grid, lo, hi, norm, perc):

  inp  = pdb.input(file)
  hier = inp.construct_hierarchy(False,False,False)
  sel  = hier.atom_selection_cache().selection('name CA')
  hier = hier.select(sel)
  xyz  = hier.models()[0].atoms().extract_xyz()
  xyz  = sum((model.atoms().extract_xyz() for model in hier.models()[1:]), xyz)
  xyz *= 1./len(hier.models())
  hier.models()[0].atoms().set_xyz(xyz)
  for model in hier.models()[1:]:
    hier.remove_model(model)

  print('Writing new PDB')
  hier.write_pdb_file(file.replace('.pdb','_pdbdist.pdb'),
                      crystal_symmetry = inp.crystal_symmetry_from_cryst1())

  grid = grid / np.nanmax(grid) * norm
  with np.errstate(invalid='ignore'):
    grid[grid < np.nanpercentile(grid, 100-perc)] = np.nan

  with open(file.replace('.pdb','_pdbdist.pml'), 'w') as cmdfile:
    for a in range(lo, hi+1):
      for b in range(lo, a):
        val = grid[a-lo,b-lo]
        if not np.isnan(val):
          print('bond resi {} and name CA, resi {} and name CA'.format(a,b),
                file=cmdfile)
          print('set_bond stick_radius, {}, resi {} or resi {}'.format(val,a,b),
                file=cmdfile)


def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().pdbdist

  if p.params.subtract:
    ref = extract_grid(p.params.subtract, 'mean', p.params.chain,
                       p.params.residue_range, p.params.states)
  else:
    ref = None

  # Extract and align grids
  grids = [extract_grid(file, p.params.mode, p.params.chain,
                        p.params.residue_range, p.params.states,
                        ref)
           for file in p.input.pdb[:4]]
  grids, (lo, hi) = align_grids(grids)

  # Combine grids
  if len(grids) == 1:
    grid = grids[0]
  elif p.params.combine == 'both':
    grid, second = grids[:2]
    second = second.T
    valid  = ~np.isnan(second)
    grid[valid] = second[valid]
  else:
    pairs = zip(*(2*[iter(grids)]))
    for n, (a, b) in enumerate(pairs):
      with np.errstate(divide='ignore', invalid='ignore'):
        merge = funcs[p.params.combine](a, b)
      if n == 0:
        grid = merge
      else:
        merge = merge.T
        valid = ~np.isnan(merge)
        grid[valid] = merge[valid]

  if p.params.full and np.isnan(grid[0,-1]):
    empty       = np.isnan(grid)
    grid[empty] = grid.T[empty]

  # Write pymol stuff
  if p.output.write_pymol_bonds:
    write_pymol_bonds(p.input.pdb[0], grid, lo, hi, p.pymol.normalize, p.pymol.percentile)
     
  # Create custom colormap
  cmaps   = plt.colormaps()
  pos     = next(m for m in cmaps if m.lower()==p.params.positive_cmap.lower())
  poscmap = plt.get_cmap(pos)
  neg     = next(m for m in cmaps if m.lower()==p.params.negative_cmap.lower())
  negcmap = plt.get_cmap(neg)
  high    = p.params.max if p.params.max is not None else np.nanmax(grid)
  low     = p.params.min if p.params.min is not None else np.nanmin(grid)
  if low < 0 < high:
    frac  = high / (high - low)
    poss  = poscmap(np.linspace(0, frac/max(frac,1-frac), int(frac * 5040)))
    negs  = negcmap(np.linspace((1-frac)/max(frac,1-frac), 0, int((1-frac) * 5040)))
    stack = np.vstack((negs, poss))
    cmap  = colors.LinearSegmentedColormap.from_list('composite', stack)
  elif low >= 0:
    cmap  = poscmap
  elif high <= 0:
    cmap  = negcmap

  # Construct label 
  label = 'C$_{\\alpha}$ '
  if   p.params.mode == 'cov' : label += 'position covariance%s ($\mathrm{\AA}^{2}$)'
  elif p.params.mode == 'cc'  : label += 'position correlation%s'
  elif p.params.mode == 'std' : label += 'distance standard deviation%s ($\mathrm{\AA}$)'
  elif p.params.mode == 'var' : label += 'distance variation%s ($\mathrm{\AA}^{2}$)'
  elif p.params.mode == 'mean': label += 'distance mean%s ($\mathrm{\AA}$)'

  if   p.params.combine == 'sub': label = label%' difference'
  elif p.params.combine == 'add': label = label%' sum'
  elif p.params.combine == 'mul': label = label[:label.find('%s')+2]%' product'
  elif p.params.combine == 'div': label = label[:label.find('%s')+2]%' ratio'
  else                          : label = label%''

  # Make plot
  fig, ax = plt.subplots(figsize=(5,4))
  im = ax.imshow(grid, origin='lower', vmin=p.params.min, vmax=p.params.max,
                 extent=(lo-.5,hi+.5,lo-.5,hi+.5), cmap=cmap)
  ax.xaxis.set_ticks_position(p.params.xlabel)
  ax.xaxis.set_label_position(p.params.xlabel)
  ax.set_frame_on(False)
  ax.set_xlabel(p.params.label)
  ax.set_ylabel(p.params.label)
  cb = fig.colorbar(im, label=label)
  fig.tight_layout()
  ax.set_yticks([t for t in ax.get_xticks() if lo <= t <= hi])
  ax.hlines(p.params.lines, lo-0.5, hi+0.5, color=p.params.lc, lw=0.5)
  ax.vlines(p.params.lines, lo-0.5, hi+0.5, color=p.params.lc, lw=0.5)
  plt.savefig(p.output.png_out, dpi=300)

  if p.params.show:
    plt.show()
