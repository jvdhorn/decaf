from __future__ import print_function, division
from iotbx import pdb
from cctbx.array_family import flex
from . import phil
import matplotlib.pyplot as plt
import numpy as np



def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().pdbdist

  print('Reading', p.input.pdb)
  hier = pdb.input(p.input.pdb).construct_hierarchy()
  sel  = hier.atom_selection_cache().selection('name CA')
  calp = hier.select(sel)
  resi = [atom.parent().parent().resseq_as_int() for atom in
          calp.models()[p.params.chain].chains()[p.params.chain].atoms()]
  lim  = p.params.limit_ca or max(resi)
  sel  = flex.size_t([resi.index(n) for n in sorted(set(resi)) if n <= lim])
  resi = [resi[i] for i in sel]
  grid = []

  for model in calp.models():
    xyz  = model.chains()[p.params.chain].atoms().extract_xyz().select(sel)
    dmat = [(xyz - coord).dot()**0.5 for coord in xyz]
    grid.append(dmat)

  grid = np.array(grid).std(axis=0)
  grid[np.triu_indices(grid.shape[0],1)] = np.nan

  lo, hi = min(resi), max(resi)
  for n in range(lo, hi):
    if n not in resi:
      grid = np.insert(grid, n-lo, np.nan, axis=0)
      grid = np.insert(grid, n-lo, np.nan, axis=1)

  fig, ax = plt.subplots(figsize=(5,4))
  im = ax.imshow(grid, origin='lower', vmin=p.params.min, vmax=p.params.max,
                 extent=(lo-.5,hi+.5,lo-.5,hi+.5), cmap=p.params.cmap)
  ax.xaxis.set_ticks_position('top')
  ax.xaxis.set_label_position('top')
  ax.set_frame_on(False)
  ax.set_xlabel('Residue number')
  ax.set_ylabel('Residue number')
  cb = fig.colorbar(im)
  cb.set_label('C$_{\\alpha}$ distance standard deviation ($\mathrm{\AA}$)')
  fig.tight_layout()
  plt.savefig(p.output.png_out, dpi=300)

  if p.params.show:
    plt.show()
