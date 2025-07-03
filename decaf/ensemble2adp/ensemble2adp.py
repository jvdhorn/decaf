from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import pdb
from . import phil
import numpy as np

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().ensemble2adp
  print('Reading', p.input.pdb)
  inp    = pdb.input(p.input.pdb)
  hier   = inp.construct_hierarchy()
  if p.input.models:
    for model in hier.models()[p.input.models:]:
      hier.remove_model(model)

  # Reduce to common atoms
  common = set()
  for model in hier.models():
    id_strs = {atom.pdb_label_columns() for atom in model.atoms()}
    if not common:
      common  = id_strs
    else:
      common &= id_strs

  hier = hier.select(
           flex.bool(atom.pdb_label_columns() in common for atom in hier.atoms())
         )
 
  # Find mean atomic positions 
  means  = hier.models()[0].atoms().extract_xyz()
  for model in hier.models()[1:]:
    means += model.atoms().extract_xyz()
  means *= 1./len(hier.models())

  print('Calculating covariances')
  adps   = []
  isos   = []
  for i_atom, mean in enumerate(means):
    M = np.array([model.atoms()[i_atom].xyz for model in hier.models()]
                 ).reshape(-1,3) - [mean]
    cov = sum((np.dot(k[:,None],k[None,:]) for k in M), np.zeros((3,3))
              ) / (M.shape[0] - 1 or 1)
    adps.append(tuple(cov[(0,1,2,0,0,1),(0,1,2,1,2,2)]))
    isos.append((M**2).sum()/M.shape[0] * (8./3) * np.pi**2)

  for model in hier.models()[1:]:
    hier.remove_model(model)

  hier.models()[0].atoms().set_xyz(means)

  for atom, adp, biso in zip(hier.models()[0].atoms(), adps, isos):
    atom.set_uij(adp)
    atom.set_b(biso)

  print('Writing new PDB')
  hier.write_pdb_file(p.input.pdb.replace('.pdb','_adp.pdb'),
                      crystal_symmetry = inp.crystal_symmetry_from_cryst1())
