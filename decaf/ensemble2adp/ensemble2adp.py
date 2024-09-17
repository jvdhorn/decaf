from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import pdb
from . import phil
import numpy as np

biso_const = (8./3) * (np.pi ** 2)

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

  if len(hier.models()) < 2:
    print('Please provide multistate PDB'); return

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
  for i_atom, atom in enumerate(hier.models()[0].atoms()):
    M = np.array([model.atoms()[i_atom].xyz for model in hier.models()]
                 ).reshape(-1,3) - [means[i_atom]]
    cov = np.dot(M.T, M) / (M.shape[0] - 1)
    atom.set_uij(tuple(cov[(0,1,2,0,0,1),(0,1,2,1,2,2)]))
    atom.set_b((M**2).sum()/M.shape[0] * biso_const)

  hier.models()[0].atoms().set_xyz(means)

  for model in hier.models()[1:]:
    hier.remove_model(model)

  print('Writing new PDB')
  hier.write_pdb_file(p.input.pdb.replace('.pdb','_adp.pdb'),
                      crystal_symmetry = inp.crystal_symmetry_from_cryst1())