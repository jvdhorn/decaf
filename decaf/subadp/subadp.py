from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import pdb
from . import phil
import numpy as np

funcs = {'sub': lambda a,b:a-b,
         'add': lambda a,b:a+b,
        }

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().subadp
  print('Reading', p.input.pdb_1)
  inp    = pdb.input(p.input.pdb_1)
  hier   = inp.construct_hierarchy(False,False,False)
  for model in hier.models()[1:]: hier.remove_model(model)
  common = {atom.pdb_label_columns() for atom in hier.atoms()}

  for file in p.input.pdb_2:
    print('Reading', file)
    # Reduce to common atoms
    second  = pdb.input(file).construct_hierarchy(False,False,False)
    for model in second.models()[1:]: second.remove_model(model)
    common &= {atom.pdb_label_columns() for atom in second.atoms()}
    hier = hier.select(
             flex.bool(atom.pdb_label_columns() in common for atom in hier.atoms())
           )
    second = second.select(
             flex.bool(atom.pdb_label_columns() in common for atom in second.atoms())
           )
    # Align atom arrays
    hier.sort_atoms_in_place()
    second.sort_atoms_in_place()
    # Apply operation
    mode = funcs[p.params.mode]
    hier.atoms().set_uij(mode(hier.atoms().extract_uij(), second.atoms().extract_uij()))
    hier.atoms().set_b(mode(hier.atoms().extract_b(), second.atoms().extract_b()))


  print('Writing new PDB')
  files    = [f.replace('.pdb','') for f in [p.input.pdb_1]+p.input.pdb_2]
  new_file = ('_%s_'%p.params.mode).join(f.strip('./').replace('/','_')
              for f in files)+'.pdb'
  hier.write_pdb_file(new_file,
                      crystal_symmetry = inp.crystal_symmetry_from_cryst1())
