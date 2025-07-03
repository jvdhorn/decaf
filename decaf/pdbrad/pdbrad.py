from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import pdb
from . import phil
import numpy as np

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().pdbrad
  print('Reading', p.input.pdb)
  inp    = pdb.input(p.input.pdb)
  hier   = inp.construct_hierarchy()
  sel    = hier.atom_selection_cache().sel_protein()
  hier   = hier.select(sel)
  xrs    = hier.extract_xray_structure()
  mean   = xrs.sites_frac().mean()
  size   = xrs.sites_frac().size()
  other  = xrs.customized_copy()
  other.scatterers().set_sites(flex.vec3_double([mean]*size))
  diff   = xrs.difference_vectors_cart(other)
  radius = diff.rms_length()
  print('Radius  :', radius)
  print('Diameter:', radius * 2)

