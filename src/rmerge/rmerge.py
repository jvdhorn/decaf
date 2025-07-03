from __future__ import division, print_function
import phil
from iotbx import mtz
from cctbx import crystal

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

def main(args, log):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().rmerge
  print('Reading', p.input.mtz_1)
  mtz_object      = mtz.object(p.input.mtz_1)
  miller_array    = mtz_object.crystals()[0].miller_set(False).array(mtz_object.get_column(
                    p.params.array).extract_values().as_double()).expand_to_p1()
  new_space_group = p.params.space_group or str(miller_array.space_group_info())
  new_symm        = crystal.symmetry(unit_cell          = miller_array.unit_cell(),
                                     space_group_symbol = new_space_group)
  print('New space group:', new_symm.space_group_info().symbol_and_number())
  array_reduced   = miller_array.customized_copy(crystal_symmetry = new_symm)
  merged          = array_reduced.merge_equivalents()
  print('R_merge:', merged.r_merge())
