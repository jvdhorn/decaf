from __future__ import division, print_function
from iotbx import mtz, ccp4_map
from cctbx import crystal, miller
from cctbx.array_family import flex
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().mtz2map
  print('Reading', p.input.mtz)
  obj     = mtz.object(p.input.mtz)
  lores   = max(p.input.resolution)
  hires   = min(p.input.resolution)
  if lores == hires: lores = 9e99
  arr     = obj.crystals()[0].miller_set(False).array(obj.get_column(p.input.lbl
            ).extract_values()).resolution_filter(lores, hires).expand_to_p1()
  data    = arr.data().as_numpy_array()
  ind     = arr.indices().as_vec3_double().as_numpy_array().astype(int)
  origin  = abs(ind).max(axis=0)
  shape   = 2 * origin + 1
  grid    = np.zeros(shape=shape, dtype=float)
  grid[:] = p.input.fill if p.input.fill is not None else min(data.min(),0)-1000
  grid[tuple((-ind + origin).T)] = data
  grid[tuple(( ind + origin).T)] = data

  flex_data = flex.double(grid.ravel())
  flex_data.reshape(flex.grid(tuple(-i-2 for i in origin),
                              tuple(i-1 for i in origin)))

  if p.output.map_out:
    label = p.output.map_out.replace('.map', '')
  else:
    label = p.input.mtz.replace('.mtz', '')

  print('Writing MAP')
  ccp4_map.write_ccp4_map(
    file_name      = label+'.map',
    unit_cell      = arr.crystal_symmetry().unit_cell().reciprocal(),
    unit_cell_grid = (1,1,1),
    space_group    = arr.crystal_symmetry().space_group(),
    map_data       = flex_data,
    labels         = flex.std_string([label+'/'+p.input.lbl]),
  )
