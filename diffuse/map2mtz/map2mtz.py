from __future__ import division, print_function
from iotbx import ccp4_map
from cctbx import crystal, miller
from cctbx.array_family import flex
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().map2mtz
  print('Reading', p.input.map)
  reader  = ccp4_map.map_reader(p.input.map)
  
  lores   = max(p.input.resolution)
  hires   = min(p.input.resolution)
  if lores == hires: lores = float('inf')

  data    = reader.map_data().as_numpy_array()
  x, y, z = data.shape
  offset  = x//2, y//2, z//2

  cell    = list(reader.unit_cell().reciprocal_parameters())
  cell[:3]= list(map(lambda x,y:x*y, cell[:3], reader.unit_cell_grid))
  symm    = crystal.symmetry(cell, 'P1')
  mset    = miller.build_set(symm, False, max_index=offset)
  indices = mset.indices().as_vec3_double().as_numpy_array().astype(int) + offset
  data    = flex.double(data[tuple(indices.T)].astype(float))
  if p.input.cutoff == None:
    mask  = data > data.min_max_mean().min
  else:
    mask  = data > p.input.cutoff
  array   = mset.array(data).select(mask).resolution_filter(lores, hires)
  mtzdset = array.as_mtz_dataset(column_root_label=p.output.lbl,column_types='J')

  if p.output.mtz_out:
    label = p.output.mtz_out.replace('.mtz', '')
  else:
    label = p.input.map.replace('.map', '')

  print('Writing MTZ')
  mtzdset.mtz_object().write(label+'.mtz')
