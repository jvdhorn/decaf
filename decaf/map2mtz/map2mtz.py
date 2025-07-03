from __future__ import division, print_function
from iotbx import ccp4_map
from cctbx import crystal, miller
from cctbx.array_family import flex
from scipy import stats
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

  if p.input.center:
    data  = np.roll(data, offset, axis=(0,1,2))

  mode    = stats.mode(data,axis=None)
  cutoff  = float((p.input.cutoff if p.input.cutoff is not None else mode.mode
            if mode.count/data.size > 0.4 else data.min() - 1))
  indices = np.mgrid[tuple(map(slice, data.shape))].reshape(3,-1).T - offset
  allowed = data.ravel() > cutoff
  indices = flex.miller_index(
              indices[allowed].ravel().view('int64,int64,int64').tolist()
            )
  data    = flex.float(data.ravel()[allowed]).as_double()
  cell    = list(reader.unit_cell().reciprocal_parameters())
  cell[:3]= list(map(lambda x,y:x*y, cell[:3], reader.unit_cell_grid))
  symm    = crystal.symmetry(cell, 'P1')

  array   = miller.set(symm, indices).array(data).resolution_filter(lores, hires)

  mtzdset = array.as_mtz_dataset(column_root_label=p.output.lbl,column_types='J')

  if p.output.mtz_out:
    label = p.output.mtz_out.replace('.mtz', '')
  else:
    label = p.input.map.replace('.map', '')

  print('Writing MTZ')
  mtzdset.mtz_object().write(label+'.mtz')
