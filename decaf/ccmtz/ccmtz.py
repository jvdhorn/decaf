from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
from . import phil
import matplotlib.pyplot as plt
import numpy as np

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().ccmtz
  hi, lo = (sorted(p.input.resolution) + [float('inf')])[:2]
  print('Reading', p.input.mtz_1)
  obj    = mtz.object(p.input.mtz_1)
  first  = obj.crystals()[0].miller_set(False).array(obj.get_column(
           p.input.lbl_1).extract_values().as_double()).resolution_filter(lo,hi)
  first = first.select(first.data() > p.input.cutoff)
  print('Reading', p.input.mtz_2)
  obj    = mtz.object(p.input.mtz_2)
  second = obj.crystals()[0].miller_set(False).array(obj.get_column(
           p.input.lbl_2).extract_values().as_double())
  second = second.select(second.data() > p.input.cutoff).resolution_filter(lo,hi)
  # Select within limits
  h,k,l  = first.indices().as_vec3_double().as_numpy_array().T
  hlim   = list(map(float, p.input.hlim.split()))
  klim   = list(map(float, p.input.klim.split()))
  llim   = list(map(float, p.input.llim.split()))
  hpos   = ( h >= min(hlim)) & ( h <= max(hlim))
  kpos   = ( k >= min(klim)) & ( k <= max(klim))
  lpos   = ( l >= min(llim)) & ( l <= max(llim))
  hneg   = (-h >= min(hlim)) & (-h <= max(hlim))
  kneg   = (-k >= min(klim)) & (-k <= max(klim))
  lneg   = (-l >= min(llim)) & (-l <= max(llim))
  sel    = flex.bool((hpos & kpos & lpos) | (hneg & kneg & lneg))
  first  = first.select(sel)
  first, second = first.common_sets(second, assert_is_similar_symmetry=False)
  print('Calculating correlation coefficients')
  bins   = first.setup_binner(n_bins=p.input.bins)
  if p.input.subtract:
    func = {'min': np.min, 'mean': np.mean}[p.input.subtract]
    first_data  = first._data.as_numpy_array()
    second_data = second._data.as_numpy_array()
    for n in bins.range_all():
      selection     = bins.array_indices(n).as_numpy_array()
      first_select  = first_data[selection]
      if first_select.size: first_data[selection] -= func(first_select)
      second_select = second_data[selection]
      if second_select.size: second_data[selection] -= func(second_select)
    first._data  = flex.double(first_data)
    second._data = flex.double(second_data)
  corr   = first.correlation(second, use_binning=True, assert_is_similar_symmetry=False)
  corr.show()
  ccall  = first.correlation(second, assert_is_similar_symmetry=False).coefficient()
  print('CC all      :', ccall)
  second = second.customized_copy(crystal_symmetry=first.crystal_symmetry())
  r1i    = first.r1_factor(first.scale(second))
  print('R-factor (I):', r1i)
  r1f    = first.f_sq_as_f().r1_factor(first.scale(second).f_sq_as_f())
  print('R-factor (F):', r1f)

  if p.input.plot: 
    x = ['$\infty$'] + ['%.2f'%corr.binner.bin_d_range(i)[1] for i in corr.binner.range_used()]
    y = [corr.data[i] for i in corr.binner.range_used()]
    xticks = [i-0.5 for i in range(len(x))]
    if len(y) >= 10:
      xsel   = [int(i//(len(y)/10)) for i in range(len(y))]+[10]
      xind   = [xsel.index(i) for i in range(11)]
      x      = [x[i] for i in xind]
      xticks = [xticks[i] for i in xind]
  
    plt.bar(range(len(y)), y, label='Binned')
    plt.plot(xticks, [ccall] * len(xticks))
    plt.plot(xticks, [ccall] * len(xticks), label='All')
    plt.xticks(xticks, x)
    plt.ylim(0, 1)
    plt.xlabel('Resolution shell limits ($\AA$)')
    plt.ylabel('Correlation coefficient')
    plt.legend()
    plt.tight_layout()
    plt.savefig('ccmtz.png', dpi=300)

    if p.input.show:
      plt.show()

