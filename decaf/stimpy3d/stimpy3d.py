from __future__ import division, print_function
from scitbx.array_family import flex
from iotbx import mtz
from scipy import ndimage, stats, optimize
from . import phil
import matplotlib.pyplot as plt
import numpy as np
import os

def nw_stats(N, mean, var, skew):

  if skew < 0: skew = 0
  Sigma = (N**2 * skew/2.)**(1./3) * var**(1./2)
  Var   = var * (1-(N/4.)**(1./3) * skew**(2./3))
  Mu    = mean - Sigma
  return (Sigma, Var, Mu)


def filter_sigma(array, sigma=3, repeat=99):

  array            = np.array(array)
  while repeat:
    std              = np.std(array)
    mean             = np.mean(array)
    min_lim, max_lim = mean - sigma * std, mean + sigma * std
    selection        = (array >= min_lim) & (array <= max_lim)
    array_select     = array[selection]
    if array_select.size < array.size:
      array   = array_select
      repeat -= 1
    else:
      array   = array_select
      repeat  = 0
  return array


def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().stimpy3d
  pre    = 'stimpy3d_' + p.input.mtz.replace('.mtz', '_' + p.input.lbl.lower())
  print('Reading', p.input.mtz)
  obj   = mtz.object(p.input.mtz)
  first = obj.crystals()[0].miller_set(False).array(obj.get_column(
          p.input.lbl).extract_values()).expand_to_p1()
  first = first.select(first.data() > p.input.cutoff)
  d_min = first.d_spacings().data().min_max_mean().min
  first.setup_binner(d_max=float('inf'), d_min=d_min, n_bins=p.input.bins)
  bins  = first.binner()
  data  = first.data().as_numpy_array()
  ind   = bins.bin_indices().as_numpy_array() - 1
  
  if p.input.write:
    drc = (p.input.directory if p.input.directory is not None else pre) + '/'
    try: os.mkdir(drc)
    except OSError: pass

  mu    = np.zeros(p.input.bins)
  sigma = np.zeros(p.input.bins)
  Ntwin = []

  print()
  for n in range(p.input.bins):
    mask            = ind == n
    data_sel        = data[mask]
    data_filt       = filter_sigma(data_sel, p.input.sigma_cutoff)
    excess          = (data_sel < data_filt.min()) | (data_sel > data_filt.max())
    mean, var, skew = data_filt.mean(), data_filt.var(), stats.skew(data_filt)
    print('Region', n, 'range:', *bins.bin_d_range(n+1))
    print('N={:.0f}; mean={:.2f}; var={:.2f}; skew={:.2f}'.format(
          data_filt.size, mean, var, skew))
    if p.input.use_fit_params:
      fit             = stats.skewnorm.fit(data_sel)
      skewnorm        = stats.skewnorm(*fit)
      mean, var, skew = skewnorm.stats(moments='mvs')
      print('Fit stats: mean={:.2f}; var={:.2f}; skew={:.2f}'.format(mean, var, skew))
    Sigma, Var, Mu  = nw_stats(p.input.N, mean, var, skew)
    print('Noisy Wilson stats: Sigma={:.2f}; Var={:.2f}; Mu={:.2f}'.format(Sigma, Var, Mu))
    muzerofunc      = lambda x: nw_stats(x, mean, var, skew)[-1]**2
    muzero          = optimize.minimize_scalar(muzerofunc, bounds=(0,9e9),
                                               method='bounded').x
    print('Mu-zero found at N={}'.format(muzero))
    if skew > 0: Ntwin.append(muzero)
    mu[n]           = max(0, Mu)
    sigma[n]        = Sigma

    print()
    if p.input.write:
      name = p.input.mtz.replace('.mtz','_region_{}_counts.txt'.format(n))
      with open(drc + name, 'w') as txt:
        for val in data_sel: print(val, end='\n', file=txt)
      plt.figure()
      y, x = plt.hist(data_filt, bins=25)[:2]
      x    = (x[:-1] + x[1:]) / 2
      if p.input.use_fit_params:
        pdf  = skewnorm.pdf(x)
        plt.plot(x, pdf * mask.sum() * (x[1] - x[0]))
      name = p.input.mtz.replace('.mtz','_region_{}_statistics.png'.format(n))
      plt.savefig(drc + name)
      plt.close()

    if p.input.remove_excess:
      data_sel_copy         = data_sel.copy()
      data_sel_copy[excess] = np.nan
      data[mask]            = data_sel_copy

  print('Average N-twin at Mu-zero over all shells:', sum(Ntwin)/len(Ntwin))
  print()

  mu         = ndimage.uniform_filter(mu, p.input.smooth_kernel)
  background = data.copy()
  signal     = data.copy()

  for n in range(p.input.bins):
    sel             = ind == n
    background[sel] = mu[n]
    signal[sel]     = sigma[n]

  processed  = data - background

  mtz_out = pre + '.mtz'

  print('Writing', mtz_out)

  keep             = ~np.isnan(data)
  first            = first.select(flex.bool(keep))
  processed_array  = first.customized_copy(data = flex.float(processed[keep]).as_double())
  background_array = first.customized_copy(data = flex.float(background[keep]).as_double())
  signal_array     = first.customized_copy(data = flex.float(signal[keep]).as_double())

  dataset = processed_array.as_mtz_dataset(column_root_label = p.input.lbl,
                                           column_types      = 'J')
  dataset.add_miller_array(background_array,
                           column_root_label = 'BACKGROUND',
                           column_types      = 'J')
  dataset.add_miller_array(signal_array,
                           column_root_label = 'SIGNAL',
                           column_types      = 'J')

  dataset.mtz_object().write(mtz_out)
