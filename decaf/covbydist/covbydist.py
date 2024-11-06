from __future__ import print_function, division
from iotbx import pdb
from matplotlib import pyplot as plt
from scipy import optimize, stats
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().covbydist

  print('Reading', p.input.pdb)
  hier = pdb.input(p.input.pdb).construct_hierarchy()
  if not p.params.include_neighbours:
    hier = hier.select(hier.atom_selection_cache().sel_chain_id('A'))
  hier = hier.select(hier.atom_selection_cache().sel_name(' CA '))
  xyz  = hier.atoms().extract_xyz().as_numpy_array().reshape(len(hier.models()),-1,3)
  mean = xyz.mean(axis=0)
  cov  = np.cov((xyz-mean).reshape(xyz.shape[0],-1).T)
  cov  = cov[0::3,0::3] + cov[1::3,1::3] + cov[2::3,2::3]
  dist = ((mean[...,None] - mean.T)**2).sum(axis=1)**.5
  cov  = cov[(dist>=p.params.min_distance) & (dist<=p.params.max_distance)]
  dist = dist[(dist>=p.params.min_distance) & (dist<=p.params.max_distance)]
  
  func = lambda r,a,b,g:a*np.exp(-r/g)+b
  fit  = optimize.curve_fit(func, dist, cov, [0,0,dist.max()/2])
  rng  = np.linspace(dist.min(), dist.max(), p.params.bins)
  
  y,x,_= stats.binned_statistic(dist, cov, bins=p.params.bins, statistic='mean')
  e,_,_= stats.binned_statistic(dist, cov, bins=p.params.bins, statistic='std')
  x    = (x[1:] + x[:-1]) / 2.
  
  plt.figure(figsize=(4,3))
  plt.errorbar(x,y,e, ls='None', fmt='.', capsize=1.5, label='Covariance spread in bin', zorder=1)
  plt.plot(rng, func(rng, *fit[0]), label='Decay fit ($\\gamma = {:.2f}$)'.format(fit[0][2]))
  plt.legend(frameon=False)
  plt.axhline(0, color='black', linewidth=0.3)
  plt.xlabel('$C_{\\alpha}$ pair distance ($\mathrm{\AA}$)')
  plt.ylabel('$C_{\\alpha}$ pair covariance ($\mathrm{\AA}^{2}$)')
  plt.tight_layout()
  plt.savefig(p.output.png_out, dpi=300)

  if p.params.show:
    plt.show()
