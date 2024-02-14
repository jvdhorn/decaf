from __future__ import division, print_function
from iotbx import mtz
from scitbx.array_family import flex
from matplotlib import pyplot as plt, colors, ticker, transforms as tr
from scipy import ndimage as ndi
import numpy as np
import warnings
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().slicemtz
  print('Reading', p.input.mtz)
  obj     = mtz.object(p.input.mtz)
  labels  = obj.column_labels()
  label   = p.input.lbl if p.input.lbl in labels else labels[3]
  hi, lo  = (sorted(p.params.resolution) + [float('inf')])[:2]
  arr     = obj.crystals()[0].miller_set(False).array(obj.get_column(label
            ).extract_values()).expand_to_p1().resolution_filter(lo, hi)
  data    = arr.data().as_numpy_array()
  ind     = arr.indices().as_vec3_double().as_numpy_array().astype(int)
  offset  = abs(ind).max(axis=0)
  shape   = 2 * offset + 1
  grid    = np.zeros(shape=shape)
  grid[:] = np.nan
  grid[tuple((-ind + offset).T)] = data
  grid[tuple(( ind + offset).T)] = data
  level   = p.params.slice.strip('hkl')
  dim     = p.params.slice.find(level)
  level   = int(level) + offset[dim]
  sel     = [slice(None)] * 3
  sel[dim]= slice(level-p.params.depth, level+1+p.params.depth)
  slab    = grid[tuple(sel)]
  name    = '_'+p.params.slice
  name   += '_d%d'%p.params.depth
  if p.params.mode != 'sum':
    name += p.params.mode
  with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=RuntimeWarning)
    func  = getattr(np, 'nan' + p.params.mode)
    layer = func(slab, axis=dim)
    nans  = np.isnan(slab).all(axis=dim)
    layer[nans] = np.nan

  # Multiply and offset
  layer   = (layer * p.params.multiply) + p.params.add

  # Apply filter
  if p.params.filter_size > 0:
    filt          = {'gaussian' : ndi.gaussian_filter,
                     'gausslap' : ndi.gaussian_laplace,
                     'uniform'  : ndi.uniform_filter,
                     'min'      : ndi.minimum_filter,
                     'max'      : ndi.maximum_filter,
                     'median'   : ndi.median_filter}[p.params.fmode]
    nans          = np.isnan(layer)
    layer[nans]   = 0.
    layer         = filt(layer, p.params.filter_size, mode='constant')
    layer[nans]   = np.nan
    if p.params.fmode == 'gaussian':
      norm          = np.ones(layer.shape)
      norm[nans]    = 0.
      norm          = filt(norm, p.params.filter_size, mode='constant')
      layer[~nans] /= norm[~nans]
    name         += '_filt%d'%p.params.filter_size
    if p.params.fmode != 'gaussian':
      name += p.params.fmode
  
  # Redefine min and max
  high    = p.params.max
  low     = p.params.min
  if high is not None:
    name += '_max%d'%high
  if low is not None:
    name += '_min%d'%low

  # Apply log scaling
  if p.params.log_intensities:
    with warnings.catch_warnings():
      warnings.simplefilter('ignore', category=RuntimeWarning)
      layer[(layer < 1) & (layer > -1)] = 0
      layer[layer >=  1] =  np.log10( layer[layer >=  1])
      layer[layer <= -1] = -np.log10(-layer[layer <= -1])
      if high is not None:
        high = np.log10(high) if high >= 1 else -np.log10(-high) if high <= -1 else 0
      if low is not None:
        low = np.log10(low) if low >= 1 else -np.log10(-low) if low <= -1 else 0
    name += '_log'

  # Identify Bragg positions
  off2d   = np.array(layer.shape) // 2
  mask    = np.full(layer.shape, np.nan)
  sc_size = p.params.sc_size
  if -1 not in sc_size:
    sc_size.pop(dim)
    b_ind = np.array(np.meshgrid(*map(range, mask.shape))).reshape(2,-1).T
    b_ind = b_ind[((b_ind-off2d) % sc_size == 0).all(axis=1)]
    mask[tuple(b_ind.T)]  = 1.
    mask[np.isnan(layer)] = np.nan

  # Transpose layer and mask
  layer = layer.T
  mask  = mask.T

  # Increase mask size
  factor  = max(0, p.params.dotsize) * 2 + 1
  bigmask = np.full(np.array(mask.shape)*factor, np.nan)
  bigmask[factor//2::factor,factor//2::factor] = mask
  mask    = bigmask

  # Overlay axes
  if p.params.axes:
    x, y  = mask.shape
    mask[x//2,:] = 1.
    mask[:,y//2] = 1.

  # Trim edges
  i, j = off2d
  clip = [-i-0.5, i+0.5, -j-0.5, j+0.5]

  if p.params.zoom is not None:
    x,y,z = (p.params.zoom+[15])[:3]
    clip  = [x-z, x+z, y-z, y+z]
    name += '_x%d_y%d'%(x,y)

  elif p.params.trim:
    notnans  = ~np.isnan(layer)
    xwhere   = np.argwhere(notnans.any(axis=0))
    ywhere   = np.argwhere(notnans.any(axis=1))
    clip[0] += xwhere.min()
    clip[1] -= layer.shape[1] - xwhere.max() - 1
    clip[2] += ywhere.min()
    clip[3] -= layer.shape[0] - ywhere.max() - 1

  else:
    name += '_full'

  # Create custom colormap
  cmaps   = plt.colormaps()
  pos     = next(m for m in cmaps if m.lower()==p.params.positive_cmap.lower())
  poscmap = plt.get_cmap(pos)
  neg     = next(m for m in cmaps if m.lower()==p.params.negative_cmap.lower())
  negcmap = plt.get_cmap(neg)
  high    = high if high is not None else np.nanmax(layer)
  low     = low if low is not None else np.nanmin(layer)
  if low < 0 < high:
    frac  = high / (high - low)
    poss  = poscmap(np.linspace(0, 1, int(frac * 5040)))
    negs  = negcmap(np.linspace(1, 0, int((1-frac) * 5040)))
    stack = np.vstack((negs, poss))
    cmap  = colors.LinearSegmentedColormap.from_list('composite', stack)
  elif low >= 0:
    cmap  = poscmap
  elif high <= 0:
    cmap  = negcmap

  # Extract transformations
  cell = arr.unit_cell().reciprocal_parameters()
  lens = list(cell[:3])
  lens.pop(dim)
  ang  = np.radians(90 - cell[3 + dim])
  asp  = lens[1]/lens[0] * np.cos(ang)

  # Make plot
  cnt  = p.params.contours
  dstr = p.params.distribution

  if p.params.save:
    fig  = plot_layer(layer, mask, low, high, cmap, ang, asp, clip, cnt, 8, False, dstr)
    fig.savefig(p.input.mtz.replace('.mtz', name+'.png'))

  else:
    fig = plot_layer(layer, mask, low, high, cmap, ang, asp, clip, cnt, 1, True, dstr)
    plt.show()


def plot_layer(layer, mask, vmin, vmax, cmap, ang, asp, clip, cnt, scale, rep, dstr):

  plt.close()
  plt.rc('font', size=scale*11.)
  figsize = 8 * scale, 6 * scale
  fig, (ax, dx, cax, pad) = plt.subplots(1,4,figsize=figsize,
    gridspec_kw={'width_ratios':[0.85, 0.04, 0.03, 0.08], 'wspace':0, 'hspace':0})
  fig.tight_layout()
  ax.set_axis_off()
  dx.set_axis_off()
  pad.set_axis_off()

  # Plot distribution
  if dstr != 'False':
    bins = 256
    log  = (dstr == 'Log')
    vals = layer[~np.isnan(layer)]
    vals[vals<vmin] = vmin
    vals[vals>vmax] = vmax
    y, x, patches = dx.hist(vals, bins=bins, orientation='horizontal',
                            range=(vmin, vmax), log=log)
    dx.hist(vals, bins=bins, orientation='horizontal', histtype='step',
            range=(vmin, vmax), ec='black', lw=.5, log=log)
    clrs = cmap(np.linspace(0, 1, bins))
    for p, c in zip(patches, clrs): p.set_facecolor(c)
    dx.margins(x=0, y=.001)
    dx.invert_xaxis()
    pos = dx.get_position()
    dx.set_position([pos.x0, pos.y0-.025, pos.width, pos.height])

  # Plot image
  j,i = np.array(layer.shape) // 2
  ext = [-i-0.5, i+0.5, -j-0.5, j+0.5]
  im  = ax.imshow(layer, cmap=cmap, vmin=vmin, vmax=vmax, aspect=asp,
                  origin='lower', extent=ext, interpolation='none')

  # Apply transformations
  trans  = ax.transData + tr.Affine2D().skew(ang,0)
  x, y   = (ext[0]+ext[1])/2., (ext[2]+ext[3])/2.
  xt, yt = trans.transform((x,y))
  xa, ya = ax.transData.transform((x,y))
  trans += tr.Affine2D().translate(xa-xt, ya-yt)

  im.set_transform(trans)

  ax.transData = trans

  ax.set_xlim((clip[0], clip[1]))
  ax.set_ylim((clip[2], clip[3]))

  # Plot contour
  if cnt and vmin <= 0 <= vmax:
    ax.contour(layer, colors='black', levels=[0], transform=trans,
               aspect=asp, origin='lower', extent=ext, linewidths=0.2)

  # Plot axes and Bragg positions
  ax.imshow(mask, cmap='binary_r', transform=trans, aspect=asp,
            origin='lower', extent=ext, interpolation='none')

  # Enforce shown values
  for _ in range(100 * rep):
    ax.imshow(layer, cmap=cmap, alpha=0, transform=trans, aspect=asp,
              origin='lower', extent=ext, interpolation='none')

  # Plot colorbar
  scf = ticker.ScalarFormatter()
  scf.set_powerlimits((-2,2))
  scf.set_useMathText(True)
  bar = fig.colorbar(im, cax=cax, format=scf)
  bar.ax.tick_params(width=2*scale, size=8.*scale, labelsize=22. * scale)
  bar.ax.yaxis.get_offset_text().set(size=22. * scale)
  bar.ax.yaxis.get_offset_text().set_position((0, 0))
  pos = bar.ax.get_position()
  bar.ax.set_position([pos.x0, pos.y0-.025, pos.width, pos.height])

  return fig
