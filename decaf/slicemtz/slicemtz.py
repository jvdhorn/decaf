from __future__ import division, print_function
from iotbx import mtz, ccp4_map
from scitbx.array_family import flex
from matplotlib import pyplot as plt, colors, ticker, transforms as tr
from scipy import ndimage as ndi, stats
import numpy as np
from . import phil

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().slicemtz
  print('Reading', p.input.mtz)
  if p.input.mtz.endswith('.mtz'):
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
    cell    = arr.unit_cell().reciprocal_parameters()
  else:
    reader  = ccp4_map.map_reader(p.input.mtz)
    grid    = reader.map_data().as_numpy_array()
    minimum = grid.min()
    cutoff  = (p.input.cutoff if p.input.cutoff is not None else minimum if
               minimum == stats.mode(grid,axis=None).mode else minimum - 1)
    grid[grid<=cutoff] = np.nan
    x, y, z = grid.shape
    offset  = x//2, y//2, z//2
    cell    = list(reader.unit_cell().parameters())
    cell[:3]= list(map(lambda x,y:x/y, cell[:3], reader.unit_cell_grid))
  if p.params.center:
    grid    = np.roll(grid, offset, axis=(0,1,2))
  level   = p.params.slice.strip('hkl')
  dim     = p.params.slice.find(level)
  level   = min(grid.shape[dim]-1, max(0, int(level) + int(offset[dim])))
  sel     = [slice(None)] * 3
  sel[dim]= slice(level-p.params.depth, level+1+p.params.depth)
  slab    = grid[tuple(sel)]
  name    = '_'+p.params.slice
  name   += '_d%d'%p.params.depth
  if p.params.mode != 'sum':
    name += p.params.mode
  with np.errstate(all='ignore'):
    func  = getattr(np, 'nan' + p.params.mode)
    layer = func(slab, axis=dim)
    nans  = np.isnan(slab).all(axis=dim)
    layer[nans] = np.nan

  # Multiply and offset
  layer   = (layer * p.params.multiply) + p.params.add
  off2d   = np.array(layer.shape) // 2

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
  if p.params.max is not None:
    high  = p.params.max
    name += '_max%d'%high
  else:
    high  = np.nanmax(layer)
  if p.params.min is not None:
    low   = p.params.min
    name += '_min%d'%low
  else:
    low   = np.nanmin(layer)

  # Apply log scaling
  if p.params.log:
    with warnings.catch_warnings():
      warnings.simplefilter('ignore', category=RuntimeWarning)
      layer[(layer < 1) & (layer > -1)] = 0
      layer[layer >=  1] =  np.log10( layer[layer >=  1])
      layer[layer <= -1] = -np.log10(-layer[layer <= -1])
      high = np.log10(high) if high >= 1 else -np.log10(-high) if high <= -1 else 0
      low = np.log10(low) if low >= 1 else -np.log10(-low) if low <= -1 else 0
    name += '_log'

  # Autoscale
  if p.params.autoscale:
    values   = layer[~np.isnan(layer)]
    while True:
      mean   = np.mean(values)
      lim    = p.params.clip * np.std(values)
      mask   = (values > mean - lim) & (values < mean + lim)
      if mask.all(): break
      else: values = values[mask]
    high = min(high, mean + lim)
    low  = max(low, mean - lim)

  # Identify Bragg positions
  mask    = np.full(layer.shape, np.nan)
  sc_size = p.params.sc_size
  if -1 not in sc_size:
    sc_size.pop(dim)
    b_ind = np.array(np.meshgrid(*map(range, mask.shape))).reshape(2,-1).T
    b_ind = b_ind[((b_ind-off2d) % sc_size == 0).all(axis=1)]
    mask[tuple(b_ind.T)]  = 1.
    mask[np.isnan(layer)] = np.nan

  # Overlay axes
  if p.params.axes:
    x, y  = mask.shape
    mask[x//2,:] = 1.
    mask[:,y//2] = 1.

  # Increase mask size
  factor  = max(0, p.params.dotsize) * 2 + 1
  bigmask = np.zeros(np.array(mask.shape)*factor)
  bigmask[factor//2::factor,factor//2::factor] = mask
  kernel  = np.zeros((factor*2-1,)*2)
  kernel[factor-1,:] += 1.
  kernel[:,factor-1] += 1.
  bigmask = ndi.convolve(bigmask, kernel, mode='constant') // 2
  bigmask[bigmask==0] = np.nan
  mask    = bigmask

  # Transpose layer and mask
  layer = layer.T
  mask  = mask.T

  # Extract transformations
  lens = list(cell[:3])
  lens.pop(dim)
  ang  = np.radians(90 - cell[3 + dim])
  asp  = lens[1]/lens[0] * np.cos(ang)

  # Zoom or trim
  i, j = off2d
  edge = [-i-0.5, i+0.5, -j-0.5, j+0.5]

  if p.params.zoom is not None:
    clip  = (p.params.zoom+[15])[:3]
    name += '_zoom'

  else:
    if p.params.trim:
      notnans  = ~np.isnan(layer)
      xwhere   = np.argwhere(notnans.any(axis=0))
      ywhere   = np.argwhere(notnans.any(axis=1))
      edge[0] += xwhere.min()
      edge[1] -= layer.shape[1] - xwhere.max() - 1
      edge[2] += ywhere.min()
      edge[3] -= layer.shape[0] - ywhere.max() - 1
    else:
      name += '_full'

    x        = edge[0] + (edge[1] - edge[0]) / 2.
    y        = edge[2] + (edge[3] - edge[2]) / 2.
    z        = max(edge[1] - edge[0], edge[3] - edge[2]) / 2.
    clip     = [x,y,z]

  # Prepare inset
  if p.params.inset is not None:
    inset = (p.params.inset+[15])[:3]
    name += '_inset'
  else:
    inset = None

  # Create custom colormap
  cmaps   = plt.colormaps()
  pos     = next(m for m in cmaps if m.lower()==p.params.positive_cmap.lower())
  poscmap = plt.get_cmap(pos)
  neg     = next(m for m in cmaps if m.lower()==p.params.negative_cmap.lower())
  negcmap = plt.get_cmap(neg)
  frac    = high / (high - low)
  if p.params.normalize_cmap:
    pfrac =    frac  / max(frac, 1-frac)
    nfrac = (1-frac) / max(frac, 1-frac)
  else:
    pfrac = nfrac = 1
  if low < 0 < high:
    poss  = poscmap(np.linspace(0, pfrac, int(frac * 5040)))
    negs  = negcmap(np.linspace(nfrac, 0, int((1-frac) * 5040)))
    stack = np.vstack((negs, poss))
    cmap  = colors.LinearSegmentedColormap.from_list('composite', stack)
  elif low >= 0:
    cmap  = poscmap
  elif high <= 0:
    cmap  = negcmap

  # Make plot
  cnt  = p.params.contours
  dstr = p.params.distribution

  if cnt:
    name += '_cnt%d'%cnt

  if p.params.save:
    fig  = plot_layer(layer, mask, low, high, cmap, ang, asp, clip, cnt, 8, dstr, inset)
    name = p.output.png_out or p.input.mtz.replace('.mtz', name+'.png')
    fig.savefig(name)

  else:
    fig = plot_layer(layer, mask, low, high, cmap, ang, asp, clip, cnt, 1, dstr, inset)
    plt.show()


def axes_to_data(ax, xy):

    return ax.transData.inverted().transform(ax.transAxes.transform(xy))


def plot_section(ax, layer, mask, vmin, vmax, asp, cmap, cnt, ang, clip):

  # Dummy plot
  j,i = np.array(layer.shape) // 2
  ext = [-i-0.5, i+0.5, -j-0.5, j+0.5]
  im  = ax.imshow(layer, cmap=cmap, alpha=0, vmin=vmin, vmax=vmax, aspect=asp,
                  origin='lower', extent=ext, interpolation='none', zorder=2)

  # Apply transformations
  trans  = ax.transData + tr.Affine2D().skew(ang,0)
  x, y   = (ext[0]+ext[1])/2., (ext[2]+ext[3])/2.
  xt, yt = trans.transform((x,y))
  xa, ya = ax.transData.transform((x,y))
  trans += tr.Affine2D().translate(xa-xt, ya-yt)

  im.set_transform(trans)

  ax.transData = trans

  # Actual plot
  if cnt:
    with np.errstate(divide='ignore', invalid='ignore'):
      copy = layer.copy()
      copy[copy<vmin]=vmin
      copy[copy>vmax]=vmax
    kw = {'cmap':cmap, 'transform':trans, 'vmin':vmin, 'vmax':vmax,
          'origin':'lower', 'extent':ext, 'levels':abs(cnt), 'zorder':0}
    if cnt < 0:
      cr = ax.contour(copy, linewidths=0.2, **kw)
    elif cnt > 0:
      cr = ax.contourf(copy, **kw)
    im = plt.cm.ScalarMappable(norm=colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
    im.set_array([])
  else:
    im = ax.imshow(layer, cmap=cmap, transform=trans, vmin=vmin, vmax=vmax,
                   aspect=asp, origin='lower', extent=ext, interpolation='none',
                   zorder=0)

  # Plot axes and Bragg positions
  ax.imshow(mask, cmap='binary_r', transform=trans, aspect=asp,
            origin='lower', extent=ext, interpolation='none', zorder=1)

  # Set plot limits
  x,y,z = clip
  if z <= 1:
    x,y = axes_to_data(ax, (x,y))
    z  *= max(i, j) * 2

  w,h = min(z,z*asp), min(z,z/asp)
  ax.set_xlim((x-w, x+w))
  ax.set_ylim((y-h, y+h))

  return im


def plot_layer(layer, mask, vmin, vmax, cmap, ang, asp, clip, cnt, scale, dstr, inset):

  plt.close()
  plt.rc('font', size=scale*11.)
  figsize = 8 * scale, 6 * scale
  fig, (ax, dx, cax, pad) = plt.subplots(1,4,figsize=figsize,
    gridspec_kw={'width_ratios':[0.85, 0.04, 0.02, 0.09], 'wspace':0, 'hspace':0})
  fig.tight_layout()
  ax.set_axis_off()
  dx.set_axis_off()
  pad.set_axis_off()

  # Plot distribution
  if dstr.lower() not in 'false0':
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

  # Plot data
  im = plot_section(ax, layer, mask, vmin, vmax, asp, cmap, cnt, ang, clip)

  # Plot inset
  if inset is not None:
    x, y = ax.transAxes.transform((0.5,0.5))
    i, j = ax.transAxes.transform((0,0))
    w, h = fig.transFigure.inverted().transform([min(x-i, y-j)] * 2)
    x, y = fig.transFigure.inverted().transform((x,y))
    ins  = fig.add_axes([x-w, y-h, w, h])
    ins.tick_params(left=False,right=False,top=False,bottom=False,labelleft=False,
                    labelright=False,labeltop=False,labelbottom=False)
    plot_section(ins, layer, mask, vmin, vmax, asp, cmap, cnt, ang, inset)
    # Draw rectangle
    corners = [axes_to_data(ins, xy) for xy in [(0,0),(0,1),(1,1),(1,0),(0,0)]]
    ax.plot(*zip(*corners), color='black', lw=0.5)

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
