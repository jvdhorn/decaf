from __future__ import division, print_function
from iotbx import mtz, ccp4_map
from scitbx.array_family import flex
from matplotlib import pyplot as plt, colors, ticker, transforms as tr
from scipy import ndimage as ndi, stats, spatial
import numpy as np
from . import phil


def log10(value):

  value = np.array(value)
  value[(-1 < value) & (value < 1)] = 0
  value[value >=  1] =  np.log10( value[value >=  1])
  value[value <= -1] = -np.log10(-value[value <= -1])

  return value


class Hull(spatial.ConvexHull):

  def mask_inside(self, pts):

    mask = np.array(
             [np.dot(eq[:-1], pts.T) + eq[-1] < 1e-12 for eq in self.equations]
           ).T.all(axis=1).ravel()

    return mask


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
    grid    = np.full(shape, np.nan)
    grid[tuple((-ind + offset).T)] = data
    grid[tuple(( ind + offset).T)] = data
    cell    = arr.unit_cell().reciprocal_parameters()
    frac    = arr.unit_cell().reciprocal().fractionalize
    orth    = arr.unit_cell().reciprocal().orthogonalize
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
    frac    = reader.unit_cell().fractionalize
    orth    = reader.unit_cell().orthogonalize
  if p.params.center:
    grid    = np.roll(grid, offset, axis=(0,1,2))
  slevel  = p.params.slice.lower().strip('hkl')
  dim     = p.params.slice.find(slevel)
  level   = min(grid.shape[dim]-1, max(0, int(slevel) + int(offset[dim])))
  label   = p.params.slice.lower().replace(slevel, str(int(level-offset[dim])))
  sel     = [slice(None)] * 3
  sel[dim]= slice(level-p.params.depth, level+1+p.params.depth)
  slab    = grid[tuple(sel)]
  name    = '_'+label
  name   += '_d%d'%p.params.depth
  if p.params.mode != 'sum':
    name += p.params.mode
  with np.errstate(all='ignore'):
    func  = getattr(np, 'nan' + p.params.mode)
    layer = func(slab, axis=dim).T
    nans  = np.isnan(slab).all(axis=dim).T
    layer[nans] = np.nan

  # Construct projection
  if p.params.projection:
    width  = 721
    if p.params.projection == 'gp':
      phi  = np.arcsin(np.linspace(-1, 1, int(width // np.pi * 2 + 1)))
    elif p.params.projection == 'eq':
      phi  = np.linspace(-1, 1, width // 2 + 1) * np.pi/2.
    elif p.params.projection == 'mercator':
      phi  = 2 * np.arctan(np.exp(np.linspace(-1, 1, width) * np.pi)) - np.pi/2.
    radius = min(orth(offset))
    xy     = 1j ** np.linspace(-2, 2, width)
    z      = np.sin(phi)
    xyz    = np.hstack((np.outer((1-z*z)**0.5, xy).view(float).reshape(-1,2),
                        z.repeat(xy.size)[...,None])) * radius
    if 'hk' in p.params.slice:
      pass
    elif 'kl' in p.params.slice:
      xyz = xyz[:,(1,2,0)]
    else:
      xyz = xyz[:,(2,0,1)]
    hkl    = frac(flex.vec3_double(xyz)).as_numpy_array().astype(int) + offset
    layer  = grid[tuple(hkl.T)].reshape(z.size, -1)
    nans   = np.isnan(layer)
    cell   = (1,1,1,90,90,90)
    del xyz, hkl

  # Multiply and offset
  layer   = (layer * p.params.multiply) + p.params.add
  off2d   = np.array(layer.shape) // 2

  # Apply filter
  if p.params.filter > 0:
    filt          = {'gaussian' : ndi.gaussian_filter,
                     'gausslap' : ndi.gaussian_laplace,
                     'uniform'  : ndi.uniform_filter,
                     'min'      : ndi.minimum_filter,
                     'max'      : ndi.maximum_filter,
                     'median'   : ndi.median_filter}[p.params.fmode]
    layer[nans]   = 0.
    layer         = filt(layer, p.params.filter, mode='constant')
    layer[nans]   = np.nan
    if p.params.fmode == 'gaussian':
      norm          = np.ones(layer.shape)
      norm[nans]    = 0.
      norm          = filt(norm, p.params.filter, mode='constant')
      layer[~nans] /= norm[~nans]
    name         += '_filt%d'%p.params.filter
    if p.params.fmode != 'gaussian':
      name += p.params.fmode
  
  # Apply log scaling
  if p.params.log:
    with np.errstate(all='ignore'):
      layer = log10(layer)
    name += '_log'

  # Autoscale
  high = np.nanmax(layer)
  low  = np.nanmin(layer)
  if p.params.autoscale:
    values = layer[~nans]
    while values.size:
      mean = np.mean(values)
      lim  = p.params.sigma * np.std(values)
      mask = (values > mean - lim) & (values < mean + lim)
      if mask.all(): break
      else         : values = values[mask]
    else:
      mean = np.nanmean(layer)
      lim  = p.params.sigma * np.nanstd(layer)
    high = min(high, mean + lim)
    low  = max(low, mean - lim)

  # Overrule high and low
  if p.params.max is not None:
    high  = log10(p.params.max) if p.params.log else p.params.max
    name += '_max%d'%high
  if p.params.min is not None:
    low   = log10(p.params.min) if p.params.log else p.params.min
    name += '_min%d'%low

  # Identify Bragg positions
  mask    = np.full(layer.shape, np.nan)
  if p.params.sc_size is not None:
    sc_size = list(p.params.sc_size)
    sc_size.pop(dim)
    b_ind = np.array(np.meshgrid(*map(range, mask.shape))).reshape(2,-1).T
    b_ind = b_ind[((b_ind-off2d) % sc_size[::-1] == 0).all(axis=1)]
    mask[tuple(b_ind.T)] = 1.

  # Overlay axes
  if p.params.axes:
    x, y  = mask.shape
    mask[:,y//2] = 1.
    mask[x//2,:] = 1.
    if p.params.projection:
      mask[:,y//4    ] = 1.
      mask[:,y//4 * 3] = 1.

  # Trim mask to convex hull
  if mask.any():
    hull  = Hull(np.array(np.where(~nans)).T)
    where = np.array(np.where(mask)).T
    mask[tuple(where[~hull.mask_inside(where)].T)] = np.nan

  # Trim
  i, j = off2d
  ext  = [-j-0.5, j+0.5, -i-0.5, i+0.5]
  clip = None

  if p.params.trim:
    notnans  = ~nans
    xwhere   = np.argwhere(notnans.any(axis=0))
    ywhere   = np.argwhere(notnans.any(axis=1))
    nwext    = ext[:]
    nwext[0] = ext[0] + xwhere.min()
    nwext[1] = ext[1] - (layer.shape[1] - xwhere.max() - 1)
    nwext[2] = ext[2] + ywhere.min()
    nwext[3] = ext[3] - (layer.shape[0] - ywhere.max() - 1)
    layer    = layer[int(nwext[2]-ext[2]):int(nwext[3]-ext[3]) or None,
                     int(nwext[0]-ext[0]):int(nwext[1]-ext[1]) or None]
    mask     = mask [int(nwext[2]-ext[2]):int(nwext[3]-ext[3]) or None,
                     int(nwext[0]-ext[0]):int(nwext[1]-ext[1]) or None]
    ext      = nwext
    clip     = []

  # Prepare zoom
  if p.params.zoom is not None:
    clip  = (p.params.zoom+[15])[:3]
    name += '_zoom'

  # Prepare inset
  if p.params.inset is not None:
    inset = (p.params.inset+[15])[:3]
    name += '_inset'
  else:
    inset = None

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

  # Extract transformations
  lens = list(cell[:3])
  lens.pop(dim)
  ang  = np.radians(90 - cell[3 + dim])
  asp  = lens[1]/lens[0] * np.cos(ang)

  # Make plot
  cnt  = p.params.contours
  dstr = p.params.distribution
  loc  = p.params.position
  fs   = p.output.figscale

  if not p.params.label:
    label = ''
  elif p.params.depth > 0:
    label = label.replace(label.strip('hkl'), 'n')

  if cnt:
    name += '_cnt%d'%cnt

  if p.params.save:
    fig  = plot_layer(layer, mask, low, high, cmap, ang, asp,
                      clip, cnt, fs, dstr, inset, ext, loc, label)
    name = p.output.png_out or p.input.mtz.replace('.mtz', name+'.png')
    fig.savefig(name)

  else:
    fig = plot_layer(layer, mask, low, high, cmap, ang, asp,
                     clip, cnt, 1, dstr, inset, ext, loc, label)
    plt.show()


def axes_to_data(ax, xy):

    return ax.transData.inverted().transform(ax.transAxes.transform(xy))


def data_to_axes(ax, xy):

    return ax.transAxes.inverted().transform(ax.transData.transform(xy))


def plot_section(ax, layer, mask, vmin, vmax, asp, cmap, cnt, ang, clip, ext):

  # Dummy plot
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
  if clip:
    x,y,z = clip
    if z <= 1:
      x,y = axes_to_data(ax, (x,y))
      z  *= max(ext[1]-ext[0], ext[3]-ext[2])
    w,h = min(z,z*asp), min(z,z/asp)
    ax.set_xlim((x-w, x+w))
    ax.set_ylim((y-h, y+h))
  else:
    j, i = layer.shape
    xcrd = np.linspace(ext[0]+0.5, ext[1]-0.5, i)
    ycrd = np.linspace(ext[2]+0.5, ext[3]-0.5, j)
    mask = ~np.isnan(layer.ravel())
    mesh = np.array(np.meshgrid(xcrd, ycrd)).reshape(2,-1).T[mask]
    axxy = data_to_axes(ax, mesh)
    x0,_ = axes_to_data(ax, (axxy[:,0].min(), 0.5))
    x1,_ = axes_to_data(ax, (axxy[:,0].max(), 0.5))
    x0,x1= x0 - 0.5, x1 + 0.5
    if clip is None: x0, x1 = min(x0, ext[0]), max(x1, ext[1])
    ax.set_xlim((x0, x1))

  return im


def plot_layer(layer, mask, vmin, vmax, cmap, ang, asp,
               clip, cnt, scale, dstr, inset, ext, loc, label):

  plt.close()
  plt.rc('font', size=5. * scale)
  figsize = 8 * scale, 6 * scale
  fig, (ax, dx, cax, pad) = plt.subplots(1,4,figsize=figsize,
    gridspec_kw={'width_ratios':[0.80, 0.04, 0.02, 0.14], 'wspace':0, 'hspace':0})
  ax.set_axis_off()
  dx.set_axis_off()
  pad.set_axis_off()
  for Axes in fig.axes:
    Axes.tick_params(left=False,right=False,top=False,bottom=False,labelleft=False,
                     labelright=False,labeltop=False,labelbottom=False)
  fig.tight_layout()
  for Axes in (dx, cax):
    pos = Axes.get_position()
    Axes.set_position([pos.x0, pos.y0, pos.width, pos.height*0.95])

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

  # Plot data
  im = plot_section(ax, layer, mask, vmin, vmax, asp, cmap, cnt, ang, clip, ext)

  # Plot inset
  if inset is not None:
    x1, y1 = ax.transAxes.transform((0.5,0.5))
    x0, y0 = ax.transAxes.transform((0,0))
    w, h   = fig.transFigure.inverted().transform([min(x1-x0, y1-y0)] * 2)
    x1, y1 = fig.transFigure.inverted().transform((x1,y1))
    x      = {'l':x1-w,'r':x1}.get(loc[1])
    y      = {'b':y1-h,'t':y1}.get(loc[0])
    ins    = fig.add_axes([x, y, w, h])
    ins.tick_params(left=False,right=False,top=False,bottom=False,labelleft=False,
                    labelright=False,labeltop=False,labelbottom=False)
    plot_section(ins, layer, mask, vmin, vmax, asp, cmap, cnt, ang, inset, ext)
    # Draw rectangle
    corners = [axes_to_data(ins, xy) for xy in [(0,0),(0,1),(1,1),(1,0),(0,0)]]
    ax.plot(*zip(*corners), color='black', lw=0.5 * scale)

  # Plot label
  if label:
    lbl = fig.add_axes([0.7,0.025,0,0])
    lbl.set_axis_off()
    lbl.text(0,1,label,fontsize=22. * scale)

  # Plot colorbar
  scf = ticker.ScalarFormatter()
  scf.set_powerlimits((-2,2))
  scf.set_useMathText(True)
  bar = fig.colorbar(im, cax=cax, format=scf)
  bar.ax.tick_params(width=2. * scale, size=8. * scale, labelsize=22. * scale)
  bar.ax.yaxis.get_offset_text().set(size=22. * scale)
  bar.ax.yaxis.get_offset_text().set_position((0, 0))

  return fig
