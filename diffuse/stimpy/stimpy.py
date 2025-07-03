from __future__ import division, print_function
from scitbx.array_family import flex
from scipy import ndimage, stats
from dxtbx.format.FormatCBFMini import FormatCBFMini as formatter
from matplotlib import pyplot as plt, colors
import os
import numpy as np
import dxtbx
from . import phil

np.random.seed(1337)
colours  = np.concatenate(( ((1,1,1),), np.random.rand(9999, 3) ))
rnd_cmap = colors.LinearSegmentedColormap.from_list('rnd_cmap', colours, N=10000)

def extract_data(frame, scale=None):

  data = frame.get_raw_data().as_numpy_array()

  detector = frame.get_detectorbase()
  
  if hasattr(detector, 'header') and 'Added=' in detector.header:
    factor  = int(float(detector.header.split('Added=')[1].split(';')[0].strip()))
    data   -= factor

  if scale is not None: data = data/scale

  return data

def identify_regions(array, diagonal=True, interpolate=True, threshold=None, min_group_size=100):
  if threshold is None:
    threshold = array.min()
  bad         = (array <= threshold)
  
  # Interpolate panel borders to extend contiguous regions
  if interpolate:
    array = interpolate_panels(array, threshold=threshold)
  
  kernels = [
   [[0,1,0],
    [1,1,1],
    [0,1,0]],
   [[1,1,1],
    [1,1,1],
    [1,1,1]],
  ]
  structure = kernels[int(diagonal)]
  regions   = np.zeros_like(array)
  counter   = 0
  values    = np.unique(array.ravel())
  for value in values:
    mask          = (array == value)
    lbl, n        = ndimage.label(mask, structure=structure)
    regions[mask] = (lbl + counter)[mask]
    counter      += n
  regions[bad]    = threshold

  # Reassign pixels from groups which are too small
  values, counts  = np.unique(regions, return_counts=True)
  too_small_list  = values[counts < min_group_size]
  too_small_mask  = np.isin(regions, too_small_list)
  too_small_mask |= bad
  ind             = ndimage.distance_transform_edt(too_small_mask,
                                                   return_distances = False,
                                                   return_indices   = True)
  regions         = regions[tuple(ind)]
  regions[bad]    = threshold

  return regions

def plot_regions(array, show=False, filename='regions.png'):
  
  plt.figure()
  plt.imshow(array, cmap=rnd_cmap)
  plt.savefig(filename, dpi=200)

def region_statistics(array):
  
  regions, counts = np.unique(array.ravel(), return_counts=True)
  counts_sorted   = sorted(zip(counts, regions))[::-1]
  for count, region in counts_sorted:
      print('Region {:5.0f}: {:8.0f} pixels'.format(region, count))
  print('Total number of regions: {}'.format(len(regions)))

def interpolate_panels(array, threshold=-1):

  mask         = (array <= threshold)
  ind          = ndimage.distance_transform_edt(mask,
                                                return_distances = False,
                                                return_indices   = True)
  interpolated = array[tuple(ind)]

  return interpolated

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

def chapman_stats(N, mean, var, skew):

  if skew < 0: skew = 0
  Sigma = (N**2 * skew/2.)**(1./3) * var**(1./2)
  Var   = var * (1-(N/4.)**(1./3) * skew**(2./3))
  Mu    = mean - Sigma
  return (Sigma, Var, Mu)

def chapman_stats_discrete(N, mean, var):

  Sigma = (N * max(var - mean, 0)) ** (1./2)
  Mu    = mean - Sigma
  return (Sigma, Mu)

def image_region_statistics(array, regions, threshold=None, save_img=True, write=True, template=None, directory='stimpy', pref='', sigma=3, N=8, discrete=False, use_fit_params=False, remove_excess=False):

  if threshold is None:
    threshold = array.min()
  bad         = (array <= threshold)

  n_regions = np.unique(regions.ravel())

  background = np.zeros_like(array)
  signal     = np.zeros_like(array)

  for region in n_regions:
    mask = (regions == region) & ~bad
    data = array[mask]
    data = filter_sigma(data, sigma) # Use only data within a certain sigma range
    excess = (array[mask] < data.min()) ^ (array[mask] > data.max()) if data.size else data
    if remove_excess: bad[mask] |= excess
    n, mean, var, skew = data.size, data.mean(), data.var(), stats.skew(data)
    print('Region {:.0f}: N={:.0f}; mean={:.2f}; var={:.2f}; skew={:.2f}'.format(region, n, mean, var, skew))

    if n and (write or use_fit_params):
      fit      = stats.skewnorm.fit(array[mask])
      skewnorm = stats.skewnorm(*fit)
      fitmean, fitvar, fitskew = skewnorm.stats(moments='mvs')
      print('Fit stats: mean={:.2f}; var={:.2f}; skew={:.2f}'.format(fitmean, fitvar, fitskew))

      plt.figure()
      hist_data = plt.hist(data, bins=20,)
      y         = hist_data[0]
      x_edges   = hist_data[1]
      x         = (x_edges[:-1] + x_edges[1:]) / 2
      pdf       = skewnorm.pdf(x)
      plt.plot(x, pdf * mask.sum() * (x[1] - x[0]))
      plt.savefig('{}/region_{:.0f}_statistics.png'.format(directory, region), dpi=200)
      plt.close()

      with open('{}/region_{:.0f}_histogram.x_y'.format(directory, region), 'w') as x_y:
        for x_val, y_val in zip(x,y): x_y.write('{} {}\n'.format(x_val, y_val))
    
      with open('{}/region_{:.0f}_counts.txt'.format(directory, region), 'w') as txt:
        contents = '\n'.join(map(str, array[mask]))
        txt.write(contents)

    else:
      fitmean, fitvar, fitskew = 0, 0, 0

    if use_fit_params:
      mean, var, skew = fitmean, fitvar, fitskew
   
    if not discrete:
      Sigma, Var, Mu = chapman_stats(N, mean, var, skew)
      print('Chapman stats: Sigma={:.2f}; Var={:.2f}; Mu={:.2f}'.format(Sigma, Var, Mu))
      print()

    else:
      Sigma, Mu = chapman_stats_discrete(N, mean, var)
      print('Chapman stats: Sigma={:.2f}; Mu={:.2f}'.format(Sigma, Mu))
      print()

    if save_img and template is not None:
      new_array       = np.zeros_like(array)
      new_array[mask] = array[mask]
      write_image(array=new_array, template=template, filename='{}/{}_region_{:.0f}.cbf'.format(directory, pref, region))

    background[(regions == region)] = max(0, Mu)
    signal[(regions == region)] = Sigma

  background[np.isnan(background)] = 0
  
  return background, signal, bad

def remove_bragg(array, polar, dilation=7, median_filter_size=9, value=None):

  array          = array.copy()
  if value is None: value = min(array.min(), 0)
  median         = ndimage.median_filter(array, median_filter_size)
  mask           = array - median > polar
  dilated        = ndimage.maximum_filter(mask, dilation)
  array[dilated] = value

  return array, dilated, value

def smooth(array, kernel=1):

  array = ndimage.uniform_filter(array, kernel)
  return array

def write_image(array, template, filename='image.img', encoding='int16'):

  array = array.copy()

  if filename.endswith('.cbf'):
    formatter.as_file(
      detector = template.get_detector(),
      beam     = template.get_beam(),
      gonio    = template.get_goniometer(),
      scan     = template.get_scan(),
      data     = flex.int32(array.astype('int32')),
      path     = filename,
    )

  elif filename.endswith('.img'):
    addfactor = abs(min(0, array.min()))
    array    += addfactor
    head_str  = 'Added={};'.format(addfactor)
    header    = template.get_detectorbase().header
    headlines = header.splitlines()
    if 'Added' in header:
      n_line            = next(n for n, line in enumerate(headlines) if 'Added' in line)
      headlines[n_line] = head_str
    else:
      n_line            = max(n for n, line in enumerate(headlines) if '=' in line)
      headlines.insert(n_line+1, head_str)
    header    = '\n'.join(headlines)
    head_len  = int(header.split('HEADER_BYTES=')[1].split(';')[0].strip())
    header    = '{:{}}'.format(header, head_len)[:head_len]

    if   encoding == 'int16': array = np.int16(array)
    elif encoding == 'int32': array = np.int32(array)
    elif encoding == 'int64': array = np.int64(array)

    with open(filename, 'w') as writer:
      writer.write(header)
      writer.write('')
      writer.write(array)

def main(args, log):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().stimpy
  pref = p.input.image.split('/')[-1].split('.')[0]
  path = p.output.directory
  if path == 'stimpy':
    path += '_' + pref

  try:
    os.mkdir(path)
  except OSError:
    pass

  if not p.output.show:
    plt.switch_backend('agg')

  polar_scale  = p.params.polar_scale
  image_scale  = p.params.image_scale
  bin_photons  = p.params.bin_photons
  
  # Prepare polar
  polar        = dxtbx.load(p.input.polar)
  polar_data   = extract_data(polar, scale=polar_scale)
  
  # Bin polar data and plot statistics
  polar_binned = polar_data // bin_photons
  region_statistics(polar_binned)
  plot_regions(polar_binned, show=p.output.show, filename='{}/{}_values.png'.format(path, pref))

  # Convert polar data into contiguous regions and plot statistics
  regions      = identify_regions(array          = polar_binned,
                                  diagonal       = p.params.diagonal,
                                  interpolate    = p.params.interpolate_panels,
                                  threshold      = p.params.threshold,
                                  min_group_size = p.params.min_group_size)
  region_statistics(regions)
  plot_regions(regions, show=p.output.show, filename='{}/{}_regions.png'.format(path, pref))

  # Show plots if desired
  if p.output.show: plt.show()
  
  # Prepare image
  image        = dxtbx.load(p.input.image)
  image_data   = extract_data(image, scale=image_scale)
  nobragg, bmask, bfill = remove_bragg(image_data,
                                polar_data,
                                dilation           = p.params.bragg_mask_dilation,
                                median_filter_size = p.params.bragg_mask_median_filter)

  # Determine background
  background, signal, badmask = image_region_statistics(nobragg,
                                         regions,
                                         threshold      = p.params.threshold,
                                         write          = p.output.write,
                                         save_img       = p.output.save_region_images,
                                         template       = image,
                                         directory      = path,
                                         pref           = pref,
                                         N              = p.params.N,
                                         sigma          = p.params.sigma,
                                         remove_excess  = p.params.remove_excess,
                                         discrete       = p.params.discrete,
                                         use_fit_params = p.params.use_fit_params)
  background   = smooth(background, p.params.smooth_background)
  corrected    = image_data - background
  badmask     ^= bmask

  if p.params.remove_bragg:
    corrected[bmask] = p.params.fill
  
  corrected[badmask] = p.params.fill
  
  if '.img' in p.input.image:
    write_image(array    = (signal),
                template = image,
                filename = '{}/{}_signal.img'.format(path, pref), encoding = p.output.encoding)
    write_image(array    = (background),
                template = image,
                filename = '{}/{}_background.img'.format(path, pref), encoding = p.output.encoding)
    write_image(array    = (corrected),
                template = image,
                filename = '{}/{}_corrected.img'.format(path, pref), encoding = p.output.encoding)
  write_image(array    = (signal),
              template = image,
              filename = '{}/{}_signal.cbf'.format(path, pref))
  write_image(array    = (background),
              template = image,
              filename = '{}/{}_background.cbf'.format(path, pref))
  write_image(array    = (corrected),
              template = image,
              filename = '{}/{}_corrected.cbf'.format(path, pref))

def run(args):

  with open('log.txt', 'w') as log:
    main(args = args, log = log)

