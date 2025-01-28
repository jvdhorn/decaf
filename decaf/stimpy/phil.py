import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  stimpy
    {
      input
        .help = "Input files"
      {
        image = None
          .type = path
          .help = 'Basic image input file'
        radial = None
          .type = path
          .help = 'Radial average image input file'
      }
      params
        .help = "Control running"
      {
        image_scale = 1.0
          .type = float
          .help = 'If input image has been scaled, it will be divided by this number'
        radial_scale = 1.0
          .type = float
          .help = 'If radial image has been scaled, it will be divided by this number'
        bin_counts = 1.0
          .type = float
          .help = 'Intensity counts within this range are considered to be the same in the radial image'
        min_group_size = 1000
          .type = int
          .help = 'Contiguous regions must have this amount of pixels, otherwise they will be redistributed'
        interpolate_panels = True
          .type = bool
          .help = 'Apply simple NN interpolation to identify contiguous regions over panel borders'
        diagonal = True
          .type = bool
          .help = 'Include diagonally touching pixels in contigous regions'
        threshold = None
          .type = float
          .help = 'All counts below the threshold are considered bad. If none: threshold is set to lowest value.'
        fill = -1000
          .type = int
          .help = 'Fill Bragg cutouts with this value'
        N = 1.0
          .type = float
          .help = 'Number of symmetry elements in the space group'
        sigma = 9e999
          .type = float
          .help = 'Disregard intensities outside this sigma range'
        remove_excess = False
          .type = bool
          .help = 'Remove data points outside sigma filter'
        remove_bragg = True
          .type = bool
          .help = 'Apply median threshold filter to remove Bragg peaks'
        mode = *continuous discrete dual
          .type = choice
          .help = 'Noisy Wilson distribution treatment'
        use_fit_params = False
          .type = bool
          .help = 'Fit a skewnorm distribution to the observed intensities and use those properties'
        smooth_background = 3
          .type = int
          .help = 'Smooth the region-based background with a kernel of this size'
        median_size = 9
          .type = int
          .help = 'Median filter mask size for identifying peaks'
        dilation_size = 5
          .type = int
          .help = 'Dilation mask size for the removal of Bragg peaks'
      }
      output
        .help = "Output parameters"
      {
        directory = stimpy
          .type = path
          .help = 'Name of directory in which to put results'
        show = False
          .type = bool
          .help = 'Show plots while running'
        encoding = int16
          .type = str
          .help = 'Bitwise encoding for img-files'
        write = False
          .type = bool
          .help = 'Write statistical details and plots'
        save_region_images = False
          .type = bool
          .help = 'Save images for each individual region'
      }
    }
    """)

  defaults = ['image','radial']
  for i, arg in enumerate(args):
    if '=' not in arg:
      args[i] = defaults.pop(0)+'='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
