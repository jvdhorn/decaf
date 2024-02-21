import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  slicemtz
    {
      input
        .help = "Input files"
      {
        mtz = None
          .type = path
          .help = 'Mtz or map input file'
        lbl = IDFF
          .type = str
          .help = 'Array label'
        cutoff = None
          .type = float
          .help = 'If input is map, only include values higher than this one'
      }
      params
        .help = "Control running"
      {
        slice = hk0
          .type = str
          .help = 'Slice to show'
        depth = 0
          .type = int
          .help = 'Additional thickness on both sides of the slice'
        multiply = 1.0
          .type = float
          .help = 'Multiply all values by this amount'
        add = 0.0
          .type = float
          .help = 'Offset all values by this amount (after multiplying)'
        mode = *sum mean median min max prod std var
          .type = choice
          .help = 'How to treat the accumulation of voxels in the depth-direction'
        min = None
          .type = float
          .help = 'Lower limit for intensity scale'
        max = None
          .type = float
          .help = 'Upper limit for intensity scale'
        resolution = 0
          .type = floats
          .help = 'Resolution limit'
        log_intensities = False
          .type = bool
          .help = 'Apply logarithmic scaling'
        axes = True
          .type = bool
          .help = 'Overlay axes'
        trim = True
          .type = bool
          .help = 'Trim edges to fill figure'
        contours = False
          .type = bool
          .help = 'Overlay zero-level contours'
        distribution = *Lin Log False
          .type = choice
          .help = 'Plot intensity distribution'
        sc_size = -1 -1 -1
          .type = ints
          .help = 'Supercell size; draw Bragg positions when given'
        zoom = None
          .type = floats
          .help = 'Zoom around index: "x, y, [number of voxels to pad]"'
        dotsize = 1
          .type = int
          .help = 'Width of the axes and Bragg position indicators (higher = thinner)'
        filter_size = 0
          .type = int
          .help = 'Before plotting, apply Gaussian filter with this sigma'
        fmode = *gaussian gausslap uniform min max median
          .type = choice
          .help = 'How to treat the accumulation of voxels in the depth-direction'
        positive_cmap = BuPu
          .type = str
          .help = 'Matplotlib colormap for positive values'
        negative_cmap = YlOrRd
          .type = str
          .help = 'Matplotlib colormap for negative values'
        save = False
          .type = bool
          .help = 'Save plot as PNG instead of showing'
      }
    }
    """)

  for i, arg in enumerate(args):
    if '=' not in arg:
      args[i] = 'mtz='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
