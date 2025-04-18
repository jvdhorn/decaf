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
        power = 1.0
          .type = float
          .help = 'Raise all values to this power'
        multiply = 1.0
          .type = float
          .help = 'Multiply all values by this amount'
        add = 0.0
          .type = float
          .help = 'Offset all values by this amount (after multiplying)'
        mode = *sum mean median min max prod std var solid
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
        projection = gp eq mercator
          .type = choice
          .help = 'Show a projection at the highest resolution'
        width = 721
          .type = int
          .help = 'Number of horizontal sampling points used in prjection mode'
        autoscale = True
          .type = bool
          .help = 'Increase visibility on very skewed data'
        sigma = 6.0
          .type = float
          .help = 'Autoscale to this many standard deviations'
        log = False
          .type = bool
          .help = 'Apply logarithmic scaling'
        axes = True
          .type = bool
          .help = 'Overlay axes'
        trim = True
          .type = bool
          .help = 'Trim edges to fill figure'
        center = False
          .type = bool
          .help = 'Put corner on origin'
        contours = 0
          .type = int
          .help = 'Plot (this many) contours instead of grid, negative for unfilled'
        distribution = *Lin Log False
          .type = choice
          .help = 'Plot intensity distribution'
        sc_size = None
          .type = ints
          .help = 'Supercell size; draw Bragg positions when given'
        zoom = None
          .type = floats
          .help = 'Zoom around index: "x, y, [number of voxels to pad]", can be fractional'
        inset = None
          .type = floats
          .help = 'Create inset with zoom around index, using zoom syntax'
        position = *bl br tl tr
          .type = choice
          .help = 'Corner in which to place the inset'
        dotsize = 1
          .type = int
          .help = 'Width of the axes and Bragg position indicators (higher = thinner)'
        filter = 0
          .type = int
          .help = 'Before plotting, apply Gaussian filter with this sigma'
        fmode = *gaussian gausslap uniform min max median
          .type = choice
          .help = 'Filter mode'
        positive_cmap = BuPu
          .type = str
          .help = 'Matplotlib colormap for positive values'
        negative_cmap = YlOrRd
          .type = str
          .help = 'Matplotlib colormap for negative values'
        normalize_cmap = True
          .type = bool
          .help = 'Match scales of positive and negative cmaps'
        overlay = black
          .type = str
          .help = 'Overlay colour'
        label = True
          .type = bool
          .help = 'Add slice label to plot'
        save = False
          .type = bool
          .help = 'Save plot as PNG instead of showing'
      }
      output
        .help = "Output files"
      {
        png_out = None
          .type = path
          .help = 'PNG output file'
        figscale = 6.0
          .type = float
          .help = 'Scalar for the output image file'
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
