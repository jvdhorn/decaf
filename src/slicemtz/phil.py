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
        mtz_1 = None
          .type = path
          .help = 'Mtz input file'
        lbl_1 = IDFF
          .type = str
          .help = 'Array label'
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
        mode = *sum mean median min max prod std var
          .type = choice
          .help = 'How to treat the accumulation of voxels in the depth-direction'
        min = None
          .type = float
          .help = 'Lower limit for intensity scale'
        max = None
          .type = float
          .help = 'Upper limit for intensity scale'
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
        show = False
          .type = bool
          .help = 'Show plot'
      }
    }
    """)

  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil.extract()
