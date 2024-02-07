import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  qdep
    {
      input
        .help = "Input files"
      {
        mtz = None
          .type = path
          .help = 'Mtz input file'
        lbl = IDFF
          .type = str
          .help = 'Array of interest'
      }
      params
        .help = "Control running"
      {
        sc_size = None
          .type = ints
          .help = 'Supercell size'
          .optional = False
        max = None
          .type = ints
          .help = 'Maximum values to consider for h, k and l'
        min = None
          .type = ints
          .help = 'Minimum values to consider for h, k and l'
        range = -1 -1 -1
          .type = ints
          .help = 'Extent of the analysis in each dimension. Default is half of the supercell.'
        resolution = 0
          .type = floats
          .help = 'Resolution cutoff (or range if two values are given)'
        strong = 0
          .type = int
          .help = 'Consider this number of strongest reflections (negative for weakest)'
        plot = False
          .type = bool
          .help = 'To plot or not to plot?'
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
