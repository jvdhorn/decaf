import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  mtzstats
    {
      input
        .help = "Input files"
      {
        mtz = None
          .type = path
          .help = 'MTZ input file'
        lbl = IDFF
          .type = str
          .help = 'Desired array label'
        cutoff = -9e99
          .type = float
          .help = 'Only consider intensities larger than this value'
        resolution = 0
          .type = floats
          .help = 'Resolution limit (or range if 2 values are provided)'
      }
      params
        .help = "Running parameters"
      {
        sg = None
          .type = str
          .help = 'Space group symbol or number'
        size = 3
          .type = ints
          .help = 'Kernel size for determining number of local maxima'
        floats = False
          .type = bool
          .help = 'Show min, max, mean, median, std and var as floats'
        bins = 1000
          .type = int
          .help = 'Number of bins for calculating Shannon entropy'
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
