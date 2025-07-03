import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  patterson
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
        resolution = 0
          .type = floats
          .help = 'Resolution limit (or range if 2 values are provided)'
        fill = 0
          .type = float
          .help = 'Fill value for missing intensities'
      }
      params
        .help = "Input files"
      {
        bins = None
          .type = int
          .help = 'If provided, subtract mean value in this many intensity bins'
        sample = 3.0
          .type = float
          .help = 'Sampling rate in 1/A - higher is finer'
        use_intensities = False
          .type = bool
          .help = 'Patterson from intensities rather than F'
        center = True
          .type = bool
          .help = 'Place origin in center of the map'
        limit = None
          .type = float
          .help = 'Limit dimensions in Angstroms'
      }
      output
        .help = "Input files"
      {
        map_out = None
          .type = path
          .help = 'Patterson output file'
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
