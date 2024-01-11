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
      }
      params
        .help = "Running parameters"
      {
        sg = None
          .type = str
          .help = 'Space group symbol or number'
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
