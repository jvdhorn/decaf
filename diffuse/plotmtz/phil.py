import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  plotmtz
    {
      input
        .help = "Input files"
      {
        mtz = None
          .type = path
          .help = 'First file'
          .multiple = True
        lbl = IDFF
          .type = str
          .help = 'First array label'
        bins = 100
          .type = int
          .help = 'Number of bins'
        show = False
          .type = bool
          .help = 'Show plot or not'
      }
      params
        .help = "Running parameters"
      {
        resolution = "9e99 0"
          .type = floats
          .help = 'Resolution cutoff (provide 2 values (in any order) for range)'
        log = True
          .type = bool
          .help = 'Logarithmic y-axis'
        legend = True
          .type = bool
          .help = 'Show legend'
        cutoff = -9e99
          .type = float
          .help = 'Only consider intensities larger than this value'
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
