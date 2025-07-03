import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  pattsize
    {
      input
        .help = "Input files"
      {
        map = None
          .type = path
          .help = 'Map input file'
        binsize = 1.
          .type = float
          .help = 'Spacing between shells'
        sigma = 2.
          .type = float
          .help = 'Sigma cutoff used to determine object size'
        full = False
          .type = bool
          .help = 'Analyze the full map, not just up to the cutoff (can take minutes)'
        show = False
          .type = bool
          .help = 'Show plot'
        scale = *linear log
          .type = choice
          .help = 'Scale y-axis log or linear'
      }
    }
    """)

  for i, arg in enumerate(args):
    if '=' not in arg:
      args[i] = 'map='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
