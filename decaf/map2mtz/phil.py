import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  map2mtz
    {
      input
        .help = "Input files"
      {
        map = None
          .type = path
          .help = 'Map input file'
        resolution = 0
          .type = floats
          .help = 'Resolution range (or max if one value provided)'
        cutoff = None
          .type = float
          .help = 'Cutoff value for valid intensities'
        center = False
          .type = bool
          .help = 'Put corner on origin'
      }
      output
        .help = "Output files"
      {
        lbl = 'IDFF'
          .type = str
          .help = 'Output file array'
        mtz_out = None
          .type = path
          .help = 'Mtz output file'
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
