import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  mtz2map
    {
      input
        .help = "Input files"
      {
        mtz = None
          .type = path
          .help = 'Mtz input file'
        lbl = 'IDFF'
          .type = str
          .help = 'Input file array'
        resolution = 0
          .type = floats
          .help = 'Resolution range (or max if one value provided)'
        fill = None
          .type = float
          .help = 'Fill value'
      }
      output
        .help = "Output files"
      {
        map_out = None
          .type = path
          .help = 'Map output file'
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
