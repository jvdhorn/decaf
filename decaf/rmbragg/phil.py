import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  rmbragg
    {
      input
        .help = "Input files"
      {
        mtz = None
          .type = path
          .help = 'Mtz input file'
        lbl = IDFF
          .type = str
          .help = 'Array of interest (used for intensity fraction)'
      }
      params
        .help = "Control running"
      {
        sc_size = None
          .type = ints
          .help = 'Supercell size'
        box = 1 1 1
          .type = ints
          .help = 'Box size, odd numbers only, 1 1 1 for just Bragg'
        hlim = '-inf +inf'
          .type = str
          .help = 'Limit reflections in h (inclusive)'
        klim = '-inf +inf'
          .type = str
          .help = 'Limit reflections in k (inclusive)'
        llim = '-inf +inf'
          .type = str
          .help = 'Limit reflections in l (inclusive)'
        fraction = None
          .type = float
          .help = 'Remove this fraction of highest intensities within box or whole map'
        keep = False
          .type = bool
          .help = 'Invert selection'
      }
      output
        .help = "Output files"
      {
        mtz_out = None
          .type = path
          .help = 'Mtz output file'
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
