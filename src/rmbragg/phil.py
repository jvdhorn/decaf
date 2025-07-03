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
        mtz_1 = None
          .type = path
          .help = 'Mtz input file'
      }
      params
        .help = "Control running"
      {
        sc_size = 5 5 10
          .type = strings
          .help = 'Supercell size'
        box = 1 1 1
          .type = strings
          .help = 'Box size, odd numbers only, 1 1 1 for just Bragg'
        keep = False
          .type = bool
          .help = 'Setting this to True inverts the selection, keeping the Bragg peaks'
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
      args[i] = 'mtz_1='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
