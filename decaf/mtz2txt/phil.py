import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  mtz2txt
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
        resolution = '9e99 0'
          .type = floats
          .help = 'Resolution range (or max if one value provided)'
      }
      output
        .help = "Output files"
      {
        txt_out = None
          .type = path
          .help = 'Txt output file'
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
