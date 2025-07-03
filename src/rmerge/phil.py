import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  rmerge
    {
      input
        .help = "Input files"
      {
        mtz_1 = None
          .type = path
          .help = 'Input mtz'
      }
      params
        .help = "Running parameters"
      {
        array = IDFF
          .type = str
          .help = 'Array of interest'
        space_group = None
          .type = str
          .help = 'Space group symbol'
      }
    }
    """)

  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil.extract()
