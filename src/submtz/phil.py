import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  submtz
    {
      input
        .help = "Input files"
      {
        mtz_1 = None
          .type = path
          .help = 'First file'
        lbl_1 = IDFF
          .type = str
          .help = 'First array label'
        mtz_2 = None
          .type = path
          .help = 'Second file'
        lbl_2 = IDFF
          .type = str
          .help = 'First array label'
      }
      params
        .help = "Running parameters"
      {
        cutoff = -9e99
          .type = float
          .help = 'Only consider intensities larger than this value'
        add = False
          .type = bool
          .help = 'Add arrays instead of subtracting'
        scale = 1
          .type = float
          .help = 'Scale factor for second map (0 = auto scale)'
      }
    }
    """)

  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil.extract()
