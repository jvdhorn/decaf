import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  filter_mtz
    {
      input
        .help = "Input files"
      {
        mtz_1 = None
          .type = path
          .help = 'Mtz input file'
        lbl_1 = IDFF
          .type = str
          .help = 'Array label'
      }
      params
        .help = "Control running"
      {
        size = 1.0
          .type = floats
          .help = 'Sigma for Gaussian, kernel for uniform; can be anisotropic'
        filter = *gaussian uniform
          .type = choice
          .help = 'Filter type'
      }
      output
        .help = "Output files"
      {
        mtz_out = None
          .type = str
          .help = 'Mtz output file'
      }
    }
    """)

  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil.extract()
