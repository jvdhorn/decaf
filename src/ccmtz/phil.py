import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  ccmtz
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
        hlim = '-inf +inf'
          .type = str
          .help = 'Limit reflections in h (inclusive)'
        klim = '-inf +inf'
          .type = str
          .help = 'Limit reflections in k (inclusive)'
        llim = '-inf +inf'
          .type = str
          .help = 'Limit reflections in l (inclusive)'
        bins = 10
          .type = int
          .help = 'Number of resolution bins'
        cutoff = -1000
          .type = float
          .help = 'Only consider intensities larger than this value'
        show = False
          .type = bool
          .help = 'Show plot'
      }
    }
    """)

  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil.extract()
