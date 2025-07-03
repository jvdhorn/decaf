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
          .help = 'First file; for R1, this should be the experimental data'
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
          .help = 'Number of resolution shells'
        subtract = min mean
          .type = choice
          .help = 'Correction for every resolution shell'
        cutoff = -9e99
          .type = float
          .help = 'Only consider intensities larger than this value'
        plot = False
          .type = bool
          .help = 'Create plot'
        show = False
          .type = bool
          .help = 'Show plot'
      }
    }
    """)

  defaults = ['mtz_1','mtz_2']
  for i, arg in enumerate(args):
    if '=' not in arg:
      args[i] = defaults.pop(0)+'='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
