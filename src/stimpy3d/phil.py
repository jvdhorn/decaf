import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  stimpy3d
    {
      input
        .help = "Input files"
      {
        mtz_1 = None
          .type = path
          .help = 'Input file'
        lbl_1 = IDFF
          .type = str
          .help = 'Array label'
        bins = 1
          .type = int
          .help = 'Number of resolution bins'
        N = 1
          .type = int
          .help = 'Symmetry multiplicity'
        smooth_kernel = 3
          .type = int
          .help = '1D kernel for smoothing the binned background'
        sigma_cutoff = 5.
          .type = float
          .help = 'Discard intensities outside this sigma range'
        cutoff = -9e99
          .type = float
          .help = 'Only consider intensities larger than this value'
        use_fit_params = False
          .type = bool
          .help = 'Use fitted moments to determine distribution stats'
        write = True
          .type = bool
          .help = 'Write additional information such as raw distributions'
        directory = None
          .type = path
          .help = 'Directory to write all files to'
        remove_excess = False
          .type = bool
          .help = 'Remove data outside of sigma filter'
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
