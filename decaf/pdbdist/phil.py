import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  pdbdist
    {
      input
        .help = "Input files"
      {
        pdb = None
          .type = path
          .help = 'PDB input file'
      }
      params
        .help = "Parameters"
      {
        show = True
          .type = bool
          .help = 'Show plot'
        cmap = inferno
          .type = str
          .help = 'Colormap'
        limit_ca = None
          .type = int
          .help = 'C-alpha cutoff'
        chain = 0
          .type = int
          .help = 'Select chain to plot'
        min = None
          .type = float
          .help = 'Min value to plot'
        max = None
          .type = float
          .help = 'Max value to plot'
      }
      output
        .help = "Output files"
      {
        png_out = pdbdist.png
          .type = path
          .help = 'Output plot png'
      }
    }
    """)

  for i, arg in enumerate(args):
    if '=' not in arg:
      args[i] = 'pdb='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
