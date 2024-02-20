import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  btrace
    {
      input
        .help = "Input files"
      {
        input = None
          .type = path
          .help = 'PDB or JSON input files'
          .multiple = True
      }
      params
        .help = "Parameters"
      {
        lines = None
          .type = floats
          .help = 'Locations for vertical lines'
        show = True
          .type = bool
          .help = 'Show plot'
        legend = False
          .type = bool
          .help = 'Plot legend'
        cmap = Set1
          .type = str
          .help = 'Colormap used for the color cycler'
        min = None
          .type = float
          .help = 'Set min y-value'
        max = None
          .type = float
          .help = 'Set max y-value'
        chain = 0
          .type = ints
          .help = 'Chains to plot'
      }
      output
        .help = "Output files"
      {
        png_out = btrace.png
          .type = path
          .help = 'Output plot png'
      }
    }
    """)

  for i, arg in enumerate(args):
    if '=' not in arg:
      args[i] = 'input='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
