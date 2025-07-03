import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  covbydist
    {
      input
        .help = "Input files"
      {
        pdb = None
          .type = path
          .help = 'Multistate PDB input file'
      }
      params
        .help = "Parameters"
      {
        show = True
          .type = bool
          .help = 'Show plot'
        include_neighbours = True
          .type = bool
          .help = 'Include neighbouring molecules in the analysis'
        min_distance = 0.0
          .type = float
          .help = 'Lowest pair distance to consider'
        max_distance = 9e999
          .type = float
          .help = 'Highest pair distance to consider'
        bins = 50
          .type = int
          .help = 'Number of bins in the plot'
      }
      output
        .help = "Output files"
      {
        png_out = covbydist.png
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
