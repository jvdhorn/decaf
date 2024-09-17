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
          .multiple = True
      }
      params
        .help = "Parameters"
      {
        show = True
          .type = bool
          .help = 'Show plot'
        positive_cmap = inferno
          .type = str
          .help = 'Matplotlib colormap for positive values'
        negative_cmap = gist_earth
          .type = str
          .help = 'Matplotlib colormap for negative values'
        residue_range = None
          .type = ints
          .help = 'Min and max residue'
        states = None
          .type = ints
          .help = 'Min and max states'
        chain = 0
          .type = int
          .help = 'Select chain to plot'
        min = None
          .type = float
          .help = 'Min value to plot'
        max = None
          .type = float
          .help = 'Max value to plot'
        mode = *cov cc std var mean
          .type = choice
          .help = 'Ensemble treatment'
        combine = *both sub div add mul
          .type = choice
          .help = 'Multiple input treatment'
        full = True
          .type = bool
          .help = 'Show full matrix if possible'
        subtract = None
          .type = path
          .help = 'Subtract intramolecular distances from this file'
        lines = None
          .type = floats
          .help = 'Draw horizontal and vertical lines at these positions'
        lc = cyan
          .type = str
          .help = 'Line colors'
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
