import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  ensemble2adp
    {
      input
        .help = "Input files"
      {
        pdb = None
          .type = path
          .help = 'PDB input file'
        models = None
          .type = int
          .help = 'Number of models to consider'
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
