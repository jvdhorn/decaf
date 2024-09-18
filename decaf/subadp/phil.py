import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  subadp
    {
      input
        .help = "Input files"
      {
        pdb_1 = None
          .type = path
          .help = 'First file'
        pdb_2 = None
          .type = path
          .help = 'Second file (can be multiple)'
          .multiple = True
      }
      params
        .help = "Parameters"
      {
        mode = *sub add
          .type = choice
          .help = 'Operator'
      }
    }
    """)

  defaults = ['pdb_1','pdb_2']
  for i, arg in enumerate(args):
    if '=' not in arg:
      defaults.append(defaults[-1])
      args[i] = defaults.pop(0)+'='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
