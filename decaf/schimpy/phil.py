import libtbx.phil

def phil_parse(args=None):
  '''
  Contains default parameters, will process commandline params and 
  changes them
  '''
  # Default parameters
  master_phil = libtbx.phil.parse("""
  schimpy
    {
      input
        .help = "Input files and parameters"
      {
        pdb_in = None
          .type = path
          .help = 'PDB file containing base model'
        tls_in = None
          .type = path
          .help = 'Files containing TLS matrices'
        tls_origin = 0 0 0
          .type = floats
          .help = 'Origin of the TLS-matrices, if not provided elsewhere'
        tls_multipliers = 1 1 1
          .type = floats
          .help = 'Multipliers for input TLS matrices - for enforcing translation or rotation'
        max_level = 1
          .type = int
          .help = 'Highest order TLS layer to analyze, -1 for all'
        min_level = 1
          .type = int
          .help = 'Lowest order TLS layer to analyze'
        reverse_levels = True
          .type = bool
          .help = 'Reverse the order levels in which motions are applied'
        sc_size = 1 1 1
          .type = ints
          .help = 'Supercell size in three directions'
        n_models = 128
          .type = ints
          .help = 'Number of models to generate, multiple possible'
        pick_extremes = None
          .type = floats
          .help = 'Do a pre-simulation, pick this many structures and fill the supercell with them'
        extremes_from_level = 1
          .type = int
          .help = 'After the first level, replace all structures by this many extremes'
        precorrelate_level = None
          .type = int
          .help = 'Perform one correlation round after picking extremes, using this level'
        randomize_after = None
          .type = ints
          .help = 'Swap ASUs randomly after simulating these levels'
        residues = '-1'
          .type = ints
          .help = 'Only randomize TLS groups if they include any of these residues, -1 for all'
        strategy = environments
          .type = str
          .help = '(DEPRECATED) Strategy to use for structure optimization'
        amplitude = 1.0
          .type = float
          .help = 'Additional amplitude for atom shifts'
        reset_b = True
          .type = bool
          .help = 'Set all B-factors to 0'
        remove_waters = True
          .type = bool
          .help = 'Remove all water molecules'
        zero_trace = False
          .type = bool
          .help = 'Set diagonal elements of screw matrix to zero'
        processes = 1
          .type = int
          .help = 'Number of parallel processes'
        seed = 0
          .type = int
          .help = 'Seed offset control (-1 for random)'
        interval = 1.0
          .type = float
          .help = 'Time delay in seconds between consecutive simulation starts'
      }
      environments
        .help = "Environment definition parameters"
      {
        rigid_only = False
          .type = bool
          .help = 'Only include interactions between helices and/or sheets'
        protein_only = True
          .type = bool
          .help = 'Only include interactions between proteins (not solvent or ligands)'
        heavy_only = False
          .type = bool
          .help = 'Only include interactions between non-hydrogen elements'
        hbond_only = False
          .type = bool
          .help = 'Only include interactions between N, O and F'
        cutoff = 3.
          .type = float
          .help = 'Distance cutoff for neighouring chains'
        stretch = 0.25
          .type = float
          .help = 'Relative distance from the contact points to the center of mass of the corresponding group'
      }
      regularization
        .help = "Parameters involving regularization"
      {
        regularize = False
          .type = bool
          .help = 'Regularize model after applying shifts'
        scope = 1
          .type = int
          .help = 'Regularize internally (0) or in environment (1)'
      }
      correlation
        .help = "Correlation parameters"
      {
        correlate = True
          .type = bool
          .help = 'Correlate shifts in supercell'
        correlate_fancy = True
          .type = bool
          .help = '(DEPRECTATED) Correlate shifts in supercell for each level'
        skip_levels = -1
          .type = ints
          .help = 'Skip these levels when applying correlations'
        sequential = 0
          .type = int
          .help = 'Correlate groups within level sequentially (-1 to reverse)'
        insideout = 0
          .type = int
          .help = 'Correlate groups from the protein interior outward (-1 to reverse)'
        use_contacts = True
          .type = bool
          .help = 'Correlate using contact points instead of vectors'
        contact_midpoint = False
          .type = bool
          .help = 'Use midpoint between groups as contact point instead of bond mean'
        contact_level = -1
          .type = int
          .help = 'Use this level to define contact points (-1 for highest)'
        allow_bad_frac = 0.0
          .type = float
          .help = 'Random threshold for allowing bad swaps'
        disallow_good_frac = 0.0
          .type = float
          .help = 'Random threshold for disallowing good swaps'
        swap_frac = 1.0
          .type = float
          .help = 'Prematurely end correlation procedure after this fraction of swaps (negative for random)'
        require_both = False
          .type = bool
          .help = 'Require the energies of both positions to go down, instead of the sum'
        uncorrelate = False
          .type = bool
          .help = 'Invert the result of the swap allowance'
        use_pbc = True
          .type = bool
          .help = 'Use periodic boundary conditions in the correlation procedure'
        energy_percentile = 0
          .type = float
          .help = 'Only swaps involving energies higher than this percentile after randomizing are considered'
        inherit = True
          .type = bool
          .help = 'Inherit contacts from highest order level; gives a finer sampling'
        correlate_after = 1
          .type = int
          .help = '(DEPRECATED) Correlate after randomizing this level (-1 for all)'
        shifts = True
          .type = bool
          .help = 'Correlate translation vectors (if not use_contacts)'
        rotvecs = True
          .type = bool
          .help = 'Correlate rotation vectors (if not use_contacts)'
        weights = 1.5
          .type = float
          .help = 'Weights are set to n_interactions to the power of this number'
        override_weights = None
          .type = floats
          .help = 'Multipliers for final weights'
        force_power = 2.
          .type = float
          .help = 'The spring force scales with the distance to this power'
      }
      params
        .help = "Control running"
      {
        resolution = 2.0
          .type = float
          .help = 'Resolution limit for structure factor calculation'
        k_sol = 0.35
          .type = float
          .help = 'Flat bulk solvent parameter k'
        b_sol = 50.0
          .type = float
          .help = 'Flat bulk solvent parameter b'
        b_cart = 0 0 0 0 0 0
          .type = floats
          .help = 'Flat bulk solvent anistropic scale matrix'
        fft_mem = 0.
          .type = float
          .help = 'Wait for this amount of memory (GB) to be available before starting the fft'
      }
      output
        .help = "output files"
      {
        write_pdb = 1
          .type = int
          .help = 'Write PDB-file containting this many structures (-1 for all)'
        write_mmcif = False
          .type = bool
          .help = 'Write mmCIF file INSTEAD OF pdb'
        write_mtz = True
          .type = bool
          .help = 'Write MTZ file?'
        single_mtz = False
          .type = bool
          .help = 'Write MTZ with phases of a single supercell'
        write_b = True
          .type = bool
          .help = 'Write rmsf B-factors?'
        pdb_out = supercell_out.pdb
          .type = path
          .help = 'Name of output file containing supercell'
        superpose = True
          .type = bool
          .help = 'Superpose structures in the pdb'
        scale_displacements = 25.
          .type = float
          .help = 'Scale factor for visual center of mass displacements'
        mtz_out = supercell_avg.mtz
          .type = path
          .help = 'Name of MTZ output file'
        b_out   = ensemble_b.json
          .type = path
          .help = 'Json file containing rmsf B-factors'
      }
    }
    """)

  defaults = ['pdb_in','tls_in']
  for i, arg in enumerate(args):
    if '=' not in arg:
      args[i] = defaults.pop(0)+'='+arg
  interpreter = master_phil.command_line_argument_interpreter()
  arguments = [interpreter.process(arg) for arg in args]
  working_phil = master_phil.fetch(sources = arguments)

  return working_phil
