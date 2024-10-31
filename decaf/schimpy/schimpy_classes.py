from __future__ import division, print_function
from scitbx import matrix
from cctbx import crystal, geometry_restraints, xray
from mmtbx.tls import analysis, tools, tls_as_xyz
from mmtbx.command_line.fmodel import fmodel_from_xray_structure_master_params
from mmtbx.utils import fmodel_from_xray_structure
from mmtbx.model import manager as model_manager_base
from mmtbx import refinement, secondary_structure
from mmtbx_tls_ext import apply_tls_shifts
from iotbx.pdb.hierarchy import af_shared_atom, atom as atom_type
from scitbx.array_family import flex
from matplotlib import pyplot as plt, gridspec
import scitbx
import iotbx
import numpy as np
import random
import csv
import json
import multiprocessing as mp
import os

DEG2RAD = np.pi/180.
DEG2RADSQ = DEG2RAD**2

def model_manager(*args, **kwargs):
  
  """
  Wrapper function for the model manager, to account for inconsistent build_grm requirement.
  """

  try:
    return model_manager_base(*args, build_grm=True, **kwargs)
  except:
    return model_manager_base(*args, **kwargs)

def get_file_extension(files):

  return files[0].split('.')[-1]

def n_from_file_name(string):

  return int(filter(str.isdigit, string.split('.')[-2]))

def msd_to_b(msd):

  return 8./3 * np.pi**2 * msd

def b_to_msd(b):

  return b / (8./3 * np.pi**2)

def tls_from_components(components, amp=1):

  t_temp = ['T11', 'T22', 'T33', 'T12', 'T13', 'T23']
  l_temp = ['L11', 'L22', 'L33', 'L12', 'L13', 'L23']
  s_temp = ['S11', 'S12', 'S13', 'S21', 'S22', 'S23', 'S31', 'S32', 'S33']
  t = [float(components[i])*amp for i in t_temp]
  l = [float(components[i])*amp for i in l_temp]
  s = [float(components[i])*amp for i in s_temp]
  return t, l, s

class TLS_Level(dict):

  def __init__(self, *args, **kwargs):

    super(dict, self).__init__(*args, **kwargs)

  def as_remark_3(self, file_name):
    
    with open(file_name, 'w') as pdb_file:
      tools.remark_3_tls(tlsos             = self.values(),
                         selection_strings = [group.selection_string
                                              for group in self.values()],
                         out               = pdb_file
      )

  def as_pymol_colors(self, file_name, pdb_hierarchy, c_only=False, palette='tab10'):

    try:
      cmap = plt.get_cmap(palette)
    except:
      cmap = plt.get_cmap('rainbow')
    mod   = len(cmap.colors) if hasattr(cmap, 'colors') else 256
    cache = pdb_hierarchy.atom_selection_cache()
    atoms = pdb_hierarchy.atoms()
    with open(file_name, 'w') as pymol_script:
      for i, group in enumerate(self.values()):
        color = '0x{:02X}{:02X}{:02X}'.format(*cmap(i%mod, bytes=True))
        selection = cache.selection(group.selection_string 
                                  + ' and element C' * c_only)
        for atom in atoms.select(selection):
           cmd = 'color {}, id {}\n'.format(color, atom.serial)
           pymol_script.write(cmd)

def tls_from_csv(input_files, max_level, min_level, origin, cache, mult, zero_trace, log=None):
  
  tls_hierarchy = {}
  mat_files = {}
  amp_files = {}
  for file_name in input_files:
    n_layer = n_from_file_name(file_name)
    if 'matrices' in file_name:
      mat_files[n_layer] = file_name
    elif 'amplitudes' in file_name:
      amp_files[n_layer] = file_name
    else:
      continue
  for n_layer in mat_files:
    if n_layer >= min_level and (n_layer <= max_level or max_level == -1):
      groups = TLS_Level()
      read_mat = open(mat_files[n_layer], 'r')
      read_amp = open(amp_files[n_layer], 'r')
      mat_csv = csv.reader(read_mat, delimiter=',')
      amp_csv = csv.reader(read_amp, delimiter=',')
      for n_line, (mat_line, amp_line) in enumerate(zip(mat_csv, amp_csv)):
        if n_line==0: mat_labels = mat_line
        else:
          components = dict(zip(mat_labels, mat_line))
          n                = int(mat_line[0])
          amp              = float(amp_line[-1])
          t,l,s            = tls_from_components(components)
          selection_string = components['label']
          selection        = cache.iselection(selection_string)
          try:
            groups[n] = TLS(t                = t,
                            l                = l,
                            s                = s,
                            amp              = amp,
                            mult             = mult,
                            zero_trace       = zero_trace,
                            origin           = origin,
                            selection_string = selection_string,
                            selection        = selection,
                            log              = log)
          except Exception as e:
            print(('Error while analyzing TLS group {} '
                   'in layer {} with selection string "{}": {}'
                   '\nUsing previous group as substitute.')
                   .format(n, n_layer, selection_string, e))
            groups[n] = groups[n-1]
      read_mat.close()
      read_amp.close()
      tls_hierarchy[n_layer] = groups
  return tls_hierarchy

def tls_from_pdb(input_files, max_level, min_level, cache, mult, zero_trace, log=None):

  tls_hierarchy = {}
  n_layer = 1
  for file_name in input_files:
    if n_layer >= min_level and (n_layer <= max_level or max_level == -1):
      groups = TLS_Level()
      pdb_inp = iotbx.pdb.input(file_name = file_name)
      pdb_hierarchy = pdb_inp.construct_hierarchy()
      tls_params = pdb_inp.extract_tls_params(pdb_hierarchy).tls_params
      for j, group in enumerate(tls_params):
        n = j + 1
        selection   = cache.selection(group.selection_string)
        try:
          groups[n] = TLS(t                = group.t,
                          l                = group.l,
                          s                = group.s,
                          origin           = group.origin,
                          mult             = mult,
                          zero_trace       = zero_trace,
                          selection_string = group.selection_string,
                          selection        = selection,
                          log              = log)
        except Exception as e:
          print(('Error while analyzing TLS group {} '
                 'in layer {} with selection string "{}": {}'
                 '\nUsing previous group as substitute.')
                 .format(n, n_layer, group.selection_string, e))
          groups[n] = groups[n-1]
      if groups:
        tls_hierarchy[n_layer] = groups
        n_layer               += 1
  return tls_hierarchy

def tls_from_json(input_file, max_level, min_level, cache, mult, zero_trace, log=None):

  tls_hierarchy = {}
  with open(input_file) as json_in:
    tls_dict = json.load(json_in)
  levels = sorted((level['level_number'],level['level_name'])
                   for level in tls_dict['level_descriptions'])
  data   = tls_dict['tls_level_data']
  for level in data.values():
    n_level = level['level_number']
    if n_level >= min_level and (n_level <= max_level or max_level == -1):
      groups = TLS_Level()
      for group in level['tls_group_data']:
        try:
          n         = group['group_number']
          modes     = group['tls_modes'][0]
          sele_str  = str(group['selection'])
          selection = cache.selection(sele_str)
          groups[n] = TLS(t                = modes['T'],
                          l                = modes['L'],
                          s                = modes['S'],
                          amp              = list(modes['amplitudes'].values())[0],
                          mult             = mult,
                          zero_trace       = zero_trace,
                          origin           = list(group['tls_origins'].values())[0],
                          selection_string = sele_str,
                          selection        = selection,
                          log              = log)
        except Exception as e:
          print(('Error while analyzing TLS group {} '
                 'in layer {} with selection string "{}": {}'
                 '\nUsing previous group as substitute.')
                 .format(n, n_level, group['selection'], e))
          groups[n] = groups[n-1]
      tls_hierarchy[n_level] = groups
  return tls_hierarchy

def write_b_factors_to_json(sequence, b_factors, file_name='ensemble_b.json'):
  
  label = file_name.replace('.json','')
  dump  = {'x'      : sequence,
           'y'      : b_factors,
           'label'  : label,
           'ls'     : '-',
           'xlabel' : 'Residue number',
           'ylabel' : 'B-factor ($\AA^{2}$)'}
  with open(file_name, 'w') as f:
    json.dump(dump, f)


def reduce_models_common(hier):

  common = set()
  for model in hier.models():
    id_strs = {atom.pdb_label_columns() for atom in model.atoms()}
    if not common:
      common  = id_strs
    else:
      common &= id_strs

  hier = hier.select(
           flex.bool(atom.pdb_label_columns() in common for atom in hier.atoms())
         )

  return hier


class Manager():

  def __init__(self, tls_in, pdb_in, tls_origin, mult, max_level, min_level,
               sc_size, zero_trace, seed_offset=0, reset_b=False,
               remove_waters=True):

    self.input_files    = tls_in
    self.pdb_in         = pdb_in
    self.origin         = tls_origin
    self.mult           = mult
    self.zero_trace     = zero_trace
    self.max_level      = max_level
    self.min_level      = min_level
    self.sc_size        = sc_size
    self.seed_offset    = seed_offset
    self.reset_b        = reset_b
    self.remove_waters  = remove_waters
    self.working_models = []

    self.process_pdb_in()
    self.build_tls_hierarchy()

  def build_tls_hierarchy(self):

    with open('tls_analysis_log.txt','w') as log:
      file_type = get_file_extension(self.input_files)
      if file_type == 'csv':
        self.tls_hierarchy = tls_from_csv(input_files = self.input_files,
                                          origin      = self.origin,
                                          mult        = self.mult,
                                          zero_trace  = self.zero_trace,
                                          cache       = self.base_cache,
                                          max_level   = self.max_level,
                                          min_level   = self.min_level,
                                          log         = log)
      elif file_type == 'pdb':
        self.tls_hierarchy = tls_from_pdb(input_files = self.input_files,
                                          cache       = self.base_cache,
                                          mult        = self.mult,
                                          zero_trace  = self.zero_trace,
                                          max_level   = self.max_level,
                                          min_level   = self.min_level,
                                          log         = log)
      elif file_type == 'json':
        self.tls_hierarchy = tls_from_json(input_file = self.input_files[0],
                                           cache      = self.base_cache,
                                           max_level  = self.max_level,
                                           min_level  = self.min_level,
                                           mult       = self.mult,
                                           zero_trace = self.zero_trace,
                                           log        = log)
      else:
        print('Unsupported file type')

  def process_pdb_in(self):

    pdb_inp             = iotbx.pdb.input(self.pdb_in)
    pdb_hierarchy       = pdb_inp.construct_hierarchy()

    # Filter out water
    selection_string    = 'not resname HOH' if self.remove_waters else 'all'
    cache               = pdb_hierarchy.atom_selection_cache()
    water_sel           = cache.selection(selection_string)
    self.base_hierarchy = pdb_hierarchy.select(water_sel)

    # Reduce to single model
    self.base_hierarchy = reduce_models_common(self.base_hierarchy)
    self.avail_coords   = [model.atoms().extract_xyz()
                           for model in self.base_hierarchy.models()]
    for model in self.base_hierarchy.models()[1:]:
      self.base_hierarchy.remove_model(model)
    self.base_hierarchy.models()[0].atoms().set_xyz(
      sum(self.avail_coords[1:], self.avail_coords[0]) * (1./len(self.avail_coords))
    )

    # Erase temperature data
    if self.reset_b:
      for atom in self.base_hierarchy.atoms():
        atom.uij_erase()
        atom.siguij_erase()
        atom.b = 0

    # Initialize stuff for later use
    self.base_cache     = self.base_hierarchy.atom_selection_cache()
    self.base_chain     = self.base_hierarchy.models()[0].chains()[0]
    self.cryst_symm     = pdb_inp.crystal_symmetry_from_cryst1()

    # Set new crystal symmetry
    cell                = self.cryst_symm.unit_cell().parameters()
    supercell           = [a*x for a, x in zip(cell, self.sc_size)]
    supercell          += cell[3:]
    self.sc_symm = crystal.symmetry(unit_cell          = supercell,
                                    space_group_symbol = 'P1')
   
  def new_model(self, n_model=None):
    
    if n_model is None:
      n_model = len(self.working_models)
    model_manager     = Model(base_hierarchy    = self.base_hierarchy,
                              cryst_symm        = self.cryst_symm,
                              sc_size           = self.sc_size,
                              sc_symm           = self.sc_symm,
                              tls_hierarchy     = self.tls_hierarchy,
                              n_model           = n_model,
                              seed_offset       = self.seed_offset)
    self.working_models.append(n_model)
    return model_manager

class Plotter():

  def __init__(self, title='deltas', subtitle='', cutoff=None, data_range=None,
                     ylabel='Counts', xlabel='Deviation from ideal ($\AA$)'):

    self.filename = '{}.eps'.format(title.lower().replace(' ','_'))
    self.cutoff = cutoff
    self.data_range = data_range
    if subtitle != '': title += '\n' + subtitle

    self.fig = plt.figure()
    if cutoff is not None:
      gs = gridspec.GridSpec(ncols=1, nrows=5, figure=self.fig)
      self.ax = self.fig.add_subplot(gs[1:,:])
      self.ax1 = self.fig.add_subplot(gs[0,:])
      self.y_top = 0
      self.y_bottom = 0
      self.ax1.set_title(title)
    else:
      self.ax = self.fig.add_subplot(111)
      self.ax.set_title(title)
    
    self.ax.set_xlabel(xlabel)
    self.ax.set_ylabel(ylabel)
    
  def add(self, deltas, label='', scale=1):
    
    n_bins = int(len(deltas)**(1./3)) or 1
    if self.data_range is not None:
      histogram = flex.histogram(data = deltas,
                                 n_slots = n_bins,
                                 data_min = self.data_range[0],
                                 data_max = self.data_range[1])
    else: histogram = flex.histogram(data = deltas, n_slots = n_bins)
    slots = histogram.slots().as_double() * scale
    x_values = histogram.slot_centers() - histogram.slot_width()/2
    
    self.ax.hist(x        = x_values,
                 bins     = n_bins,
                 weights  = slots,
                 label    = str(label),
                 histtype = 'step')
    
    if self.cutoff is not None:
      self.ax1.hist(x        = x_values,
                    bins     = n_bins,
                    weights  = slots,
                    label    = str(label),
                    histtype = 'step')
      self.y_bottom = max(self.y_bottom,
                          slots[histogram.get_i_slot(self.cutoff)+1],
                          slots[histogram.get_i_slot(-self.cutoff)-1])
      self.y_top = min(self.y_top or self.ax1.dataLim.y1,
                       slots[histogram.get_i_slot(self.cutoff)],
                       slots[histogram.get_i_slot(-self.cutoff)])

  def plot(self):
    
    if self.cutoff is not None:
      self.ax.set_ylim(top = self.y_bottom * 1.1)
      self.ax1.set_ylim(bottom = self.y_top * 0.9)
      self.ax.spines['top'].set_visible(False)
      self.ax1.spines['bottom'].set_visible(False)
      self.ax1.set_xticks([])

      # Format y-axis cutoff lines 
      d = .015 
      kwargs = dict(transform=self.ax.transAxes, color='k', clip_on=False)
      self.ax.plot((-d,+d),(1,1), **kwargs)
      self.ax.plot((1-d,1+d),(1,1), **kwargs)
      kwargs.update(transform=self.ax1.transAxes)
      self.ax1.plot((-d,d),(0,0), **kwargs) 
      self.ax1.plot((1-d,1+d),(0,0), **kwargs)
     
    if filter(None,self.ax.get_legend_handles_labels()): self.ax.legend()
    self.fig.savefig(self.filename)

def chain_id_generator(start=0, seq=None):
  
  def base_seq(n, base): return base_seq(n//base-1, base) + [n%base] if n>=0 else []
  if seq is None: seq='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'
  while True: yield ''.join(seq[i] for i in base_seq(start, len(seq))); start+=1

def int_generator(start=0):

  while True: yield start; start+=1

def pick_from(collection, n=2):

  collection = tuple(collection)

  if n == 0:
    yield ()
  else:
    for i, item in enumerate(collection):
      for others in pick_from(collection[i+1:], n=n-1):
        yield (item,) + others


class Model():
  
  def __init__(self,
               base_hierarchy,
               cryst_symm,
               sc_size,
               sc_symm,
               n_model,
               tls_hierarchy,
               seed_offset = 0):
      
    self.cryst_symm        = cryst_symm
    self.sc_size           = sc_size
    self.sc_symm           = sc_symm
    self.tls_hierarchy     = tls_hierarchy
    self.working_chains    = []
    self.base_hierarchy    = base_hierarchy
    self.seed_offset       = seed_offset
    self.n_model           = n_model
    self.n_chains          = 0
    self.msd_collected     = 0
    self.f_collected       = 0
    self.shifts_collected  = 0

    # Extract conversion matrices

    self.orthomat = sc_symm.unit_cell().orthogonalization_matrix()
    self.fractmat = sc_symm.unit_cell().fractionalization_matrix()

  def build_model(self,
                  method='random',
                  avail_coords=[],
                  amp=1.0,
                  reverse=False,
                  residues={-1},
                  correlate=True,
                  correlate_fancy=True,
                  correlate_sequential=0,
                  correlate_insideout=0,
                  include_shifts=True,
                  include_rotvecs=True,
                  weight_exp=2.,
                  spring_weights=2.,
                  use_contacts=True,
                  correlate_after=-1,
                  rigid_only=False,
                  protein_only=True,
                  heavy_only=False,
                  hbond_only=False,
                  dist_cutoff=3.,
                  stretch=1.,
                  use_com_midpoint=False,
                  inherit_contacts=True,
                  contact_level=-1,
                  allow_bad_frac=0.,
                  disallow_good_frac=0.,
                  pairs_frac=1.,
                  skip=(),
                  override=None,
                  adopt=None,
                  pick_extremes=None,
                  extremes_from_level=1,
                  precorrelate_level=1,
                  randomize_after=[],
                  require_both=False,
                  uncorrelate=False,
                  use_pbc=True,
                  percentile=0,
                  max_level=-1,
                  max_chains=9e999,
                  regularize=False,
                  scope=1):

    # Store correlation parameters
    self.weight_exp      = weight_exp
    self.spring_weights  = spring_weights
    self.seq             = correlate_sequential
    self.iso             = correlate_insideout
    self.reverse         = reverse
    self.abf             = allow_bad_frac
    self.dgf             = disallow_good_frac
    self.require_both    = require_both
    self.uncorrelate     = uncorrelate
    self.use_pbc         = use_pbc
    self.percentile      = percentile
    self.stretch         = stretch
    self.pairs_frac      = pairs_frac
    self.skip            = skip
    self.override        = override
    self.randomize_after = randomize_after or []

    # Setup global paramters
    self.max_chains        = max_chains
    self.max_level         = max_level
    self.working_hierarchy = self.base_hierarchy.deep_copy()
    working_model          = self.working_hierarchy.models()[0]
    working_chains         = list(working_model.chains())
    self.n_atoms_in_chain  = working_model.atoms().size()
    n_asu_iter             = int_generator()
    self.chain_id_iter     = chain_id_generator(len(working_chains))
    symm_ops               = self.cryst_symm.space_group()
    size_h, size_k, size_l = self.sc_size
    translations           = ((h,k,l) for h in range(size_h)
                                      for k in range(size_k)
                                      for l in range(size_l))

    try:
      restraints_manager     = model_manager(
                                 model_input      = None,
                                 pdb_hierarchy    = self.base_hierarchy,
                                 crystal_symmetry = self.sc_symm,
                                 log              = open(os.devnull,'w')
                               ).get_restraints_manager()
    except:
      restraints_manager = None

    # Build chains based on symmetry
    for translation in translations:
      for symm_op in symm_ops:
        atoms = af_shared_atom()
        for chain in working_chains:
          new_chain    = chain.detached_copy()
          new_chain.id = next(self.chain_id_iter)
          working_model.append_chain(new_chain)
          atoms.extend(new_chain.atoms())
        self.n_chains   += 1

        chain_manager = Chain(tls_hierarchy      = self.tls_hierarchy,
                              atoms              = atoms,
                              translation        = translation,
                              symm_op            = symm_op,
                              parent_cell        = self.cryst_symm.unit_cell(),
                              n_chain            = next(n_asu_iter),
                              restraints_manager = restraints_manager,
                              residues           = set(residues))
        self.working_chains.append(chain_manager)
   
    # Clean up
    for chain in working_chains: working_model.remove_chain(chain)
    self.working_hierarchy.atoms().reset_serial()

    # Fill supercell with chains
    self.all_move_in_position()
    
    # Write TLS group colors for PyMOL
    if self.n_model == 0:
      self.write_tls_colors()

    # Initialize environments
    if correlate or precorrelate_level or regularize or self.n_model == 0:
      self.init_environments(dist_cutoff=dist_cutoff,
                             rigid_only=rigid_only,
                             protein_only=protein_only,
                             heavy_only=heavy_only,
                             hbond_only=hbond_only,
                             use_com_midpoint=use_com_midpoint,
                             inherit_contacts=inherit_contacts,
                             need_restraints=regularize,
                             contact_level=contact_level)

    # Replace starting models by some extremes
    if pick_extremes:
      self.seed(self.seed_offset)
      for chain in self.working_chains:
        chain.set_random(n_level_select=extremes_from_level, amp=amp)
      self.seed()
      self.replace_with_extremes(n=pick_extremes)
      if precorrelate_level:
        swap_shifts = (extremes_from_level == precorrelate_level)
        self.correlate_level(precorrelate_level, swap_shifts=swap_shifts)
      for chain in self.working_chains:
        chain.clear_shifts()

    # Randomize starting models from multistate input
    else:
      self.seed()
      for chain in self.working_chains:
        coords = avail_coords[np.random.randint(len(avail_coords))]
        chain.set_coordinates(chain.apply_operations(coords))

    # Introduce shifts
    if method == None or len(self.tls_hierarchy) == 0 or self.max_level == 0:
      pass
    elif method == 'analysis':
      self.init_proxies()
      self.random_analysis(amp=amp)
    elif method == 'piet':
      self.all_set_random(amp=amp, regularize_each_level=True, environments=True, reverse=True)
    else:
      if correlate and correlate_fancy:
        self.shift_and_correlate(amp=amp)
      else:
        self.all_set_random(amp=amp, max_level=correlate_after, reverse=reverse)
        if correlate:
          self.correlate(include_shifts=include_shifts, include_rotvecs=include_rotvecs,
                         weight_exp=weight_exp, spring_weights=spring_weights,
                         use_contacts=use_contacts, abf=allow_bad_frac,
                         dgf=disallow_good_frac)
        if correlate_after > -1: self.all_set_random(amp=amp, min_level=correlate_after+1)
      if regularize:
        special = 'nonbonded' if method == 'nb' else None
        self.regularize(environments=scope, special=special)

    # Conformation overrides
    if adopt:
      pairs = zip(*2*[iter(adopt)])
      for source, target in pairs:
        self.working_chains[target].adopt_from(self.working_chains[source])

  def rank_chains_by_rmsd(self):

    curr = 0
    rank = 0
    todo = set(chain.n_chain for chain in self.working_chains) - {curr}
    self.working_chains[curr].rank = rank

    while todo:

      coord_a = self.working_chains[curr].base_coordinates()

      def key(other):
        coord_b = self.working_chains[other].base_coordinates()
        return (coord_b - coord_a).dot().as_numpy_array().sum()

      rank += 1
      curr  = min(todo, key=key)
      todo -= {curr}
      self.working_chains[curr].rank = rank

  def write_all_environments(self, smooth=False):

    if hasattr(self, 'environments'):
      if smooth:
        self.rank_chains_by_rmsd()
      master_hier = self.environments.template.deep_copy()
      master_hier.atoms().reset_i_seq()
      calphas     = master_hier.get_peptide_c_alpha_selection()
      for model in master_hier.models():
        master_hier.remove_model(model)
      for chain in sorted(self.working_chains, key=lambda chain:chain.rank):
          xyz  = self.environments.get_environment_sites_cart(chain.n_chain)
          xyz  = chain.apply_reverse_operations(xyz)
          hier = self.environments.template.deep_copy()
          hier.atoms().set_xyz(xyz)
          hier = hier.select(calphas)
          for chain in hier.chains():
            chain.atoms().reset_i_seq()
          for model in hier.models():
            master_hier.append_model(model.detached_copy())
      master_hier.write_pdb_file('all_environments.pdb',
                                 crystal_symmetry=self.cryst_symm)

  def write_com_displacements(self, scale=25.):

    orig          = flex.vec3_double([self.base_hierarchy.atoms().extract_xyz().mean()])
    hier          = type(self.working_hierarchy)()
    model         = type(self.working_hierarchy.models()[0])()
    atom          = type(self.working_hierarchy.atoms()[0])
    chain_id_iter = chain_id_generator()
    serial        = 1
    diffs         = list()

    hier.append_model(model)

    chain    = type(next(self.working_hierarchy.chains()))()
    chain.id = next(chain_id_iter)
    res_iter = chain_id_generator()
    model.append_chain(chain)
    resgr         = type(next(self.working_hierarchy.residue_groups()))()
    atmgr         = type(next(self.working_hierarchy.atom_groups()))()
    atmgr.resname = next(res_iter)
    chain.append_residue_group(resgr)
    resgr.append_atom_group(atmgr)

    for chn in self.working_chains:
      old_mean    = chn.apply_operations(orig)[0]
      new_mean    = chn.coordinates().mean()
      amp_mean    = tuple((n-o) * scale + o for o, n in zip(old_mean, new_mean))
      old         = atom()
      old.xyz     = old_mean
      old.name    = 'OLD'
      old.element = 'X'
      old.serial  = '{:5d}'.format(serial)[-5:]
      atmgr.append_atom(old)
      new         = atom()
      new.xyz     = amp_mean
      new.name    = 'NEW'
      new.element = 'X'
      new.serial  = '{:5d}'.format(serial)[-5:]
      atmgr.append_atom(new)
      serial     += 1
      diffs.append(tuple(n-o for o, n in zip(old_mean, new_mean)))

    hier.write_pdb_file('com_displacements.pdb', crystal_symmetry = self.sc_symm)

    diffs = np.array(diffs)
    ref_x, ref_y, ref_z = abs(diffs).max(axis=0)
    col_vecs = (diffs / (ref_x, ref_y, ref_z) * 128 + 128).astype(int)
    col_vecs[col_vecs < 0]   = 0
    col_vecs[col_vecs > 255] = 255

    ort  = self.sc_symm.unit_cell().orthogonalize
    dist = min(((flex.vec3_double((ort(corner),)) - ort((.5,.5,.5))).dot() ** .5)[0]
               for corner in ((0,0,0),(0,0,1),(0,1,0),(0,1,1))) * 2 * 0.95

    with open('col_displacements.pml','w') as out:
      for n, (a, b, c) in enumerate(col_vecs):
        col = '0x{:02x}{:02x}{:02x}'.format(a, b, c)
        out.write(('dist _d{:d}, v. and id {:d}, v. and id {:d}, {:.2f}\n'
                  ).format(n+1,n+1,n+1, dist))
        out.write('color {}, visible and id {:d}\n'.format(col, n+1))
        out.write('color {}, _d{:d}\n'.format(col, n+1))
      out.write('group d, _d*\n')
      out.write('hide labels\n')

    # Calculate displacement dot products
    nchains = list(range(self.n_chains))
    zero    = flex.vec3_double((self.working_chains[0].init_coordinates.mean(),))
    xrs     = xray.structure(crystal_symmetry = self.sc_symm)
    ort     = xrs.unit_cell().orthogonalize
    for n in nchains:
      mean = self.working_chains[n].apply_operations(zero)[0]
      xrs.add_scatterer(
        xray.scatterer(
          label     = 'DUM',
          occupancy = 0,
          site      = xrs.unit_cell().fractionalize(mean)
        )
      )

    if self.use_pbc:
      trans   = np.mgrid[-1:2,-1:2,-1:2].reshape(3,-1).T
    else:
      trans   = np.array([[0,0,0]])
    frac    = xrs.sites_frac()
    cart    = xrs.sites_cart()
    opts    = [ort(flex.vec3_double(trans+frac[n])) for n in nchains]
    dists   = dict()
    for first in nchains:
      first_cart = cart[first]
      for second in nchains[first:]:
        vecs = opts[second] - first_cart
        dots = vecs.dot()
        best = list(dots).index(dots.min_max_mean().min)
        dists[first,second] = dists[second,first] = vecs[best]

    dot_products = list()
    for first in nchains:
      first_diff = diffs[first]
      for second in nchains[first:]:
        second_diff = diffs[second]
        dot         = np.dot(first_diff, second_diff)
        dist        = dists[first,second]
        dot_products.append(dist + (dot,))

    frc      = self.cryst_symm.unit_cell().fractionalize
    with open('dot_products_by_distance.out','w') as cart,\
         open('dot_products_by_distance_frac.out', 'w') as frac:
      for a, b, c, d in dot_products:
        cart.write('{:8f} {:8f}\n'.format((a*a+b*b+c*c)**0.5, d))
        frac.write('{:8f} {:8f} {:8f} {:8f}\n'.format(*frc((a,b,c))+(d,)))

  def rmsd_per_asu(self):

    rmsds    = list()

    # Calculate rmsds
    for nasu in range(self.n_chains):
      msd_byatom   = self.calc_msd(slc=slice(nasu, nasu+1))
      rmsd         = (sum(msd_byatom) / msd_byatom.size()) ** 0.5
      rmsds.append(rmsd)

    dists = self.get_distances('asu')

    for n,d in enumerate(dists):
      print()
      for i in sorted(d, key=d.__getitem__):
        print('{:5d}{:5d}{:7.2f}{:6.3f}'.format(n,i,d[i], rmsds[i]))

  def rmsd_per_unit_cell(self):

    symm_ops = len(self.cryst_symm.space_group())
    rmsds    = list()

    # Calculate rmsds
    for ncell in range(self.n_chains // symm_ops):
      start        = ncell * symm_ops
      msd_byatom   = self.calc_msd(slc=slice(start, start+symm_ops))
      rmsd         = (sum(msd_byatom) / msd_byatom.size() / symm_ops) ** 0.5
      rmsds.append(rmsd)

    dists = self.get_distances('cell')

    for n,d in enumerate(dists):
      print()
      for i in sorted(d, key=d.__getitem__):
        print('{:4d}{:4d}{:7.2f}{:6.3f}'.format(n,i,d[i], rmsds[i]))

  def chains_as_models(self):

    self.all_move_back()
    model          = self.working_hierarchy.models()[0]
    chains_per_asu = len(tuple(self.working_hierarchy.chains())) // self.n_chains

    for n_chain, chain in enumerate(self.working_hierarchy.chains()):
      if chain.is_protein():
        if n_chain % chains_per_asu == 0:
          chain_id_iter = chain_id_generator()
          new_model     = type(model)()
          self.working_hierarchy.append_model(new_model)
        chain.id  = next(chain_id_iter)
        new_model.append_chain(chain.detached_copy())
        new_model.atoms().reset_serial()

    self.working_hierarchy.remove_model(model)
    self.sc_symm = self.cryst_symm

  def add_coms_to_hierarchy(self):

    model = self.working_hierarchy.models()[0]
    chain = type(model.chains()[0])()
    resgr = type(model.chains()[0].residue_groups()[0])()
    atmgr = type(model.chains()[0].residue_groups()[0].atom_groups()[0])()
    atom = type(self.working_hierarchy.atoms()[0])
    for working_chain in self.working_chains:
      new_atom         = atom()
      new_atom.xyz     = working_chain.coordinates().mean()
      new_atom.name    = 'DUM'
      new_atom.element = 'X'
      atmgr.append_atom(new_atom)
    atmgr.resname = 'COM'
    chain.id      = '99'
    resgr.append_atom_group(atmgr)
    chain.append_residue_group(resgr)
    model.append_chain(chain)
    self.working_hierarchy.atoms().reset_serial()

  def close(self, keep_hierarchy = False, force = True):
    
    if not force and self.n_model == 0: return

    # Clean up attributes (for pickling)
    if not keep_hierarchy: 
      if hasattr(self, 'working_hierarchy'):
        del self.working_hierarchy
      if hasattr(self, 'sc_symm'):
        del self.sc_symm
    if hasattr(self, 'tls_hierarchy'):
      del self.tls_hierarchy
    if hasattr(self, 'working_chains'):
      del self.working_chains[:]
    if hasattr(self, 'chain_id_iter'):
      del self.chain_id_iter
    if hasattr(self, 'environments'):
      del self.environments
    if hasattr(self, 'intermolecular_proxies'):
      del self.intermolecular_proxies

  def coordinates(self,fractional=False):

    coordinates = self.working_hierarchy.atoms().extract_xyz()
    if fractional: return coordinates * self.fractmat
    else         : return coordinates

  def seed(self, seed=None):

    if seed is None and self.seed_offset != 1:
      seed = self.n_model + self.seed_offset
      
    np.random.seed(seed)
    random.seed(seed)

  def all_move_in_position(self):
    
    for chain in self.working_chains:
      chain.move_in_position()

  def all_move_back(self):
    
    for chain in self.working_chains:
      chain.move_back()

  def all_revert(self):

    for chain in self.working_chains:
      chain.revert()

  def all_set_random(self,
                     amp                       = 1.0,
                     levels                    = None,
                     min_level                 = 0,
                     max_level                 = -1,
                     reverse                   = False,
                     regularize_each_level     = False,
                     environments              = False):
   
    levels = levels or sorted(self.tls_hierarchy.keys(), reverse=reverse)
    levels = [lvl for lvl in levels if lvl <= self.max_level or self.max_level < 0]
    if max_level < 0: max_level = max(levels or [0])
    for n_level in levels:
      if min_level <= n_level <= max_level:
        for chain in self.working_chains:
          if chain.n_chain < self.max_chains:
            chain.set_random(n_level_select=n_level, amp=amp)
        if regularize_each_level:
          print('Regularizing level {}'.format(n_level))
          self.regularize(environments=environments)

  def all_coms(self, orig=False):

    centers_of_mass = []
    for chain in self.working_chains:
      if orig:
        coordinates = chain.reverted_coordinates()
      else:
        coordinates = chain.coordinates()
      centers_of_mass.append(coordinates.mean())
    return flex.vec3_double(centers_of_mass)

  def all_rotvecs(self):

    rotvecs = []
    for chain in self.working_chains:
        rotvecs.append(chain.rotvec())
    return flex.vec3_double(rotvecs)
  
  def all_cancel_rotations(self):

    for chain in self.working_chains:
      chain.cancel_rotations()

  def optimize_vectors(self, points, shifts=None, rotvecs=None, weight_exp=2,
                       allow_bad_frac=0., disallow_good_frac=0.):

    assert any((shifts, rotvecs))
    if shifts == None: shifts = rotvecs
    if rotvecs == None: rotvecs = shifts
    shifts = shifts.as_numpy_array()
    rotvecs = rotvecs.as_numpy_array()
    L,D = shifts.shape
    if hasattr(self, 'environments'):
      indices = np.array(self.environments.indices)
      weights = np.array(list(map(len, self.environments.dstlist))) ** weight_exp
    else:
      points = (points * self.fractmat).as_numpy_array() % 1
      N   = min(12, L-1)
      R   = (0,1,-1)
      stack = np.row_stack([points+(i,j,k) for i in R for j in R for k in R])
      stack = (flex.vec3_double(stack) * self.orthomat).as_numpy_array().reshape(-1,L,D)
      indices = np.array([np.square(point-stack).sum(axis=2).min(axis=0).argsort()[1:N+1]
                          for point in stack[0]])
      weights = np.ones(N)
    weights = weights.reshape(1,-1,1)/weights.sum()
    pairs = np.random.permutation(np.array(np.triu_indices(L,1)).T)
    order = np.arange(L)
    for pair in pairs:
      i,j = pair
      selection = indices[pair]
      nb_shifts = (shifts[selection] * weights).sum(axis=1)
      nb_rotvecs = (rotvecs[selection] * weights).sum(axis=1)
      old_shift_i, old_shift_j = shifts[pair]
      old_rotvec_i, old_rotvec_j = rotvecs[pair]
      new_shift_i, new_rotvec_i = (self.working_chains[i].symm_rot_double * (
                                    self.working_chains[j].symm_rot_inverse 
                                  * flex.vec3_double([old_shift_j, old_rotvec_j])
                                  )).as_numpy_array()
      new_shift_j, new_rotvec_j = (self.working_chains[j].symm_rot_double * (
                                    self.working_chains[i].symm_rot_inverse 
                                  * flex.vec3_double([old_shift_i, old_rotvec_i])
                                  )).as_numpy_array()
      old_shifts = np.square(nb_shifts - (old_shift_i, old_shift_j)).sum()
      old_rotvecs = np.square(nb_rotvecs - (old_rotvec_i, old_rotvec_j)).sum()
      new_shifts = np.square(nb_shifts - (new_shift_i, new_shift_j)).sum()
      new_rotvecs = np.square(nb_rotvecs - (new_rotvec_i, new_rotvec_j)).sum()

      allow = new_shifts < old_shifts and new_rotvecs < old_rotvecs
      rnd   = np.random.random()

      if (allow and rnd >= disallow_good_frac) or (not allow and rnd < allow_bad_frac):
        shifts[pair]  = (new_shift_i, new_shift_j)
        rotvecs[pair] = (new_rotvec_i, new_rotvec_j)
        order[pair]   = order[pair[::-1]]

    return order

  def optimize_contacts(self, weight_exp=2, spring_weights=2,
                        allow_bad_frac=0., disallow_good_frac=0.):

    assert hasattr(self, 'environments')
    contact = np.array(self.environments.contact)
    indices = np.array(self.environments.indices)
    inverse = np.array(self.environments.inverse)
    translt = np.array(self.environments.translt)
    weights = np.array(list(map(len, self.environments.dstlist))) ** weight_exp
    all_contacts = []
    for chain in self.working_chains:
      all_contacts.append(
        chain.apply_operations(chain.rotmat().as_mat3() * flex.vec3_double(contact))
      + chain.trans()
      )
    all_contacts = np.array(all_contacts)

    L, D = indices.shape
    exp = spring_weights / 2.
    weights = weights.reshape(-1,1)/weights.sum()
    pairs = np.random.permutation(np.array(np.triu_indices(L,1)).T)
    order = np.arange(L)
    for pair in pairs:
      i,j = pair
      ref_i = all_contacts[indices[i],inverse] + translt[i]
      ref_j = all_contacts[indices[j],inverse] + translt[j]
      old_i = all_contacts[i]
      old_j = all_contacts[j]
      new_i = self.working_chains[i].apply_operations(
                self.working_chains[j].apply_reverse_operations(flex.vec3_double(old_j))
              ).as_numpy_array()
      new_j = self.working_chains[j].apply_operations(
                self.working_chains[i].apply_reverse_operations(flex.vec3_double(old_i))
              ).as_numpy_array()
      if exp == 1.:
        old_spread = ((np.square(ref_i - old_i) + np.square(ref_j - old_j)) * weights).sum()
        new_spread = ((np.square(ref_i - new_i) + np.square(ref_j - new_j)) * weights).sum()
      else:
        old_spread = ((np.square(ref_i - old_i).sum(axis=1)**exp
                     + np.square(ref_j - old_j).sum(axis=1)**exp) * weights[...,0]).sum()
        new_spread = ((np.square(ref_i - new_i).sum(axis=1)**exp
                     + np.square(ref_j - new_j).sum(axis=1)**exp) * weights[...,0]).sum()

      allow = new_spread < old_spread
      rnd   = np.random.random()

      if (allow and rnd >= disallow_good_frac) or (not allow and rnd < allow_bad_frac):
        all_contacts[i] = new_i
        all_contacts[j] = new_j
        order[pair]     = order[pair[::-1]]

    return order

  def swap_coordinates(self, first, second):

    chain_a = self.working_chains[first]
    chain_b = self.working_chains[second]
    chain_a.move_back()
    chain_b.move_back()
    coord_a = chain_a.coordinates()
    coord_b = chain_b.coordinates()
    chain_a.set_coordinates(coord_b)
    chain_b.set_coordinates(coord_a)
    chain_a.move_in_position()
    chain_b.move_in_position()

  def correlate(self, include_shifts=True, include_rotvecs=True,
                weight_exp=2., spring_weights=2., use_contacts=True,
                abf=0., dgf=0.):

    if use_contacts:
      order      = self.optimize_contacts(weight_exp=weight_exp, spring_weights=spring_weights,
                                          allow_bad_frac=abf, disallow_good_frac=dgf)
    else:
      points     = self.all_coms(orig=True)
      shifts     = points - self.all_coms() if include_shifts else None
      rotvecs    = self.all_rotvecs() if include_rotvecs else None
      order      = self.optimize_vectors(points,
                                         shifts=shifts,
                                         rotvecs=rotvecs,
                                         weight_exp=weight_exp,
                                         allow_bad_frac=abf,
                                         disallow_good_frac=dgf)
   
    done = set()
    x = 0
    while len(done) < len(order):
      done.add(x)
      y = order[x]
      if y in done: x = (set(order) - done or {-1}).pop()
      else: self.swap_coordinates(x,y); x=y

  def optimize_level(self, n_level, weight_exp=2., spring_weights=2., seq=0, iso=0,
                     allow_bad_frac=0., disallow_good_frac=0., require_both=False,
                     uncorrelate=False, percentile=0, stretch=1., pairs_frac=1., 
                     override=None, swap_shifts=False):

    contact = np.array(self.environments.get_contact(n_level))
    conshft = np.array(self.environments.get_conshft(n_level))
    conlens = np.array(self.environments.get_conlens(n_level))
    indices = {i: np.array(array) for i, array in
               self.environments.get_indices(n_level).items() if any(array)}
    translt = {i: np.array(array) for i, array in
               self.environments.get_translt(n_level).items()}
    weights = {i: np.array(list(map(len,array))) ** weight_exp for i, array in
               self.environments.get_dstlist(n_level).items()}
    weights = {i: array.reshape(-1,1)/array.sum() for i, array in weights.items()}
    grpsort = self.environments.get_grpsort(n_level)

    exp     = spring_weights / 2.
    chains  = self.working_chains

    if override is not None and n_level == 1:
      weights = {i: np.array((override+[1]*len(w))[:len(w)])[...,None] * w
                 for i,w in weights.items()}

    for chain in chains:
      for n_group, mod in chain.mod_hierarchy[n_level].items():
        com  = tuple(mod.get_com())
        cont = contact[n_group][chain.n_chain]
        contact[n_group][chain.n_chain]  = (cont - com) * stretch + com

    def group_energies(g):
      idxs = indices[g]
      trns = translt[g]
      clen = conlens[g]
      wght = weights[g]
      engy = []
      for i in range(len(chains)):
        ref_i = contact[tuple(idxs[i].T)] + trns[i]
        old_i = contact[g,i,:clen]
        engy.append((np.square(ref_i - old_i) * wght).sum())
      return engy

    percentiles = {g: np.percentile(group_energies(g), percentile) for g in indices}

    order   = {i: np.arange(len(chains)) for i in indices}

    active  = [i for i in indices if chains[0].get_modifier(n_level,i).is_active]
    pairs   = np.fromiter((j for p in zip(*np.triu_indices(len(chains),1))
                           for i in active for j in (i,)+p), dtype=int).reshape(-1,3)
    np.random.shuffle(pairs)
    if pairs_frac < 0:
      pairs_frac = np.random.random() * abs(pairs_frac)
    pairs   = pairs[:int(len(pairs) * pairs_frac)]

    if seq:
      pairs = pairs[pairs[:,0].argsort(kind='stable')[::seq]]

    if iso:
      pairs = pairs[np.argsort(map(grpsort.index,pairs[:,0]), kind='stable')[::iso]]

    i_mask = j_mask = Ellipsis

    for g,i,j in pairs:

      idxs = indices[g]
      trns = translt[g]
      clen = conlens[g]
      wght = weights[g]

      ref_i = contact[tuple(idxs[i].T)] + trns[i]
      ref_j = contact[tuple(idxs[j].T)] + trns[j]

      old_i = contact[g,i,:clen]
      old_j = contact[g,j,:clen]

      if swap_shifts:
        shift = flex.vec3_double(conshft[g,j,:clen] - conshft[g,i,:clen])
        new_i = old_i + (chains[i].symm_rot_double *        shift )
        new_j = old_j + (chains[j].symm_rot_double * (-1. * shift))
      else:
        new_i = chains[i].apply_operations(
                  chains[j].apply_reverse_operations(flex.vec3_double(old_j))
                ).as_numpy_array()
        new_j = chains[j].apply_operations(
                  chains[i].apply_reverse_operations(flex.vec3_double(old_i))
                ).as_numpy_array()

      if not self.use_pbc:
        i_mask = ~trns[i].any(axis=1)
        j_mask = ~trns[j].any(axis=1)

      if exp == 1.:
        old_i_energy = (np.square(ref_i - old_i) * wght)[i_mask].sum()
        new_i_energy = (np.square(ref_i - new_i) * wght)[i_mask].sum()
        old_j_energy = (np.square(ref_j - old_j) * wght)[j_mask].sum()
        new_j_energy = (np.square(ref_j - new_j) * wght)[j_mask].sum()
      else:
        old_i_energy = (np.square(ref_i - old_i).sum(axis=1)**exp * wght[...,0])[i_mask].sum()
        new_i_energy = (np.square(ref_i - new_i).sum(axis=1)**exp * wght[...,0])[i_mask].sum()
        old_j_energy = (np.square(ref_j - old_j).sum(axis=1)**exp * wght[...,0])[j_mask].sum()
        new_j_energy = (np.square(ref_j - new_j).sum(axis=1)**exp * wght[...,0])[j_mask].sum()

      if require_both:
        allow = (new_i_energy < old_i_energy) and (new_j_energy < old_j_energy)
      else:
        allow = (new_i_energy + new_j_energy) < (old_i_energy + old_j_energy)

      limit = percentiles[g]
      allow = allow and old_i_energy >= limit and old_j_energy >= limit

      rnd   = np.random.random()
      allow = (allow and rnd >= disallow_good_frac) or (not allow and rnd < allow_bad_frac)

      if uncorrelate:
        allow = not allow

      if allow:
        contact[g,i,:clen] = new_i
        contact[g,j,:clen] = new_j
        conshft[g,[i,j]]   = conshft[g,[j,i]]
        order[g][[i,j]]    = order[g][[j,i]]

    if n_level == 1:
      ranks = np.argsort(group_energies(min(indices)))
      for chain, rank in zip(self.working_chains, ranks):
        chain.rank = rank

    return order

  def swap_shifts(self, n_level, n_group, first, second):

    mod_a    = self.working_chains[first].get_modifier(n_level, n_group)
    mod_b    = self.working_chains[second].get_modifier(n_level, n_group)
    shifts_a = mod_a.get_shifts()
    shifts_b = mod_b.get_shifts()
    mod_a.set_shifts(shifts_b)
    mod_b.set_shifts(shifts_a)

  def swap_positions(self, n_level, n_group, first, second):

    mod_a    = self.working_chains[first].get_modifier(n_level, n_group)
    mod_b    = self.working_chains[second].get_modifier(n_level, n_group)
    coords_a = mod_a.base_coordinates()
    coords_b = mod_b.base_coordinates()
    mod_a.set_coordinates(coords_b)
    mod_b.set_coordinates(coords_a)

  def correlate_level(self, n_level, swap_shifts=True):

    orders = self.optimize_level(n_level, swap_shifts=swap_shifts, weight_exp=self.weight_exp,
                                 spring_weights=self.spring_weights, seq=self.seq, iso=self.iso,
                                 allow_bad_frac=self.abf, disallow_good_frac=self.dgf,
                                 require_both=self.require_both, uncorrelate=self.uncorrelate,
                                 percentile=self.percentile, stretch=self.stretch,
                                 pairs_frac=self.pairs_frac, override=self.override)

    for n_group, order in orders.items():
      done = set()
      x = 0
      while len(done) < len(order):
        done.add(x)
        y = order[x]
        if y in done: x = (set(order) - done or {-1}).pop()
        else:
          if swap_shifts: self.swap_shifts(n_level, n_group, x, y)
          else: self.swap_positions(n_level, n_group, x, y)
          x=y

  def shift_and_correlate(self, amp=1.0):

    levels = sorted(self.tls_hierarchy.keys())
    levels = [lvl for lvl in levels if lvl <= self.max_level or self.max_level < 0]
    if self.reverse: levels = levels[::-1]
    for n_level in levels:
      for chain in self.working_chains:
        if chain.n_chain < self.max_chains:
          chain.set_random(n_level_select=n_level, amp=amp)
      if n_level not in self.skip:
        self.correlate_level(n_level)
      if n_level in self.randomize_after:
        self.randomly_swap_all_chains()

  def randomly_swap_all_chains(self):

    for i in reversed(range(self.n_chains)):
      j = int(np.random.random() * (i+1))
      self.swap_chains(i, j)

  def swap_chains(self, first, second):

    first  = self.working_chains[first]
    second = self.working_chains[second]

    for n_level, level in first.mod_hierarchy.items():
      for n_group, modifier in level.items():
        first_shifts = modifier.get_shifts()
        modifier.set_shifts(second.get_modifier(n_level, n_group).get_shifts())
        second.get_modifier(n_level, n_group).set_shifts(first_shifts)

  def replace_with_extremes(self, n=()):

    if len(n) == 1 and n[0]%1 == 0:
      pick  = int(n[0])
      odds  = np.ones(max(2,pick))
      odds /= odds.sum()
    else:
      pick  = len(n)
      odds  = np.array((list(n) + [1.-n[0]])[:max(2,pick)])
      odds /= odds.sum()

    if pick in (0, self.n_chains): return

    extremes   = tuple(self.find_extremes(pick))
    if pick < 2:
      extremes = extremes + (None,)
    others     = sorted(set(self.working_chains) - set(extremes))

    for chain in others:
      other = extremes[np.random.choice(len(extremes), p=odds)]
      chain.adopt_from(other)

    for chain in sorted(set(extremes) - {None}):
      other = others[np.random.choice(len(others))]
      chain.adopt_from(other)

  def find_extremes(self, n=2):

    def get_coords(chain, cache=dict()):
      result = cache.get(chain.n_chain)
      if result is None:
        result = chain.base_coordinates()
        cache[chain.n_chain] = result
      return result

    def get_rmsd(a, b, cache=dict()):
      result = cache.get((a.n_chain,b.n_chain))
      if result is None:
        result = flex.sum((get_coords(a) - get_coords(b)).dot()) ** 0.5
        cache[a.n_chain,b.n_chain] = result
      return result

    def key(chain):
      return flex.sum((get_coords(chain)-chain.init_coordinates).dot()) ** 0.5

    candidates = sorted(self.working_chains, key=key)[::-1][:10*n]

    if n == 1:
      return candidates[:1]

    best_score = -1
    best_combo = ()
    for combo in pick_from(candidates, n):
      score = 0
      for a, b in pick_from(combo, 2):
        score += get_rmsd(a, b)
      if score > best_score:
        best_score = score
        best_combo = combo

    return best_combo

  def random_analysis(self, amp=1.0):

    # Start randomization and analysis by level of motion
    intermolplot = Plotter('Intermolecular distance deviations',
                           'Between {} structures'.format(self.n_chains),
                           data_range=(-3,3))
    intramolplot = Plotter('Bond distance deviations',
                           'Average over {} structures'.format(self.n_chains),
                           cutoff=0.005,
                           data_range=(-1.5,1.5))
    maxdeltaplot = Plotter('Maximum absolute bond distance deviations',
                           'Out of {} structures'.format(self.n_chains),
                           ylabel='Structures',
                           xlabel='Maximum deviation from ideal ($\AA$)',
                           data_range=(0,3))
    for n_level in sorted(self.tls_hierarchy):
      for chain in self.working_chains: chain.set_random(n_level, amp=amp)
      intermolplot.add(self.intermolecular_proxies.deltas(self.coordinates()),
                       'Level {}'.format(n_level))
      for chain in self.working_chains:
        if chain.n_chain == 0: bonded_deltas = chain.bonded_deltas()
        else: bonded_deltas.extend(chain.bonded_deltas())
      intramolplot.add(bonded_deltas,
                       'Level {}'.format(n_level),scale=(1./self.n_chains))
      max_bonded_deltas = flex.double(max(abs(chain.bonded_deltas()))
                                        for chain in self.working_chains)
      maxdeltaplot.add(max_bonded_deltas,
                       'Level {}'.format(n_level))
    intermolplot.plot()
    intramolplot.plot()
    maxdeltaplot.plot()

  def all_is_in_position(self):

    return all(chain.is_in_position for chain in self.working_chains)

  def none_is_in_position(self):

    return not any(chain.is_in_position for chain in self.working_chains)

  def xray_structure(self):

    xray_structure = self.working_hierarchy.extract_xray_structure(
                     crystal_symmetry = self.sc_symm)
    return xray_structure

  def calc_f_model(self,
                   high_resolution = 2.,
                   k_sol           = 0.0,
                   b_sol           = 0.0,
                   b_cart          = 6 * [0.0]):

    self.close(keep_hierarchy = True, force = False)

    xrs        = self.xray_structure()
    params     = fmodel_from_xray_structure_master_params.extract()

    # Set params
    params.high_resolution = high_resolution
    params.fmodel.k_sol    = k_sol
    params.fmodel.b_sol    = b_sol
    params.fmodel.b_cart   = b_cart

    miller_set = xrs.build_miller_set(anomalous_flag=False, d_min=high_resolution)
    f_obs      = miller_set.array(flex.double(miller_set.size()))

    f_model    = fmodel_from_xray_structure(xray_structure = xrs,
                                            f_obs          = f_obs,
                                            params         = params,)

    return f_model.f_model

  def calc_msd(self, slc=slice(None,None)):

    for n, chain in enumerate(self.working_chains[slc]):
      if n == 0:
        sum_sd  = chain.calc_sd().deep_copy()
      else:
        sum_sd += chain.calc_sd()
    return sum_sd

  def write_tls_colors(self):

    for n_level, level in self.tls_hierarchy.items():
      level.as_pymol_colors(file_name     = 'color_tls_level_{}.pml'.format(n_level),
                            pdb_hierarchy = self.base_hierarchy)

  def write_pdb(self, pdb_out=None):

    if not pdb_out: pdb_out = 'supercell_out_{}.pdb'.format(self.n_model)
    self.working_hierarchy.write_pdb_file(pdb_out,
                                          crystal_symmetry = self.sc_symm)

  def init_proxies(self, intermolecular_cutoff=4., write_pymol_dashes=False):

    xrs = self.xray_structure()

    # Setup intermolecular proxies
    n_atoms = self.n_atoms_in_chain
    pair_asu_table = xrs.pair_asu_table(distance_cutoff = intermolecular_cutoff)
    for i_seq, pair_asu_dict in enumerate(pair_asu_table.table()):
      i_chain = i_seq//n_atoms
      for j_seq in pair_asu_dict:
        if j_seq//n_atoms == i_chain: pair_asu_dict.erase(j_seq)
    self.intermolecular_proxies = geometry_restraints.bond_sorted_asu_proxies(
      pair_asu_table = pair_asu_table)
    if write_pymol_dashes:
      with open('dashes_intermolecular.pml', 'w') as dashes:
        dashes.write(
          self.intermolecular_proxies.simple.as_pymol_dashes(self.working_hierarchy)
        )

  def init_environments(self, dist_cutoff=3., rigid_only=False,
                        protein_only=True, heavy_only=False, hbond_only=False,
                        use_com_midpoint=False, inherit_contacts=True,
                        contact_level=-1, need_restraints=False):

    self.environments = Environments(chains            = self.working_chains,
                                     crystal_symmetry  = self.sc_symm,
                                     base_symmetry     = self.cryst_symm,
                                     base_hierarchy    = self.base_hierarchy,
                                     write_environment = not self.n_model,
                                     dist_cutoff       = dist_cutoff,
                                     rigid_only        = rigid_only,
                                     protein_only      = protein_only,
                                     heavy_only        = heavy_only,
                                     hbond_only        = hbond_only,
                                     use_com_midpoint  = use_com_midpoint,
                                     inherit_contacts  = inherit_contacts,
                                     contact_level     = contact_level,
                                     need_restraints   = need_restraints,
                        )

  def regularize(self, environments=False, special=None):

    shifts_total = None
    for chain in self.working_chains:
      if environments:
        shifts = chain.regularize(environments = self.environments,
                                  keep_shifts  = True,
                                  special      = special)
      else:
        shifts = chain.regularize(keep_shifts  = True,
                                  special      = special)
      if shifts_total is None:
        shifts_total  = shifts.dot()
      else:
        shifts_total += shifts.dot()
    self.shifts_total = shifts_total


class Environments():

  def __init__(self,
               chains,
               crystal_symmetry,
               base_symmetry,
               base_hierarchy,
               dist_cutoff       = 3.,
               protein_only      = True,
               rigid_only        = False,
               heavy_only        = False,
               hbond_only        = False,
               use_com_midpoint  = False,
               inherit_contacts  = True,
               contact_level     = -1,
               write_environment = False,
               need_restraints   = False):

    self.chains   = chains
    self.crystal_symmetry = crystal_symmetry
    self.base_hierarchy = base_hierarchy
    self.fractmat = crystal_symmetry.unit_cell().fractionalization_matrix()
    self.orthomat = crystal_symmetry.unit_cell().orthogonalization_matrix()

    # Replace chains with mean position
    means    = self.fractmat * flex.vec3_double(
                chain.coordinates().mean() for chain in chains
               )

    # Find atom farthest away from mean and
    # set search radius to twice this distance
    xyz      = chains[0].coordinates()
    radius   = 2 * np.sqrt(max((xyz-xyz.mean()).dot()))

    # Create dummy xray.structure object with means
    xrs = xray.structure(crystal_symmetry = crystal_symmetry)
    for mean in means:
      xrs.add_scatterer(
        xray.scatterer(
          label     = 'DUM',
          occupancy = 0,
          site      = mean
        )
      )
    # Create table listing all nearest neighbour chains (inc. symm. mates)
    pair_sym_table = xrs.pair_asu_table(
                       distance_cutoff = radius
                     ).extract_pair_sym_table(
                       skip_j_seq_less_than_i_seq       = False,
                       all_interactions_from_inside_asu = True
                     )

    # Get reference coordination vectors from original crystal symmetry
    selection      = flex.bool(self.base_hierarchy.atoms().size(), True)
    if protein_only:
      selection   &= base_hierarchy.atom_selection_cache().sel_protein()
    if rigid_only:
      ss           = secondary_structure.manager(self.base_hierarchy,
                                                   log = open(os.devnull,'w'))
      selection   &= ss.helix_selection() | ss.beta_selection()
    if heavy_only:
      selection   &= flex.bool(not atom.element_is_hydrogen()
                               for atom in self.base_hierarchy.atoms())
    if hbond_only:
      selection   &= flex.bool(atom.element.strip() in {'N', 'O', 'F'}
                               for atom in self.base_hierarchy.atoms())
    actual_ids     = [i for i, b in enumerate(selection) if b]
    base_xrs       = self.base_hierarchy.extract_xray_structure(
                      crystal_symmetry = base_symmetry
                     )
    base_sym_table = base_xrs.select(selection).pair_asu_table(
                       distance_cutoff = dist_cutoff
                     ).extract_pair_sym_table(
                       skip_j_seq_less_than_i_seq       = False,
                       all_interactions_from_inside_asu = True
                     )
    sym_edge_list  = base_sym_table.symmetry_edge_list()
    sites_frac     = base_xrs.sites_frac()
    cell           = base_symmetry.unit_cell()
    raw_symmetries = [pair[2] for pair in sym_edge_list]
    symmetries     = list(set(raw_symmetries))
    mean           = flex.vec3_double((sites_frac.mean(),))
    ref_vectors    = flex.vec3_double(
                       (cell.orthogonalization_matrix() * 
                         (symm.r().as_double() * mean + symm.t().as_double() - mean)
                       )[0]
                     for symm in symmetries)
    ref_round      = ref_vectors.round(1) # To avoid rounding errors while sorting
    ref_sorted     = list(zip(*sorted(zip(ref_round, ref_vectors))))[1]
    sym_sorted     = list(zip(*sorted(zip(ref_round, symmetries))))[1]
    
    # Sort and filter environments by mean-mean vector
    self.sorted_environments = {}
    for i_seq, pair_sym_dict in enumerate(pair_sym_table):
      vectors    = []
      for j_seq, symm_ops in pair_sym_dict.items():
        for symm_op in symm_ops:
          first  = means[i_seq:i_seq+1]
          second = means[j_seq:j_seq+1] + symm_op.t().as_double()
          vector = (self.orthomat * (self.chains[i_seq].symm_rot_inverse * (second-first))
                   ).round(1)[0]
          # Filter by reference vectors
          if min((ref_round - vector).dot()) <= 1.:
            vectors.append((vector,(j_seq, symm_op)))
      environment = list(zip(*sorted(vectors)))[1] if vectors else tuple()
      self.sorted_environments[i_seq] = environment

    environ = [value for key, value in sorted(self.sorted_environments.items())]

    # Prepare connectivity tables
    self.contact = 0.5 * flex.vec3_double(ref_sorted) + base_xrs.sites_cart().mean()
    self.indices = [list(zip(*env))[0] for env in environ]
    self.translt = [[self.crystal_symmetry.unit_cell().orthogonalize(conn[1].t().as_double())
                     for conn in env] for env in environ]
    self.inverse = [next(l for l,(j,k) in enumerate(environ[n]) if (j,k.t().as_double()) ==
                    (0,(-sym.t()).as_double())) for n,sym in environ[0]]
    self.comdist = [(x*x+y*y+z*z)**0.5 for x,y,z in ref_sorted]
    self.dstlist = [[] for _ in symmetries]
    for i, j, sym in sym_edge_list:
        dist = flex.sqrt(
                 (
                   cell.orthogonalize(
                     (sym.r().as_double() * cell.fractionalize(xyz[j:j+1]))
                     + sym.t().as_double()
                   ) - xyz[i]
                 ).dot()
               )[0]
        self.dstlist[sym_sorted.index(sym)].append(dist)

    # Prepare complex connectivity
    self.lvl_indices = dict()
    self.lvl_translt = dict()
    self.lvl_conlens = dict()
    self.lvl_dstlist = dict()
    self.lvl_grpsort = dict()

    if inherit_contacts and chains[0].mod_hierarchy:
      group_supers = dict()
      group_lookup = dict()
      prev_level = None
      lvl_groupcoms = dict()
      mod_hierarchy = chains[0].mod_hierarchy
 
      for n_level, level in sorted(mod_hierarchy.items()):
        lvl_groupcoms[n_level] = group_coms = dict()
        for n_group, modifier in level.items():
          curr = (n_level, n_group)
          for atom in modifier.working_atoms():
            if curr not in group_supers:
              prev = (prev_level, group_lookup.get(atom.i_seq))
              group_supers[curr] = [curr] + group_supers.get(prev,[])
            group_lookup[atom.i_seq] = n_group
          group_coms[n_group] = cell.fractionalize(modifier.get_com())
        prev_level = n_level
      group_subs = {a:sorted(b for b, subs in group_supers.items() if a in subs[1:2])
                    for a in group_supers}
      max_level = n_level

      contacting_groups = {n:dict() for n in level.keys()}
      for i, pair_sym_dict in enumerate(base_sym_table):
        i = actual_ids[i]
        own_group = group_lookup.get(i)
        for j, sym_ops in pair_sym_dict.items():
          j = actual_ids[j]
          other_group = group_lookup.get(j)
          for sym in sym_ops:
            if (own_group and other_group) \
             and not (own_group == other_group and sym.is_unit_mx()):
              a = cell.orthogonalize(
                   (sym.r().as_double() * cell.fractionalize(xyz[j:j+1]))
                    + sym.t().as_double()
                  )
              b = xyz[i]
              dist = (flex.sqrt((a-b).dot())[0], ((a+b)*0.5)[0])
              if (other_group, sym) not in contacting_groups.get(own_group):
                contacting_groups[own_group][(other_group, sym)] = []
              contacting_groups[own_group][(other_group, sym)].append(dist)
      umx = type(sym)()
      for i, groups in contacting_groups.items():
        entry = [(j,sym,list(zip(*dist))) for (j,sym),dist in groups.items()]
        entry.sort(key=lambda x:(((umx,)+sym_sorted).index(x[1]),x[0]))
        contacting_groups[i] = entry
      contacts = dict()
      for i, entry in contacting_groups.items():
        com = group_coms[i]
        if use_com_midpoint:
          contacts[i] = flex.vec3_double([
            cell.orthogonalize(
              (((sym.r().as_double() * flex.vec3_double((group_coms[j],)))
              + sym.t().as_double() - com) * 0.5 + com)[0]
            ) for j,sym,dist in entry
          ])
        else:
          contacts[i] = flex.vec3_double([
            flex.vec3_double(dist[1]).mean() for j,sym,dist in entry
          ])
      lvl_unique = dict()
      lvl_select = dict()
      lvl_select[n_level] = dict()
      for n_group, vec in sorted(contacts.items()):
        last = sum(map(len, lvl_select[n_level].values()))
        lvl_select[n_level][n_group] = flex.size_t_range(last, last+vec.size())
      if contacts:
        contacts = flex.vec3_double(
                     sum(map(tuple,list(zip(*sorted(contacts.items())))[1]),())
                   )
      else:
        contacts = flex.vec3_double()
      for n_level, level in sorted(mod_hierarchy.items(), reverse=True):
        if n_level not in lvl_select:
          lvl_select[n_level] = dict()
          for n_group in level:
            subs = [lvl_select[lvl][grp] for lvl,grp in group_subs[(n_level,n_group)]]
            lvl_select[n_level][n_group] = flex.size_t(sum(map(tuple, subs),()))
        lvl_unique[n_level] = dict()
        for n_group, selection in lvl_select[n_level].items():
          vectors = contacts.select(selection)
          unique = flex.bool(sum((vectors-vec).dot()<1)==1 for vec in vectors)
          lvl_unique[n_level][n_group] = unique
 
      for chain in self.chains:
        atoms = af_shared_atom(atom_type() for _ in contacts)
        atoms.set_xyz(chain.apply_operations(contacts))
        for n_level, groups in lvl_select.items():
          for n_group, sel in groups.items():
            mod = chain.get_modifier(n_level, n_group)
            unq = lvl_unique[n_level][n_group]
            for atom, is_contact in zip(atoms.select(sel), unq):
              mod.add_atom(atom, contact = is_contact)

      if contacting_groups:
        contacting_groups = sum(map(tuple, 
                              list(zip(*sorted(contacting_groups.items())))[1]
                            ),())
      else:
        contacting_groups = tuple()
 
      for n_level, level in mod_hierarchy.items():
        self.lvl_conlens[n_level] = dict()
        self.lvl_dstlist[n_level] = dict()
        inverse    = dict()
        sc_indices = dict()
        sc_translt = dict()
        for n_group in level:
          inverse[n_group] = []
          sel = lvl_select[n_level][n_group]
          unq = lvl_unique[n_level][n_group]
          con = contacts.select(sel).select(unq)
          cg  = [contacting_groups[i] for i,u in zip(sel,unq) if u]
          bad = []
          for n, (j,sym,dist) in enumerate(cg):
            j     = next(g for l,g in group_supers[(max_level, j)] if l == n_level)
            sel_j = lvl_select[n_level][j]
            unq_j = lvl_unique[n_level][j]
            con_j = contacts.select(sel_j).select(unq_j)
            j_inverse = (cell.orthogonalize(
                         (sym.r().as_double() * cell.fractionalize(con_j))\
                        + sym.t().as_double()
                        ) - con[n]).dot()
            nearest = (j_inverse == min(j_inverse)).iselection()[0] if j_inverse else -1
            inverse[n_group].append(nearest)
            if not j_inverse or min(j_inverse) > 1:
              bad.append(n)
 
          sc_indices[n_group] = []
          sc_translt[n_group] = []
          for n, chain in enumerate(self.chains):
            indices = []
            translt = []
            for (j, sym, dist), r in zip(cg, inverse[n_group]):
              j   = next(g for l,g in group_supers[(max_level, j)] if l == n_level)
              idx = ((umx,)+sym_sorted).index(sym)
              indices.append((j, ((n,)+self.indices[n])[idx], r))
              translt.append(([(0.,0.,0.)]+self.translt[n])[idx])
            sc_indices[n_group].append(indices)
            sc_translt[n_group].append(translt)
          self.lvl_conlens[n_level][n_group] = len(cg)
          self.lvl_dstlist[n_level][n_group] = [dist[0][:1 if n in bad else None]
                                                for n, (j, sym, dist) in enumerate(cg)]

        self.lvl_indices[n_level] = sc_indices
        self.lvl_translt[n_level] = sc_translt
        group_coms = lvl_groupcoms[n_level]
        grpsort = lambda i: (flex.vec3_double((cell.orthogonalize(group_coms[i]),))
                            - xyz.mean()).dot()[0]
        self.lvl_grpsort[n_level] = sorted(group_coms, key=grpsort)

    else:
      group_supers = dict()
      group_lookup = dict()
      prev_level = None
      mod_hierarchy = chains[0].mod_hierarchy
      for n_level, level in sorted(mod_hierarchy.items()):
        group_coms = dict()
        for n_group, modifier in level.items():
          curr = (n_level, n_group)
          for atom in modifier.working_atoms():
            if curr not in group_supers:
              prev = (prev_level, group_lookup.get(atom.i_seq))
              group_supers[curr] = [curr] + group_supers.get(prev,[])
            group_lookup[atom.i_seq] = n_group
          group_coms[n_group] = cell.fractionalize(modifier.get_com())
        contacting_groups = {n:dict() for n in level.keys()}
        for i, pair_sym_dict in enumerate(base_sym_table):
          own_group = group_lookup.get(i)
          for j, sym_ops in pair_sym_dict.items():
            other_group = group_lookup.get(j)
            for sym in sym_ops:
              if (own_group and other_group) \
               and not (own_group == other_group and sym.is_unit_mx()):
                a = cell.orthogonalize(
                      (sym.r().as_double() * cell.fractionalize(xyz[j:j+1]))
                      + sym.t().as_double()
                    )
                b = xyz[i]
                dist = (flex.sqrt((a-b).dot())[0], ((a+b)*0.5)[0])
                if (other_group, sym) not in contacting_groups.get(own_group):
                  contacting_groups[own_group][(other_group, sym)] = []
                contacting_groups[own_group][(other_group, sym)].append(dist)
        umx = type(sym)()
        for i, groups in contacting_groups.items():
          entry = [(j,sym,list(zip(*dist))) for (j,sym),dist in groups.items()]
          entry.sort(key=lambda x:(((umx,)+sym_sorted).index(x[1]),x[0]))
          contacting_groups[i] = entry
        contacts = dict()
        for i, entry in contacting_groups.items():
          com = group_coms[i]
          if use_com_midpoint:
            contacts[i] = flex.vec3_double([
              cell.orthogonalize(
                (((sym.r().as_double() * flex.vec3_double((group_coms[j],)))
                + sym.t().as_double() - com) * 0.5 + com)[0]
              ) for j,sym,dist in entry
            ])
          else:
            contacts[i] = flex.vec3_double([
              flex.vec3_double(dist[1]).mean() for j,sym,dist in entry
            ])
        contact_subs = dict()
        for i, conts in contacts.items():
          contact_subs[i] = {
            n_other: [
              sorted(other,
#                key = lambda j:(flex.vec3_double((other[j].get_com(),))-cont).dot()[0]
                key = lambda j:min((other[j].working_atoms().extract_xyz()-cont).dot())
              )[0] for cont in conts
            ] for n_other, other in mod_hierarchy.items() if n_other > n_level
          }
        inverse = dict()
        for i, entry in contacting_groups.items():
          inverse[i] = []
          for n, (j, sym, dist) in enumerate(entry):
            j_inverse = cell.orthogonalize(
                         (sym.r().as_double() * cell.fractionalize(contacts[j]))\
                        + sym.t().as_double()
                        ) - contacts[i][n]
            inverse[i].append((j_inverse.dot() < 1.).iselection()[0])
        for chain in self.chains:
          for i, entry in contacts.items():
            atoms        = af_shared_atom(atom_type() for _ in entry)
            atoms.set_xyz(chain.apply_operations(entry))
            for sup in group_supers[(n_level, i)]:
              is_contacts = (n_level, i) == sup
              chain.get_modifier(*sup).add_atoms(atoms, contacts = is_contacts)
            for lvl, sub in contact_subs[i].items():
              for atom, other in zip(atoms, sub):
                chain.get_modifier(lvl, other).working_atoms().append(atom)
        sc_indices = dict()
        sc_translt = dict()
        for i, entry in sorted(contacting_groups.items()):
          sc_indices[i] = []
          sc_translt[i] = []
          for n, chain in enumerate(self.chains):
            indices = []
            translt = []
            for (j, sym, dist), r in zip(contacting_groups[i], inverse[i]):
              idx = ((umx,)+sym_sorted).index(sym)
              indices.append((j, ((n,)+self.indices[n])[idx], r))
              translt.append(([(0.,0.,0.)]+self.translt[n])[idx])
            sc_indices[i].append(indices)
            sc_translt[i].append(translt)
        
        self.lvl_indices[n_level] = sc_indices
        self.lvl_translt[n_level] = sc_translt
        self.lvl_conlens[n_level] = {i: len(entry)
                                     for i,entry in contacting_groups.items()}
        self.lvl_dstlist[n_level] = {i: [dist[0] for j, sym, dist in entry]
                                     for i,entry in contacting_groups.items()}
        grpsort = lambda i: (flex.vec3_double((cell.orthogonalize(group_coms[i]),))
                            - xyz.mean()).dot()[0]
        self.lvl_grpsort[n_level] = sorted(group_coms, key=grpsort)
 
        prev_level = n_level

    # Create template hierarchy
    hierarchy     = base_hierarchy.deep_copy()
    model         = hierarchy.models()[0]
    chains        = model.chains()
    chain_id_iter = chain_id_generator(start=1)
    i_seq, env    = list(self.sorted_environments.items())[0]
    chain_ids     = []
    for chain in chains: chain.id = 'A'
    for asu in env:
      chain_id = next(chain_id_iter)
      chain_ids.append(chain_id)
      for chain in chains:
        new_chain    = chain.detached_copy()
        new_chain.id = chain_id
        model.append_chain(new_chain)
    hierarchy.atoms().reset_serial()
    hierarchy.atoms().set_xyz(
      self.get_environment_sites_cart(i_seq)
    )
    self.template = hierarchy.deep_copy()

    # Increase cell size to ensure vacuum
    new_cell = tuple(ax * 10 for ax in cell.parameters()[:3]) + cell.parameters()[3:]
    cs_environment = self.crystal_symmetry.customized_copy(
                       unit_cell = new_cell
                     )

    # Get restraints manager
    if need_restraints:
      self.restraints_manager     = model_manager(
                                    model_input      = None,
                                    pdb_hierarchy    = hierarchy,
                                    crystal_symmetry = cs_environment,
                                    log              = open(os.devnull,'w')
                                    ).get_restraints_manager()
    else:
      self.restraints_manager = None
    
    # Write files with mean positions, environments and intermolecular dashes
    if write_environment:
      with open('environment.pdb', 'w') as writer:
        write = lambda string='': writer.write('REMARK   3 {}\n'.format(string))
        write()
        write('ENVIRONMENT DETAILS.')
        write('  DISTANCE CUTOFF: {:.2f} A'.format(dist_cutoff))
        write('  ATOMS INCLUDED: {}'.format('RIGID AND '*rigid_only + ('HBOND ONLY'
              if hbond_only else 'HEAVY ONLY' if heavy_only else 'PROTEIN ONLY' if
              protein_only else 'ALL')))
        write('  CHAIN A NEIGHBOURS: {}'.format(len(chain_ids)))
        write()
        for n, chain_id in enumerate(chain_ids):
          write(' CHAIN '+chain_id)
          write('  SHORT CONTACTS WITH CHAIN A: {}'.format(len(self.dstlist[n])))
          write('  COM-COM DISTANCE TO CHAIN A: {:.2f}'.format(self.comdist[n]))
        write()

      hier_copy = hierarchy.deep_copy()
      model     = hier_copy.models()[0]

      for n_level, level in self.chains[0].mod_hierarchy.items():
        chain    = type(model.chains()[0])()
        chain.id = next(chain_id_iter)
        model.append_chain(chain)
        res_iter = chain_id_generator(start=0)
        int_iter = int_generator(start=1)
        for n_mod, mod in level.items():
          resgr = type(model.chains()[0].residue_groups()[0])()
          atmgr = type(model.chains()[0].residue_groups()[0].atom_groups()[0])()
          resgr.resseq  = '{:4d}'.format(next(int_iter))
          atmgr.resname = next(res_iter)
          chain.append_residue_group(resgr)
          resgr.append_atom_group(atmgr)
          for atom in mod.get_contacts():
            atom = atom.detached_copy()
            atom.name    = 'DUM'
            atom.element = 'X'
            atmgr.append_atom(atom)


      for chain in hier_copy.chains():
        chain.atoms().reset_serial()

      hier_copy.write_pdb_file('environment.pdb',
                               crystal_symmetry=base_symmetry,
                               open_append=True)

  def get_contact(self, n_level):

      nan = float('nan')
      mods = self.chains[0].mod_hierarchy[n_level]
      shape = (max(mods)+1,
               len(self.chains),
               max(mod.get_contacts().extract_xyz().size() for mod in mods.values()),
               3)
      result = [ [[(nan,)*shape[3]] * shape[2] for _ in range(shape[1])]
                 for _ in range(shape[0]) ]

      for n, chain in enumerate(self.chains):
        modifiers = chain.mod_hierarchy[n_level]
        for i, mod in modifiers.items():
          xyz = list(mod.get_contacts().extract_xyz())
          result[i][n][:len(xyz)] = xyz

      return result

  def get_conshft(self, n_level):

      result = self.get_contact(n_level)

      for n, chain in enumerate(self.chains):
        modifiers = chain.mod_hierarchy[n_level]
        for i, mod in modifiers.items():
          xyz = list(mod.get_contact_shifts())
          result[i][n][:len(xyz)] = xyz

      return result

  def get_indices(self, n_level):

      return self.lvl_indices[n_level]

  def get_translt(self, n_level):

      return self.lvl_translt[n_level]

  def get_conlens(self, n_level):

      n_contacts = self.lvl_conlens[n_level]
      return [n_contacts.get(n, 0) for n in range(max(n_contacts)+1)]

  def get_dstlist(self, n_level):

      return self.lvl_dstlist[n_level]
    
  def get_grpsort(self, n_level):

      return self.lvl_grpsort[n_level]
    
  def get_environment_sites_cart(self, i_seq):

    sites_cart = self.chains[i_seq].coordinates()
    for j_seq, symm_op in self.sorted_environments[i_seq]:
      trans_j = self.crystal_symmetry.unit_cell().orthogonalize(
                  symm_op.t().as_double()
                )
      sites_cart.extend(self.chains[j_seq].coordinates() + trans_j)
    return sites_cart
     
  def regularize_chain_in_environment(self, chain, special=None):

    sites_cart = lbfgs(
      sites_cart                  = self.get_environment_sites_cart(chain.n_chain),
      geometry_restraints_manager = self.restraints_manager.geometry,
      special                     = special,
    )
    chain.set_coordinates(sites_cart[:chain.atoms.size()])

  def regularize_all(self):

    for chain in self.chains:
      self.regularize_chain_in_environment(chain)

class Chain():

  def __init__(self,
               tls_hierarchy,
               atoms,
               translation,
               symm_op,
               parent_cell,
               n_chain            = None,
               restraints_manager = None,
               residues           = {-1}):
   
    self.tls_hierarchy      = tls_hierarchy
    self.n_chain            = n_chain
    self.rank               = n_chain
    self.atoms              = atoms
    self.init_coordinates   = self.atoms.extract_xyz()
    self.restraints_manager = restraints_manager
    self.residues           = residues
    
    # define positioning parameters
    self.is_in_position    = False
    self.translation       = translation
    self.symm_op           = symm_op
    self.parent_cell       = parent_cell
    self.sc_trans          = self.parent_cell.orthogonalize(self.translation)
    self.symm_trans_double = self.symm_op.t().as_double()
    self.symm_trans        = self.parent_cell.orthogonalize(self.symm_trans_double)
    self.symm_rot_double   = self.symm_op.r().as_double()
    self.symm_rot_inverse  = self.symm_op.r().inverse().as_double()

    self.build_chain()

  def build_chain(self):
    
    self.mod_hierarchy = {}
    for n_level, level in self.tls_hierarchy.items():
      self.mod_hierarchy[n_level] = {}
      for n_group, tls_obj in level.items():
        atoms     = self.working_atoms().select(tls_obj.selection)
        modifier  = Modifier(working_atoms   = atoms,
                             tls_obj         = tls_obj,
                             parent          = self)
        self.mod_hierarchy[n_level][n_group] = modifier
        modifier.is_active = (-1 in self.residues or 
                              bool(self.residues.intersection(modifier.residues)))

  def get_modifier(self, n_layer, n_group):

    return self.mod_hierarchy[n_layer][n_group]

  def clear_shifts(self):

    for n_level, level in self.mod_hierarchy.items():
      for n_group, modifier in level.items():
        if hasattr(modifier, 'shifts'):
          modifier.shifts *= 0.

  def adopt_from(self, other):

    for n_level, level in self.mod_hierarchy.items():
      for n_group, modifier in level.items():
        if other is None:
          modifier.set_shifts(modifier.get_shifts() * 0.)
        else:
          modifier.set_shifts(other.get_modifier(n_level, n_group).get_shifts())

  def reverted_coordinates(self):

    return self.apply_operations(self.init_coordinates)

  def coordinates(self):

    return self.working_atoms().extract_xyz()

  def bonded_deltas(self):

    return self.restraints_manager.geometry.pair_proxies().bond_proxies.deltas(
             self.coordinates()
           )

  def calc_sd(self):

    return (self.base_coordinates() - self.init_coordinates).dot()

  def trans(self):

    trans = (self.base_coordinates() - self.init_coordinates).mean()
    if self.is_in_position:
      trans = tuple((self.symm_rot_double * flex.vec3_double((trans,)))[0])
    return trans

  def rotmat(self):

    base = self.base_coordinates()
    init = self.init_coordinates
    sample = flex.random_selection(base.size(), 3)
    base_sample = matrix.sqr((base.select(sample) - base.mean()).as_double()).transpose()
    init_sample = matrix.sqr((init.select(sample) - init.mean()).as_double()).transpose()
    rotmat = (base_sample * init_sample.inverse())
    return rotmat

  def rotvec(self):
    
    rotmat = self.rotmat()
    ang, ax = rotmat.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle()
    rotvec = tuple(ang * ax)
    if self.is_in_position:
      rotvec = tuple((self.symm_rot_double * flex.vec3_double((rotvec,)))[0])
    return rotvec

  def base_coordinates(self):

    coordinates = self.coordinates()
    if self.is_in_position: coordinates = self.apply_reverse_operations(coordinates)
    return coordinates

  def set_coordinates(self, coordinates):

    self.working_atoms().set_xyz(coordinates)

  def revert(self):

    self.set_coordinates(self.init_coordinates)
    if self.is_in_position:
      self.is_in_position = False
      self.move_in_position()

  def working_atoms(self):

    return self.atoms

  def set_random(self, amp=1., n_level_select=-1):

    for n_level, level in self.mod_hierarchy.items():
      if n_level_select in (-1, n_level):
        for modifier in level.values():
          if modifier.is_active:
            modifier.set_random(amp=amp)

  def apply_operations(self, coordinates):
    
    new_coordinates = (self.symm_rot_double * coordinates
                      + self.symm_trans + self.sc_trans)
    return new_coordinates

  def apply_reverse_operations(self, coordinates):

    new_coordinates = (self.symm_rot_inverse * 
                      (coordinates - self.symm_trans - self.sc_trans))
    return new_coordinates

  def cancel_rotations(self):

    was_in_position = self.is_in_position
    self.move_back()
    com = self.coordinates().mean()
    old_com = self.init_coordinates.mean()
    shift = [a-b for a,b in zip(com, old_com)]
    self.set_coordinates(self.init_coordinates + shift)
    if was_in_position: self.move_in_position()

  def move_in_position(self):
    
    if not self.is_in_position:
      coordinates         = self.coordinates()
      new_coordinates     = self.apply_operations(coordinates)
      self.is_in_position = True
      self.set_coordinates(new_coordinates)

  def move_back(self):

    if self.is_in_position:
      coordinates         = self.coordinates()
      new_coordinates     = self.apply_reverse_operations(coordinates)
      self.is_in_position = False
      self.set_coordinates(new_coordinates)

  def regularize(self, environments=None, keep_shifts=False, special=None):

    if keep_shifts:
      init_coordinates = self.base_coordinates()

    if environments is not None:
      environments.regularize_chain_in_environment(self, special)
    else:
      sites_cart = lbfgs(
                     sites_cart                  = self.coordinates(),
                     geometry_restraints_manager = self.restraints_manager.geometry,
                     special                     = special,
                   )
      self.set_coordinates(sites_cart)
    
    if keep_shifts:
      return self.base_coordinates() - init_coordinates

def lbfgs(sites_cart,
          geometry_restraints_manager,
          special                     = None,
          geometry_restraints_flags   = None,
          lbfgs_termination_params    = None):

  if geometry_restraints_flags is None:
    if special is None:
      geometry_restriants_flags = geometry_restraints.flags.flags(
        bond        = True,
        nonbonded   = True,
        angle       = True,
        dihedral    = True,
        chirality   = True,
        planarity   = True,
        parallelity = True
      )
    elif special == 'nonbonded':
      geometry_restriants_flags = geometry_restraints.flags.flags(
        bond        = False,
        nonbonded   = True,
        angle       = False,
        dihedral    = False,
        chirality   = False,
        planarity   = False,
        parallelity = False
      )
  if lbfgs_termination_params is None:
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = 25
    )
  refinement.geometry_minimization.lbfgs(
    sites_cart                         = sites_cart,
    geometry_restraints_manager        = geometry_restraints_manager,
    geometry_restraints_flags          = geometry_restraints_flags,
    lbfgs_termination_params           = lbfgs_termination_params,
    correct_special_position_tolerance = False
  )
  return sites_cart

class Modifier():

  def __init__(self, tls_obj, working_atoms, parent):

    self.atoms  = working_atoms
    self.origin = tls_obj.origin
    self.r      = tls_obj.r
    self.eps    = tls_obj.eps
    self.amp    = tls_obj.amp
    self.parent = parent
    self.shifts = None
    self.size   = working_atoms.size()
    self.residues = set(atom.fetch_labels().resseq_as_int() for atom in self.atoms)
    self.contacts = flex.size_t()
    self.is_active = True

  def set_random(self, amp=1.):

    coordinates       = self.base_coordinates()
    sites_cart        = coordinates - self.origin
    d0                = matrix.col(tls_as_xyz.step_i__get_dxdydz(self.r, None)) * self.amp * amp
    d_r_M_V           = tls_as_xyz.formula_49(self.r, None) * self.amp * amp
    new_coordinates   = apply_tls_shifts(
      sites_cart      = sites_cart,
      R_ML_transposed = self.r.R_ML.transpose(),
      R_ML            = self.r.R_ML,
      d0              = d0,
      d_r_M_V         = d_r_M_V,
      s               = matrix.col((self.r.sx,self.r.sy,self.r.sz)),
      wy_lx           = self.r.wy_lx,
      wz_lx           = self.r.wz_lx,
      wz_ly           = self.r.wz_ly,
      wx_ly           = self.r.wx_ly,
      wx_lz           = self.r.wx_lz,
      wy_lz           = self.r.wy_lz,
      origin          = self.origin)
    self.shifts       = new_coordinates - coordinates
    self.set_coordinates(new_coordinates)

  def working_atoms(self):

    return self.atoms

  def get_shifts(self):

    return self.shifts or flex.vec3_double(self.atoms.size())

  def set_shifts(self, shifts):

    coordinates     = self.base_coordinates()
    new_coordinates = coordinates + shifts
    self.shifts     = shifts
    self.set_coordinates(new_coordinates)
 
  def base_coordinates(self):
 
    coordinates = self.working_atoms().extract_xyz()
    if self.parent.is_in_position:
      coordinates = self.parent.apply_reverse_operations(coordinates)
    if self.shifts:
      coordinates = coordinates - self.shifts
    return coordinates
  
  def set_coordinates(self, coordinates):
    
    if self.parent.is_in_position:
      coordinates = self.parent.apply_operations(coordinates)
    self.working_atoms().set_xyz(coordinates)

  def add_atoms(self, atoms, contacts = False):

    self.working_atoms().extend(atoms)
    if contacts:
      self.contacts.extend(
        flex.size_t_range(self.atoms.size() - atoms.size(), self.atoms.size())
      )

  def add_atom(self, atom, contact = False):

    self.working_atoms().append(atom)
    if contact:
      self.contacts.append(self.working_atoms().size()-1)

  def get_contacts(self):

    return self.working_atoms().select(self.contacts)

  def get_contact_shifts(self):

    return self.get_shifts().select(self.contacts)

  def get_com(self):

    return self.working_atoms()[:self.size].extract_xyz().mean()

class Collector():

  def __init__(self, base_hierarchy, cryst_symm, miller_set=None):

    self.base_hierarchy   = base_hierarchy
    self.cryst_symm       = cryst_symm
    self.cache            = base_hierarchy.atom_selection_cache()

    if miller_set: self.init_f_collector(miller_set)
    self.init_msd_collector(self.cache)
    self.init_shifts_collector(self.cache)
    self.hierarchy        = None

  def init_f_collector(self, miller_set):

    array_size            = miller_set.size()
    self.f_collector      = mp.Array('d', 2 * array_size)
    self.i_collector      = mp.Array('d', array_size)
    self.miller_set       = miller_set
    self.f_count          = mp.Value('i', 0)

  def init_msd_collector(self, cache):

    self.msd_collector    = mp.Array('d', cache.n_seq)
    self.msd_count        = mp.Value('i', 0)
    self.msd_total        = mp.Value('i', 0)

  def init_shifts_collector(self, cache):

    self.shifts_collector = mp.Array('d', cache.n_seq)
    self.shifts_count     = mp.Value('i', 0)

  def f_collect(self, f_model):

    with self.f_count.get_lock():
      f_buffer              = np.frombuffer(self.f_collector.get_obj())
      i_buffer              = np.frombuffer(self.i_collector.get_obj())
      f_model               = f_model.data().as_numpy_array()
      
      f_buffer             += f_model.view('double')
      f_model              *= f_model.conjugate()
      i_buffer             += f_model.real
      self.f_count.value   += 1

    return self.f_count.value

  def msd_collect(self, msd, count=1):

    with self.msd_count.get_lock():
      msd_buffer            = np.frombuffer(self.msd_collector.get_obj())
      msd_buffer           += msd
      self.msd_count.value += 1
      self.msd_total.value += count

    return self.msd_count.value

  def shifts_collect(self, shifts, count=1):

    with self.shifts_count.get_lock():
      shifts_buffer            = np.frombuffer(self.shifts_collector.get_obj())
      shifts_buffer           += shifts
      self.shifts_count.value += count

    return self.shifts_count.value

  def pdb_collect(self, result):

    if self.hierarchy is None:
      self.hierarchy = result.working_hierarchy
      self.sc_symm   = result.sc_symm
    else:
      for model in result.working_hierarchy.models():
        self.hierarchy.append_model(model.detached_copy())

  def write_b(self, b_out='ensemble_b.json', ca_only=True):
   
    with self.msd_count.get_lock():
      if self.msd_count.value:
        # Set selection and atom sequence based on C-alpha boolean, filter out conformers
        if ca_only:
          selection = self.cache.selection('name CA and (altid " " or altid A)')
        else:
          selection = self.cache.selection('(altid " " or altid A)')
    
        msd = flex.double(self.msd_collector).select(selection)
        msd = 1./self.msd_total.value * msd
      
        sequence  = map(int, self.cache.resid_list.select(selection))
        # Convert to B-factors
        b_factors = msd_to_b(msd)
        # Write json file for plotting
        write_b_factors_to_json(sequence  = tuple(sequence),
                                b_factors = tuple(b_factors),
                                file_name = b_out)

  def write_shifts(self,
                   shifts_out = 'regularization_shifts.pdb',
                   tls_colors = 'tls_colors.pml'):

    with self.shifts_count.get_lock():
      if self.shifts_count.value:
        hierarchy = self.base_hierarchy.deep_copy()
        shifts_total = flex.double(self.shifts_collector)
        hierarchy.atoms().set_b(msd_to_b(shifts_total/self.shifts_count.value))
        hierarchy.write_pdb_file(shifts_out,
                                 crystal_symmetry = self.cryst_symm,
                                 open_append      = False)

  def write_mtz(self, mtz_out='supercell_avg.mtz'):

    with self.f_count.get_lock():
      if self.f_count.value:
        f_sq_avg   = flex.double(self.i_collector.get_obj()) / self.f_count.value
        f_complex  = np.frombuffer(self.f_collector.get_obj()).view('complex').copy()
        f_avg_sq   = flex.double(
                       (f_complex*f_complex.conjugate()).real/self.f_count.value**2
                     )
        
        # I_total  = <F^2>
        itot_array = self.miller_set.array(f_sq_avg)
        # I_Bragg  = <F>^2
        ibrg_array = itot_array.customized_copy(data = f_avg_sq)
        # I_diff   = <F^2> - <F>^2
        idff_array = itot_array.customized_copy(data = f_sq_avg - f_avg_sq)
       
        dataset    = itot_array.as_mtz_dataset(column_root_label = 'ITOT',
                                               column_types      = 'J')
        dataset.add_miller_array(ibrg_array,
                                 column_root_label = 'IBRG',
                                 column_types      = 'J')
        dataset.add_miller_array(idff_array,
                                 column_root_label = 'IDFF',
                                 column_types      = 'J')
       
        dataset.mtz_object().write(mtz_out)

  def write_pdb(self, pdb_out='supercell_out.pdb', mmcif=False):

    if self.hierarchy:
      if not mmcif:
        self.hierarchy.write_pdb_file(pdb_out,
                                      crystal_symmetry = self.sc_symm)
      else:
        self.hierarchy.write_mmcif_file(pdb_out.replace('.pdb','.cif'),
                                        crystal_symmetry = self.sc_symm)

class TLS():

  def __init__(self, t, l, s, origin, selection_string, selection, amp=1.,
               zero_trace=False, mult=(1.,1.,1.), eps=1.e-6, log=None):

    if zero_trace: s[0]=s[4]=s[8]=0

    self.log               = log
    self.eps               = eps
    self.T                 = matrix.sym(sym_mat3=t) * mult[0]
    self.L                 = matrix.sym(sym_mat3=l)*DEG2RADSQ * mult[1]
    self.S                 = matrix.sqr(s)*DEG2RAD * mult[2]
    self.amp               = np.sqrt(amp)
    self.origin            = origin
    self.selection_string  = selection_string
    self.selection         = selection

    try:
      self.r = analysis.run(T              = self.T,
                            L              = self.L,
                            S              = self.S,
                            implementation = 'c++',
                            log            = self.log).self_check()
    except:
      self.r = analysis.run(T              = self.T,
                            L              = self.L,
                            S              = self.S,
                            log            = self.log).self_check()
  
