from __future__ import division, print_function
from cctbx import miller
import iotbx
import sys
import glob
from . import phil
from . import schimpy_classes as tls
import multiprocessing as mp
import time
import gc
import psutil
import traceback

def run(args):

  scope       = phil.phil_parse(args = args)
  if not args: scope.show(attributes_level=2); return
  p           = scope.extract().schimpy
  if p.input.tls_in:
    tls_in    = [file_name for expression in p.input.tls_in.split()
                           for file_name  in glob.glob(expression)]
  else:
    tls_in    = [p.input.pdb_in]

  print('Initializing')
  max_level  = max(p.input.max_level, p.correlation.contact_level)
  manager    = tls.Manager(tls_in      = tls_in,
                           pdb_in      = p.input.pdb_in,
                           tls_origin  = p.input.tls_origin,
                           mult        = p.input.tls_multipliers,
                           zero_trace  = p.input.zero_trace,
                           max_level   = max_level,
                           min_level   = p.input.min_level,
                           sc_size     = p.input.sc_size,
                           seed_offset = p.input.seed,
                           reset_b     = p.input.reset_b)

  miller_set = miller.build_set(crystal_symmetry = manager.sc_symm,
                                d_min            = p.params.resolution,
                                anomalous_flag   = False)

  collector  =    tls.Collector(miller_set       = miller_set,
                                base_hierarchy   = manager.base_hierarchy,
                                cryst_symm       = manager.cryst_symm)

  pool       =          mp.Pool(processes        = p.input.processes,
                                initializer      = initializer,
                                initargs         = (collector, p, manager, mp.Value('d')),
                                maxtasksperchild = 1, # Force memory cleanup after every simulation
                                )
  
  try:
    for result in pool.imap_unordered(func     = model_worker,
                                      iterable = range(max(p.input.n_models))):
      
      n_model = result.n_model
      if hasattr(result, 'working_hierarchy'):
        print('Collecting pdb {}'.format(n_model))
        collector.pdb_collect(result)
      print("Model {} finished".format(n_model))

      # Write stuff if appropriate
      if p.output.write_b and result.msd_collected in p.input.n_models:
        print('Writing ensemble B {}'.format(result.msd_collected))
        file_name = file_number(p.output.b_out, result.msd_collected)
        collector.write_b(file_name)
      if p.output.write_mtz and result.f_collected in p.input.n_models:
        print('Writing MTZ {}'.format(result.f_collected))
        file_name = file_number(p.output.mtz_out, result.f_collected)
        collector.write_mtz(file_name)
      if collector.shifts_count.value and result.shifts_collected in p.input.n_models:
        print('Writing regularization shifts {}'.format(result.shifts_collected))
        file_name = file_number('regularization_shifts.pdb', result.shifts_collected)
        collector.write_shifts(file_name)

    pool.close()
    pool.join()

  except Exception as e:
    print('Oopsie, something went wrong!')
    print(''.join(traceback.format_exception(*sys.exc_info())))

    pool.terminate()
    pool.join()

    if p.output.write_b:
      print('Writing ensemble B')
      collector.write_b(p.output.b_out)
    if p.output.write_mtz:
      print('Writing MTZ')
      collector.write_mtz(p.output.mtz_out)
    if collector.shifts_count.value:
      print('Writing regularization shifts')
      collector.write_shifts()
  
  finally:
    if p.output.write_pdb:
      as_mmcif = p.output.write_mmcif and not p.output.superpose
      print('Writing {}'.format(['PDB','mmCIF'][as_mmcif]))
      collector.write_pdb(p.output.pdb_out, mmcif=as_mmcif)
    n_models = max([collector.f_count.value, collector.msd_count.value])
    print('Processed {} models, exiting'.format(n_models))

def file_number(name, number):
  file_split    = name.split('.')
  file_split[0] = file_split[0]+'_{}'.format(number)
  file_name     = '.'.join(file_split)
  return file_name

def initializer(collect, params, man, sttime):

  global collector, p, manager, start_time
  collector, p, manager, start_time = collect, params, man, sttime

def model_worker(n_model):

  keep = p.output.write_pdb > n_model or p.output.write_pdb < 0

  try:
    # Enforce start intervals
    with start_time.get_lock():
      time.sleep(max(0.1, start_time.value - time.time()))
      start_time.value = time.time() + p.input.interval
    # Start building
    print('Building model {}'.format(n_model))
    model = manager.new_model(n_model)
    model.build_model(method=p.input.strategy,
                      avail_coords=manager.avail_coords,
                      amp=p.input.amplitude,
                      reverse=p.input.reverse_levels,
                      residues=p.input.residues,
                      correlate=p.correlation.correlate,
                      correlate_fancy=p.correlation.correlate_fancy,
                      correlate_sequential=p.correlation.sequential,
                      correlate_insideout=p.correlation.insideout,
                      include_shifts=p.correlation.shifts,
                      include_rotvecs=p.correlation.rotvecs,
                      weight_exp=p.correlation.weights,
                      spring_weights=p.correlation.force_power,
                      use_contacts=p.correlation.use_contacts,
                      correlate_after=p.correlation.correlate_after,
                      use_com_midpoint=p.correlation.contact_midpoint,
                      inherit_contacts=p.correlation.inherit,
                      contact_level=p.correlation.contact_level,
                      allow_bad_frac=p.correlation.allow_bad_frac,
                      disallow_good_frac=p.correlation.disallow_good_frac,
                      pairs_frac=p.correlation.swap_frac,
                      require_both=p.correlation.require_both,
                      uncorrelate=p.correlation.uncorrelate,
                      percentile=p.correlation.energy_percentile,
                      skip=p.correlation.skip_levels,
                      override=p.correlation.override_weights,
                      max_level=p.input.max_level,
                      rigid_only=p.environments.rigid_only,
                      protein_only=p.environments.protein_only,
                      heavy_only=p.environments.heavy_only,
                      hbond_only=p.environments.hbond_only,
                      dist_cutoff=p.environments.cutoff,
                      stretch=p.environments.stretch,
                      regularize=p.regularization.regularize,
                      scope=p.regularization.scope,
                      )
 
    if p.output.write_b:
      print('Calculating mean square displacements of {}'.format(n_model))
      msd = model.calc_msd()
      model.msd_collected = collector.msd_collect(msd, model.n_chains)
    
    if p.output.write_mtz:

      while (psutil.virtual_memory().available / 1024 ** 3 < p.params.fft_mem):
        time.sleep(5)

      print('Calculating f_model {}'.format(n_model))
      f_model = model.calc_f_model(high_resolution = p.params.resolution,
                                   k_sol           = p.params.k_sol,
                                   b_sol           = p.params.b_sol,
                                   b_cart          = list(map(float,p.params.b_cart)))
      print('Collecting f_model {}'.format(n_model))
      model.f_collected = collector.f_collect(f_model)

    if hasattr(model, 'shifts_total'):
      print('Collecting shifts {}'.format(n_model))
      model.shifts_collected = collector.shifts_collect(model.shifts_total, model.n_chains)

    if n_model == 0 and p.input.max_level > 0:
      model.write_com_displacements()
      model.write_all_environments()

    if keep: model.add_coms_to_hierarchy()
    if p.output.superpose: model.chains_as_models()

    model.close(keep_hierarchy = keep)

    gc.collect()

    return model

  except KeyboardInterrupt as e:
    return e
  except Exception as e:
    print('Oopsie, something went wrong!')
    print(''.join(traceback.format_exception(*sys.exc_info())))
    return e
  finally:
    model.close(keep_hierarchy = keep)
    return model

