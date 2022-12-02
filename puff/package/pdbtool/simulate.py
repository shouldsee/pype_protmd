#!/usr/bin/env python

import os
import shutil
import copy 

import pdbstruct
import util
import asa

import namd
import amber
import gromacs


backbone_atoms = [ \
    "OXT" "H", "H1", "H2", "H3","HN", "N", "C", "O", "CA", "HA"]

    
def make_asa_pdb(
    in_pdb, out_pdb, backbone_atoms=['N','C','O'], k=1.0):
  mol = pdbstruct.Molecule(in_pdb)
  indices = []
  pdbstruct.add_radii(mol.atoms())
  asas = asa.calculate_asa(mol.atoms(), 1.4)
  for i, atom in enumerate(mol.atoms()):
    if asas[i] > 9.0 or atom.type in backbone_atoms:
      atom.bfactor = k
    else:
      atom.bfactor = 0.0
  mol.write_pdb(out_pdb)


def make_backbone_pdb(
    in_pdb, out_pdb, backbone_atoms=['N','C','O'], k=1.0):
  mol = pdbstruct.Molecule(in_pdb)
  indices = []
  for atom in mol.atoms():
    atom.bfactor = 0.0
  for i, atom in enumerate(mol.atoms()):
    if atom.type in backbone_atoms:
      atom.bfactor = k
    else:
      atom.bfactor = 0.0
  mol.write_pdb(out_pdb)


# adaptor routines to individual simulation modules
# 2 types of generic constraints are provided
#    1. constraint="backbone" - backbone atom constrained
#    2. constraint="surface" - backbone and surface constrained
# These can be made directly by hand by using a
# constraint_pdb, following the convention of NAMD. 

def make_init_files(force_field, raw_pdb, name, constrain=False):
  "Returns topology and coordinate fnames"
  if force_field in ['NAMD']:
    top = os.path.abspath(name + '.psf')
    crds = os.path.abspath(name + '.pdb')
    if not os.path.isfile(top):
      namd.make_namd_from_pdb(raw_pdb, top, crds)
  elif force_field in ['AMBER']:
    top = os.path.abspath(name + '.top')
    crds = os.path.abspath(name + '.crd')
    if not os.path.isfile(top):
      amber.convert_pdb_to_amber(raw_pdb, top, crds)
  elif force_field in ['GROMACS']:
    top = os.path.abspath(name + '.top')
    crds = os.path.abspath(name + '.gro')
    if not os.path.isfile(top):
      gromacs.prep_solvation(raw_pdb, name, True, constrain)
  else:
    raise ValueError, "Couldn't recognize force-field " + force_field
  return top, crds
    

def make_constraint_pdb(top, crds, constraint, constraint_pdb):
  "Builds constraint pdb file"
  if '.psf' in top:  
    if 'surface' in constraint:
      temp_pdb = util.temp_fname('.pdb')
      temp2_pdb = util.temp_fname('.pdb')
      namd.convert_namd_pdb_to_pdb(crds, temp_pdb)
      make_asa_pdb(temp_pdb, temp2_pdb, backbone_atoms, 20.0)
      namd.convert_pdb_to_namd_pdb(temp2_pdb, constraint_pdb)
      os.remove(temp_pdb)
      os.remove(temp2_pdb)
    elif 'backbone' in constraint:
      temp_pdb = util.temp_fname('.pdb')
      temp2_pdb = util.temp_fname('.pdb')
      namd.convert_namd_pdb_to_pdb(crds, temp_pdb)
      make_backbone_pdb(temp_pdb, temp2_pdb, backbone_atoms, 20.0)
      namd.convert_pdb_to_namd_pdb(temp2_pdb, constraint_pdb)
      os.remove(temp_pdb)
      os.remove(temp2_pdb)
  if '.top' in top:
    temp_pdb = util.temp_fname('.pdb')
    amber.convert_amber_to_pdb(top, crds, temp_pdb)
    if 'surface' in constraint:
      make_asa_pdb(temp_pdb, constraint_pdb, backbone_atoms)
    elif 'backbone' in constraint:
      make_backbone_pdb(temp_pdb, constraint_pdb, backbone_atoms)
    if constraint_pdb:
      shutil.copy(crds, constraint_pdb.replace('.pdb', '.crd'))
    os.remove(temp_pdb)


def fetch_parms(
    top, crds, force_field, constraint_pdb, simulation_type, name):
  if force_field in ['NAMD']:
    md_package = 'namd'
  elif force_field in ['AMBER']:
    md_package = 'amber'
  elif force_field in ['GROMACS']:
    md_package = 'gromacs'

  source_parms = eval('%s.%s_parms' % (md_package, simulation_type))
  parms = source_parms.copy()

  parms['force_field'] = force_field
  parms['topology'] = top
  parms['input_crds'] = crds
  parms['output_name'] = name
  parms['cutoff'] = 12.0
  if constraint_pdb:
    parms['constraint_pdb'] = constraint_pdb
    parms['restraint'] = constraint_pdb

  return parms


def test_files(*args):
  for a in args:
    if not os.path.isfile(a):
      return False
  return True
  
  
def get_restart_files(name):
  top = name + '.top'
  gro = name + '.gro'
  psf = name + '.psf'
  coor = name + '.coor'
  crd = name + '.crd'
  rst = name + '.rst'
  if test_files(top, gro):
    return gromacs.get_restart_files(name)
  elif test_files(psf, coor):
    return namd.get_restart_files(name)
  elif test_files(top, crd) or test_files(top, rst):
    return amber.get_restart_files(name)
  else:
    raise IOError, "can't find restart files %s" % name
  

def merge_simulations(md_name, sim_dirs):
  if not sim_dirs:
    return
  sample_amber_top = '%s/%s.top' % (sim_dirs[0], md_name)
  sample_gro_crds = '%s/%s.gro' % (sim_dirs[0], md_name)
  if os.path.isfile(sample_gro_crds):
    gromacs.merge_simulations(md_name, sim_dirs)
  elif os.path.isfile(sample_amber_top):
    amber.merge_simulations(md_name, sim_dirs)
  else:
    namd.merge_simulations(md_name, sim_dirs)

  # calculate time spent in simulations
  fnames = ['%s/%s.time' % (pulse, md_name) 
            for pulse in sim_dirs
            if os.path.isfile(pulse)]
  vals = [open(f).read().split()[0] for f in fnames]
  time = sum([float(val) for val in vals])
  open(md_name+'.pulse.time', 'w').write(util.elapsed_time_str(time))


    
def write_soup_to_restart_files(soup, name, force_field):
  if force_field in ['NAMD']:
    return namd.write_soup_to_restart_files(soup, name)
  elif force_field in ['AMBER']:
    return amber.write_soup_to_restart_files(soup, name)
  elif force_field in ['GROMACS']:
    return gromacs.write_soup_to_restart_files(soup, name)


def soup_from_restart_files(top, crds, vels, skip_solvent=True):
  if crds.endswith('.gro'):
    return gromacs.soup_from_restart_files(
        top, crds, vels, skip_solvent)
  elif top.endswith('.top'):
    return amber.soup_from_restart_files(top, crds, vels)
  elif top.endswith('.psf'):
    return namd.soup_from_restart_files(top, crds, vels)


def run_parms(parms):
  if parms['force_field'] in 'NAMD':
    namd.run(parms)
  elif parms['force_field'] in 'AMBER':
    amber.run(parms)
  elif parms['force_field'] in 'GROMACS':
    gromacs.run(parms)
  

def minimize(top, crds, force_field, name, constraint_pdb=""):
  parms = fetch_parms(
      top, crds, force_field, constraint_pdb, 'minimization', name)
  run_parms(parms)


def langevin(
    top, crds, force_field, vels, n_step, temp, name, 
    n_step_per_snapshot=50, constraint_pdb=""):
  parms = fetch_parms(
      top, crds, force_field, constraint_pdb, 
      'langevin_thermometer', name)
  parms['input_vels'] = vels
  parms['n_step_dynamics'] = n_step
  parms['n_step_per_snapshot'] = n_step_per_snapshot
  parms['temp_thermometer'] = "%.1f" % temp
  parms['temp_initial'] = str(temp)
  parms['gamma_ln'] = 5.0
  run_parms(parms)


def constant(
    top, crds, force_field, vels, n_step, name, 
    n_step_per_snapshot=50, constraint_pdb=""):
  parms = fetch_parms(top, crds, force_field, constraint_pdb,
                      'constant_energy', name)
  parms['input_vels'] = vels
  parms['n_step_dynamics'] = n_step
  parms['cutoff'] = 12.0
  parms['n_step_per_snapshot'] = n_step_per_snapshot
  run_parms(parms)


def preequilibrate(in_name, out_name, temperature):
  if os.path.isfile(in_name + '.gro'):
    gromacs.preequilibrate(in_name, out_name, temperature)
  elif os.path.isfile(in_name + '.top'):
    amber.preequilibrate(in_name, out_name, temperature)
  elif os.path.isfile(in_name + '.psf'):
    namd.preequilibrate(in_name, out_name, temperature)


def continue_if_file(fname, msg):
  if not os.path.isfile(fname):
    print "Running:", msg
    return True
  print "Skipping:", msg
  return False


def pulse(
    in_name, md_name, force_field, n_step, pulse_fn, 
    n_step_per_pulse=100, constraint_pdb=""):
  """
  Takes as argument, a first order function:
    def pulse_fn(soup):
  that updates the velocities at the beginnig of every pulse
  """
  timer = util.Timer()

  top, crds, vels = get_restart_files(in_name)

  parms = fetch_parms(top, crds, force_field, constraint_pdb, 
                      'constant_energy', md_name)
  parms['input_md_name'] = in_name
  parms['input_vels'] = vels
  parms['n_step_dynamics'] = n_step
  parms['n_step_per_snapshot'] = n_step_per_pulse // 2
  parms['n_step_per_pulse'] = n_step_per_pulse

  config = 'config'
  if not os.path.isfile(config):
    util.write_dict(parms, config)

  if 'GROMACS' in force_field:
    trj = md_name + '.trr'
  elif 'AMBER' in force_field: 
    trj = md_name + '.trj'
  elif 'NAMD' in force_field:
    trj = md_name + '.dcd'
  else:
    raise ValueError, "unrecognized force-filed", force_field

  config = md_name + ".config"

  n_pulse = parms['n_step_dynamics'] / n_step_per_pulse
  n_step_list = [n_step_per_pulse for i in range(n_pulse)]
  n_excess_step = parms['n_step_dynamics'] % n_step_per_pulse
  if n_excess_step > 0:
    n_pulse += 1
    n_step_list.append(n_excess_step)
  
  if os.path.isfile(trj):
    print "simulation already run."
    return

  save_dir = os.getcwd()

  pulse_parms = copy.deepcopy(parms)
  pulse_parms['topology'] = os.path.abspath(parms['topology'])
  if 'thermometer_type' in pulse_parms:
    del pulse_parms['thermometer_type']
  in_md_name = parms['input_md_name']
  pulse_parms['input_md_name'] = os.path.abspath(in_md_name)

  pulses = ["pulse%d" % i for i in range(n_pulse)]

  for pulse, n_step in zip(pulses, n_step_list):
    os.chdir(save_dir)
    util.goto_dir(pulse)

    try:
      test_top, test_crds, test_vels = \
        get_restart_files(md_name)
      md_exists = True
    except:
      md_exists = False

    if md_exists:
      print "Skipping:", pulse
    else:
      print "Running:", pulse
      pulse_parms['n_step_dynamics'] = n_step

      pulse_in_top, pulse_in_crds, pulse_in_vels = \
        get_restart_files(pulse_parms['input_md_name'])

      soup = soup_from_restart_files(
          pulse_in_top, pulse_in_crds, pulse_in_vels)
      pulse_fn(soup)
      crds, vels = write_soup_to_restart_files(
          soup, md_name + '.pulse.in', force_field)
      pulse_parms['input_crds'] = crds
      pulse_parms['input_vels'] = vels

      run_parms(pulse_parms)

    pulse_parms['input_md_name'] = os.path.abspath(md_name)

  os.chdir(save_dir)
  open(md_name+'.time', 'w').write(timer.str()+'\n')
  
  merge_simulations(md_name, pulses)

  for pulse in pulses: 
    os.system('rm -rf %s' % pulse)

  util.write_dict(parms, config)


    
    
