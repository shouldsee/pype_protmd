#!/usr/bin/env python

"""
This module is designed to perform MD simulations that
systematically apply a local residue-based perturbation to a
protein using a variety of different MD packages.
"""

import os
import shutil
import glob
import time
import random
import math

import util
import pdbstruct
import vector3d
import simulate


#################################################
# Statistical heating routines
#################################################
# Velocity units: angstrom/ps
# Energy units: Da*(angstrom/ps)^2
# Mass units: Da
################################################


backbone_atoms = \
    ["OXT" "H", "H1", "H2", "H3","HN", "N", "C", "O", "CA", "HA"]


def maxwell_velocity(temp, mass):
  "temp in K, mass in a.m.u., output in angstroms/ps"
  factor = 8314.47148 * temp / mass # in (m/s)^2
  convert_to_ang_per_ps = 0.01 # from m/s to angstroms/ps
  return random.gauss(0, math.sqrt(factor)) * convert_to_ang_per_ps


def mean_energy(temp, n_degree_of_freedom):
  boltzmann = 8314.47148  # Da (m/s)^2 K^-1
  convert_to_ang_per_ps = 1.0E-4 # (m/s)^2 to (angstroms/ps)^2
  return n_degree_of_freedom * 0.5 * boltzmann * \
         temp * convert_to_ang_per_ps # Da (angstroms/ps)^2


def random_energy(temp, n_degree_of_freedom):
  average = mean_energy(temp, n_degree_of_freedom);
  std_dev = math.sqrt(average)
  return random.gauss(average, std_dev)

  
def kinetic_energy(atoms):
  en = 0.0
  for a in atoms:
    vel = a.vel.length()
    en += 0.5 * a.mass * vel * vel
  return en


def gas_randomize(atoms, temp):
  "Heats residues uses gas approximation: vel: angstroms/ps"
  for atom in atoms:
    atom.vel.x = maxwell_velocity(temp, atom.mass)
    atom.vel.y = maxwell_velocity(temp, atom.mass)
    atom.vel.z = maxwell_velocity(temp, atom.mass)


def anderson_velocity_scale(atoms, temp, n_degree_of_freedom):
  "Scales velocity of atoms to energy at temp. Vel: angstroms/ps"
  target_energy = mean_energy(temp, n_degree_of_freedom)
  kin = kinetic_energy(atoms)
  if vector3d.is_near_zero(kin):
    gas_randomize(atoms, temp)
  else:
    scaling_factor = math.sqrt(target_energy / kin)
    for atom in atoms:
      atom.vel.scale(scaling_factor)


def gas_heat_sidechain(soup, i, heating_temp, backbone_atoms):
  atoms = [atom for atom in soup.residue(i).atoms() 
           if atom.type not in backbone_atoms]
  gas_randomize(atoms, heating_temp)


def gas_fn(i, heating_temp, backbone_atoms):
  return lambda soup: gas_heat_sidechain(
      soup, i, heating_temp, backbone_atoms)
                                                   
                                                   
#################################################
# Rotational motion 
#################################################
# Inertia units: Da*angstrom^2
# Rotational velocity: radians/ps
################################################


def moment_of_inertia(atom, axis, anchor):
  "Units in Da*angstroms^^2"
  r = atom.pos - anchor
  r_perp = r.perpendicular_vec(axis)
  r_len = r_perp.length()
  return atom.mass * r_len * r_len


def rotational_velocity(atom, axis, anchor):
  r = atom.pos - anchor
  r_perp = r.perpendicular_vec(axis)
  vel_perp = atom.vel.perpendicular_vec(axis)
  vel_tang = vel_perp.perpendicular_vec(r_perp)
  pos_ref = vector3d.CrossProductVec(axis, r_perp)
  if vector3d.dot(vel_tang, pos_ref) < 0.0:
    sign = -1.0
  else:
    sign = 1.0
  if vector3d.is_near_zero(r_perp.length()):
    result = 0.0
  else:
    result = sign * vel_tang.length() / r_perp.length()
  return result
  

def radial_velocity(atom, axis, anchor):
  r = atom.pos - anchor
  r_perp = r.perpendicular_vec(axis)
  if vector3d.is_near_zero(r_perp.length()):
    return vector3d.Vector3d(0, 0, 0)
  return atom.vel.parallel_vec(r_perp)
  
  
def total_moment_of_inertia(atoms, axis, anchor):
  moments = [moment_of_inertia(atom, axis, anchor)
             for atom in atoms]
  return sum(moments)
    

def weighted_rotational_velocity(atoms, axis, anchor):
  moments = [ \
      moment_of_inertia(atom, axis, anchor) for atom in atoms]
  total_moment = sum(moments)
  weights = [moment / total_moment for moment in moments]
  rot_vels = [ \
      rotational_velocity(atom, axis, anchor) for atom in atoms]
  weighted_rot_vels = [ \
      rot_vel*weight for rot_vel, weight in zip(rot_vels, weights)]
  return sum(weighted_rot_vels)


def add_rotational_velocity(atoms, rot_vel, axis, anchor):
  for atom in atoms:
    r_perp = (atom.pos - anchor).perpendicular_vec(axis)
    v_tang_dir = vector3d.CrossProductVec(axis, r_perp)
    v_tang_dir_len = v_tang_dir.length()
    if vector3d.is_near_zero(v_tang_dir_len):
      v_tang = vector3d.Vector3d(0., 0., 0.)
    else:
      v_new_len = rot_vel * r_perp.length()
      v_tang = v_tang_dir.scaled_vec(v_new_len / v_tang_dir_len)
    atom.vel += v_tang
  
  
def get_axis_anchor(res, i):
  chi_topology = pdbstruct.get_res_chi_topology(res.type)

  p = [res.atom(atom_type).pos for atom_type in chi_topology[i]]
  axis = p[2] - p[1]
  anchor = p[2]
  return axis, anchor
  
    
def atoms_affected_by_chi(atoms, i):
  return [atom for atom in atoms 
          if pdbstruct.get_atom_sidechain_nesting(atom.type) >= i]


def get_rot_vel_chi(res, i):
  axis, anchor = get_axis_anchor(res, i)    
  atoms = atoms_affected_by_chi(res.atoms(), i)
  return weighted_rotational_velocity(atoms, axis, anchor)


def get_random_chi_rot_vel(res, i, temp):
  axis, anchor = get_axis_anchor(res, i)
  atoms = atoms_affected_by_chi(res.atoms(), i)
  moment = total_moment_of_inertia(atoms, axis, anchor)
  n_atom = len(atoms)
  energy = random_energy(temp, 3*n_atom)
  return math.sqrt(2 * energy / moment)


def add_rot_vel_to_chi(res, i, target_rot_vel):
  axis, anchor = get_axis_anchor(res, i)
  atoms = atoms_affected_by_chi(res.atoms(), i)
  add_rotational_velocity(atoms, target_rot_vel, axis, anchor)
  
  
def randomize_clean_high_chi_with_attractor(
    res, temp, mean_chi_values, 
    max_delta_chi=vector3d.DEG2RAD*45.0):
  "Chi-pulses the res at temp. Vel: angstroms per picosecond."

  n_chi = pdbstruct.get_n_chi(res)
  rot_vels = [get_rot_vel_chi(res, i) for i in range(n_chi)]

  # clear all velocities
  for atom in res.atoms():
    atom.vel.set(0.0, 0.0, 0.0)

  for i_chi in reversed(range(n_chi)):
    # if chi angle has exceeded bounds, change direction
    # else keep same direction
    chi = pdbstruct.calculate_chi(res, i_chi)
    delta_chi = vector3d.normalize_angle(chi - mean_chi_values[i])
    if delta_chi > max_delta_chi:
      sign = -1.0
    elif delta_chi < -max_delta_chi:
      sign = 1.0
    else:
      if rot_vels[i_chi] > 0.0:
        sign = +1.0
      else:
        sign = -1.0

    target_rot_vel = sign * get_random_chi_rot_vel(res, i_chi, temp)
    add_rot_vel_to_chi(res, i_chi, target_rot_vel)

  anderson_velocity_scale(res.atoms(), temp, 3*len(res.atoms()))


def rip_fn(soup, i, heating_temp):
  max_delta_chi = 60.0 * vector3d.DEG2RAD
  n_chi = pdbstruct.get_n_chi(soup.residue(i))
  mean_chis = [\
      pdbstruct.calculate_chi(soup.residue(i), j) 
      for j in range(n_chi)]
  func = randomize_clean_high_chi_with_attractor
  pulse_fn = lambda soup: func(
      soup.residue(i), heating_temp, mean_chis, max_delta_chi)
  return pulse_fn


# local residue-based perturbation


def wait_for_short_time():
  time.sleep(2.0 + random.randint(1, 100)/10.0)


def wait_for_file_to_be_created(fname):
  is_checking = False
  while not os.path.isfile(fname):
    if not is_checking:
      print 'Waiting for %s...' % fname
      is_checking = True
    wait_for_short_time()
  wait_for_short_time()


def minimize_pdb_or_wait(parms):
  top_dir = parms['top_dir']
  min_dir = os.path.join(top_dir, 'min')

  min_start_fname = min_dir + '/min.started'
  min_finish_fname = min_dir + '/min.finished'
  if os.path.isfile(min_finish_fname):
    return
  if os.path.isfile(min_start_fname):
    wait_for_file_to_be_created(min_finish_fname)
    return
  if not os.path.isdir(min_dir):
    os.makedirs(min_dir)
  open(min_start_fname, 'w').write('done')

  save_dir = os.getcwd()

  util.goto_dir(top_dir)

  target_raw_pdb = os.path.basename(parms['raw_pdb'])
  if not os.path.isfile(target_raw_pdb):
    shutil.copy(parms['raw_pdb'], '.')
  raw_pdb = os.path.basename(parms['raw_pdb'])
  raw_pdb = os.path.abspath(raw_pdb)
  if 'in_top' in parms and 'in_crds' in parms:
    top = os.path.abspath('sim.top')
    shutil.copy(parms['in_top'], top)
    crds = os.path.abspath('sim.crd')
    shutil.copy(parms['in_crds'], crds)
  else:
    top, crds = simulate.make_init_files(
        parms['force_field'], 
        os.path.basename(parms['raw_pdb']), 'sim')
  
  util.goto_dir('min')
  force_field = parms['force_field']
  simulate.minimize(top, crds, force_field, 'min')
  os.chdir(save_dir)

  open(min_finish_fname, 'w').write('done')

  
def equilibrate_or_wait(parms):
  temp_back = parms['temp_back']
  top_dir = parms['top_dir']
  equil_dir = '%s/equil/%dk' % (top_dir, temp_back)
  out_md_name = equil_dir + '/md'
  
  equil_start_fname = equil_dir + '/md.started'
  equil_finish_fname = equil_dir + '/md.finished'
  save_dir = os.getcwd()
  if os.path.isfile(equil_finish_fname):
    return out_md_name
  if os.path.isfile(equil_start_fname):
    wait_for_file_to_be_created(equil_finish_fname)
    return out_md_name
  util.goto_dir(equil_dir)
  open(equil_start_fname, 'w').write('done')

  in_md_name = os.path.join(top_dir, 'min', 'min')
  simulate.preequilibrate(in_md_name, 'md', temp_back)

  open(equil_finish_fname, 'w').write('done')

  util.goto_dir(save_dir)
  return out_md_name 


def heat_residue(config):
  done = config + '.done'
  if os.path.isfile(done):
    return

  parms = eval(open(config, 'r').read())
  residue_dir = parms['residue_dir']
  top_dir = parms['top_dir']
  rel_residue_dir = os.path.relpath(
      residue_dir, top_dir)

  skip_dir = False
  if os.path.isdir(residue_dir):
    skip_dir = True
  else:
    try:
      os.makedirs(residue_dir)
    except:
      skip_dir = True
  if skip_dir:
    print "Skipping:", rel_residue_dir
    return

  print "Running:", rel_residue_dir

  minimize_pdb_or_wait(parms)
  in_md_name = equilibrate_or_wait(parms)
  
  save_dir = os.getcwd()
  os.chdir(residue_dir)

  temp_perturb = parms['temp_perturb']
  thermometer_type = parms['thermometer_type']
  i = parms['residue']

  if thermometer_type == "rip":
    top, crds, vels = simulate.get_restart_files(
        in_md_name)
    soup = simulate.soup_from_restart_files(
        top, crds, vels)
    pulse_fn = rip_fn(soup, i, temp_perturb)
  elif thermometer_type == "gas":
    pulse_fn = gas_fn(i, temp_perturb, backbone_atoms)

  force_field = parms['force_field']
  n_step = parms['n_step']
  n_step_per_pulse = parms['n_step_per_pulse']
  simulate.pulse(
      in_md_name, 'md', force_field, n_step, 
      pulse_fn, n_step_per_pulse)

  os.chdir(save_dir)
  open(done, 'w').write('done')

  wait_for_short_time()



def rotatable_residues(pdb):
  soup = pdbstruct.Soup(pdb)
  residues = \
      [i for i, res in enumerate(soup.residues()) \
       if pdbstruct.get_res_chi_topology(res.type) \
         and res.type != 'PRO']
  return residues  
  
  
def create_residue_configs(parms):
  residues = parms['residues']
  raw_pdb = os.path.abspath(parms['raw_pdb'])
  if not residues:
    residues = rotatable_residues(raw_pdb)
  group_dirs = []
  top_dir = os.path.abspath(parms['sim_dir'])
  print "RIP directory:", os.path.abspath(top_dir)
  for (perturb_type, temp_back, temp_perturb, n_step) \
      in parms['heating_cases']:
    group_dir = \
        '%s/%dk-rip-%dk' % (top_dir, temp_back, temp_perturb)
    group_dirs.append(group_dir)
    if not os.path.isdir(group_dir):
      os.makedirs(group_dir)
    for r in residues:
      config = '%s/res%d.config' % (group_dir, r+1)
      rel_config = os.path.relpath(
          config, group_dir)
      if not os.path.isfile(config):
        new_parms = {}
        new_parms['constraint'] = parms['constraint']
        new_parms['n_step'] = n_step
        new_parms['n_step_per_pulse'] = 100
        new_parms['raw_pdb'] = raw_pdb
        if 'in_crds' in parms and 'in_top' in parms:
          new_parms['in_top'] = os.path.abspath(parms['in_top'])
          new_parms['in_crds'] = os.path.abspath(parms['in_crds'])
        new_parms['force_field'] = parms['force_field']
        new_parms['residue'] = r
        new_parms['top_dir'] = top_dir
        new_parms['temp_back'] = temp_back
        new_parms['temp_perturb'] = temp_perturb
        new_parms['residue_dir'] = '%s/%d' % (group_dir, r+1)
        new_parms['residue'] = r
        new_parms['thermometer_type'] = perturb_type
        util.write_dict(new_parms, config)
        print "Running: creating", rel_config
      else:
        print "Skipping: creating", rel_config
  return group_dirs 
  

def process_group_config(config):
  done = config + '.done'
  if os.path.isfile(done):
    print config, "already processed."
    return
  parms = util.read_dict(config)
  for group_dir in create_residue_configs(parms):
    os.chdir(group_dir)
    residue_configs = [
        os.path.abspath(f) 
        for f in glob.glob('*.config')
        if not os.path.isfile(f + '.done')]
    save_dir = os.getcwd()
    for residue_config in residue_configs:
      os.chdir(save_dir)
      fail = residue_config + '.fail'
      if not os.path.isfile(fail):
        try:
          heat_residue(residue_config)
        except:
          open(fail, 'w').write('fail')
  

usage = """
usage: perturb.py working_dir

- working_dir: directory contains .config files to process

Sample .config file:
====================
{ 'raw_pdb':            '2gch.pdb',
  'sim_dir':            '2gch/gbsa',
  'force_field':        'AMBER',       # or 'GROMACS' or 'AMBER'
  'constraint':         '',            # 'backbone' or 'surface'
  'n_step_per_pulse':   100,
  'residues':           [],
  'heating_cases':      [('rip', 10, 26, 5000), 
                         ('gas', 10, 300, 5000),
                         ('rip', 300, 300, 10000)],
}
"""


if __name__ == "__main__":  
  import sys
  if len(sys.argv) < 2:
    print usage
    sys.exit(1)
  configs = sys.argv[1:]
  save_dir = os.getcwd()
  for full_config in configs:
    os.chdir(save_dir)
    config_dir = os.path.dirname(full_config)
    if config_dir:
      os.chdir(config_dir)
    config = os.path.basename(full_config)
    parms = util.read_dict(config)
    print parms
    if 'residues' in parms:
      process_group_config(config)
    elif 'residue' in parms:
      heat_residue(config)
    else:
      print "Can't parse %s file" % config


  

