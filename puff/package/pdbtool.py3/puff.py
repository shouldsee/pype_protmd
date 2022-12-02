#!/usr/bin/env python
# encoding: utf-8

import os
import shutil
import glob
import time
import random
import math
import re

import numpy

import util
import pdbstruct
import vector3d
import simulate
import perturb


convert_Nm_to_cal = 0.000238846
avogadro = 6.02214179E23
work_in_pNAng_to_kCalPerMol = \
    1E-12*1E-10*convert_Nm_to_cal*avogadro
momentum_in_daAngPerPs_to_force_in_pN = 1.66E-1
work_in_DaAngsPerPsSq_to_pNAng = 1.66E-1
md_step_in_ps = 0.001 


# Strategy objects to implement an "apply" function that
# modifies a protein Soup object. Strategy objects are 
# used so that states can be stored and processed between
# pulses.

def average_vel(atoms):
  momentum = vector3d.Vector3d()
  mass = 0.0
  for a in atoms:
    momentum += a.vel.scaled_vec(a.mass)
    mass += a.mass
  return momentum.scaled_vec(1.0/mass)
  

def add_vel_to_atoms(atoms, vel_diff):
  for a in atoms:
    a.vel_last = a.vel
    a.vel += vel_diff
    a.work_delta = \
        vector3d.dot(vel_diff, a.vel) * a.mass * md_step_in_ps * \
        work_in_DaAngsPerPsSq_to_pNAng
    

def get_atoms_of_residues(soup, residues):
  atoms = []
  for i in residues:
    atoms.extend(soup.residue(i).atoms())
  return atoms
  
  
class PushApartByAcc(object):
  """
  The simplest type of force, a constant force. The
  acceleration is just taken from the target_val
  """
  def __init__(
      self, domain1, domain2, target_val, dt=0.1,
      temp=None, is_backbone_only=False, force_fname='md.acc'):
    self.domain1 = domain1
    self.domain2 = domain2
    self.dt = dt
    self.target_val = target_val
    self.temp = temp
    self.force_fname = os.path.abspath(force_fname)
    self.is_backbone_only = is_backbone_only
    
  def change_velocities(self):
    diff_vel = self.target_val * self.dt
    diff_axis_vel2to1 = self.axis2to1.scaled_vec(diff_vel)
    self.vel_diff = vector3d.dot(diff_axis_vel2to1, self.axis2to1)
    self.vel = vector3d.dot(self.axis_vel2to1, self.axis2to1)
    diff_axis_vel2to1.scale(0.5)
    add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
    add_vel_to_atoms(self.atoms2, -diff_axis_vel2to1)

  def print_str(self):
    separation = self.disp2to1.length()
    force = self.vel_diff*self.mass*momentum_in_daAngPerPs_to_force_in_pN
    new_vel2to1 = \
        average_vel(self.atoms1) - average_vel(self.atoms2)
    new_vel = vector3d.dot(new_vel2to1, self.axis2to1)
    displacement = new_vel*md_step_in_ps
    work = force*displacement*work_in_pNAng_to_kCalPerMol
    return '%f %f %f %f' % (self.vel_diff, self.vel, separation, work, work2)
    
  def apply(self, soup):
    self.soup = soup
    if self.temp:
      atoms = self.soup.atoms()
      perturb.anderson_velocity_scale(
          atoms, self.temp, 3*len(atoms))
    self.atoms1 = get_atoms_of_residues(soup, self.domain1)
    self.atoms2 = get_atoms_of_residues(soup, self.domain2)
    if self.is_backbone_only:
      atoms1 = \
          [a for a in self.atoms1 if a.type in perturb.backbone_atoms]
      atoms2 = \
          [a for a in self.atoms2 if a.type in perturb.backbone_atoms]
    self.mass = 0.0
    for i in self.domain1 + self.domain2:
      for a in self.soup.residue(i).atoms():
        self.mass += a.mass
    self.disp2to1 = pdbstruct.get_center(self.atoms1) \
                     - pdbstruct.get_center(self.atoms2)
    self.axis2to1 = self.disp2to1.normal_vec()
    self.vel2to1 = average_vel(self.atoms1) \
                    - average_vel(self.atoms2)
    self.axis_vel2to1 = self.vel2to1.parallel_vec(self.axis2to1)
    self.change_velocities()
    
    f = open(self.force_fname, 'a')
    f.write(self.print_str() + '\n')
    f.close()


class PushApartByVel(PushApartByAcc):
  def print_str(self):
    separation = self.disp2to1.length()
    force = self.vel_diff*self.mass*momentum_in_daAngPerPs_to_force_in_pN
    displacement = self.target_val*md_step_in_ps
    work = force*displacement*work_in_pNAng_to_kCalPerMol
    work2 = 0.0
    for i in self.domain1 + self.domain2:
      for a in self.soup.residue(i).atoms():
        work2 += a.work_delta
    work2 *= work_in_pNAng_to_kCalPerMol
    return '%f %f %f %f %f %f %f' % \
       (self.vel_diff, self.vel, separation, self.target_val, work, self.mass, work2)

  def change_velocities(self):
    target_axis_vel2to1 = self.axis2to1.scaled_vec(self.target_val)
    diff_axis_vel2to1 = target_axis_vel2to1 - self.axis_vel2to1
    self.vel_diff = vector3d.dot(diff_axis_vel2to1, self.axis2to1)
    self.vel = vector3d.dot(self.axis_vel2to1, self.axis2to1)
    diff_axis_vel2to1.scale(0.5)
    add_vel_to_atoms(self.atoms1, diff_axis_vel2to1)
    add_vel_to_atoms(self.atoms2, -diff_axis_vel2to1)
  
  
class PushApartByVelAndRip(PushApartByVel):
  def __init__(
      self, domain1, domain2, target_val, 
      i_rip, temp_rip, init_pdb, 
      temp=None, is_backbone_only=False, force_fname='md.acc'):
    PushApartByVel.__init__(
        self, domain1, domain2, target_val,
        temp, is_backbone_only, force_fname)
    self.temp_rip = temp_rip
    self.i_rip = i_rip
    init_soup = pdbstruct.Soup(init_pdb)
    self.rip_fn = perturb.rip_fn(init_soup, i_rip, temp_rip)

  def change_velocities(self):
    PushApartByVel.change_velocities(self)
    self.rip_fn(self.soup)


class PushApartByMovingConstraint(PushApartByVel):

  def __init__(
      self, domain1, domain2, target_val, dt=0.1,
      temp=None, is_backbone_only=False, force_fname='md.acc'):
    PushApartByVel.__init__(
        self, domain1, domain2, target_val,
        temp, is_backbone_only, force_fname)
    self.time = 0.0
    self.disp_init = 0.0
    self.dt = dt
    self.constraint_vel = target_val

  def change_velocities(self):
    disp = self.disp2to1.length()
    if self.time == 0.0:
      self.disp_init = disp
    disp_target = self.time*self.constraint_vel + self.disp_init
    dist_delta = disp_target - disp
    self.target_val = dist_delta / self.dt
    PushApartByVel.change_velocities(self)
    self.time += self.dt  

  def print_str(self):
    s = PushApartByVel.print_str(self)
    s += ' %f' % (self.target_val)
    return s


class RotateLoop(PushApartByAcc):
  def __init__(
      self, domain, target_val, temp=None, 
      is_backbone_only=False, force_fname='md.acc'):
    self.domain = domain
    self.target_val = target_val
    self.temp = temp
    self.force_fname = os.path.abspath(force_fname)
    self.center_init = None
    self.is_backbone_only = is_backbone_only

  def apply(self, soup):
    self.soup = soup
    if self.temp:
      atoms = self.soup.atoms()
      perturb.anderson_velocity_scale(
          atoms, self.temp, 3*len(atoms))
    pivot1 = self.soup.residue(self.domain[0]).atom("CA")
    pivot2 = self.soup.residue(self.domain[-1]).atom("CA")
    rot_center = pivot1.pos
    axis = (pivot2.pos - pivot1.pos).normal_vec()
    atoms = get_atoms_of_residues(soup, self.domain)
    if self.is_backbone_only:
      atoms = [
          a for a in atoms 
          if a.type in perturb.backbone_atoms]
    rot_vel = perturb.weighted_rotational_velocity(
        atoms, axis, rot_center)
    diff_rot_vel = self.target_val - rot_vel
    perturb.add_rotational_velocity(
        atoms, diff_rot_vel, axis, rot_center)
    center = pdbstruct.get_center(atoms)
    if self.center_init is None:
      self.center_init = center
    disp = vector3d.pos_distance(center, self.center_init)
    f = open(self.force_fname, 'a')
    f.write("%f %f %f %f\n" % \
        (disp, rot_vel, diff_rot_vel, self.target_val))
    f.close()

def convert_to_deg(a):
  b = a*vector3d.RAD2DEG
  if b<0:
    return b+360
  else:
    return b

class RotateDna(object):
  
  def __init__(
      self, axis_pair, strand1, strand2, protein, 
      target_rot_vel, force_fname='md.acc'):
    self.axis_pair = axis_pair
    self.strand1 = strand1
    self.strand2 = strand2
    self.protein = protein
    self.target_rot_vel = target_rot_vel
    self.force_fname = os.path.abspath(force_fname)
    
  def update_file(self, s):
    f = open(self.force_fname, 'a')
    f.write(s + '\n')
    f.close()
    
  def apply(self, soup):    
    i, j = self.axis_pair
    a1 = soup.residue(i).atom("CA")
    a2 = soup.residue(j).atom("CA")
    p1 = a1.pos
    p2 = a2.pos
    axis = p2 - p1

    spokes = [soup.residue(i).atom("C3'") for i in self.strand1]
    disp_vecs = [s.pos - p1 for s in spokes]
    radial_vecs = [d.perpendicular_vec(axis) 
                   for d in disp_vecs]
    r0 = radial_vecs[0]
    angs = [vector3d.vec_dihedral(r0, axis, r) 
            for r in radial_vecs]
    angs = [vector3d.RAD2DEG*a for a in angs]

    n_turn = 0
    cum_angs = []
    last_ang = angs[0]
    for spoke, ang in zip(spokes, angs):
      if last_ang > 0 and ang < 0:
        n_turn += 1
      cum_ang = ang + n_turn*360
      cum_angs.append(cum_ang)
      last_ang = ang

    total = cum_angs[-1]
    offset1 = (cum_angs[-1] % 360)/2. 
    offset2 = total - offset1

    mol = pdbstruct.Molecule()
    for atom, ang in zip(spokes, cum_angs):
      if ang >= offset1 and ang <= offset2:
        mol.insert_atom(atom)
    center = pdbstruct.get_center(mol.atoms())
    
    def get_rot_vel(residues):
      atoms = []
      for i in residues:
        atoms.extend(soup.residue(i).atoms())
      rot_vel = perturb.weighted_rotational_velocity(
        atoms, axis, center)
      return rot_vel*vector3d.RAD2DEG

    strands = self.strand1 + self.strand2

    pre_w_strands = get_rot_vel(strands)
    pre_w_protein = get_rot_vel(self.protein)

    for i in strands:
      for a in soup.residue(i).atoms():
        rot_vel = perturb.rotational_velocity(a, axis, center)
        rad_vel = perturb.radial_velocity(a, axis, center)
        diff_rot_vel = self.target_rot_vel - rot_vel
        perturb.add_rotational_velocity(
            [a], diff_rot_vel, axis, center)
        a.vel = a.vel - rad_vel;

    for i in self.protein:
      for a in soup.residue(i).atoms():
        rot_vel = perturb.rotational_velocity(a, axis, center)
        perturb.add_rotational_velocity(
            [a], -rot_vel, axis, center)

    post_w_strands = get_rot_vel(strands)
    post_w_protein = get_rot_vel(self.protein)

    s = "%f %f %f %f " % \
        (pre_w_strands, pre_w_protein,
         post_w_strands, post_w_protein)
    self.update_file(s)
    

class TwistByRotVel(PushApartByAcc):
  def __init__(
      self, twist_domain, anchor_domain,
      target_val, temp=None, is_backbone_only=False, force_fname='md.acc'):
    self.twist_domain = twist_domain
    self.anchor_domain = anchor_domain
    self.target_val = target_val
    self.temp = temp
    self.force_fname = os.path.abspath(force_fname)
    self.is_backbone_only = is_backbone_only

  def apply(self, soup):
    self.soup = soup
    if self.temp:
      atoms = self.soup.atoms()
      perturb.anderson_velocity_scale(
          atoms, self.temp, 3*len(atoms))
    twist_atoms = get_atoms_of_residues(
        soup, self.twist_domain)
    if self.is_backbone_only:
      twist_atoms = [
          a for a in twist_atoms 
          if a.type in perturb.backbone_atoms]
    rot_center = pdbstruct.get_center(twist_atoms)
    anchor_atoms = get_atoms_of_residues(
        soup, self.anchor_domain)
    anchor_center = pdbstruct.get_center(anchor_atoms)
    axis = anchor_center - rot_center
    rot_vel = perturb.weighted_rotational_velocity(
        twist_atoms, axis, rot_center)
    diff_rot_vel = self.target_val - rot_vel
    perturb.add_rotational_velocity(
        twist_atoms, diff_rot_vel, axis, rot_center)
    f = open(self.force_fname, 'a')
    f.write("%f %f %f\n" % \
        (rot_vel, diff_rot_vel, self.target_val))
    f.close()


class PushApartByRotVel(PushApartByAcc):

  def __init__(
      self, domain1, domain2, i_pivot, 
      target_val, temp=None, is_backbone_only=False,
      force_fname='md.acc'):
    PushApartByAcc.__init__(
        self, domain1, domain2, 
        target_val, temp, is_backbone_only, force_fname)
    self.i_pivot = i_pivot

  def change_velocities(self):
    center1 = pdbstruct.get_center(self.atoms1)
    center2 = pdbstruct.get_center(self.atoms2)
    rot_center = self.soup.residue(self.i_pivot).atom("CA").pos
    rot_disp1 = center1 - rot_center
    rot_disp2 = center2 - rot_center
    rot_axis = vector3d.CrossProductVec(rot_disp1, rot_disp2)
    # positive means rotating closer to atoms2
    rot_vel1 = perturb.weighted_rotational_velocity(
        self.atoms1, rot_axis, rot_center)
    # positive means rotating away from atoms1
    rot_vel2 = perturb.weighted_rotational_velocity(
        self.atoms2, rot_axis, rot_center)
    rot_vel2to1 = rot_vel1 - rot_vel2
    diff_rot_vel2to1 = -self.target_val - rot_vel2to1
    diff_rot_vel2to1 *= 0.5
    perturb.add_rotational_velocity(
        self.atoms1, diff_rot_vel2to1, rot_axis, rot_center)
    perturb.add_rotational_velocity(
        self.atoms2, -diff_rot_vel2to1, rot_axis, rot_center)
    self.vel_diff = 2.0 * diff_rot_vel2to1
    self.vel = rot_vel2to1

  
def get_puff_info(md_dir):
  "Returns times, dists, vels, forces, work of a PUFF trajectory"

  config = md_dir + '/config'
  print(os.getcwd(), config)
  parms = eval(open(config, 'r').read())
  s = parms['pulse_strategy']
  if 'n_step_per_pulse' in parms:
    dt = 0.001*parms['n_step_per_pulse']
  else:
    dt = 0.1
  p = re.search(r'\(.+\)', s)
  parms = eval(p.group())
  residues = parms[0] + parms[1]
  target_vel = float(parms[2])
  print("get_puff_info vel", target_vel)

  md = md_dir + '/md'
  top, crds, vels = simulate.get_restart_files(md)
  soup = simulate.soup_from_restart_files(top, crds, vels)
  mass = 0.0
  for i in residues:
    for a in soup.residue(i).atoms():
      mass += a.mass
  print("get_puff_info mass", mass)

  data = util.read_array(md_dir + '/md.acc')
  times = [dt*i for i in range(data.shape[1])]
  vel_diffs = data[0,:]
  forces = vel_diffs*mass*momentum_in_daAngPerPs_to_force_in_pN
  pre_vels = data[1,:]
  dists = data[2,:]

  if data.shape[0] == 4:
    post_vels = data[3,:]
  else:
    post_vels = numpy.ones((len(pre_vels)), float)
    post_vels *= target_vel

  work = []
  area = 0.0
  for f, v in zip(forces, post_vels):
    displacement = md_step_in_ps*v
    area += f*md_step_in_ps*v*work_in_pNAng_to_kCalPerMol
    work.append(area)
  
  label = "%.2f Å/ps" % target_vel

  return label, times, dists, pre_vels, post_vels, forces, work

#############################################

nuc_atom_types = { \
  'G': """C1* O1P O6 C5* 1H2 O4* H3* 
        2H2 C8 C2 C6 C5 C4 1H5* O2P P 
        H1* 2H2* C2* N1 N2 N3 C4* N7 O5* H8 
        2H5* N9 C3* H1 H4* O3* 1H2*
        """,
  'C': """C1* O1P C5* 1H2* O4* 1H4 H3* C3*
         2H4 O2 C2 C6 C5 C4 1H5* O2P P H1*
          2H2* C2* N1 N3 N4 C4* O5* 2H5*
           H6 H4* H5 O3*
        """,
  'A': """C1* O1P C5* 1H2* O4* 1H6 H3* C3* 
        2H6 C8 C4* C2 C6 C5 C4 1H5* O2P P 
        H1* 2H2* C2* N1 N3 N6 N7 O5* H8 
        2H5* N9 H2 H4* O3*
        """,
  'T': """C1* H3T O1P O4 C5* 1H2* O4* 1H7 H3* 
        C3* 2H7 O2 C2 C7 C6 C5 C4 1H5* 3H7 
        O2P P H1* 2H2* C2* N1 N3 C4* O5* 2H5*
         H3 H6 H4* O3*
         """ 
}


def extract_sidechain_atoms(atoms):
  result = [a for a in atoms.split() if '*' not in a]
  return ' '.join(result)


sidechain_atom_types = {}
for key in nuc_atom_types:
  atoms = nuc_atom_types[key]
  sidechain_atom_types[key] = extract_sidechain_atoms(atoms)
  

def dna_sidechain(residue):
  nuc = residue.type.strip()
  nuc_sc_atom_types = sidechain_atom_types[nuc]
  return [a for a in residue.atoms() 
          if a.type in nuc_sc_atom_types]

    
def get_i_residue(soup, tag):
  chain_id, res_num = tag.split(":")
  res_num = int(res_num)
  for i, r in enumerate(soup.residues()):
    if r.chain_id == chain_id and r.num == res_num:
      return i
  return None


def get_basepair(soup, label1, label2):
  i = get_i_residue(soup, label1)
  j = get_i_residue(soup, label2)
  return (i, j)


def get_dna_duplex_parameters(soup, base_pair1, base_pair2):
  res1a, res1b = [soup.residue(k) for k in base_pair1]
  center1 = pdbstruct.get_center(
      dna_sidechain(res1a) + dna_sidechain(res1b))

  res2a, res2b = [soup.residue(k) for k in base_pair2]
  center2 = pdbstruct.get_center(
      dna_sidechain(res2a) + dna_sidechain(res2b))

  axis = center2 - center1

  backbone_pos1 = res1a.atom("C3*").pos
  r1 = backbone_pos1 - center1
  r1_perp = r1.perpendicular_vec(axis)
  r1_on_axis = backbone_pos1 - r1_perp

  backbone_pos2 = res2a.atom("C3*").pos
  r2 = backbone_pos2 - center2
  r2_perp = r2.perpendicular_vec(axis)
  r2_on_axis = backbone_pos2 - r2_perp

  separation = r2_on_axis - r1_on_axis
  dihedral = vector3d.vec_dihedral(r2_perp, axis, r1_perp)

  return axis, center1, center2, separation, dihedral


class RotateDuplex(PushApartByAcc):
  def __init__(
      self, duplex, basepair1, basepair2, target_val, temp=None, 
      is_backbone_only=False, force_fname='md.acc'):
    self.duplex = duplex
    self.basepair1 = basepair1
    self.basepair2 = basepair2
    self.target_val = target_val
    self.temp = temp
    self.force_fname = os.path.abspath(force_fname)

  def apply(self, soup):
    self.soup = soup
    if self.temp:
      atoms = self.soup.atoms()
      perturb.anderson_velocity_scale(
          atoms, self.temp, 3*len(atoms))
    axis, center1, center2, separation, dihedral = \
      get_dna_duplex_parameters(
          self.soup, self.basepair1, self.basepair2) 
    translation = vector3d.Translation(separation)
    rotation = vector3d.Rotation(axis, dihedral, center1)
    transform = rotation*translation
    ratio = separation.length()/dihedral

    atoms = get_atoms_of_residues(soup, self.duplex)

    vel_axis = average_vel(atoms).parallel_vec(axis) 
    vel_target = axis.normal_vec().scaled_vec(self.target_val)
    add_vel_to_atoms(atoms, vel_target - vel_axis)
    
    rot_vel = perturb.weighted_rotational_velocity(
        atoms, axis, center1)
    rot_vel_target = self.target_val/ratio
    rot_vel_diff = rot_vel_target - rot_vel
    perturb.add_rotational_velocity(
        atoms, rot_vel_diff, axis, center1)

    # f = open(self.force_fname, 'a')
    # f.write("%f %f %f %f\n" % \
    #     (disp, rot_vel, diff_rot_vel, self.target_val))
    # f.close()


class AngleConstraint(PushApartByAcc):
  def __init__(self, constraints, temp=None):
    self.constraints = constraints
    self.temp = temp

  def apply(self, soup):
    self.soup = soup
    if self.temp:
      atoms = self.soup.atoms()
      perturb.anderson_velocity_scale(
          atoms, self.temp, 3*len(atoms))
    for i, j, target_angle in self.constraints:
      center = soup.residue(i).atom("CA")
      pivot = soup.residue(i).atom("C")
      target = soup.residue(j).atom("CA")
      arm1 = pivot.pos - center.pos
      arm2 = target.pos - center.pos
      axis = vector3d.CrossProductVec(arm1, arm2)
      angle = vector3d.vec_dihedral(arm1, axis, arm2)
      angle_diff = target_angle - angle
      if abs(angle_diff) < 1*vector3d.DEG2RAD:
        sign = angle_diff/abs(angle_diff)
        atoms = soup.residue(j).atoms()
        seq_sep = abs(j-i)
        if seq_sep > 10:
          mag = 10
        else:
          mag = 20-2*i+10
        perturb.add_rotational_velocity(
            atoms, mag*sign, axis, center)

#############################################


# processes the .puff-config files

def preequilibrate(
    top, crds, force_field="AMBER", temperature=300, n_step=1000):
  save_dir = os.getcwd()
  util.goto_dir('min')
  simulate.minimize(top, crds, force_field, 'min')
  util.goto_dir('../equil')
  top, crds, vels = simulate.get_restart_files('../min/min')
  simulate.langevin(
      top, crds, force_field, vels, n_step, temperature, 'md')
  os.chdir(save_dir)
    

def make_abs_fnames_in_dict(parms):
  keywords = ['pdb', 'top', 'md', 'crd', 'dir', 'top']
  for key in parms:
    if isinstance(parms[key], str):
      is_fname = False
      for keyword in keywords:
        if keyword in key:
          is_fname = True
          break
      if is_fname:
        parms[key] = os.path.abspath(parms[key])


def setup_run_or_wait_equilibrate(parms):
  save_dir = os.getcwd()
  util.goto_dir(parms['top_dir'])
  equil_finish_fname = os.path.abspath('equil_finished')
  rel_dir = os.path.relpath('equil', save_dir)
  msg = "equilibration"
  if simulate.continue_if_file(equil_finish_fname, msg):
    equil_start_fname = os.path.abspath('equil_started')
    if os.path.isfile(equil_start_fname):
      print("Waiting for other node: equilibrate", rel_dir)
      perturb.wait_for_file_to_be_created(equil_finish_fname)
    else:
      print("Running: equilibrate", rel_dir)
      open(equil_start_fname, 'w').write('done')
      if 'in_top' in parms and 'in_crds' in parms:
        top = os.path.abspath('sim.top')
        shutil.copy(parms['in_top'], top)
        crds = os.path.abspath('sim.crd')
        shutil.copy(parms['in_crds'], crds)
      else:
        top, crds = simulate.make_init_files(
            parms['force_field'], parms['in_pdb'], 'sim')
      if 'n_equilibrate_step' in parms:
        n_equilibrate_step = parms['n_equilibrate_step']
      else:
        n_equilibrate_step = 1000
      preequilibrate(
          top, crds, parms['force_field'],
          parms['temp_back'], n_equilibrate_step)
      open(equil_finish_fname, 'w').write('done')
  os.chdir(save_dir)

  
def process_config(config):
  done = os.path.abspath(config + '.done')
  rel_config = os.path.relpath(config, os.getcwd())
  msg = "PUFF simulation defined in " + rel_config
  if not simulate.continue_if_file(done, msg):
    return
    
  fail = os.path.abspath(config + '.fail')
  if os.path.isfile(fail):
    print("Skipping: previously failed", rel_config)
    return
  
  save_dir = os.getcwd()
  parms = eval(open(config).read())
  make_abs_fnames_in_dict(parms)

  if 'in_md' in parms and parms['in_md']:
    in_md = parms['in_md']
  else:
    setup_run_or_wait_equilibrate(parms)
    in_md = parms['top_dir'] + '/equil/md'
  
  util.goto_dir(parms['sim_dir'])
  this_config = 'config'
  if not os.path.isfile(this_config):
    util.write_dict(parms, this_config)
  print("Working directory:", os.path.relpath(os.getcwd(), save_dir))
  print("Strategy:", parms['pulse_strategy'])
  push_strategy = eval(parms['pulse_strategy'])

  try:
    simulate.pulse(
        in_md, 'md', parms['force_field'], parms['n_step'], 
        lambda soup: push_strategy.apply(soup),
        parms['n_step_per_pulse'])
    open(done, 'w').write('done')
  except:
    open(fail, 'w').write('fail')

  os.chdir(save_dir)


# helper routines to make config files 
# that process_config can read

def make_push_apart_config(
    top_dir, pull_dir, pdb, domain1, domain2, vel, n_step, 
    temp=None, temp_back=300, is_backbone_only=False, in_md=None,
    n_equilibrate_step=100000, n_step_per_pulse=100,
    force_field='AMBER', in_top=None,
    in_crds=None):
  parms = {}
  parms['top_dir'] = os.path.abspath(top_dir)
  parms['force_field'] = force_field
  parms['n_step'] = n_step
  parms['temp_back'] = temp_back
  parms['n_step_per_pulse'] = n_step_per_pulse

  name = 'v%.3f_t%.0f' % (vel, n_step/1000.0)
  sim_dir = pull_dir + '/' + name
  parms['sim_dir'] = os.path.abspath(sim_dir)
  parms['pulse_strategy'] = \
      'PushApartByVel(%s, %s, %f, %s, %s)' % \
          (str(domain1), str(domain2), vel, 
           str(temp), str(is_backbone_only))

  if in_md:
    parms['in_md'] = in_md
  else:
    parms['in_pdb'] = os.path.abspath(pdb)
    parms['n_equilibrate_step'] = n_equilibrate_step
    if in_top and in_crds:
      parms['in_top'] = os.path.abspath(in_top)
      parms['in_crds'] = os.path.abspath(in_crds)

  config = sim_dir.replace('/', '_') + '.puff-config'
  util.write_dict(parms, config)


usage = """ 
  puff.py - reads files of .puff-config to start a PUFF force
            simulation

  Example of config files we can process here:

  { 
    'in_pdb'             : '/home/bosco/puff/hairpin.pdb', 
    'n_equilibrate_step' : 5000, 
    'n_step'             : 100000,
    'n_step_per_pulse'   : 100,
    'force_field'        : 'AMBER', 
    'temp_back'          : 300,
    'pulse_strategy'     : 'PushApartByVel([0, 1, 2, 3, 4, 5, 6], [11, 12, 13, 14, 15, 16, 17], 0.400000, 300)', 
    'sim_dir'            : '/home/bosco/puff/hairpin/v0.40_t100/300k-rip-300k', 
    'top_dir'            : '/home/bosco/puff/hairpin' 
  }
"""


if __name__ == '__main__':
  import sys
  if len(sys.argv) < 2:
    print(usage)
    sys.exit(1)
  working_dir = os.path.abspath(sys.argv[1])
  os.chdir(working_dir)
  configs = glob.glob('*.puff-config')
  for full_config in configs:
    os.chdir(working_dir)
    config = os.path.basename(full_config)
    config_dir = os.path.dirname(full_config)
    if config_dir:
      os.chdir(config_dir)
    process_config(config)

