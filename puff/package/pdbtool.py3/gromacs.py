#!/usr/bin/env python
# encoding: utf-8
"""
gromacs.py

Created by Bosco Ho on 2009-10-13.
Copyright (c) 2009 __MyCompanyName__. All rights reserved.
"""

import os
import shutil
import glob
import xdrlib

import util
import pdbstruct

pdb2gmx = "pdb2gmx"
trjconv = "trjconv"
grompp = "grompp"
editconf = "editconf"
mdrun = "mdrun"
genion = "genion"
genbox = "genbox"

if 'talyn' in util.run_with_output('hostname'):
  mdrun = "mpiexec -np 12 /home/bosco/bin/gromacs-4.0.7/bin/mdrun"
  pdb2gmx = "/home/bosco/bin/gromacs-4.0.7/bin/pdb2gmx"
  trjconv = "/home/bosco/bin/gromacs-4.0.7/bin/trjconv"
  grompp = "/home/bosco/bin/gromacs-4.0.7/bin/grompp"
  editconf = "/home/bosco/bin/gromacs-4.0.7/bin/editconf"
  genion = "/home/bosco/bin/gromacs-4.0.7/bin/genion"
  genbox = "/home/bosco/bin/gromacs-4.0.7/bin/genbox"

# positions are in nanometers
# velocities are in nanometers/picosecond

n_dim = 3

class TrrReader:
  def __init__(self, fname):
    self.fname = fname
    self.file = open(self.fname, 'r')

    self.u = xdrlib.Unpacker(self.file.read(200))
    self.magic = self.u.unpack_int()
    self.version = self.u.unpack_string()
    self.size_ir = self.u.unpack_int()
    self.size_e = self.u.unpack_int()
    self.size_box = self.u.unpack_int()
    self.size_vir = self.u.unpack_int()
    self.size_pres = self.u.unpack_int()
    self.size_top = self.u.unpack_int()
    self.size_sym = self.u.unpack_int()
    self.size_x = self.u.unpack_int()
    self.size_v = self.u.unpack_int()
    self.size_f = self.u.unpack_int()
    self.n_atom = self.u.unpack_int()
    self.step = self.u.unpack_int()
    self.nre = self.u.unpack_int()
    self.t = self.u.unpack_float()
    self.lamb = self.u.unpack_float()

    self.calc_precision()
    self.offset = self.u.get_position()
    self.calc_size_frame()

    self.file.seek(0, 2)
    size_file = self.file.tell()
    self.n_frame = size_file / self.size_frame

  def calc_size_frame(self):
    n_vec = 0
    if self.size_box: n_vec += n_dim
    if self.size_vir: n_vec += n_dim
    if self.size_pres: n_vec += n_dim
    if self.size_x: n_vec += self.n_atom
    if self.size_v: n_vec += self.n_atom
    if self.size_f: n_vec += self.n_atom
    self.size_frame = (n_vec*n_dim*self.precision + self.offset)

  def calc_precision(self):
    "Returns 4 for single precision, and 8 for double precision"
    if self.size_box:
      self.precision = self.size_box/n_dim/n_dim
    elif self.size_x:
      self.precision = self.size_x/n_dim
    elif self.size_v:
      self.precision = self.size_v/n_dim
    elif self.size_f:
      self.precision = self.size_f/n_dim
    else:
      raise ValueError("Cannot determine precision")
    if self.precision not in [4, 8]:
      raise ValueError("Precision not single or double!")
    
  def next_3_reals(self):
    if self.precision == 4:
      return [self.u.unpack_float() for i in range(3)]
    if self.precision == 8:
      return [self.u.unpack_double() for i in range(3)]

  def __len__(self):
    return self.n_frame
    
  def __repr__(self):
    return "< Gromacs TRR Coord file %s with %d frames of %d atoms >" % \
             (self.trj, self.n_frame, self.n_atom)

  def __getitem__(self, i_frame):
    if i_frame < - 1*self.n_frame or i_frame >= self.n_frame :
      raise IndexError
    elif i_frame < 0 :
      return self.__getitem__(self.n_frame + i_frame)
    self.i_frame = i_frame
    self.file.seek(self.offset + i_frame*self.size_frame)
    self.u = xdrlib.Unpacker(self.file.read(self.size_frame))
    box, positions, velocities, forces = None, None, None, None
    if self.size_box:
      box = [self.next_3_reals() for i in range(n_dim)]
    if self.size_vir:
      dummy = [self.next_3_reals() for i in range(n_dim)]
    if self.size_pres:
      dummy = [self.next_3_reals() for i in range(n_dim)]
    if self.size_x:
      positions = [self.next_3_reals() for i in range(self.n_atom)]
    if self.size_v:
      velocities = [self.next_3_reals() for i in range(self.n_atom)]
    if self.size_f:
      forces = [self.next_3_reals() for i in range(self.n_atom)]
    return box, positions, velocities, forces


class Trajectory(object):
  def __init__(self, name):
    self.top = name + '.top'
    self.gro = name + '.gro'
    self.trr = name + '.trr'
    self.trr_reader = TrrReader(self.trr)
    self.n_frame = self.trr_reader.n_frame
    self.soup = SoupFromGromacs(self.top, self.gro)
    self.atoms = self.soup.atoms()
    self.load_frame(0)

  def load_frame(self, i_frame):
    if i_frame < -self.n_frame or i_frame >= self.n_frame :
      raise IndexError
    if i_frame < 0 :
      return self.load_frame(self.n_frame + i_frame)
    box, positions, velocities, forces = self.trr_reader[i_frame]
    self.i_frame = i_frame
    n = len(self.atoms)
    for i in range(n):
      self.atoms[i].pos.set(
          positions[i][0]*10,
          positions[i][1]*10,
          positions[i][2]*10)
      self.atoms[i].vel.set(
          velocities[i][0]*10,
          velocities[i][1]*10,
          velocities[i][2]*10)


def read_top(top):
  lines = open(top).readlines()
  masses = []
  total_qtot = 0.0
  skip = True
  top_dir = os.path.dirname(top)
  for l in lines:
    if skip:
      if 'chain topologies' in l:
        skip = False
      continue
    if l.startswith("#include"):
      itp = l.split()[1][1:-1]
      itp = os.path.join(top_dir, itp)
      if os.path.isfile(itp):
        inc_masses, inc_qtot = read_top(itp)
        masses.extend(inc_masses)
        total_qtot += inc_qtot
    if l.startswith(";"):
      break
  skip = True
  qtot = None
  for l in lines:
    if skip:
      if '[ atoms ]' in l:
        skip = False
      continue
    if l.startswith('['):
      break
    if l.startswith(";"):
      continue    
    if not l.strip():
      continue
    words = l.split()
    n = int(words[0])
    res_num = int(words[2])
    res_type = words[3]
    mass = float(words[7])
    qtot = float(words[10])
    masses.append(mass)
  if qtot is not None:
    total_qtot += qtot
  return masses, total_qtot
  

def read_gro(gro):
  positions, velocities, names, nums = [], [], [], []
  lines = open(gro).readlines()
  for l in lines[2:-1]:
    words = l.split()
    if len(l) <  43:
      continue
    positions.append([10.*float(w) for w in l[20:44].split()])
    if len(l) < 63:
      continue
    velocities.append([10.*float(w) for w in l[44:68].split()])
    names.append(l[9:15].strip())
    nums.append(int(l[15:20]))
  box = [float(w) for w in lines[-1].split()]
  return positions, velocities, box, names


def AtomFromGroLine(line):
  """
  Returns an Atom object from an atom line in a pdb file.
  Converts from nm/ps to angstroms/ps
  """
  atom = pdbstruct.Atom()
  atom.res_num = int(line[0:5])
  atom.res_type = line[5:8]
  atom.type = line[10:15].strip(" ")
  if atom.res_type == "ILE" and atom.type == "CD":
    atom.type = "CD1"
  element = ''
  for c in line[12:15]:
    if not c.isdigit() and c != " ":
      element += c
  if element[:2] in pdbstruct.two_char_elements:
    atom.element = element[:2]
  else:
    atom.element = element[0]
  atom.num = int(line[15:20])
  x = 10.0*float(line[20:28])
  y = 10.0*float(line[28:36])
  z = 10.0*float(line[36:44])
  atom.pos.set(x, y, z)
  if len(line) > 62:
    x = 10.0*float(line[44:52])
    y = 10.0*float(line[52:60])
    z = 10.0*float(line[60:68])
    atom.vel.set(x, y, z)
  return atom


def SoupFromGromacs(top, gro, skip_solvent=True):
  atoms = []
  temp_list = []
  timer = util.Timer()
  lines = open(gro, 'r').readlines()
  remaining_text = ""
  n_remaining_text = 0
  for i_line, line in enumerate(lines[2:-1]):
    atom = AtomFromGroLine(line)
    if skip_solvent and atom.res_type == "SOL":
      remaining_text = "".join(lines[i_line+2:-1])
      n_remaining_text = len(lines[i_line+2:-1])
      break
    atoms.append(atom)
  box = [float(w) for w in lines[-1].split()]
  masses, q_tot = read_top(top)
  for a, mass in zip(atoms, masses):
    a.mass = mass
  soup = pdbstruct.Polymer()
  curr_res_num = -1
  for a in atoms:
    if curr_res_num != a.res_num:
      res = pdbstruct.Residue(
          a.res_type, a.chain_id, a.res_num)
      soup.append_residue_no_renum(res.copy())
      curr_res_num = a.res_num
    soup.insert_atom(-1, a)
  soup.box = box
  soup.remaining_text = remaining_text
  soup.n_remaining_text = n_remaining_text
  return soup


def get_restart_files(name):
  top = os.path.abspath(name + '.top')
  crds = os.path.abspath(name + '.gro')
  vels = ''
  return top, crds, vels


def soup_from_restart_files(top, crds, vels, skip_solvent=True):
  return SoupFromGromacs(top, crds, skip_solvent)


def write_soup_to_gro(soup, gro):
  f = open(gro, 'w')
  f.write("Generated by gromacs.py\n")
  atoms = soup.atoms()
  n_atom = len(atoms) + soup.n_remaining_text
  f.write(str(n_atom) + '\n')
  for a in atoms:
    if a.res_type == "ILE" and a.type == "CD1":
      a.type = "CD"
    # GRO doesn't care about actual numbering,
    # will wrap when the columns are fill
    res_num = a.res_num % 100000
    atom_num = a.num % 100000
    s = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n" % \
        (res_num, a.res_type, a.type, atom_num, 
         a.pos.x*0.1, a.pos.y*0.1, a.pos.z*0.1,
         a.vel.x*0.1, a.vel.y*0.1, a.vel.z*0.1)
    f.write(s)
  if soup.remaining_text:
    f.write(soup.remaining_text)
  f.write("%10.5f%10.5f%10.5f\n" % (soup.box[0], soup.box[1], soup.box[2]))
  f.close()
  

def write_soup_to_restart_files(soup, name):
  write_soup_to_gro(soup, name + '.gro')
  return name + '.gro', ''
  

def convert_to_pdb(gro, pdb):
  util.run_with_output_file(
      '%s -f %s -o %s' % (editconf, gro, pdb), pdb)
  # fix the naming conflict with ILE CD
  lines = open(pdb, 'r')
  new_lines = []
  for l in lines:
    res = l[17:20]
    atom = l[13:16]
    if res == "ILE" and atom == "CD ":
      new_l = l[:13] + "CD1" + l[16:]
    else:
      new_l = l
    new_lines.append(new_l)
  open(pdb, 'w').write(''.join(new_lines))


def delete_backup_files(tag):
  for f in util.re_glob('*', '^#' + tag):
    os.remove(f)
    

def check_file(f):
  if not os.path.isfile(f):
    raise IOError("%f not found" % f)

ions_mdp = """
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep   ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0  ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01    ; Energy step size
nsteps      = 50000   ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist     = 1       ; Frequency to update the neighbor list and long range forces
ns_type     = grid    ; Method to determine neighbor list (simple, grid)
rlist       = 1.0     ; Cut-off for making neighbor list (short range forces)
coulombtype = PME     ; Treatment of long range electrostatic interactions
rcoulomb    = 1.0     ; Short-range electrostatic cut-off
rvdw        = 1.0     ; Short-range Van der Waals cut-off
pbc         = xyz     ; Periodic Boundary Conditions (xyz/no)
"""
def neutralize_with_salt(in_top, in_gro, out_name):
  masses, qtot = read_top(in_top)
  counter_ion_charge = -int(qtot)

  if counter_ion_charge == 0:
    shutil.copy(in_gro, out_name + '.gro')
    return

  in_mdp = out_name + '.in.mdp'
  open(in_mdp, 'w').write(ions_mdp)
  top = out_name + '.top'
  if in_top != top:
    shutil.copy(in_top, top)
  tpr = out_name + '.salt.tpr'
  out_mdp = out_name + '.mdp'
  util.run_with_output_file(
       '%s -f %s -po %s -c %s -p %s -o %s' \
            % (grompp, in_mdp, out_mdp, in_gro, top, tpr),
        tpr)
  check_file(tpr)

  gro = out_name + '.gro'
  log = out_name + '.genion.log'
  in_f = out_name + '.genion.in'
  open(in_f, 'w').write('12') # pick solvent
  charge_str = " -pname Na -nname Cl "
  if counter_ion_charge > 0:
    charge_str += " -np %d " % counter_ion_charge
  else:
    charge_str += " -nn %d " % abs(counter_ion_charge)
  util.run_with_output_file(
       '%s -g %s -s %s -o %s -p %s %s  < %s' \
           % (genion, log, tpr, gro, top, charge_str, in_f),
       gro + '.salt')
  check_file(gro)


def name_mangle(in_pdb, out_pdb):
  s = pdbstruct.Polymer(in_pdb)
  new_res = {'HIS':'HID', 'CYS':'CYM', 'LYS':'LYN'}
  for r in s.residues():
    if r.type in new_res:
      r.set_type(new_res[r.type])
  dnas = ['DA', 'DT', 'DC', 'DG']
  for c in s.chains():
    for (i, suffix) in [(0, '5'), 
                        (c.n_residue()-1, '3')]:
      r = c.residue(i)
      if r.type.strip() in dnas:
        new_r_type = r.type.strip() + suffix
        r.set_type(new_r_type)
  for r in s.residues():
    r.change_atom_type('OP1', 'O1P')
    r.change_atom_type('OP2', 'O2P')
  s.write_pdb(out_pdb)
  

def prep_solvation(pdb, name, is_neutralize=True, is_constrain=False):
  save_dir = os.getcwd()
  solv_dir = os.path.join(os.path.dirname(name), 'solv')
  full_pdb = os.path.abspath(pdb)
  util.goto_dir(solv_dir)
  pdb_copy = os.path.basename(pdb)
  if not os.path.isfile(pdb_copy):
    shutil.copy(full_pdb, pdb_copy)
  raw_gro = name + '.raw.gro'
  top = name + '.top'
  itp = name + '_posre.itp'
  cmd = '%s -ignh -ff amber99 -missing -f %s -o %s -p %s -i %s' \
      % (pdb2gmx, pdb, raw_gro, top, itp)
  util.run_with_output_file(cmd, '%s.pdb2gmx' % name)
  check_file(raw_gro)

  box_gro = name + '.box.gro'
  dist_to_edge = 1.0 # in nm 
  util.run_with_output_file(
      '%s -f %s -o %s -c -d %f -bt cubic' \
          % (editconf, raw_gro, box_gro, dist_to_edge),
      '%s.box' % name)
  check_file(box_gro)

  solvated_gro = name + '.solvated.gro'
  util.run_with_output_file(
      '%s -cp %s -cs spc216.gro -o %s -p %s' \
          % (genbox, box_gro, solvated_gro, top),
       '%s.solvated' % name)
  check_file(solvated_gro)

  if is_neutralize:
    neutralize_with_salt(top, solvated_gro, name)
  gro = name + '.gro'
  convert_to_pdb(gro, name + '.pdb')
  fnames = util.re_glob(
      '*', 
      os.path.basename(name) + \
      r'[^\.]*\.(pdb|itp|gro|mdp|top)$')
  for fname in fnames:
    shutil.copy(fname, save_dir)
  delete_backup_files(name)
  os.chdir(save_dir)


minimization_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_name' : 'min', 
  'force_field': 'GROMACS',
  'constraint_pdb': '',
  'n_step_minimization' : 100, 
} 

min_mdp = """
; template .mdp file used as input into grompp to generate energy minimization for mdrun

; Parameters describing what to do, when to stop and what to save
integrator  = steep    ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0   ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01     ; Energy step size
nsteps      = $n_step_minimization    ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist     = 1        ; Frequency to update the neighbor list and long range forces
ns_type     = grid     ; Method to determine neighbor list (simple, grid)
rlist       = 1.0      ; Cut-off for making neighbor list (short range forces)
coulombtype = PME      ; Treatment of long range electrostatic interactions
rcoulomb    = 1.0      ; Short-range electrostatic cut-off
rvdw        = 1.0      ; Short-range Van der Waals cut-off
pbc         = xyz      ; 3-dimensional periodic boundary conditions (xyz|no)
"""

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_name' : 'md', 
  'force_field': 'GROMACS',
  'solvent_state': 2,
  'surface_area': 1,
  'constraint_pdb': '',
  'n_step_per_snapshot' : 50, 
  'n_step_dynamics' : 1000, 
  'temp_initial': 0.0, # ignored if it is 0.0
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.gro', 
  'output_name' : 'md', 
  'force_field': 'GROMACS',
  'cutoff': 16.0,  # non-bonded cutoff
  'constraint_pdb': '',
  'random_seed' : 2342, 
  'temp_thermometer' : 300.0, 
  'temp_initial': 0.0, # ignored if it is 0.0
  'n_step_per_snapshot' : 50, 
  'n_step_dynamics' : 1000, 
} 

dynamics_mdp = """
; Run parameters
integrator      = md            ; leap-frog integrator
nsteps          = $n_step_dynamics  ; time = n_step_dynamics*dt
dt              = 0.001         ; time-step in fs
; Output control
nstxout         = $n_step_per_snapshot  ; save coordinates every 0.05 ps
nstvout         = $n_step_per_snapshot  ; save velocities every 0.05 ps
nstenergy       = $n_step_per_snapshot  ; save energies every 0.05 ps
nstlog          = $n_step_per_snapshot  ; update log file every 0.05 ps
; Bond parameters
continuation    = yes           ; first dynamics run
constraints     = hbonds        ; bonds from heavy atom to H, constrained
; constraints     = all-bonds     ; all bonds (even heavy atom-H bonds) constrained
; constraint_algorithm = lincs    ; holonomic constraints 
; lincs_iter      = 1             ; accuracy of LINCS
; lincs_order     = 4             ; also related to accuracy
; Neighborsearching
ns_type         = grid          ; search neighboring grid cels
nstlist         = 5             ; 10 fs
rlist           = 1.0           ; short-range neighborlist cutoff (in nm)
rcoulomb        = 1.0           ; short-range electrostatic cutoff (in nm)
rvdw            = 1.0           ; short-range van der Waals cutoff (in nm)
; Periodic boundary conditions
pbc             = xyz           ; 3-dimensional periodic boundary conditions (xyz|no)
; Electrostatics
; coulombtype     = PME           ; Particle Mesh Ewald for long-range electrostatics
; pme_order       = 4             ; cubic interpolation
; fourierspacing  = 0.16          ; grid spacing for FFT
; Pressure coupling is on
pcoupl          = Parrinello-Rahman   ; Pressure coupling on in NPT
pcoupltype      = isotropic     ; uniform scaling of box vectors
tau_p           = 2.0           ; time constant, in ps
ref_p           = 1.0           ; reference pressure, in bar
compressibility = 4.5e-5        ; isothermal compressibility of water, bar^-1
; Dispersion correction
DispCorr        = EnerPres      ; account for cut-off vdW scheme
"""

temp_mdp = """
; Temperature coupling is on
tcoupl          = V-rescale     ; modified Berendsen thermostat
tc-grps         = Protein Non-Protein   ; two coupling groups - more accurate
tau_t           = 0.1   0.1     ; time constant, in ps
ref_t           = %s  %s   ; reference temperature, one for each group, in K
"""

vel_mdp = """
; Velocity generation
gen_vel          = yes       ; assign velocities from Maxwell distribution
gen_temp         = %s     ; temperature for Maxwell distribution
gen_seed         = -1        ; generate a random seed
"""

constraint_mdp = """
define           =  -DPOSRES
"""

def make_mdp(parms):
  """
  Using GROMACS, we will always solvate, use periodic 
  boundaries and Particle Ewald Mesh.
  """
  mdp = ""
  if 'n_step_minimization' in parms:
    cp = min_mdp
    mdp = util.replace_dict(cp, parms)  
  else:
    cp = dynamics_mdp
    mdp += util.replace_dict(cp, parms)
    if 'temp_thermometer' in parms:
      t = parms['temp_thermometer']
      mdp += temp_mdp % (t, t)
    else:
      mdp += "; Temperature coupling is off\n"
      mdp += "tcoupl          = no\n"
    if 'temp_initial' in parms and parms['temp_initial'] > 0.0:
      t = parms['temp_initial']
      mdp += vel_mdp % t
    else:        
      mdp += "; Velocity generation\n"
      mdp += "gen_vel         = no            ; Velocity generation is off \n"
    mdp = "title           = Template for constant temperature/pressure\n" +\
          mdp
  if 'restraint' in parms and parms['restraint']:
    mdp += "define            = -DPOSRES  ; position restrain the protein\n"
  return mdp
  
  
def replace_include_file(chain, r1, r2):
  lines = open(chain, 'r').readlines()
  new_lines = []
  for l in lines:
    if l.startswith('#include'):
      new_lines.append(l.replace(r1, r2))
    else:
      new_lines.append(l)
  open(chain, 'w').write(''.join(new_lines))


def run(parms):
  name = parms['output_name']
  in_mdp = name + '.in.mdp'
  open(in_mdp, 'w').write(make_mdp(parms))
  mdp = name + '.mdp'

  # Copies across topology files and does
  # some name mangling so that all new topology
  # files point to each other correctly
  top = name + '.top'
  in_top = parms['topology']
  shutil.copy(in_top, top)
  in_name = os.path.basename(in_top).replace('.top', '')
  in_dir = os.path.dirname(in_top)
  file_tag = "%s/%s_*itp" % (in_dir, in_name)
  new_files = [top]
  for f in glob.glob(file_tag):
    new_f = os.path.basename(f)
    new_f = new_f.replace(in_name, name)
    shutil.copy(f, new_f)
    new_files.append(new_f)
  for f in new_files:
    replace_include_file(f, in_name + "_", name + "_")

  in_gro = name + '.in.gro'
  shutil.copy(parms['input_crds'], in_gro)
  tpr = name + '.tpr'
  util.run_with_output_file(
       '%s -f %s -po %s -c %s -p %s -o %s' \
            % (grompp, in_mdp, mdp, in_gro, top, tpr),
        tpr)

  stopwatch = util.Timer()
  util.run_with_output_file(
       '%s -v -deffnm %s' % (mdrun, name),
       name + '.mdrun')
  stopwatch.stop()
  open(name + '.time', 'w').write(stopwatch.str())

  delete_backup_files(name)


def merge_simulations(name, pulses):
  save_dir = os.getcwd()
  trr_fname = name + '.trr'
  trrs = [os.path.join(p, trr_fname) for p in pulses]
  trrs = [t for t in trrs if os.path.isfile(t)]
  f = open(trr_fname, 'w')
  for trr in trrs:
    trr = TrrReader(trr)
    for i_frame in range(trr.n_frame-1):
      trr.file.seek(trr.size_frame*i_frame)
      txt = trr.file.read(trr.size_frame)
      f.write(txt)
  f.close()
  for pulse in reversed(pulses):
    gro = "%s/%s.gro" % (pulse, name)
    if os.path.isfile(gro):
      for ext in ['.top', '.itp', '.tpr', '.mdp', '.in.mdp', '.gro']:
        os.system('cp %s/%s*%s .' % (pulse, name, ext))
      os.chdir(save_dir)
      return
  os.chdir(save_dir)
    
    
def preequilibrate(in_name, out_name, temperature):
  top, gro, dummy = get_restart_files(in_name)
  restraint_name = 'restraint'
  parms = langevin_thermometer_parms.copy()
  parms['restraint'] = True
  parms['topology'] = top
  parms['input_crds'] = gro
  parms['output_name'] = restraint_name
  parms['temp_thermometer'] = temperature
  parms['temp_initial'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 5000
  run(parms)  

  top, gro, dummy = get_restart_files(restraint_name)
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = gro
  parms['output_name'] = out_name
  parms['temp_thermometer'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 5000
  run(parms)  

