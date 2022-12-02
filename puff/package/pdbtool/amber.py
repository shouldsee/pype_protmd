#!/usr/bin/env python

"""
This modules interaces with the amber molecular dynamics package,
which should already be installed and available on the path. The
routines are used to access the topology files and get the frames
of a trajectory. The interface is mainly through pdb files and
strings.
"""


import os
import copy
import shutil
import util
import pdbstruct
import string
import vector3d
import pdbtext

ambpdb = 'ambpdb'
tleap = 'tleap'
sander = 'sander'

if 'talyn' in util.run_with_output('hostname'):
  ambpdb = '/home/bosco/bin/amber11/exe/ambpdb'
  tleap = '/home/bosco/bin/amber11/exe/tleap'
  sander = '/home/bosco/bin/amber11/exe/sander'

module_dir = os.path.dirname(__file__) 
if module_dir is "": 
  module_dir = "."


def convert_amber_to_pdb(top, crd, pdb):
  """
  Converts an amber .crd (or .rst) file to .pdb.
  Beware that ambpdb craps out if the input filenames
  are too long!
  """
  # Since ambpdb can't handle long filenames
  # if the top and crd files are not in the current directory,
  # will copy into temporary files in the current directory  
  temp_fnames = []
  temp_pdb = util.temp_fname('.pdb')
  temp_fnames.append(temp_pdb)
  temp_top = top
  if os.path.dirname(top):
    temp_top = util.temp_fname('.top')
    shutil.copy(top, temp_top)
    temp_fnames.append(temp_top)
  temp_crd = crd
  if os.path.dirname(crd):
    root, ext = os.path.splitext(crd)
    temp_crd = util.temp_fname(ext)
    shutil.copy(crd, temp_crd)
    temp_fnames.append(temp_crd)
  os.system("%s -bres -p %s < %s > %s 2> /dev/null" \
             % (ambpdb, temp_top, temp_crd, temp_pdb))
  # clean-up the amber format to conform with standard
  # pdb format with HETATM fields and chain id's
  lines = open(temp_pdb, 'r').readlines()
  i_chain = 0
  is_water = False
  new_lines = []
  for line in lines:
    is_water = False
    if line.startswith('ATOM'):
      is_water = line[17:20] == 'WAT'
    if line.startswith('ATOM'): # normal protein chain
      if is_water:
        new_lines.append("HETATM" + line[6:17] + "HOH" + line[20:])
      else: # normal protein chain
        chain_id = string.ascii_uppercase[i_chain]
        new_lines.append(line[:21] + chain_id + line[22:])
    elif line.startswith('TER'):
      if not is_water:
        new_lines.append(line)
    else:
      new_lines.append(line)
    if not is_water and line.startswith('TER'):
      i_chain += 1
  txt = ''.join(new_lines)
  if not txt.strip():
    raise ValueError, "Failed to make PDB from %s" % top
  else:
    open(pdb, 'w').write(txt)
  for fname in temp_fnames:
    os.remove(fname)


def is_top_has_box_dimension(top):
  attributes = get_attribute_list(top, 'POINTERS', 8)
  return int(attributes[27]) > 0
  

def convert_pdb_to_amber(pdb, top, crd, solvent_buffer=0.0): 
  """Convert a .pdb file into amber .top and .crd file."""

  script = "# generate %s and %s files from %s\n" % (top, crd, pdb)
  script += "\n"
  script += "# load in amber force field\n"
  script += "source leaprc.ff99SB\n"
  script += "\n"
  script += "# use AMBER6 GB radii for igb=1, gbparm=2\n"
  script += "set default PBradii mbondi\n"
 
  # strip pdb and then load into soup
  bare_pdb = util.temp_fname('.pdb')
  txt = open(pdb, 'r').read()
  txt = pdbtext.strip_other_nmr_models(txt)
  txt = pdbtext.strip_lines(txt, lambda l: l.startswith('HETATM'))
  txt = pdbtext.strip_lines(txt, lambda l: l.startswith('ANISOU'))
  txt = pdbtext.strip_lines(txt, lambda l: l.startswith('CONECT'))
  txt = pdbtext.strip_lines(txt, lambda l: l.startswith('MASTER'))
  txt = pdbtext.strip_alternative_atoms(txt)
  # txt = pdbtext.strip_hydrogens(txt)
  txt = pdbtext.renumber_residues(txt)
  open(bare_pdb, 'w').write(txt)
  soup = pdbstruct.Soup(bare_pdb)
  os.remove(bare_pdb)
  
  script += "\n"
  residues = [r.type for r in soup.residues()]

  if 'PHD' in residues:
    script += "\n"
    script += "# non-standard amino acid PHD"
    script += "source leaprc.gaff\n"
    script += "loadAmberPrep %s/parms/phd.prepin\n" % module_dir
    script += "loadAmberParams %s/parms/phd.frcmod\n" % module_dir

  in_pdb = pdb.replace('.pdb', '.tleap.pdb')
  soup.write_pdb(in_pdb)

  script += "\n"
  script += " # insert protein\n"
  script += "pdb = loadpdb %s\n" % in_pdb
  # script += "check pdb\n"

  script += "\n"
  script += " # disulfide bonds\n"
  n = len(soup.residues())
  for i in range(n):
    for j in range(i+1, n):
      if soup.residue(i).type in 'CYS' and soup.residue(j).type in 'CYS':
        p1 = soup.residue(i).atom('SG').pos
        p2 = soup.residue(j).atom('SG').pos
        if vector3d.pos_distance(p1, p2) < 3.0:
          soup.residue(i).type = 'CYX'
          for atom in soup.residue(i).atoms():
            atom.res_type = 'CYX'
          soup.residue(j).type = 'CYX'
          for atom in soup.residue(j).atoms():
            atom.res_type = 'CYX'
          script += "bond pdb.%d.SG pdb.%d.SG\n" % (i+1, j+1)

  soup.write_pdb(in_pdb)

  if solvent_buffer > 0.0:
    script += "\n"
    script += " # add explicit waters\n"
    script += "solvateBox pdb TIP3PBOX %f iso\n" % solvent_buffer
    script += "\n"

  script += "\n"
  script += " # write files\n"
  script += "saveAmberParm pdb %s %s\n" % (top, crd)
  script += "\n"
  script += "quit\n"
  script += "\n"

  name = crd.replace('.crd', '')
  tleap_in = name + ".tleap.in"
  open(tleap_in, "w").write(script)
  tleap_out = name + ".tleap.out"
  os.system("%s -f %s > %s" % (tleap, tleap_in, tleap_out))
  if 'FATAL' in open(tleap_out).read():
    raise ValueError, "Failed to make topology file %s" % top

  convert_amber_to_pdb(top, crd, name + '.pdb')

  os.remove("leap.log")
  

def get_attribute_list(top, attribute, length):
  """Gets a list of attribute (of length) from .top file."""

  f = open(top, "rU")

  # skip lines until tag is found
  tag = "%%FLAG %s" % attribute
  while True:
    line = f.readline()
    if len(line) == 0:
      raise ValueError, "Can't find %s in file." % tag
    if line.startswith(tag):
      break

  # skip format line
  f.readline()

  # read string containing property
  s = ""
  while True:
    line = f.readline()
    if len(line) == 0 or line.startswith("%"):
      break
    s += line.replace("\n", "")

  f.close()
  
  return [s[i:i+length] for i in xrange(0, len(s)-1, length)]
    

def save_crds(fname, crds):
  "Saves coordinates to an amber restart file."

  f = open(fname, "w")

  f.write("ACE".ljust(80) + "\n")
  f.write("%5d  0.0000000E+00\n" % (len(crds) // 3))

  p = ["%12.7f" % x for x in crds]

  n_val_per_line = 6
  r = len(p) % n_val_per_line
  if r > 0: 
    p.extend([""] * (n_val_per_line - r))

  for i in xrange(0, len(p), n_val_per_line):
    f.write("".join(p[i:i + n_val_per_line]) + "\n")

  f.close()  


def convert_crd_to_trj_frame(crd):
  vals = [float(word) for word in util.words_in_file(crd)[1:]]
  lines = []
  line = ''
  for i in range(0, len(vals)):
    line += "%8.3f" % vals[i]
    if (i % 10) == 9:
      lines.append(line)
      line = ''
  if line:
    lines.append(line)
  return '\n'.join(lines) + '\n'


def load_rst_into_soup(soup, rst):
  f = open(rst, "r")
  f.readline()

  n_atom = int(f.readline().split()[0])
  n_crd = n_atom * 3
  n_line = n_crd / 6
  if n_crd % 6 > 0:
    n_line += 1

  line_list = [f.readline()[:-1] for i in range(0, n_line)]
  s = "".join(line_list)
  vals = [float(s[i:i+12]) for i in xrange(0, len(s), 12)]
  if len(vals) <> n_crd:
    raise ValueError, "Improper number of coordinates in rst file."
    
  for i, atom in enumerate(sorted(soup.atoms(), pdbstruct.cmp_atom)):
    atom.pos.set(vals[i*3], vals[i*3+1], vals[i*3+2])

  line_list = [f.readline()[:-1] for i in range(0, n_line)]
  s = "".join(line_list)
  vals = [float(s[i:i+12]) for i in xrange(0, len(s), 12)]
  if len(vals) <> n_crd:
    raise ValueError, "Improper number of coordinates in rst file."

  convert_amber_vel_to_angstroms_per_ps = 20.455
  for i, atom in enumerate(sorted(soup.atoms(), pdbstruct.cmp_atom)):
    atom.vel.set(vals[i*3], vals[i*3+1], vals[i*3+2])
    atom.vel.scale(convert_amber_vel_to_angstroms_per_ps)

  f.close()


def get_restart_files(name):
  top = os.path.abspath(name + '.top')
  crds = os.path.abspath(name + '.rst')
  if not os.path.isfile(crds):
    crds = os.path.abspath(name + '.crd')
  vels = ''
  return top, crds, vels

  
def SoupFromAmberTop(top, crds_or_rst):
  "Make a protein from an AMBER structure"

  temp_pdb = util.temp_fname('.pdb')
  convert_amber_to_pdb(top, crds_or_rst, temp_pdb)
  soup = pdbstruct.Soup(temp_pdb)
  if crds_or_rst.endswith('.rst'):
    load_rst_into_soup(soup, crds_or_rst)

  soup.box_dimension_str = None
  if is_top_has_box_dimension(top):
    lines = open(crds_or_rst, "r").readlines()
    i = -1
    while lines[i].strip() == "":
      i -= 1
    soup.box_dimension_str = lines[i][:-1]

  masses = get_attribute_list(top, 'MASS', 16)
  for m, a in zip(masses, soup.atoms()):
    a.mass = float(m)

  os.remove(temp_pdb)
  return soup


def soup_from_restart_files(top, crds, vels):
  return SoupFromAmberTop(top, crds)


def write_soup_to_rst(soup, rst):
  f = open(rst, "w")

  f.write(" ".ljust(80) + "\n")
  f.write("%5d  0.0000000E+00\n" % len(soup.atoms()))

  i = 0
  for atom in sorted(soup.atoms(), pdbstruct.cmp_atom):
    f.write("%12.7f%12.7f%12.7f" % (atom.pos.x, atom.pos.y, atom.pos.z))
    i += 1
    if i % 2 == 0:
      f.write("\n")
      i = 0
  if len(soup.atoms()) % 2 <> 0:
    f.write("\n")

  i = 0
  convert_to_amber_vel = 1.0 / 20.455
  for atom in sorted(soup.atoms(), pdbstruct.cmp_atom):
    vx = atom.vel.x * convert_to_amber_vel
    vy = atom.vel.y * convert_to_amber_vel
    vz = atom.vel.z * convert_to_amber_vel
    f.write("%12.7f%12.7f%12.7f" % (vx, vy, vz))
    i += 1
    if i % 2 == 0:
      f.write("\n")
  if len(soup.atoms()) % 2 <> 0:
    f.write("\n")

  if soup.box_dimension_str is not None:
    f.write(soup.box_dimension_str + '\n')
    
  f.close()
  
  
def write_soup_to_restart_files(soup, name):
  write_soup_to_rst(soup, name + '.rst')
  return name + '.rst', ''
  

minimization_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.crd', 
  'output_name' : 'min', 
  'force_field': 'GBSA',
  'solvent_state': 2,
  'surface_area': 1,
  'cutoff': 16.0,  # non-bonded cutoff
  'periodic': 0,  # none: 0, constant-volume: 1, constant-pressure: 2
  'constraint_pdb': '',
  'n_step_minimization' : 100, 
  'n_step_steepest_descent' : 10, 
  'shake_on': False,
} 

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.crd', 
  'output_name' : 'md', 
  'force_field': 'GBSA',
  'solvent_state': 2,
  'surface_area': 1,
  'constraint_pdb': '',
  'cutoff': 16.0,  # non-bonded cutoff
  'periodic': 0,  # none: 0, constant-volume: 1, constant-pressure: 2
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'shake_on': False,
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds' : 'in.crd', 
  'output_name' : 'md', 
  'force_field': 'GBSA',
  'solvent_state': 2,
  'surface_area': 1,
  'cutoff': 16.0,  # non-bonded cutoff
  'periodic': 0,  # none: 0, constant-volume: 1, constant-pressure: 2
  'constraint_pdb': '',
  'random_seed' : 2342, 
  'thermometer_type': 3, # contant-energy: 0, berendsen: 1, andersen: 2, langevin: 3
  'gamma_ln': 5,   
  'temp_thermometer' : 300.0, 
  'temp_initial': 0.0, # ignored if it is 0.0
  'pressure': 0, # no-pressure-scaling: 0, isotropic-pressure-scaling: 1
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'n_step_per_thermostat' : 100, 
  'shake_on': False,
} 


def make_sander_input_file(parms):
  lines = []
  lines.append("amber.py generated sander input file") 
  lines.append("&cntrl") 
  if parms['shake_on']:
    lines.append("  ntf = 2, ntb = $periodic, cut = $cutoff,") 
  else:
    lines.append("  ntf = 1, ntb = $periodic, cut = $cutoff,") 
  lines.append("  igb = $solvent_state, gbsa = $surface_area, nrespa = 1,") 
  if parms['shake_on']:
    lines.append("  ntc = 2, tol = 1.0d-8,") 
  else:
    lines.append("  ntc = 1") 
  if parms['constraint_pdb']:
    lines.append("  ntr = 1")
  if 'n_step_minimization' in parms:
    lines.append("  imin = 1, nmropt = 0,")
    lines.append("  ntmin = 1, maxcyc = $n_step_minimization, ncyc = $n_step_steepest_descent,") 
  elif 'n_step_dynamics' in parms:
    if parms['input_crds'].endswith('.rst'):
      lines.append("  ntx = 5, irest = 1,")
    else:
      lines.append("  ntx = 1,")
    lines.append("  ntpr = $n_step_per_snapshot, ntave = 0, ntwr = 500, iwrap = 0, ioutfm = 0, "
                         "ntwx = $n_step_per_snapshot, ntwv = $n_step_per_snapshot, "
                         "ntwe = $n_step_per_snapshot, ")
    lines.append("  nstlim = $n_step_dynamics, nscm = 50, dt = 0.001,") 
    if parms['periodic'] == 2:
      lines.append("  ntp = 1,") 
    if 'thermometer_type' in parms:
      lines.append("  ntt = $thermometer_type, ig = $random_seed, " + \
                   "temp0 = $temp_thermometer, tempi = $temp_initial, vlimit = 0.0, ")
      if parms['thermometer_type'] == 2:
        lines.append("  vrand = $steps_per_thermostat, ")
      elif parms['thermometer_type'] == 3:
        lines.append("  gamma_ln = $gamma_ln, ")
  else:
    raise "Can't parse parameters to run "
  lines.append("/") 
  if parms['constraint_pdb']:
    lines.append("Restrained atoms from %s" % parms['constraint_pdb'])
    restraint_weight = 20.0
    lines.append(str(restraint_weight))
    lines.append("FIND")
    lines.append("* * S *")
    lines.append("* * B *")
    lines.append("* * 3 *")
    lines.append("* * E *")
    lines.append("SEARCH")
    mol = pdbstruct.Molecule()
    mol.read_pdb(parms['constraint_pdb'])
    for i, atom in enumerate(mol.atoms()):
      if atom.bfactor > 0.0:
        lines.append("ATOM %d %d" % (i+1, i+1))
    lines.append("END")
    lines.append("END")
  lines.append("")
  return util.replace_dict("\n".join(lines), parms)  


def run(in_parms):
  "Runs AMBER with in_parms"

  parms = copy.deepcopy(in_parms)
  name = parms['output_name']
  config = name + ".config"
  if util.is_same_dict_in_file(in_parms, config):
    print "simulation already run."
    return

  input_top = parms['topology'] 
  new_top = name + '.top'
  shutil.copy(input_top, new_top)

  input_crd = parms['input_crds']
  if input_crd.endswith('.crd'): 
    new_crd = name + '.in.crd'
  else:
    new_crd = name + '.in.rst'
  shutil.copy(input_crd, new_crd)
  
  if 'n_step_minimization' in parms:
    rst = name + ".crd"
  else:
    rst = name + ".rst"
  trj = name + ".trj"
  vel = name + ".vel"
  ene = name + ".ene"
  inf = name + ".inf"
  sander_out = name + ".sander.out"
  sander_in = name + ".sander.in"

  open(sander_in, "w").write(make_sander_input_file(parms))
  cmd = "%s -O -i %s -o %s -p %s -c %s -r %s -x %s -v %s -e %s -inf %s" \
          % (sander, sander_in, sander_out, new_top, new_crd, rst, trj, vel, ene, inf)
  if parms['constraint_pdb']:
      cmd += " -ref %s" % parms['constraint_pdb'].replace('.pdb', '.crd')
  util.run_with_output_file(cmd, name)


def read_time_blocks(fname):
  "Returns a list of dictionaries containing energy values from sander out file"

  allowed_keys = ['EKtot', 'Etot', 'TEMP(K)', 
                  'TIME(PS)', 'EPtot', 'NSTEP']

  def extract_to_dict(line, ext_dict):
    words = line.split()
    for i in range(len(words)):
      if words[i] == "=":
        key, value = words[i-1], words[i+1]
        if key in allowed_keys:
          ext_dict[key] = value

  blocks = []
  block_dict = {}
  is_results = False
  for line in open(fname):
    if is_results:
      if 'A V E R A G E S' in line:
        is_results = False
      elif '----' in line:
        if len(block_dict.keys()) > 0:
          blocks.append(block_dict.copy())
      else:
        extract_to_dict(line, block_dict)    
    else:
      if '4' in line and 'RESULTS' in line:
        is_results = True
        i = -1
  return blocks


class TrjReader:
  """
  Class that reads successive sets of coordinates from Amber trajectory (or
  velocity) files. Note: for coordinate files, the box dimensions (3 floats)
  are included at the end of each frame, so this is taken care of in the code
  to read the frame, in order to calculate the correct frame size. However,
  the box dimensions are *not* included in velocity coordinate files.
  """
  
  def __init__(self, top, trj):
    """
    trj: name of amber trajectory (or velocity) file containing sets of
         coordinates 
    top: name of amber top file of trajectory
    """
    
    self.trj = trj
    self.top = top

    if '.trj' in trj:
      self.is_skip_box_dims = is_top_has_box_dimension(self.top)
    else:
      self.is_skip_box_dims = False

    if self.trj.split(".")[-1].strip().lower() == "gz":
      self._file_open = gzip.GzipFile
    else:
      self._file_open = open
    
    self.file = self._file_open(self.trj, "r")

    # skip header line
    self.pos_start_frame = len(self.file.readline())

    # calculate the size of each frame
    self.n_atom = len(get_attribute_list(self.top, 'ATOM_NAME', 4))
    if self.n_atom == 0:
      raise "No names found in .top file"

    n_line = (3 * self.n_atom) / 10
    if (3 * self.n_atom) % 10 > 0: 
      n_line += 1

    if self.is_skip_box_dims:
      n_line += 1

    self.size_frame = 0
    for i in range(0, n_line):
      self.size_frame += len(self.file.readline())
      
    # calculate n_frame
    self.file.seek(0, 2)
    pos_eof = self.file.tell()
    self.n_frame = int((pos_eof - self.pos_start_frame) / self.size_frame)

  def __getitem__(self, i):
    "Gets the list of floats for frame i."

    # Find place in file
    if i < - 1*self.n_frame or i >= self.n_frame :
      raise IndexError
    elif i < 0 :
      return self.__getitem__(self.n_frame + i)
    else :
      self.file.seek(self.pos_start_frame + i*(self.size_frame))

    # read frame as a string
    s = self.file.read(self.size_frame).replace("\n", "")

    # convert string to float
    vals = [float(s[i:i+8]) for i in xrange(0, len(s), 8)]
    
    if self.is_skip_box_dims:
      vals = vals[:-3]

    if len(vals) <> 3 * self.n_atom:
      raise ValueError, "Improper number of coordinates in frame."
    
    return vals

  def __len__(self):
    return self.n_frame
    
  def __del__(self):
    self.file.close()

  def __repr__(self):
    return "< Amber Coord file %s with %d frames of %d atoms >" % \
             (self.trj, self.n_frame, self.n_atom)
    

class Trajectory:
  """
  Class to hold chains of a pdb that changes according to a dcd
  
  Data:
    soup: collection of chains that get updated
    n_frame: number of frames in trajectory
  
  Methods:
    load_frame: loads the i'th frame from dcdectory
  """
  def __init__(self, name):
    self.topology = name + '.top'
    coor_trj_fname = name + '.trj'
    self.coor_traj = TrjReader(self.topology, coor_trj_fname)
    vel_trj_fname = name + '.vel'
    if os.path.isfile(vel_trj_fname):
      self.vel_traj = TrjReader(self.topology, vel_trj_fname)
    else:
      self.vel_traj = None
    save_crds('temp.crd', self.coor_traj[0])
    self.soup = SoupFromAmberTop(self.topology, 'temp.crd')
    os.remove('temp.crd')
    self.n_frame = len(self.coor_traj)
    self.i_frame = 0
    self.load_frame(self.i_frame)

  def load_frame(self, i):
    if i < -self.n_frame or i >= self.n_frame :
      raise IndexError
    if i < 0 :
      return self.load_frame(self.n_frame + i)
    crds = self.coor_traj[i]
    for j, a in enumerate(self.soup.atoms()):
      k = 3*j
      a.pos.set(crds[k], crds[k+1], crds[k+2])
    if self.vel_traj is not None:
      vels = self.vel_traj[i]
      for j, a in enumerate(self.soup.atoms()):
        k = 3*j
        a.vel.set(vels[k], vels[k+1], vels[k+2])
    self.i_frame = i


def merge_trajectories(top, trajs, out_traj):
  """
  Given a list of traj filenames (trajs), merges them into one complete
  trajectory (out_traj) using top to work out the number of atoms, and
  hence the size of the frame of the trajectory.
  """
  
  trj_reader = TrjReader(top, trajs[0])
  pos_start_frame = trj_reader.pos_start_frame  
  size_frame = trj_reader.size_frame
  del trj_reader

  shutil.copy(trajs[0], out_traj)

  merge_traj_file = open(out_traj, "ab+")
  for traj in trajs[1:]:
    traj_file = open(traj, "rb")
    traj_file.seek(-1, 2)
    eof = traj_file.tell()
    traj_file.seek(pos_start_frame)
    while traj_file.tell() < eof:
      merge_traj_file.write(traj_file.read(size_frame)) 
    traj_file.close()
  merge_traj_file.close()


def merge_simulations(name, pulses):
  """
  Given a list of directories with partial trajectories in each directory
  with the same name for the md, will splice them together into one uber
  simulation.
  """
  
  shutil.copy(os.path.join(pulses[0], name + '.sander.in'), name + '.sander.in')
  shutil.copy(os.path.join(pulses[0], name + '.top'), name + '.top')
  shutil.copy(os.path.join(pulses[-1], name + '.rst'), name + '.rst')

  # merge energies of pulses into one energy file
  f = open(name + '.energy', 'w')
  f.write('[\n')
  n_step = 0
  time = 0.0
  for pulse in pulses:
    energy_fname = os.path.join(pulse, name + '.energy')
    if os.path.isfile(energy_fname):
      blocks = eval(open(energy_fname).read())
    else:
      sander_out = os.path.join(pulse, name + '.sander.out')
      blocks = read_time_blocks(sander_out)
      for block in blocks:
        block_n_step = int(block['NSTEP'])
        block_time = float(block['TIME(PS)'])
        block['NSTEP'] = str(block_n_step + n_step)
        block['TIME(PS)'] = str(block_time + time)
        f.write(str(block) + ',\n')
    n_step = int(blocks[-1]['NSTEP'])
    time = float(blocks[-1]['TIME(PS)'])
  f.write(']\n')
  f.close()
    
  trajs = [os.path.join(pulse, name + '.trj') for pulse in pulses]
  merge_trajectories(name + '.top', trajs, name + '.trj')

  vels = [os.path.join(pulse, name + '.vel') for pulse in pulses]
  merge_trajectories(name + '.top', vels, name + '.vel')
  

def get_energy(top, crd):
  "Returns potential energy of top and crd by running sander"
  
  top = os.path.abspath(top)
  crd = os.path.abspath(crd)

  util.goto_dir('energy-temp')

  parms = minimization_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'energy'
  parms['n_step_minimization'] = 0
  parms['n_step_steepest_descent'] = 0

  run(parms)

  lines = open('energy.sander.out').readlines()

  energies = {}
  is_results = False
  for i, line in enumerate(lines):
    if not is_results:
      if '4' in line and 'RESULTS' in line:
        is_results = True
    else:
      if 'NSTEP' in line and 'ENERGY' in line:
        util.goto_dir('..')
        os.system('rm -rf energy-temp')
        return float(lines[i+1].split()[1])

  raise "Can't find energies"


def preequilibrate(in_name, out_name, temperature):
  """
  The problem with equilibration is that since at low
  temperatures the simulations are run at constant energy, we
  need to prequilibrate at constant energy as well. With
  implicit solvation, this is an unstable process as the system
  quickly rises to 20 K regardless of the starting conformation.
  We need to let the system find a better minimum *at constant
  energy*.
  """

  top, crd, dummy = get_restart_files(in_name)
  md_dirs = ['heat1', 'const2', 'heat3']

  util.goto_dir(md_dirs[0])
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'md'
  parms['temp_thermometer'] = temperature
  parms['temp_initial'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)
  in_top, in_crds, in_vels = get_restart_files('md')

  util.goto_dir(os.path.join('..', md_dirs[1]))
  parms = constant_energy_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'md'
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 10000
  run(parms)
  in_top, in_crds, in_vels = get_restart_files('md')

  util.goto_dir(os.path.join('..', md_dirs[2]))
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = crd
  parms['output_name'] = 'md'
  parms['temp_thermometer'] = temperature
  parms['temp_initial'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)
  
  util.goto_dir('..')
  merge_simulations('md', md_dirs)


