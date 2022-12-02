#!/usr/bin/env python

import os
import sys
import copy
import struct
import shutil
import util
import pdbstruct.molecule as molecule
import pdbstruct.polymer as polymer
import pdbstruct.protein as protein
from pdbstruct.soup import Soup
import pdbtext
import pdbstruct

vmd = '/Users/bosco/Applications/VMD\ 1.8.7.app/Contents/vmd/vmd_MACOSXX86'
psfgen = 'psfgen'
namd = 'namd2'


module_dir = os.path.dirname(__file__) 
if module_dir is "": 
  module_dir = "."


def convert_namd_pdb_to_pdb(namd_pdb, standard_pdb):
  "doesn't handle hydrogens so far"
  
  lines = open(namd_pdb, "r").readlines()
  f = open(standard_pdb, "w")
  for line in lines:
    if line.startswith("ATOM"):
      res_type = line[17:21]
      atom_type = line[12:16]
      
      new_atom_type = ""

      if atom_type[0] == "H":
        new_atom_type = atom_type[-1] + atom_type[:-1]
      if atom_type[1:] == "OT1":
        new_atom_type = " O  "
      if atom_type[1:] == "OT2":
        new_atom_type = " OXT"
      if atom_type[1:] == "OH2":
        new_atom_type = " O  "
      if res_type == "ILE " and atom_type == " CD ":
        new_atom_type = " CD1"
      if new_atom_type:
        line = line[:12] + new_atom_type + line[16:]

      if res_type == "HSE ":
        line = line[:17] + "HIS " + line[21:]

      if res_type == "TIP3":
        line = "HETATM" + line[6:17] + "HOH " + line[21:]
        if atom_type == " OH2":        
          line = line[:12] + " O  " + line[16:]
          
      f.write(line)
  f.close()
  

def convert_pdb_to_namd_pdb(standard_pdb, namd_pdb):
  soup = Soup(standard_pdb)
  for res in soup.residues():
    if res.type == "ILE" and res.has_atom('CD1'):
      res.atom('CD1').type = 'CD'
    if res.has_atom('OXT'):
      res.atom('OXT').type = 'OT2'
      if res.has_atom('O'):
        res.atom('O').type = 'OT1'
    for atom in res.atoms():
      if atom.type[0].isdigit() and atom.type[1] == "H":
        atom.type = atom.type[1:] + atom.type[0]
  soup.write_pdb(namd_pdb)


def load_masses_to_soup(psf, soup):
  lines = open(psf, "r").readlines()
  n_atom = 0
  for i, line in enumerate(lines):
    if "NATOM" in line:
      n_atom = int(line.split()[0])
      start = i+1
      break
  if n_atom == 0:
    return
  masses = [float(line.split()[7]) 
            for line in lines[start:start+n_atom]]
  for a, m in zip(soup.atoms(), masses):
    a.mass = m 
    

def SoupFromNamdPsf(psf, in_coor_pdb, in_vel_pdb):
  soup = Soup(in_coor_pdb)
  load_masses_to_soup(psf, soup)
  vel_soup = Soup(in_vel_pdb)
  for atom, vel_atom in zip(soup.atoms(), vel_soup.atoms()):
    v = vel_atom.pos
    atom.vel.set(v.x, v.y, v.z)
  return soup


def get_restart_files(name):
  top = os.path.abspath(name + '.psf')
  crds = os.path.abspath(name + '.coor')
  vels = os.path.abspath(name + '.vel')
  return top, crds, vels


def soup_from_restart_files(top, crds, vels):
  return SoupFromNamdPsf(top, crds, vels)  
  
  
def write_soup_to_vel_pdb(in_soup, out_vel_pdb):
  soup = copy.deepcopy(in_soup)
  for atom in soup.atoms():
    atom.pos.set(atom.vel.x, atom.vel.y, atom.vel.z)
  soup.write_pdb(out_vel_pdb)


def write_soup_to_restart_files(soup, name):
  write_soup_to_vel_pdb(soup, name + '.vel')
  soup.write_pdb(name + '.coor')
  return name + '.coor', name + '.vel'


def clean_pdb(pdb, out_pdb):
  def is_non_water_hetatm(line):
    if not line.startswith("HETATM"):
      return False
    if line[17:20] != 'HOH':
      return True
    return False
  txt = open(pdb, 'r').read()
  txt = pdbtext.strip_other_nmr_models(txt)
  txt = pdbtext.strip_lines(txt, is_non_water_hetatm)
  txt = pdbtext.strip_lines(
      txt, lambda l: l.startswith('ANISOU'))
  txt = pdbtext.strip_lines(
      txt, lambda l: l.startswith('CONECT'))
  txt = pdbtext.strip_lines(
      txt, lambda l: l.startswith('MASTER'))
  txt = pdbtext.strip_alternative_atoms(txt)
  txt = pdbtext.strip_hydrogens(txt)
  txt = pdbtext.renumber_residues(txt)
  f = open(out_pdb, 'w')
  f.write(txt)
  f.close()


solvate_tcl = """
# Solvate
# Set minimum padding
set pad 13

# Run solvate with automatic padding option
package require solvate
resetpsf
solvate $in_psf $in_pdb -o $name.solv -rotate -rotinc 5 -t $pad

# Find the periodic box size.
mol load psf $name.solv.psf
set mol1 [molinfo top get id]
animate read pdb $name.solv.pdb $mol1

set sel1 [atomselect $mol1 all]
set cent1 [measure center $sel1]
set size1 [measure minmax $sel1]
mol delete $mol1

set sizemin [lindex $size1 0]
set sizemax [lindex $size1 1]

set xsize [expr [lindex $sizemax 0]-[lindex $sizemin 0]]
set ysize [expr [lindex $sizemax 1]-[lindex $sizemin 1]]
set zsize [expr [lindex $sizemax 2]-[lindex $sizemin 2]]

# Find the padding in each direction to make the box cubic
set lmax $xsize
if {$ysize > $lmax} {
	set lmax $ysize
}
if {$zsize > $lmax} {
	set lmax $zsize
}

# I like my boxe size to be a nice even number
set maxsize [expr int(ceil($lmax))]
if {$maxsize%2 != 0} { set maxsize [expr $maxsize +1] }

# Calculate additional padding
set xpad [expr {$pad+0.5*($maxsize-$xsize)}]
set ypad [expr {$pad+0.5*($maxsize-$ysize)}]
set zpad [expr {$pad+0.5*($maxsize-$zsize)}]

puts ":Padding: $xpad $ypad $zpad"
puts ":Box size: $lmax"

# Adjust padding for nonzero center of mass. These are used to manually set the padding in each direction (both + and -)
set xplus [expr $xpad - [lindex $cent1 0]]
set xmin [expr $xpad + [lindex $cent1 0]]
set yplus [expr $ypad - [lindex $cent1 1]]
set ymin [expr $ypad + [lindex $cent1 1]]
set zplus [expr $zpad - [lindex $cent1 2]]
set zmin [expr $zpad + [lindex $cent1 2]]

# Rerun solvate on the original structure using calculated padding to make the box cubic
resetpsf
solvate $name.psf $name.pdb -o $name.solv -rotate -rotinc 5 -x $xmin +x $xplus -y $ymin +y $yplus -z $zmin +z $zplus

# Check that it worked
mol delete all
mol load psf $name.solv.psf
set mol1 [molinfo top get id]
animate read pdb $name.solv.pdb $mol1

set sel1 [atomselect $mol1 all]
set cent1 [measure center $sel1]
set size1 [measure minmax $sel1]

set sizemin [lindex $size1 0]
set sizemax [lindex $size1 1]

set xsize [expr [lindex $sizemax 0]-[lindex $sizemin 0]]
set ysize [expr [lindex $sizemax 1]-[lindex $sizemin 1]]
set zsize [expr [lindex $sizemax 2]-[lindex $sizemin 2]]

puts ":Final size: $xsize, $ysize, $zsize"
puts ":Length: $sizemax"
puts ":Center: $cent1"
puts ":Size: $size1"

# Ionize
package require autoionize
# Find original charge
set all [atomselect $mol1 all]
set q [measure sumweights $all weight charge]

# Determine the number and type of ions to use
set natom [expr round($q)]
if {$natom < 0} {
	set natom [expr abs($natom)]
	autoionize -psf $name.solv.psf -pdb $name.solv.pdb -nna $natom -ncl 0 -o $name
} elseif {$natom > 0} {
	autoionize -psf $name.solv.psf -pdb $name.solv.pdb -nna 0 -ncl $natom -o $name
} elseif {$natom == 0} {
	exec cp $name.solv.psf $name.psf
	exec cp $name.solv.pdb $name.pdb
}

# Check that it worked
mol delete all
mol load psf $name.psf
set mol1 [molinfo top get id]
animate read pdb $name.pdb $mol1
set all [atomselect $mol1 all]
puts ":Old charge: $q, New charge: [measure sumweights $all weight charge]"
mol delete all
"""

def solvate_psf(in_psf, in_pdb, out_name):
  parms = {
      'in_psf': in_psf,
      'in_pdb': in_pdb,
      'name': out_name,
  }
  tcl = util.replace_dict(solvate_tcl, parms)
  tcl_name = out_name + '.tcl'
  open(tcl_name, 'w').write(tcl)
  util.run_with_output_file(
      '%s -dispdev text -eofexit < %s'  % (vmd, tcl_name),
      out_name + '.solvate')
  log_name = out_name + '.solvate.log'


# Routines to generate psf files
  
_psfgen_template = """
package require psfgen 
topology $module_dir/$topology
alias residue HIS HSE 
alias atom ILE CD1 CD 
# insert protein
# insert water
guesscoord 
writepdb $out_pdb 
writepsf $psf 
"""

_chain_psfgen_template = """
segment $chain_id { 
 pdb $chain_pdb 
} 
coordpdb $chain_pdb $chain_id
"""

_water_psfgen_template = """
pdbalias residue HOH TIP3
pdbalias residue WAT TIP3
segment wat { 
 auto none
 pdb $water_pdb 
} 
pdbalias atom HOH O OH2 
pdbalias atom WAT O OH2 
coordpdb $water_pdb wat
"""

def make_namd_from_pdb(pdb, psf, out_pdb): 
  """
  Creates charmm .pdb and .psf file. Can only make NAMD
  generate CHARMM topology files, not OPLS ones - that
  requires XPLOR but I don't care. Still, if OPLS
  topology files are provided - can still run.
  """
  name = psf.replace('.psf', '')
  in_pdb = name + '.clean.pdb' 
  clean_pdb(pdb, in_pdb)
  soup = Soup(in_pdb)
  chains = soup.chains()
  water_names = ['WAT', "TIP", "TIP3", "HOH"]
  waters = [chain for chain in chains 
            if chain.n_residue() > 0 and \
              chain.residue(0).type in water_names]
  chains = [chain for chain in chains 
            if chain.n_residue() > 0 and \
              chain.residue(0).type not in water_names]
  
  raw_psf = name + '.raw.psf'
  raw_pdb = name + '.raw.pdb'
  replace = { 'pdb': in_pdb, 
              'out_pdb': raw_pdb, 
              'psf': raw_psf,
              'module_dir': module_dir,
              'topology': 'parms/charmm22.topology'
             }
  template = _psfgen_template
  if len(waters) > 0:
    template = template.replace(
        "# insert water", _water_psfgen_template)
    water_pdb = name + ".waters.pdb"
    water_soup = Soup()
    for water in waters:
      water_soup.insert_chain(water)
    water_soup.write_pdb(water_pdb)
    replace['water_pdb'] = water_pdb
  chains_template = ""
  for i, chain in enumerate(chains):
    id = 'ch%d' % i
    pdb = name + '.' + id + '.pdb'
    chain_soup = Soup()
    chain_soup.append_chain(chain)
    chain_soup.write_pdb(pdb)
    chain_replace = { 'chain_id': id, 'chain_pdb': pdb }
    chain_template = _chain_psfgen_template
    chains_template += util.replace_dict(
        chain_template, chain_replace)
  template = template.replace(
      "# insert protein", chains_template)     
  template = util.replace_dict(template, replace)

  psfgen_in = name + ".psfgen.in"
  open(psfgen_in, "w").write(template)
  util.run_with_output_file(
    "%s %s" % (psfgen, psfgen_in), name + '.psfgen')

  solvate_psf(raw_psf, raw_pdb, 'sim')

  os.remove(in_pdb)


def SoupFromPsf(psf):
  "read a psf and parse into protein chains."

  f = open(psf, "r")
  line = f.readline()
  while 'NATOM' not in line:
    line = f.readline()
  n_atom = int(line[0:8])  
  lines = [f.readline() for i in range(0, n_atom)]
  f.close
  
  # work out the chains in the psf file
  chain_info_list = []
  chain_char = ""
  res_num = -1
  chain_type = ""
  for i, line in enumerate(lines):
    word_list = line.split()
    new_chain_type = word_list[1]
    new_res_num = int(word_list[2])
    if chain_type == "":
      chain_type = new_chain_type
      res_num = new_res_num
      chain_info_list.append([chain_type, i, i+1])
    elif chain_type is "WATERS" and \
        new_chain_type is "WATERS" and \
        res_num <> new_res_num:
      res_num = new_res_num
      chain_info_list.append([chain_type, i, i+1])
    elif chain_type <> new_chain_type:
      chain_type = new_chain_type
      chain_info_list.append([chain_type, i, i+1])
    else:
      chain_info_list[-1][2] = i+1

  # make a list of chains using the objects in protein
  chains = []
  for name, first, last in chain_info_list:
    if "WATERS" in name:
      chain = polymer.Polymer()
    else:
      chain = protein.Protein()
  
    curr_res_num = -1
    for i, line in enumerate(lines[first:last]):
        word_list = line.split()
        atom = molecule.Atom()
        atom.res_num = int(word_list[2])
        atom.res_type = word_list[3][0:3]
        atom.type = word_list[4]
        atom.element = atom.type[0]
        atom.type = atom.type.strip(" ")
        atom.num = i + first
        try:
          atom.mass = float(word_list[7])
        except:
          atom.mass = 0.0

        if curr_res_num != atom.res_num:
          res = polymer.Residue(
              atom.res_type, atom.chain_id, atom.res_num, '')
          chain.append_residue_no_renum(res)
          curr_res_num = atom.res_num
        
        chain.insert_atom(chain.n_residue()-1, atom)
    chains.append(chain)
    
  soup = Soup()
  for chain in chains:
    soup.append_chain(chain)

  # load masses in [a.m.u.]
  masses = [float(line[47:60]) for line in open(psf, "r")
            if "CHAIN" in line]
  for a, m in zip(soup.atoms(), masses):
    a.mass = m 

  return soup


# Routines to prepare namd minimization/molecular-dynamics runs

minimization_parms = { 
  'topology' : 'in.psf', 
  'input_crds' : 'in.pdb', 
  'output_name' : 'min', 
  'force_field': 'NAMD', 
  'constraint_pdb': ''
} 

constant_energy_parms = { 
  'topology' : 'in.top', 
  'input_crds': 'in.coor',
  'input_vels': 'in.vel',
  'output_name' : 'md', 
  'random_seed' : 2342, 
  'force_field': 'NAMD', 
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'constraint_pdb': '',
} 

langevin_thermometer_parms = { 
  'topology' : 'in.top', 
  'input_crds': 'in.coor',
  # use either vel or temp_background to set initial velocities
  'input_vels': '',
  'temp_thermometer' : 300.0, 
  'force_field': 'NAMD', 
  'output_name' : 'md', 
  'random_seed' : 2342, 
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'constraint_pdb': '',
} 

couple_parms = { 
  'topology' : 'in.top', 
  'pdb' : 'in.crd', 
  'input_crds': 'in.coor',
  'input_vels': 'in.vel',
  'force_field': 'NAMD',
  'output_name' : 'couple', 
  'random_seed' : 2342, 
  'couple_pdb': 'couple.pdb',
  'temp_couple' : 300.0, 
  'n_step_per_snapshot' : 5, 
  'n_step_dynamics' : 1000, 
  'constraint_pdb': '',
} 


_constraint_template="""
# constraint
constraints on
consExp 2
consKFile $constraint_pdb
consKCol B
consRef $constraint_pdb
"""

md_template = """# dcdectory output
dcdFile $output_name.dcd
dcdFreq $n_step_per_snapshot
velDcdFile $output_name.vel.dcd
velDcdFreq $n_step_per_snapshot

# molecular dynamics timesteps
numSteps $n_step_dynamics
timeStep 1  # fs
firstTimeStep 0
stepsPerCycle 20
nonBondedFreq 1
fullElectFrequency 2
"""

force_field_template = """
# forcefield parameters
$psf_type
parameters $module_dir/$parameter
exclude	scaled1-4
1-4scaling 1.0
dielectric 1.0

# cutoffs
switching on
switchDist 10.0
cutoff 12.0
pairListDist 14.0

# SHAKE
rigidBonds water
molly on

# Periodic Boundary Conditions

wrapWater           on
wrapAll             on
wrapNearest         on

# PME (for full-system periodic electrostatics)
PME                 yes
PMEGridSpacing      1.0
"""

box_template = """
# Periodic Box Definition
cellBasisVector1    $len_x     0.0        0.0
cellBasisVector2    0.0        $len_y     0.0
cellBasisVector3    0.0        0.0        $len_z
cellOrigin          $x_origin  $y_origin  $z_origin
"""
def periodic_box_template(pdb):
  template = box_template
  p = pdbstruct.Polymer(pdb)
  atoms = p.atoms()
  parms = {}
  for v in ['x', 'y', 'z']:
    vals = [a.pos.__dict__[v] for a in atoms]
    v_min, v_max = min(vals), max(vals)
    parms["len_"+v] = v_max - v_min + 0.5
    parms[v+"_origin"] = sum(vals)/float(len(vals))
  return util.replace_dict(template, parms)

  
pressure_template = """
# Constant Pressure Control (variable volume)
useGroupPressure      no ;# needed for rigidBonds
useFlexibleCell       no
useConstantArea       no

langevinPiston        on
langevinPistonTarget  1.01325 ;#  in bar -> 1 atm
langevinPistonPeriod  50.0
langevinPistonDecay   25.0
langevinPistonTemp    $temperature
"""
  
def make_namd_input_file(parms):
  lines = []

  name = parms['output_name']
  lines.append("# initial structure and coordinates")
  lines.append("structure $topology")
  lines.append("coordinates $input_crds")
  lines.append("")
  
  lines.append("# output")
  lines.append("outputName $output_name")
  lines.append("binaryOutput no")
  lines.append("")

  shutil.copy(parms['topology'], name + '.psf')
  if parms['constraint_pdb']:
    shutil.copy(parms['constraint_pdb'], name + '.constraint')

  lines.extend(force_field_template.splitlines())

  if 'constraint_pdb' in parms:
    if parms['constraint_pdb']:
      lines.extend(_constraint_template.splitlines())

  if 'n_step_dynamics' in parms:
    lines.append("# initial temperature/velocities")
    if 'input_vels' in parms and parms['input_vels'] <> "":
      lines.append("velocities $input_vels")
    elif 'temp_init' in parms and parms['temp_init'] > 0.0:
      lines.append("temperature $temp_init")
      lines.append("seed $random_seed")
    else:
      raise IndexError, \
          "No initial velocity information for dynamics run"

  if 'n_step_dynamics' in parms:
    lines.extend(md_template.splitlines())
  else:
    lines.append("# minimization parameters")  
    lines.append("minimization on")

  xsc = parms['topology'].replace('.psf', '.xsc')
  if os.path.isfile(xsc):
    new_xsc = name + '.in.xsc'
    shutil.copy(xsc, new_xsc)
    lines.append("extendedSystem %s" % new_xsc)
  else:
    box_template = periodic_box_template(parms['input_crds'])
    lines.extend(box_template.splitlines())

  if 'temp_thermometer' in parms:
    lines.append("# temperature coupling")
    lines.append("langevin on")
    lines.append("langevinHydrogen on")
    if 'temp_couple' in parms:
      lines.append("langevinTemp $temp_couple")
      lines.append("langevinFile $couple_pdb")
      lines.append("langevinCol B")
    else:
      lines.append("langevinTemp $temp_thermometer")
      lines.append("langevinDamping 10.0")
    txt = util.replace_dict(
        pressure_template, 
        {'temperature': parms['temp_thermometer']})
    lines.extend(txt.splitlines())


  return util.replace_dict("\n".join(lines), parms)  


def run(in_parms):
  """
  Read parms and creates the appropriate NAMD input files
  for simulation
  """
  name = in_parms['output_name']
  namd_in = name + ".in"
  config = name + ".config"

  parms = copy.deepcopy(in_parms)
  if util.is_same_dict_in_file(parms, config):
    print "simulation already run."
    return

  parms['parameter'] = 'parms/charmm22.parameter'
  parms['psf_type'] =  'paraTypeCharmm on'

  parms['module_dir'] = module_dir

  open(namd_in, "w").write(make_namd_input_file(parms))
  
  stopwatch = util.Timer()
  util.run_with_output_file(
      "%s %s" % (namd, namd_in), name + '.namd')
  stopwatch.stop()
  open(name + '.time', 'w').write(stopwatch.str())

  for line in open(name + '.namd.log', "r"):
    if 'ERROR' in line and 'Info' not in line: 
      raise "NAMD failure: '%s'" % line.replace("\n", "")

  util.write_dict(in_parms, config)


def pulse(in_parms, n_step_per_pulse, reset_vel_func):

  config = in_parms['output_name'] + ".config"
  if util.is_same_dict_in_file(in_parms, config):
    print "simulation already run."
    return

  name = in_parms['output_name']
  shutil.copy(in_parms['topology'], name + '.psf')
  if in_parms['constraint_pdb']:
    shutil.copy(in_parms['constraint_pdb'], name + '.constraint')

  n_pulse = in_parms['n_step_dynamics'] // n_step_per_pulse
  n_step_list = [n_step_per_pulse for i in range(n_pulse)]
  n_excess_step = in_parms['n_step_dynamics'] % n_step_per_pulse
  if n_excess_step > 0:
    n_pulse += 1
    n_step_list.append(n_excess_step)

  pulse_in_coor = util.insert_path(in_parms['input_crds'], '..')
  pulse_in_vel = util.insert_path(in_parms['input_vels'], '..')

  stopwatch = util.Timer()

  pulses = ["pulse%d" % i for i in range(n_pulse)]
  for pulse, n_step in zip(pulses, n_step_list):
    util.goto_dir(pulse)

    pulse_parms = copy.deepcopy(in_parms)
    pulse_parms['topology'] = util.insert_path(in_parms['topology'], '..') 
    pulse_parms['n_step_dynamics'] = n_step

    if 'constraint_pdb' in pulse_parms and pulse_parms['constraint_pdb']: 
      pulse_parms['constraint_pdb'] = \
          util.insert_path(in_parms['constraint_pdb'], '..')

    pulse_parms['input_crds'] = pulse_in_coor
    pulse_parms['input_vels'] = pulse_in_vel
    soup = SoupFromRestartFiles(pulse_parms['topology'], 
                              pulse_parms['input_crds'],
                              pulse_parms['input_vels'])
    reset_vel_func(soup)

    pulse_parms['input_vels'] = name + ".in.vel"
    write_soup_to_vel(soup, pulse_parms['input_vels'])

    run(pulse_parms)

    pulse_in_coor = '../%s/%s.coor' % (pulse, name)
    pulse_in_vel = '../%s/%s.vel' % (pulse, name)
    os.chdir('..')
    
  stopwatch.stop()
  open(name + '.time', 'w').write(stopwatch.str())

  merge_simulations(name, pulses)
  for pulse in pulses: os.system('rm -rf %s' % pulse)

  util.write_dict(in_parms, config)


def merge_simulations(name, pulses):
  for ext in ['.psf', '.coor', '.vel']:
    fname = '%s%s' % (name, ext)
    shutil.copy('%s/%s' % (pulses[-1], fname), fname)
  trajs = [os.path.join(pulse, name + '.dcd') for pulse in pulses]
  merge_trajectories(name + '.psf', trajs, name + '.dcd')
  vels = [os.path.join(pulse, name + '.vel.dcd') for pulse in pulses]
  merge_trajectories(name + '.psf', vels, name + '.vel.dcd')
  

def check_dcd_byte_order(dcd):
  if sys.byteorder in 'big':
    option = '-B'
  elif sys.byteorder in 'little':
    option = '-L'
  else:
    raise "Couldn't find o.s. byte order %s" % sys.byteorder
  util.run_with_output('flipdcd %s %s' % (option, dcd))


# Routines to handle the CHARMM dcd trajectories
class DcdReader:
  """
  Read frames from a DCD file in terms of a list (3xn) of floats.

  Data: 
    fname, n_fixed_atom, remarks, 
    free_atom_indices, posFirstFrame, 
    size_frame, n_frame, n_atom

  Methods:
    d = DCD(fname) - open fname and read header
    len(d) - return number of frames
    d[10] - return a frame as a list (3xn) of floats.
  """

  def __init__(self, fname):
    self.fname = fname

    check_dcd_byte_order(self.fname)    
    
    self._file = open(fname, 'rb')

    # Read header
    self._read_fmt_vals('i4c')
    # if self._read_fmt_vals('i4c') <> (84, 'C', 'O', 'R', 'D') :
    #  raise "DCD format error 1"

    self._read_fmt_val('i') # nSet

    self.iStart = self._read_fmt_val('i')
    self.nStepSave = self._read_fmt_val('i')

    # skip some 
    self._read_fmt_val('5i')
    self.n_fixed_atom = self._read_fmt_val('i')
    self.timeStep = self._read_fmt_val('d')
    self._read_fmt_vals('9i')
    if self._read_fmt_val('i') <> 84 :
      raise "DCD format error 2"
    size = self._read_fmt_val('i')
    if (size - 4) % 80 <> 0 :
      raise "DCD format error 3"

    self.remarks = []
    n_line = self._read_fmt_val('i')
    for i in range(0, n_line):
      s = "".join(self._read_fmt_vals('80c'))
      self.remarks.append(s.strip())

    if self._read_fmt_val('i') <> size :
      raise "DCD format error 4"
    if self._read_fmt_val('i') <> 4 :
      raise "DCD format error 5"
    self.n_atom = self._read_fmt_val('i')
    if self._read_fmt_val('i') <> 4 :
      raise "DCD format error 6"

    if self.n_fixed_atom > 0:
      fmt = "%di" % (self.n_atom - self.n_fixed_atom)
      size = struct.calcsize(fmt)
      if self._read_fmt_val('i') <> size:
        raise "DCD format error 7"
      self.free_atom_indices = self._read_fmt_val(fmt)
      for i in range(len(self.free_atom_indices)):
        self.free_atom_indices[i] -= 1
      if self._read_fmt_val('i') <> size :
        raise "DCD format error 8"
    else:
      self.free_atom_indices = ()

    self.pos_start_frame = self._file.tell()
    size_first_frame = struct.calcsize('%df6i' % (3 * self.n_atom))
    self.posFirstFrame = self._file.tell() + size_first_frame

    n_free_atom = self.n_atom - self.n_fixed_atom
    self.size_frame = struct.calcsize('%df6i' % (3 * n_free_atom))

    if self.n_fixed_atom > 0:
      self._first_frame_x_vals = [0.0 for i in range(self.n_atom)]
      self._first_frame_y_vals = [0.0 for i in range(self.n_atom)]
      self._first_frame_z_vals = [0.0 for i in range(self.n_atom)]
      crds_fmt = '%df' % self.n_atom
      size = struct.calcsize(crds_fmt)

      if size <> self._read_fmt_val('i'):
        raise "DCD format error 9"
      self._first_frame_x_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 10"

      if size <> self._read_fmt_val('i'):
        raise "DCD format error 11"
      self._first_frame_y_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 12"

      if size <> self._read_fmt_val('i'):
        raise "DCD format error 13"
      self._first_frame_z_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 14"

    self._file.seek(0, 2)
    self.n_frame = int((self._file.tell()-self.posFirstFrame)/self.size_frame) + 1

  def _read_fmt_vals(self, fmt):
    return struct.unpack(fmt, self._file.read(struct.calcsize(fmt)))

  def _read_fmt_val(self, fmt):
    return self._read_fmt_vals(fmt)[0]
    
  def __getitem__(self, i):
    if i < - 1*self.n_frame or i >= self.n_frame :
      raise IndexError
    if i < 0 :
      return self.__getitem__(self.n_frame + i)

    if self.n_fixed_atom == 0 :
      self._file.seek(self.pos_start_frame + i*self.size_frame)

      crds_fmt = "%df" % self.n_atom
      size = struct.calcsize(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 9"
      x_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 10"
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 11"
      y_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 12"
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 13"
      z_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 14"

      frame = []
      for x, y, z in zip(x_vals, y_vals, z_vals):
        frame.extend((x, y, z))

    else:
      # first frame cached
      if i == 0:
        return copy.copy(self._firstFrame)

      self._file.seek(self.pos_start_frame + (i-1)*self.size_frame)

      crds_fmt = '%df' % len(self.free_atom_indices)
      size = struct.calcsize(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 9"
      free_x_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 10"
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 11"
      free_y_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 12"
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 13"
      free_z_vals = self._read_fmt_vals(crds_fmt)
      if size <> self._read_fmt_val('i'):
        raise "DCD format error 14"

      x_vals = copy.copy(self._first_frame_x_vals)
      y_vals = copy.copy(self._first_frame_y_vals)
      z_vals = copy.copy(self._first_frame_z_vals)
      k = 0
      for index in self.free_atom_indices:
        x_vals[index] = free_x_vals[k]
        y_vals[index] = free_y_vals[k]
        z_vals[index] = free_z_vals[k]
        k += 1

      frame = []
      for x, y, z in zip(x_vals, y_vals, z_vals):
        frame.extend((x, y, z))

    return frame

  def __len__(self):
    return self.n_frame

  def __del__(self):
    self._file.close()

  def __repr__(self):
    return "< DCD %s with %d frames of %d atoms (%d fixed) >" % \
             (self.fname, self.n_frame, self.n_atom, self.n_fixed_atom)


class Trajectory:
  def __init__(self, name):
    self.topology = name + '.psf'
    coor_fname = name + '.dcd'    
    self.coor_traj = DcdReader(coor_fname)
    vel_fname = name + '.vel.dcd'
    if os.path.isfile(vel_fname):
      self.vel_traj = DcdReader(vel_fname)
    else:
      self.vel_traj = None
    self.soup = SoupFromNamdPsf(self.topology, coor_fname, vel_fname)
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


def merge_trajectories(psf, dcds, out_dcd):
  dcd_reader = DcdReader(dcds[0])
  pos_start_frame = dcd_reader.pos_start_frame  
  size_frame = dcd_reader.size_frame
  del dcd_reader

  shutil.copy(dcds[0], out_dcd)

  merge_dcd_file = open(out_dcd, "ab+")
  for dcd in dcds[1:]:
    dcd_file = open(dcd, "rb")
    dcd_file.seek(-1, 2)
    eof = dcd_file.tell()

    dcd_file.seek(pos_start_frame)
    while dcd_file.tell() < eof:
      merge_dcd_file.write(dcd_file.read(size_frame)) 

    dcd_file.close()
  merge_dcd_file.close()


def preequilibrate(in_name, out_name, temperature):
  top, coor, vel = get_restart_files(in_name)
  restraint_name = 'restraint'
  parms = langevin_thermometer_parms.copy()
  parms['restraint'] = True
  parms['topology'] = top
  parms['input_crds'] = coor
  parms['input_vels'] = vel
  parms['output_name'] = restraint_name
  parms['temp_thermometer'] = temperature
  parms['temp_init'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)  

  top, coor, vel = get_restart_files(restraint_name)
  parms = langevin_thermometer_parms.copy()
  parms['topology'] = top
  parms['input_crds'] = coor
  parms['input_vels'] = vel
  parms['output_name'] = out_name
  parms['temp_thermometer'] = temperature
  parms['n_step_per_snapshot'] = 50
  parms['n_step_dynamics'] = 1000
  run(parms)  

