#!/usr/bin/env python

pymol = "/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL"

import copy
import os
import util
import pdbstruct
import vector3d
import asa


def get_res_name(res):
  return '%s:%d' % (res.chain_id, res.num)


def split_resname(resname):
  "Returns (chain_id, res_num)"
  words = resname.split(":")
  if len(words) == 2:
    return (words[0], int(words[1]))
  else:
    return (' ', int(words[0]))

  
def find_ca_of_resname(atoms, resname):
  chain_id, res_num = split_resname(resname)
  for atom in atoms:
    if chain_id == atom.chain_id and res_num == atom.res_num:
      if "CA" == atom.type:
        return atom
  raise IndexError, "Can't find atom %s" % resname


def pymol_id_from_resname(resname):
  chain_id, res_num = split_resname(resname)
  if chain_id == " ":
    return "resi %d" % res_num
  else:
    return "(chain %s and resi %d)" % (chain_id, res_num)
  
    
def get_pdb_transform(pdb, center_res, top_res):
  """
  Returns a transformation matrix that centers pdb to 
  center_res on the z-axis and moves top_res above center_res
  on the y-axis
  """
  soup = pdbstruct.Soup(pdb)
  atoms = soup.atoms()
  soup_center = pdbstruct.get_center(atoms)
  translation = vector3d.Translation(-soup_center)
  soup.transform(translation)
  result = translation

  center_atom = find_ca_of_resname(soup.atoms(), center_res)
  view = vector3d.Vector3d(0, 0, 1)
  axis = vector3d.CrossProductVec(view, center_atom.pos)
  angle = vector3d.vec_dihedral(view, axis, center_atom.pos)
  rotation = vector3d.RotationAtOrigin(axis, angle)
  soup.transform(rotation)
  result = rotation * translation

  top_atom = find_ca_of_resname(soup.atoms(), top_res)
  top_dir = vector3d.Vector3d(0, 1, 0)
  axis = view.copy()
  angle = vector3d.vec_dihedral(top_dir, axis, top_atom.pos)
  rotation2 = vector3d.RotationAtOrigin(axis, angle)
  result = rotation2*result
  
  del soup
  return result


def rescale_soup_bfactors(
    in_soup, temp_pdb, lower_bfactor, upper_bfactor):
  """
  Returns max_bfactor after rescale (needed for worm
  calculation)
  """
  soup = in_soup.copy()
  bfactors = []
  bfactors = [a.bfactor for a in soup.atoms()]
  # cut-off max_values
  if upper_bfactor:
    for j in range(len(bfactors)):
      if bfactors[j] > upper_bfactor:
        bfactors[j] = upper_bfactor
    # add fake water atom with upper bound just in case no 
    # atom with upper bound exists, will delete later
    dummy_atom = pdbstruct.Atom()
    dummy_atom.pos = soup.atoms()[0].pos.copy()
    dummy_atom.type = "O"
    dummy_atom.bfactor = upper_bfactor
    dummy_res = pdbstruct.Residue('XXX', '', 9999)
    dummy_res.insert_atom(dummy_atom)
    dummy_chain = pdbstruct.Polymer()
    dummy_chain.append_residue_no_renum(dummy_res)
    soup.append_chain(dummy_chain)
  # cut-off below min_val to zero
  if lower_bfactor:
    for j in range(len(bfactors)):
      if bfactors[j] < lower_bfactor:
        bfactors[j] = 0.0
  for a, bfactor in zip(soup.atoms(), bfactors):
    a.bfactor = bfactor
  soup.write_pdb(temp_pdb)
  return max(bfactors)


def make_new_bfactor_pdbs(pdbs, lower_bfactor, upper_bfactor):
  "Returns list of new_pdbs, and max_bfactor"
  max_bfactor = 0
  new_pdbs = []
  for pdb in pdbs:
    soup = pdbstruct.Soup(pdb)
    new_pdb = util.fname_variant(pdb)
    this_max_bfactor = rescale_soup_bfactors(
        soup, new_pdb, lower_bfactor, upper_bfactor)
    if this_max_bfactor > max_bfactor:
      max_bfactor = this_max_bfactor
    new_pdbs.append(new_pdb)
    del soup
  return new_pdbs, max_bfactor


def rescale_polar_soup_bfactors(
    in_soup, temp_pdb, lower_bfactor, upper_bfactor):
  """
  Returns max_bfactor after rescale (needed for worm
  calculation)
  """
  soup = copy.deepcopy(in_soup)
  bfactors = []
  bfactors = [a.bfactor for a in soup.atoms()]
  if upper_bfactor is None:
    upper_bfactor = max(bfactors)
  # cut-off max_values
  if upper_bfactor:
    for j in range(len(bfactors)):
      if bfactors[j] > upper_bfactor:
        bfactors[j] = upper_bfactor
      if bfactors[j] < -upper_bfactor:
        bfactors[j] = -upper_bfactor
    # add fake water atom with upper bound just in case no 
    # atom with upper bound exists, will delete later
    dummy_chain = pdbstruct.Polymer()
    dummy_atom = pdbstruct.Atom()
    dummy_atom.pos = soup.atoms()[0].pos.copy()
    dummy_atom.type = "O"
    dummy_atom.bfactor = upper_bfactor
    dummy_res = pdbstruct.Residue('XXX', '', 9999)
    dummy_res.insert_atom(dummy_atom)
    dummy_chain.append_residue_no_renum(dummy_res)
    dummy_atom = pdbstruct.Atom()
    dummy_atom.pos = soup.atoms()[0].pos.copy()
    dummy_atom.type = "O"
    dummy_atom.bfactor = -upper_bfactor
    dummy_res = pdbstruct.Residue('XXX', '', 9998)
    dummy_res.insert_atom(dummy_atom)
    dummy_chain.append_residue_no_renum(dummy_res)
    soup.append_chain(dummy_chain)
  # cut-off below min_val to zero
  if lower_bfactor:
    for j in range(len(bfactors)):
      if -lower_bfactor < bfactors[j] < lower_bfactor:
        bfactors[j] = 0.0
  for a, bfactor in zip(soup.atoms(), bfactors):
    a.bfactor = bfactor
  soup.write_pdb(temp_pdb)
  return max(bfactors)


def make_new_polar_bfactor_pdbs(pdbs, lower_bfactor, upper_bfactor):
  "Returns list of new_pdbs, and max_bfactor"
  max_bfactor = 0
  new_pdbs = []
  for pdb in pdbs:
    soup = pdbstruct.Soup(pdb)
    new_pdb = util.fname_variant(pdb)
    this_max_bfactor = rescale_polar_soup_bfactors(
        soup, new_pdb, lower_bfactor, upper_bfactor)
    if this_max_bfactor > max_bfactor:
      max_bfactor = this_max_bfactor
    new_pdbs.append(new_pdb)
    del soup
  return new_pdbs, max_bfactor


def bgcolor_script(bg_color='black'):
  return  "cmd.bg_color('%s')\n" % bg_color

  
def load_pdbs_script(pdbs, bg_color='black'):
  "Returns pymol script, name of pdbs"
  script = ""
  for pdb in pdbs:
    script += "load %s\n" % pdb
  script += "hide everything\n"
  return script
  

def color_chain_script(pdbs):
  names = [os.path.basename(p).replace('.pdb', '') 
           for p in pdbs]
  colors = ['util.color_chains("(%s and elem c)")\n' % n 
            for n in names]
  return ''.join(colors)


def color_bfactor_script():
  module_dir = os.path.dirname(__file__)
  color_py = os.path.join(module_dir, "color_b.py")
  script = ""
  script += "run " + color_py + "\n"
  script += "color_b all, gradient=wr\n"
  script += "remove resn XXX\n"
  return script

def color_polar_bfactor_script():
  module_dir = os.path.dirname(__file__)
  color_py = os.path.join(module_dir, "color_b.py")
  script = ""
  script += """cmd.spectrum("b", 'blue_white_red', selection="all");\n"""
  script += "remove resn XXX\n"
  return script


def cartoon_script():
  script = ""
  script += "set cartoon_flat_sheets, 0\n"
  script += "set cartoon_loop_radius, 0.4\n"
  script += "set cartoon_tube_radius, 0.4\n"
  script += "cartoon auto\n"
  script += "show cartoon\n"
  return script


def putty_script(scale_max=4.0):
  script = ""
  script += "set cartoon_flat_sheets, 0\n"
  script += "set cartoon_putty_scale_max, %f\n" % scale_max
  script += "set cartoon_putty_radius, 0.4\n"
  script += "set cartoon_putty_scale_power, 1\n"
  script += "cartoon putty\n"
  script += "show cartoon\n"
  return script


def sticks_above_bfactor_script(lower_bfactor):
  script = ""
  script += "select hot, b > %f or b < -%f\n" % \
             (lower_bfactor, lower_bfactor)
  script += "select cold, b < %f and b > -%f\n" % \
             (lower_bfactor, lower_bfactor)
  script += "show sticks, hot\n"
  script += "hide sticks, cold\n"
  script += "deselect\n"
  return script
    
    
def ligands_as_sticks_script(pdbs, color=""):
  script = ""
  for pdb in pdbs:
    name = os.path.basename(pdb).replace('.pdb', '')
    soup = pdbstruct.Soup(pdb)
    for res in soup.residues():
      if res.type not in pdbstruct.res_name_to_char:
        if res.type not in "HOH":
          chain_id_script = ""
          if res.chain_id.strip():
            chain_id_script = "and chain %s" % res.chain_id
          script += \
              "show sticks, %s %s and resn %s and resi %d\n" \
                % (name, chain_id_script, res.type, res.num)
          if color:
            script += \
              "color %s, %s %s and resn %s and resi %d\n" \
                % (color, name, chain_id_script, res.type, res.num)
  script += "show nonbonded\n"
  return script


def highlight_res_script(highlight_res):
  script = ""
  script += "select highlight, %s\n" % \
              pymol_id_from_resname(highlight_res)
  script += "show sticks, highlight\n"
  script += "color green, highlight\n"
  script += "deselect\n"
  return script


def hide_backbone_sticks_script():
  script = ""
  script += "hide sticks, hydro\n"
  script += "select bb, name c+o+n+h+oxt\n"
  script += "hide sticks, bb\n"
  script += "select nuc, resn A+U+T+C+G+A3+U3+T3+C3+G3+A5+U5+T5+C5+G5+DA+DT+DC+DG\n"
  script += "select nuc_bb, name P+O1P+O2P+OP1+Op2+O3'+C3'+C2'+C1'+O4'+C4'+C5'+O5'\n"
  script += "hide cartoon, nuc and not nuc_bb\n"
  script += "deselect\n"
  return script


def run_pymol_script(script, width=500, height=500):
  temp_pml = util.temp_fname('.pml')
  f = open(temp_pml, "w")
  f.write(script)
  f.close()
  cmd = pymol + " -W %d -H %d" % (width, height)
  is_quit = False
  for l in script.splitlines():
    if l.startswith('quit'):
      is_quit = True
      break
  if is_quit:
    cmd += ' -iqx '
  else:
    cmd += ' -q '
  cmd += temp_pml
  util.run_with_output(cmd)
  os.remove(temp_pml)


def transformed_soup_from_pdb(
    pdb, center_res=None, top_res=None, 
    width=None, height=None, frame_residues=None):
  soup = pdbstruct.Soup(pdb)
  if center_res and top_res:
    transform = get_pdb_transform(pdb, center_res, top_res)
    soup.transform(transform)
  if frame_residues:
    resnames = [pymol_id_from_resname(r) for r in frame_residues]
    soup.frame_pymol_script = "zoom (%s)\n" % ' or '.join(resnames)
  if width: soup.width = width
  if height: soup.height = height
  return soup
  

def get_scale_max(max_bfactor, upper_bfactor):
  scale_max = 4.0
  if max_bfactor is not None and upper_bfactor is not None:
    if max_bfactor < upper_bfactor:
      scale_max *= max_bfactor / upper_bfactor
      if scale_max < 1:
        scale_max = 1
  return scale_max

  
def bfactor_script(pdb, lower_bfactor=None, upper_bfactor=None, 
                   max_bfactor=None, is_putty=False):
  "Returns script that displays bfactors of pdb" 
  script = load_pdbs_script([pdb])
  script += color_bfactor_script()
  if is_putty:
    script += putty_script(get_scale_max(
         max_bfactor, upper_bfactor))
  else:
    script += cartoon_script()
    script += "cartoon tube\n"
  if not is_putty:
    if lower_bfactor is not None:
      script += sticks_above_bfactor_script(lower_bfactor)
    else:
      script += "show sticks\n"
  script += ligands_as_sticks_script([pdb])
  script += hide_backbone_sticks_script()
  return script


def filterd_file(fname, filter_fn):
  lines = open(fname, 'r').readlines()
  new_lines = filter(filter_fn, lines)
  new_fname = util.fname_variant(fname)
  open(new_fname, 'w').write(''.join(new_lines))
  return new_fname
  
  
def soup_to_bfactor_png(
    soup, png, bfactors, lower_bfactor=None, upper_bfactor=None,
    highlight_res=None, is_putty=False):
  temp_pdb = util.temp_fname('.pdb')
  soup.load_residue_bfactors(bfactors)
  max_bfactor = rescale_soup_bfactors(soup, temp_pdb, 
                lower_bfactor, upper_bfactor)
  # delete solvent waters
  temp2_pdb = filterd_file(
      temp_pdb, lambda l: l[17:20] not in ['SOL', 'TIP', 'NA+'])
  script = bgcolor_script('white')
  script += bfactor_script(
       temp2_pdb, lower_bfactor, upper_bfactor, max_bfactor, is_putty)
  if highlight_res is not None:
    script += highlight_res_script(highlight_res)
    script += hide_backbone_sticks_script()
  if 'frame_pymol_script' in soup.__dict__:
    script += soup.frame_pymol_script
  script += "clip far, -20\n"
  script += "save %s\n" % png
  script += "quit"
  width, height = 480, 480
  if 'width' in soup.__dict__:
    width = soup.width
  if 'height' in soup.__dict__:
    height = soup.height
  run_pymol_script(script, width, height)
  os.remove(temp_pdb)
  os.remove(temp2_pdb)


def make_pdb_png(
    png, pdbs, bgcolor="white", center_res=None, top_res=None,
    highlight_res=None, is_sticks=True, is_putty=False,
    width=480, height=480):
  temp_pdbs = []
  if center_res and top_res:
    transform = get_pdb_transform(pdbs[0], center_res, top_res)
    for i in range(len(pdbs)):
      soup = pdbstruct.Soup(pdbs[i])
      soup.transform(transform)
      new_pdb = util.fname_variant(pdbs[i])
      soup.write_pdb(new_pdb)
      temp_pdbs.append(new_pdb)
      pdbs[i] = new_pdb
      del soup
  if 'transparent' in bgcolor:
    script = 'set opaque_background, off\n'
  else: 
    script = bgcolor_script(bgcolor)
  script += load_pdbs_script(pdbs)
  script += color_chain_script(pdbs)
  if is_putty:
    script += putty_script(get_scale_max(
        max_bfactor, upper_bfactor))
  else:
    script += cartoon_script()
  if not is_sticks:
    script += "hide sticks\n"
  else:
    script += "show sticks\n"
  script += ligands_as_sticks_script(pdbs)
  if highlight_res:
    script += highlight_res_script(highlight_res)
  script += hide_backbone_sticks_script()
  # script += "clip far, 5\n"
  script += "save %s\n" % png
  script += "quit"

  run_pymol_script(script, width, height)

  for pdb in temp_pdbs:
    os.remove(pdb)


def make_pair_pdb_png(soup, i, j, png):
  bfactors = [0.0 for k in range(soup.n_residue())]
  bfactors[i] = 1.0
  bfactors[j] = 1.0
  soup_to_bfactor_png(soup, png, bfactors, 0.5, 1.0)


def convert_soup_to_ala(soup):
  for chain in soup.chains():
    if isinstance(chain, pdbstruct.Protein):
      for i in range(chain.n_residue()):
        chain.mutate(i, "ALA")
  

usage = """

  Copyright (c) 2007 Bosco Ho

  Starts pymol using various complex options. Uses the
  color_b.py module from Robert L. Campbell to do two-tone
  coloring.

  Usage: show.py [options] pdb ...
  
  options: -pxsbq -o out_file -c res -t res -h res -l lower -u upper

  -s:   turn-off sticks 
  -e:   peptide display mode (sphere for bb, sticks for sc)
  -x:   use blue-white-red coloring
  -g:   background color
  -h:   highlight residue with sticks and green
  -b:   color with b-factors
  -p:   putty mode with b-factors
  -l:   lower bound to display residues with b-factors
  -u:   upper bound saturation for b-factor so 
        colors can be scaled properly
  -f:   show all sticks
  -o:   output_file
  -q:   quit after rendering
  pdb:  name of pdb files, multiple pdb files can be entered.
"""


if __name__=="__main__":

  import sys, os, getopt, util, shutil
  

  opts, args = getopt.getopt(sys.argv[1:], "sqbpxewo:g:b:c:t:u:l:h:")

  if len(args) == 0:
    print usage
    sys.exit(1)

  lower_bfactor = None
  upper_bfactor = None
  highlight_res = None
  center_res = None
  top_res = None
  is_bfactor = False
  is_putty = False
  is_sticks = True
  is_peptide = False
  is_polar = False
  background_color = 'black'
  out_file = None
  is_quit = False 
  is_kill_water = False
  
  for opt, a in opts:
    if '-s' in opt:
      is_sticks = False
    if '-g' in opt:
      background_color = a
    if '-e' in opt:
      is_peptide = True
    if '-p' in opt:
      is_bfactor = True
      is_putty = True
      is_sticks = False
    if '-l' in opt:
      lower_bfactor = float(a) 
    if '-u' in opt:
      upper_bfactor = float(a) 
    if '-b' in opt:
      is_bfactor = True
    if '-x' in opt:
      is_polar = True
    if '-h' in opt:
      highlight_res = a
    if '-c' in opt:
      center_res = a
    if '-t' in opt:
      top_res = a
    if '-o' in opt:
      out_file = a
    if '-q' in opt:
      is_quit = True
    if '-w' in opt:
      is_kill_water = True
    
  if lower_bfactor is not None or upper_bfactor is not None:
    if not is_polar and not is_bfactor:
      is_bfactor = True
      
  temp_pdbs = []

  # check if pdb files actually exist
  pdbs = []
  for arg in args:
    pdb = os.path.abspath(arg)
    if os.path.isfile(pdb):
      if not pdb.endswith('.pdb'):
        renamed_pdb = pdb + '.pdb'
        shutil.copy(pdb, renamed_pdb)
        pdb = renamed_pdb
        temp_pdbs.append(renamed_pdb)
      pdbs.append(pdb)
    else:
      pdb = pdb + '.pdb'
      if os.path.isfile(pdb):
        pdbs.append(pdb)
      else:
        print("Can't find %s" % arg)
        sys.exit(1)

  if center_res and top_res:
    transform = get_pdb_transform(pdbs[0], center_res, top_res)
    for i in range(len(pdbs)):
      soup = pdbstruct.Soup(pdbs[i])
      soup.transform(transform)
      new_pdb = util.fname_variant(pdbs[i])
      soup.write_pdb(new_pdb)
      temp_pdbs.append(new_pdb)
      pdbs[i] = new_pdb
      del soup

  if is_kill_water:
    for i in range(len(pdbs)):
      lines = open(pdbs[i], 'r').readlines()
      new_lines = []
      for l in lines:
        if l[17:20] not in ['TIP', 'SOL']:
          new_lines.append(l)
      new_pdb = util.fname_variant(pdbs[i])
      open(new_pdb, 'w').write(''.join(new_lines))
      pdbs[i] = new_pdb
      temp_pdbs.append(new_pdb)
  
  max_bfactor = None
  if is_bfactor:
    rescaled_pdbs, max_bfactor = make_new_bfactor_pdbs(
        pdbs, lower_bfactor, upper_bfactor)
    temp_pdbs.extend(rescaled_pdbs)
    pdbs = rescaled_pdbs
  elif is_polar:
    rescaled_pdbs, max_bfactor = make_new_polar_bfactor_pdbs(
        pdbs, lower_bfactor, upper_bfactor)
    temp_pdbs.extend(rescaled_pdbs)
    pdbs = rescaled_pdbs
    
    
  script = bgcolor_script(background_color)
  script += load_pdbs_script(pdbs)
  if is_bfactor:
    script += color_bfactor_script()
  elif is_polar:
    script += color_polar_bfactor_script()
  else:
    script += color_chain_script(pdbs)
  if is_putty:
    script += putty_script(get_scale_max(
        max_bfactor, upper_bfactor))
  else:
    script += cartoon_script()

  if lower_bfactor is not None:
    script += sticks_above_bfactor_script(lower_bfactor)
  else:
    script += "show sticks\n"
  script += ligands_as_sticks_script(pdbs)
  if highlight_res is not None:
    script += highlight_res_script(highlight_res)
    
  if is_peptide:
    script += """
      select bb, name ca+n+h+o+c+oxt+h1+h2+h3+ch3+hh31+hh32+hh33+3hh3+2hh3+1hh3
      show sphere, bb
      hide stick, bb
      util.cbaw bb
      select sc, not bb and not hydro
      show stick, sc
      hide sphere, sc
      util.cbag sc
      set sphere_quality, 2
    """
  script += hide_backbone_sticks_script()

  if not is_sticks:
    script += "hide sticks\n"

  if out_file:
    script += "save %s\n" % out_file

  if is_quit:
    script += "quit\n"

  run_pymol_script(script, 800, 800)
  
  for pdb in temp_pdbs:
    os.remove(pdb)


