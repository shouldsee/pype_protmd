#!/usr/bin/env python

import os, sys, util

vmd = r'/Applications/VMD\ 1.8.7.app/Contents/Resources/VMD.app/Contents/MacOS/VMD'

def has_files(*fnames):
  for fname in fnames:
    if not os.path.isfile(fname):
      return False
  return True


s = """
# turn on lights 0 and 1
light 0 on
light 1 on
light 2 off
light 3 off

# position the stage and axes
axes location off
stage location off

# position and turn on menus
menu main on
menu graphics on
menu main move 900 60
menu graphics move 900 340

display projection orthographic
display depthcue on
display cuestart 1.75
display cueend 3.00
display cuemode Linear

mol new "$top" type $ext_top
mol delrep 0 0
mol addfile "$trj" type $ext_trj
mol selection {backbone}
mol representation NewCartoon 0.40 20.00 2.50 0
mol addrep top
mol selection {(not (protein and backbone or water) or name CA) and not name "[0-9]?H.*"}
mol representation CPK 1.4 2 8 6
mol addrep top
"""


if len(sys.argv) < 2:
  print "Usage: vmdtraj.py md"
  sys.exit(1)
    
in_md = sys.argv[1]

parms = {
    'ext_top': 'parm7',
    'ext_trj': 'crd',
    'top': in_md + '.top', 
    'trj': in_md + '.trj' ,
    }
if not has_files(parms['top'], parms['trj']):
  parms = {
      'ext_top': 'gro',
      'ext_trj': 'trr',
      'top': in_md + '.gro', 
      'trj': in_md + '.trr' ,
      }
  if not has_files(parms['top'], parms['trj']):
    print "Can't recognize trajectories with basename", in_md

cmd_fname = util.temp_fname('.vmd')
script = util.replace_dict(s, parms)
open(cmd_fname, 'w').write(script)
os.system(vmd + ' -e ' + cmd_fname)
os.remove(cmd_fname)
