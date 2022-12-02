#!/usr/bin/env python
import sys, os, getopt, util, analysis

pymol = "/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL"

usage = """
Shows an amber molecular dynamics trajectory with pymol using various
display options. Defaults to cartoon display with stick sidechains and
hydrogen atoms removed.

Usage: pytraj.py [-f] [-d] md [res]

-f:   show all atoms, including hydrogens
-d:   show dot density of the atoms
-e:   show in peptide mode
md:   name of trajectory, assumes md.top and md.trj
res:  highlight this residue in sequential numbering, starting from 1
"""

opts, args = getopt.getopt(sys.argv[1:], "dfe")
if len(args) == 0:
  print(usage)
  sys.exit(1)

def is_opt(opt, opts):
  return opt in [a for a,b in opts]
  
script = ""

md = os.path.abspath(args[0])
path, name = os.path.split(md)
if not os.path.isfile(md + '.top') or \
   not os.path.isfile(md + '.trj'):
  print("Can't find amber trajectory files.")
  sys.exit(1)
script += "load %s.top\n" % md
script += "load %s.trj\n" % md

script += "hide\n"

script += "select backbone, name c+o+n+h+oxt\n"
script += "deselect\n"

if is_opt("-f", opts):
  script += "show sticks\n"
else:
  script += "show sticks\n"
  script += "hide sticks, backbone\n"
  script += "hide sticks, hydro\n"
  script += "set cartoon_flat_sheets, 0\n"
  script += "set cartoon_loop_radius, 0.4\n"
  script += "show cartoon\n"

if is_opt("-d", opts):
  script += "show dots\n"
  script += "set dot_density, 2\n"
  if not is_opt("-f", opts):
    script += "hide dots, backbone\n"

if is_opt("-e", opts):
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

# choose grey for carbon atoms
script += 'util.cba(29,"%s")\n' % name
# script += "cmd.bg_color('white')\n"

# calculate secondary structure for all frames in trajectory
n_frame = analysis.open_trajectory(md).n_frame
ss_lines = ['dss all, %s\n' % (i+1) for i in range(n_frame)]
script += ''.join(ss_lines)

# color selected residue
if len(args) > 1:
  script += "select x, resi %s\n" % args[1]
  script += "util.cba(26, 'x')\n" # 154 = pink carbons
  script += "center x\n"
  script += "hide sticks, hydro or not x\n"
  script += "select neighbors, (byres center around 5)\n"
  script += "show sticks, neighbors and not backbone and not hydro\n"
  script += "deselect\n"

temp = util.temp_fname('.pml')
open(temp, "w").write(script)

os.system(pymol + " " + temp + "  > /dev/null")

if os.path.isfile(temp):
  os.remove(temp)

