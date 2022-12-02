#!/usr/bin/env python

import os
import glob
import shutil

import util
import analysis
import pdbstruct


def find_pdb_file(search_dir):
  os.chdir(search_dir)
  pdbs = util.re_glob('*pdb', r'^\d\w\w\w.*[.]pdb$')
  pdbs = [p for p in pdbs if 'hollow' not in p]
  if len(pdbs):
    return pdbs[0]
  else:
    return None
  

def make_rmsd_data_files(
    sim_dir, out_dir, i_ps, min_val, max_val):
  save_dir = os.getcwd()

  print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
  print(" ")
  print(" Analyzing RIP perturbations in:")
  print("  ", sim_dir)
  print(" Averaged values in the", str(i_ps)+"th ps")
  print(" Results stored in:")
  print("  ", out_dir)
  print()
  print()

  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)
    
  ref_pdb = out_dir+'/sim.pdb'
  shutil.copy(sim_dir + '/../sim.pdb', ref_pdb)
  raw_pdb = find_pdb_file('%s/..' % sim_dir)
  if raw_pdb:
    shutil.copy(raw_pdb, out_dir+'/raw.pdb')

  # calculate ca_rmsd.ave in each pertrubation directory
  os.chdir(sim_dir)
  res_dirs = [d for d in glob.glob('*') if d.isdigit()]
  for res_dir in res_dirs:
    os.chdir(res_dir)
    if os.path.isdir('pulse0'):
      print("Incomplete heating in", os.getcwd())
    elif not os.path.isfile('ca_rmsd') or \
         not os.path.isfile('ca_rmsd.ave'):
      print("Extracting Ca-RMSD from trajectory of RIP on residue", end=' ')
      print(res_dir)
      analysis.process_md('md')
    os.chdir(sim_dir)
  
  # extract the i_ps averages and make the matrix of fluctuations
  res_data = util.get_numbered_dir_data('ca_rmsd.ave', i_ps)
  
  # calculate the background limits
  if max_val is None:
    max_val = max(res_data.flatten())
  if min_val is None:
    flat = res_data.flatten()
    min_val = flat.mean() + 2.0*flat.std()

  ref_soup = pdbstruct.Soup(ref_pdb)

  util.goto_dir(out_dir)

  # write the ca_rmsd map for the given i_ps time period
  util.write_array('map', res_data)

  # sum the amount of ca_rmsd for each perturbation
  y_vals = util.count_n_vals_per_column(res_data, min_val)
  util.write_column('strength', y_vals)
  ref_soup.load_residue_bfactors(y_vals)
  ref_soup.write_pdb('strength.pdb')

  # sum the amount of ca_rmsd susceptibility for flexibility
  y_vals = util.count_n_vals_per_column(
      res_data.transpose(), min_val)
  util.write_column('flexibility', y_vals)
  ref_soup.load_residue_bfactors(y_vals)
  ref_soup.write_pdb('flexibility.pdb')

  # for residues that gets knocked around, calculate average
  # amount of deviation
  n_residue = res_data.shape[0]
  y_vals = []
  for j in range(n_residue):
    fluctuations = []
    for i in range(n_residue):
      v = res_data[i,j]
      if abs(i-j) > 3 and v > min_val:
        fluctuations.append(v-min_val)
    if fluctuations:
      mean = sum(fluctuations)/float(len(fluctuations))
    else:
      mean = 0.0
    y_vals.append(mean)
  util.write_column('fluctuation', y_vals)
  ref_soup.load_residue_bfactors(y_vals)
  ref_soup.write_pdb('fluctuation.pdb')
  
  # calculate multiplier to produce a
  # readable bfactor in pdb files
  median = 1.0
  multiplier = 1.0
  while median > max_val:
    multiplier *= 10.0
    median = 1.0 / multiplier
  
  # make pdb files showing pathways
  os.chdir(sim_dir)
  for res_dir in res_dirs:
    if os.path.isdir('%s/pulse0' % d):
      print("Incomplete heating for residue %d" % d)
      continue
    out_pdb = '%s/%s.pdb' % (out_dir, res_dir)
    if os.path.isfile(out_pdb):
      continue
    i = int(res_dir) - 1
    res_vals = res_data[i,:].copy()
    res_vals *= multiplier
    md_name = '%s/md' % res_dir
    print("Extracting last snapshot from trajectory of RIP on residue", end=' ') 
    print(res_dir)
    try:
      analysis.make_pdb_from_traj(md_name, i_ps, out_pdb, res_vals)
    except:
      print("Problem making", out_pdb)

  # write limits for later processing
  try:
    temp_back, dummy, temp_rip = sim_dir.split('/')[-1].split('-')
  except:
    temp_back = "can't find"
    temp_rip = "can't find"
  parms = {
      'min_val': min_val,
      'max_val': max_val,
      'multiplier': multiplier,
      'i_ps': i_ps,
      'temp_back': temp_back,
      'temp_rip': temp_rip
  }
  words = sim_dir.split('/')
  if len(words) > 3:
    parms['name'] = '/'.join(words[-3:])
  else:
    parms['name'] = sim_dir
  if raw_pdb:
    parms['pdb'] = os.path.basename(raw_pdb)[:4].upper()
  open(out_dir+'/parms', 'w').write(repr(parms))

  os.chdir(save_dir)


def make_pngs(
    data_dir, center=None, top=None, width=480, height=480, 
    frame_residues=[]):

  import plot
  import showpdb
  import pylab
  import Image

  def make_profile_pdb_img(
      png, soup, y_vals, xlabel='residue', ylabel='', 
      min_val=None, max_val=None, width=600, is_putty=False):

    if min_val is None:
      min_val = 0
    if max_val is None:
      max_val = max(y_vals)*1.1
    
    first_num=1
  
    temp_png = 'temp-1.png'
    pylab.clf()
    pylab.rc('axes', labelsize=12)
    pylab.rc('xtick', labelsize=12)
    pylab.rc('ytick', labelsize=12)
    pylab.rc('legend', fontsize=12)

    pylab.rc('xtick', direction='out')
    pylab.rc('ytick', direction='out')
    pylab.axes([0.125, 0.27, 0.3754, 0.08])
    plot.pylab_bins(y_vals, first_num, max_val)
    pylab.xlabel(xlabel)
    pylab.xlim([0, len(y_vals)])
    if ylabel: pylab.ylabel(ylabel)
    pylab.savefig(temp_png, dpi=600)

    im1 = Image.open(temp_png)
    (l, u, r, b) = im1.getbbox()
    im2 = im1.crop((l, int(b*.6), int(r*.6), int(b*.82)))
    im3 = plot.im_to_width(im2, width)

    out_png2 = png.replace('.png', '-pymol.png')
    n_res = soup.n_residue()
    y_vals = list(y_vals)
    while len(y_vals) < n_res:
      y_vals.append(0.0)
    showpdb.soup_to_bfactor_png(
        soup, out_png2, y_vals, min_val, max_val, is_putty=is_putty)  
  
    im4 = Image.open(out_png2)
    im5 = plot.im_to_width(im4, int(width*.75))
    
    new_h = im3.size[1] + im5.size[1] + 30
    im6 = Image.new("RGB", (width, new_h), "white")
    im6.paste(im5, (int(0.5*.25*width), 0))
    im6.paste(im3, (0, im5.size[1]))
    im6.save(png)

    os.remove(temp_png)

  os.chdir(data_dir)
  
  res_data = util.read_array('map')
  
  parms = util.read_dict('parms')
  min_val = parms['min_val']
  max_val = parms['max_val']
  multiplier = parms['multiplier']

  png = 'map.png' 
  if not os.path.isfile(png):
    plot.intensity_map_png(
        res_data, png, "", min_val, max_val, 
        "perturbing residue", "responding residue", 
        1, is_diagonal=True)
    print(png)

  ref_soup = showpdb.transformed_soup_from_pdb(
      out_dir + '/sim.pdb', center, top, width, height, 
      frame_residues)

  png = 'strength.png'
  if not os.path.isfile(png):
    data = util.read_column('strength')
    lo = data.mean() + data.std()
    make_profile_pdb_img(
        png, ref_soup, data, 'perturbing residue',
        'strength', lo, 30, 1000)
    print(png)

  png = 'flexibility.png'
  if not os.path.isfile(png):
    data = util.read_column('flexibility')
    make_profile_pdb_img(
        png, ref_soup, data, 'responding residue',
        'flexibility', 0, 15, 1000, True)
    print(png)

  png = 'fluctuation.png' 
  if not os.path.isfile(png):
    data = util.read_column('fluctuation')
    make_profile_pdb_img(
        png, ref_soup, data, 'responding residue',
        'av-fluct', None, 10, 1000, True)
    print(png)

  pdbs = util.re_glob('*.pdb', r'^\d+[.]pdb')
  for pdb in pdbs:
    png = pdb.replace('.pdb', '.png')
    if os.path.isfile(png):
      continue
    soup = showpdb.transformed_soup_from_pdb(
      pdb, center, top, width, height, frame_residues)
    res = soup.residue(int(pdb[:-4])-1)
    highlight = "%s:%d" % (res.chain_id, res.num)
    res_vals = [r.atom("CA").bfactor for r in soup.residues() if r.has_atom("CA")]
    showpdb.soup_to_bfactor_png(
        soup, png, res_vals, 0,
        15*multiplier, highlight, True)
    print(png)


template = """
<html>
<head>
  <title>${attributes['title']}</title>
  <link rel="stylesheet" type="text/css" 
        media="screen" 
        href="../../../../../doc/style.css" />
</head>

<body>
  %if 'header_html' in attributes:
    ${attributes['header_html']}
  %endif

  <div style="text-align:center; width:100%;
              font-weigth: normal;
              font-size: 2.5em;   
              font-family: Times, serif;
              letter-spacing:.1em; 
              padding-top:2em; padding-bottom:0em;">
    ${attributes['header']}
  </div>

  <div style="text-align:center; width:100%;
              font-weigth: normal;
              font-size:.80em;   
              padding-bottom:8em;">
    ${attributes['info']}
  </div>

  <center>
  <table align="center" cellspacing="0" 

         style="text-align:center; margin:0">
    <style>
      td {
         vertical-align: top; 
         text-align: center; 
         font-size: 1.00em; 
         padding-bottom: 5em;
      }
    </style>
    <tr>
      <td>
        column = C&alpha;-RMSD response to a RIP simulation<br>
        <img src='map.png' style='height:400px'>
      </td>
      <td>
        susceptibility of response to RIP of neighbors<br>
        <img src='flexibility.png' style='height:400px'>
      </td>
    </tr>
    <tr>
      <td>
        no. of responding residues for each RIP <br>
        <img src='strength.png' style='height:400px'>
      </td>
      <td>
        averaged C&alpha;-RMSD fluctuations above background<br>
        <img src='fluctuation.png' style='height:400px'>
      </td>
    </tr>
  </table>
  </center>
  
  <br>
  <br>
  
  <div style="text-align:center; font-size:1em; padding-bottom:1em;">
    Last snapshot of trajectory where C&alpha;-RMSD is 
    indicated by redness and thickness of chain.
  </div>
  <div style="padding:4em;">
    % for (name, png) in attributes['residue_pngs']:
      <div style="float:left; color: #AAA;
                  font-size: 0.75em; 
                  border:1px dashed #CCC; 
                  border-top:none; border-left:none; 
                  padding:5px">
        ${name}<br>
        <img src="${png}" style="width:200px">
      </div>
    % endfor
    <br clear=all>
  </div>
</body></html>
"""
def make_index_html(
    data_dir, header_html=None):
  import plot
  import showpdb

  attributes = {}

  parms = util.read_dict(data_dir + '/parms')

  if header_html:
    attributes['header_html'] = header_html
  
  title = "Large local perturbations of RIP in " + parms['name']
  attributes['title'] = title
  header = 'Response to RIP in ' + parms['name']
  attributes['header'] = header

  info = 'deviations averaged over the %dth ps of simulation' \
            % parms['i_ps']
  info += "<br>"
  info += "background temperature = " + parms['temp_back']
  info += "<br>"
  info += "rotational temperature = " + parms['temp_rip']
  info += "<br>"
  info += "significant response defined by C&alpha;-RMSD deviation "
  info += "> %.f &Aring;" % parms['min_val']
  if 'pdb' in parms:
    info += "<br>"
    info += "PDB accession code: " + parms['pdb']
  attributes['info'] = info

  save_dir = os.getcwd()
  util.goto_dir(data_dir)
  pngs = util.re_glob('*png', r'^\d+[.]png$')
  residues = [int(png[:-4]) for png in pngs]
  residues.sort()
  os.chdir(save_dir)
  
  sim_soup = pdbstruct.Soup('%s/sim.pdb' % data_dir)

  raw_pdb = data_dir + '/raw.pdb'
  if os.path.isfile(raw_pdb):
    raw_soup = pdbstruct.Soup(raw_pdb)
  else:
    raw_soup = None

  attributes['residue_pngs'] = []
  for k, i_res in enumerate(residues):
    png = str(i_res) + '.png' 
    res = sim_soup.residue(i_res - 1)
    name = "%s-%d" % (res.type, res.num)
    if raw_soup:
      res = raw_soup.residue(i_res - 1)
      name += " [%d]" % (res.num)
    attributes['residue_pngs'].append((name, png))

  plot.make_html(data_dir + '/index.html', template, attributes)


usage = """

  Copyright (c) 2007 Bosco Ho

  Calculates the flexibilities of a series of RIP perturbations

  Usage: flexibility.py [-h -c res -t res] sim_dir i_ps
  
  sim_dir:   directory that holds individual RIP perturbations
  i_ps:      the i'th ps that is averaged for the flexibilities

  -h:   also generates images and .html file (requires
        matplotlib, numpy, pymol, and the python imaging library)
  -c:   residue to center
  -t:   residue on top
"""

if __name__ == "__main__":
  import sys
  import getopt

  opts, args = getopt.getopt(sys.argv[1:], "hc:t:")
  if len(args)<2:
    print(usage)
    sys.exit(1)
  sim_dir = os.path.abspath(args[0])
  if not os.path.isdir(sim_dir):
    print(sim_dir, "doesn't exits")
    sys.exit(1)
  i_ps = int(args[1])
  
  is_html = False
  center_res = None
  top_res = None
  for opt, a in opts:
    if '-h' in opt:
      is_html = True
    if '-c' in opt:
      center_res = a
    if '-t' in opt:
      top_res = a

  out_dir = sim_dir + '/ca_rmsd.%dps' % i_ps
  if not os.path.isdir(out_dir):
    min_val, max_val = 6, 10
    make_rmsd_data_files(
        sim_dir, out_dir, i_ps, min_val, max_val)
  else:
    print("Trajectories already processed")
  if is_html:
    make_pngs(out_dir, center_res, top_res, 600, 600, [])
    make_index_html(out_dir)
