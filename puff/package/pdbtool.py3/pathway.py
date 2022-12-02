#!/usr/bin/env python
# encoding: utf-8

import os
import optparse
import shutil
import numpy
import util
import pdbstruct
import vector3d
import analysis
import math
import asa


def get_cutoff_from_map(map_data):
  def get_avg_std(vals):
    n = float(len(vals))
    avg = sum(vals)/n
    diff_sqs = [(v-avg)**2 for v in vals]
    var = sum(diff_sqs)/n
    std = math.sqrt(var)
    return avg, std

  def get_cutoff(avg, std):
    return avg + 2.0*std

  data = [d for d in map_data.flatten() if d > 0.0]
  return get_cutoff(*get_avg_std(data))


def n_above_threshold(vals, low):
  return sum([1 for v in vals if v >= low])


def n_pair_in_map(in_map, low):
  return n_above_threshold(in_map.flatten(), low)


def n_pair_in_symm_map(sym_map, low):
  return n_pair_in_map(sym_map, low) // 2


def block_2d_array(data, min_val):
  result = data.copy()
  n_x, n_y = numpy.shape(data)
  for i in range(n_x):
    for j in range(n_y):
      if result[i,j] >= min_val:
        result[i,j] = 1
      else:
        result[i,j] = 0
  return result


def count_n_vals_per_column(res_data, min_val):
  data = block_2d_array(res_data, min_val)
  n_bin = data.shape[0]
  for i in range(n_bin):
    data[i,i] = 0
  result = numpy.zeros(n_bin, int)
  for i in range(n_bin):
    result[i] = sum(data[i,:])
  return result


def make_pathway_pdbs(
   pathway_map, max_val, in_dir, out_dir, i_ps):
  """
  Returns the multiplier for bfactors written to PDB files
  """
  median = 1.0
  multiplier = 1.0
  while median > max_val:
    multiplier *= 10.0
    median = 1.0 / multiplier
  for i in range(pathway_map.shape[0]):
    i_res = i+1
    res_dir = in_dir + '/%d' % i_res
    if os.path.isdir(res_dir):
      pdb = '%s/%d.pdb' % (out_dir, i_res)
      if not os.path.isfile(pdb):
        analysis.make_pdb_from_traj(
          '%s/md' % res_dir, i_ps, pdb,
          multiplier*pathway_map[i,:].copy())
        print("Made", os.path.basename(pdb))
  return multiplier
  
  
def is_connected(i, j, soup, cutoff=3.5):
  if i == j:
    return False
  min_dist = 1000.0
  for atom_i in soup.residue(i).atoms():
    for atom_j in soup.residue(j).atoms():
      dist = vector3d.pos_distance(atom_i.pos, atom_j.pos)
      if dist < min_dist:
        min_dist = dist
  return min_dist < cutoff


backbone = ['CA', 'HA', 'N', 'H', 'O', 'C']
def is_sidechain_connected(i, j, soup, cutoff=3.5):
  if abs(i-j) <= 2:
    return False
  min_dist = 1000.0
  sidechain_atoms_i = [a for a in soup.residue(i).atoms() 
                       if a.type not in backbone]
  for atom_i in sidechain_atoms_i:
    for atom_j in soup.residue(j).atoms():
      dist = vector3d.pos_distance(atom_i.pos, atom_j.pos)
      if dist < min_dist:
        min_dist = dist
  return min_dist < cutoff


def merge_residues(residues1, residues2):
  for i in residues2:
    if i not in residues1:
      residues1.append(i)
  
  
def get_connected_residues(i, residues, soup):
  connected_residues = [i]
  neighbors = [j for j in residues if is_connected(i, j, soup)]
  connected_residues.extend(neighbors)
  non_neighbors = [j for j in residues if j not in neighbors]
  if non_neighbors:
    for j in neighbors:
      connected_to_j = \
          get_connected_residues(j, non_neighbors, soup)
      merge_residues(connected_residues, connected_to_j)
  return connected_residues


def indices_above_threshold(data, i, min_val):
  n = data.shape[0]
  return [j for j in range(n) if data[i,j] >= min_val and i != j]


def neighbors_above_threshold(i, data, min_val, soup):
  residues = indices_above_threshold(data, i, min_val)
  return [j for j in residues if is_sidechain_connected(i, j, soup)]

  
def filter_disconnected_residues(i, data, min_val, soup):
  n_res = data.shape[0]
  residues = indices_above_threshold(data, i, min_val)
  neighbors = [j for j in residues 
               if is_sidechain_connected(i, j, soup)]
  non_neighbors = [j for j in residues if j not in neighbors]
  connected_residues = [i]
  for j in neighbors:
    connected_to_j = \
        get_connected_residues(j, non_neighbors, soup)
    merge_residues(connected_residues, connected_to_j)
  for j in range(n_res):
    if j not in connected_residues:
      data[i, j] = 0.0    


def filter_indirectly_coupled(i, data, min_val, soup):
  n_res = data.shape[0]
  residues = indices_above_threshold(data, i, min_val)
  neighbors = [j for j in residues 
               if is_sidechain_connected(i, j, soup)]
  for j in range(n_res):
    if j not in neighbors:
      data[i, j] = 0.0    


def filter_short_pathways(i, data, min_val):
  n_res = data.shape[0]
  if n_above_threshold(data[i,:], min_val) < 3:
    for j in range(n_res):
      data[i, j] = 0.0


def filter_immediate_neighbors(i, data, min_val):
  n_res = data.shape[0]
  test_j = [j for j in range(i-1, i+1+1) 
            if j>=0 and j<n_res and i!=j]
  for j in test_j:
    data[i,j] = 0.0


def pathway_and_coupling_from_heatflow(
    heatflow_map, min_val, out_dir, pathway_fname, coupling_fname):
  if os.path.isfile(pathway_fname) \
      and os.path.isfile(coupling_fname):
    pathway_map = util.read_array(pathway_fname)
    coupling_map = util.read_array(coupling_fname)
  else:
    pathway_map = heatflow_map.copy()
    n_res = pathway_map.shape[0]
    coupling_map = numpy.zeros((n_res, n_res), float)
    for i in range(n_res):
      pathway_pdb = '%s/%d.pdb' % (out_dir, i+1)
      if not os.path.isfile(pathway_pdb):
        for j in range(n_res):
          pathway_map[i,j] = 0.0
      else:
        pathway_soup = pdbstruct.Soup(pathway_pdb)
        for j in range(n_res):
          if pathway_map[i,j] <= min_val:
            pathway_map[i,j] = 0.0
        filter_immediate_neighbors(i, pathway_map, min_val)
        filter_disconnected_residues(
            i, pathway_map, min_val, pathway_soup)
        filter_short_pathways(i, pathway_map, min_val)
        for j in neighbors_above_threshold(
            i, pathway_map, min_val, pathway_soup):
          coupling_map[i,j] = 1.0
        del pathway_soup
    util.write_array(pathway_fname, pathway_map)
    util.write_array(coupling_fname, coupling_map)
  return pathway_map, coupling_map


def get_map_from_res_dirs(in_dir, data_fname, i_ps, save_fname):
  if os.path.isfile(save_fname):
    data = util.read_array(save_fname)
  else:
    save_dir = os.getcwd()
    util.goto_dir(in_dir)
    res_dirs = util.re_glob('*', r'^\d+$')
    for res_dir in res_dirs:
      res_data_fname = res_dir + '/' + data_fname
      if not os.path.isfile(res_data_fname):
        os.chdir(res_dir)
        print("Extracting variables from", res_dir)
        analysis.process_md('md')
        os.chdir('..')
    print(os.getcwd())
    data = util.get_numbered_dir_data(data_fname, i_ps)
    util.goto_dir(save_dir)
    util.write_array(save_fname, data)
  return data
  

def get_filtered_map_and_save(data, filter_fn, fname):
  if os.path.isfile(fname):
    return util.read_array(fname)
  else:
    filtered_map = util.filtered_array(data, filter_fn)
    util.write_array(fname, filtered_map)
    return filtered_map


def find_first_pdb_file(search_dir):
  os.chdir(search_dir)
  pdbs = util.re_glob('*pdb', r'^\d\w\w\w.*[.]pdb$')
  pdbs = [p for p in pdbs if 'hollow' not in p]
  if len(pdbs):
    return pdbs[0]
  else:
    return None


def is_backbone_hbonded(soup, i, j):
  if i==j:
    return False

  def get_hbond_atoms_i(i):
    result = []
    n_res = len(soup.residues())
    atom_list = [(i, "O"), (i, "H"),
                 (i-1, "O"), (i-1, "H"), (i-2, "O"),
                 (i+1, "H"), (i+1, "O"), (i+2, "H")]
    for j, atom_type in atom_list:
      if j>=0 and j<n_res:
        if soup.residue(j).has_atom(atom_type):
          result.append(soup.residue(j).atom(atom_type))
    return result

  hbond_atoms_i = get_hbond_atoms_i(i)
  hbond_atoms_j = get_hbond_atoms_i(j)
  n_hbond = 0
  for atom_i in hbond_atoms_i:
    for atom_j in hbond_atoms_j:
      if atom_i.type != atom_j.type:
         if vector3d.pos_distance(atom_i.pos, atom_j.pos) < 3.5:
            n_hbond += 1
  return n_hbond >= 2


def generate_datafiles(
    in_dir, data_fname, pdb, out_dir, i_ps):
  if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

  sim_pdb = out_dir + '/sim.pdb'
  if not os.path.isfile(sim_pdb):
    shutil.copy(pdb, out_dir + '/sim.pdb')
    raw_pdb = find_first_pdb_file(os.path.dirname(pdb))
    if raw_pdb:
      shutil.copy(raw_pdb, out_dir + '/raw.pdb')
  raw_pdb = out_dir + '/raw.pdb'
  if not os.path.isfile(raw_pdb):
    raw_pdb = None
  soup = pdbstruct.Soup(sim_pdb)
  n_res = soup.n_residue()

  is_buried = asa.get_is_buried_state(
        soup, pdb.replace('.pdb', '.buried'))

  heatflow_map_fname = out_dir + '/heatflow_map'
  if os.path.isfile(heatflow_map_fname):
    heatflow_map = util.read_array(heatflow_map_fname)
  else:
    heatflow_map = get_map_from_res_dirs(
        in_dir, data_fname, i_ps, out_dir+'/heatflow_map')
    util.write_array(heatflow_map_fname, heatflow_map)
  min_val = get_cutoff_from_map(heatflow_map)
  max_val = min_val*2.0
  multiplier = \
      make_pathway_pdbs(
          heatflow_map, max_val, in_dir, out_dir, i_ps)

  pathway_map, coupling_map = \
      pathway_and_coupling_from_heatflow(
          heatflow_map, min_val, out_dir, 
          out_dir+'/pathway_map', 
          out_dir + '/coupling_map')
  buried_coupling_map = get_filtered_map_and_save(
      coupling_map, 
      lambda i,j: not is_buried[i] or not is_buried[j],
      out_dir + '/buried_coupling_map')  
  tertiary_coupling_map = get_filtered_map_and_save(
      coupling_map, 
      lambda i,j: is_backbone_hbonded(soup, i, j) \
                    or abs(i-j)<=3,
      out_dir + '/tertiary_coupling_map')  
  buried_tertiary_coupling_map = get_filtered_map_and_save(
      buried_coupling_map, 
      lambda i,j: is_backbone_hbonded(soup, i, j) \
                    or abs(i-j)<=3,
      out_dir + '/buried_tertiary_coupling_map')

  for i in range(n_res):
    partners = indices_above_threshold(coupling_map, i, 0.5)
    if partners:
      i_res = i+1
      data = out_dir + '/%s-kin.ave' % i_res
      if not os.path.isfile(data):
        res_dir = in_dir + '/%d' % i_res
        raw_data = '%s/kin.ave' % res_dir
        shutil.copy(raw_data, data)

  parms = {
     'min_val': min_val,
     'max_val': max_val,
     'multiplier': multiplier,
     'title':  "%s (%.3f-%.3f)" % (data_fname, min_val, max_val)
  }
  util.write_dict(parms, out_dir + '/parms')
  
  
template = """
<html>
<head><title>${attributes['name']}</title></head>
<body>

  <div style="font-size:2em; padding:2em; text-align:left">
    ${attributes['name']}
  </div>

  % for pngs in attributes['top_part']:
    <table>
    % for png in pngs:
      <td style="width:350px; padding-bottom:1.5em; text-align:center">
        ${png.replace('_', ' ').replace('.png', '').replace(' map', '')}
        <br>
        <img src="${png}" style="width:350px">
      </td>
    % endfor
    </table>
  % endfor

  % for (text, i, couplings) in attributes['couplings']:
  <table style="width:100%; border-top:1px solid #CCC">
    <td style="width:260px; vertical-align:top;">
      <span style="font-size:1.5em;">${text}</span>
      <br />
      <img src="${i}-kin-slice.png">
    </td>
    <td style="width:260px; vertical-align:top;">
      <img src="${i}-pathway.png">
    </td>
    <td style="border-left:1px solid #CCC">
      % for coupling in couplings:
      <img src="${coupling}">
      % endfor
    </td>
  </table>
  % endfor

</body></html>
"""


def generate_images_and_html(name, out_dir, center, top, frame_residues):

  import plot
  import showpdb

  @plot.skip_if_png_exists
  @plot.success_on_png
  def make_rip_slice_png(
      png, data_fname, i_ps_list, max_y, min_y, gap_y, highlight_x):
    data = util.read_array(data_fname)
    y_vals_list = [data[i_ps,:] for i_ps in i_ps_list]
    n_ps, n_res = data.shape
    x_vals = list(range(1, n_res+1))
    legend_labels = ["%d ps" % i for i in i_ps_list]
    # make xticks and xtick_labels
    first_num = 1
    locs = [first_num]
    tens = [ \
        i for i in range(first_num, n_res+first_num) 
        if (i % 10) == 0]
    locs.extend(tens)
    xticks = [loc-first_num+1 for loc in locs]
    xtick_labels = []
    for i in locs:
      if i % 50 == 0 or i == first_num:
        label = '%d' % i
      else:
        label = ''
      xtick_labels.append(label)
    xlabel = 'responding residue'
    ylabel = os.path.basename(data_fname).replace('.ave', '')
    if 'ca_rmsd' in ylabel:
      ylabel = "Ca RMSD [Å]"
    elif 'rmsd' in ylabel:
      ylabel = "RMSD [Å]"
    elif 'kin' in ylabel:
      ylabel = "kin [kcal/mol]"
    plot.make_trace_png(
        png, x_vals, y_vals_list, xlabel, ylabel, 
        min_y, max_y, gap_y, legend_labels,
        xticks, xtick_labels)
    plot.resize_png(png, 250)


  @plot.skip_if_png_exists
  @plot.success_on_png
  def make_imgs_from_map(fname, soup):
    filtered_map = util.read_array(fname)
    plot.intensity_map_png(
        filtered_map, fname+'.png', 
        "%d couplings" % n_pair_in_map(filtered_map, 0.5),
         0.5, 1.0,
        'perturbing residue', 'responding residue', 
        is_legend=False)
    png = fname+'_profile.png'
    display_map = plot.symmetric_array(filtered_map)
    y_vals = count_n_vals_per_column(display_map, 0.5)
    soup.load_residue_bfactors(y_vals)
    soup.write_pdb(png.replace('.png', '.pdb'))
    showpdb.soup_to_bfactor_png(soup, png, y_vals,  0.5, 1)

  
  @plot.skip_if_png_exists
  @plot.success_on_png
  def plot_map_png_from_file(
      in_fname, png, title, min_val, max_val,
      xlabel, ylabel, z_gap, is_diagonal=False, is_legend=True):
    data = util.read_array(in_fname)
    plot.intensity_map_png(
        data, png, title, min_val, max_val,
        xlabel, ylabel, z_gap, 
        is_diagonal=is_diagonal, is_legend=is_legend)


  @plot.skip_if_png_exists
  @plot.success_on_png
  def make_pathway_img(
      data, i, min_val, max_val, multiplier, soup, out_png):
    b_factors = []
    max_x = multiplier*max_val
    for val in data[i,:]:
      if val > max_val:
        b_factors.append(max_x)
      elif val < min_val:
        b_factors.append(0.0)
      else:
        b_factors.append(multiplier*val)
    showpdb.soup_to_bfactor_png(
        soup, out_png, b_factors, 
        min_val*multiplier,
        max_val*multiplier, 
        showpdb.get_res_name(soup.residue(i)))
    plot.resize_png(out_png, 250)


  @plot.skip_if_png_exists
  @plot.success_on_png
  def make_coupling_img(soup, i, j, png):
    bfactors = [0.0 for k in range(soup.n_residue())]
    bfactors[i] = 2.0
    bfactors[j] = 1.0
    showpdb.soup_to_bfactor_png(
        soup, png, bfactors, 0.5, 2.0, 
        showpdb.get_res_name(soup.residue(i)))
    plot.resize_png(png, 250)


  parms = util.read_dict(out_dir + '/parms')
  
  plot_map_png_from_file(
      out_dir + '/heatflow_map', out_dir + '/heatflow_map.png', 
      parms['title'], parms['min_val'], parms['max_val'],
      'perturbing residue', 'responding residue', 0.01,
      is_diagonal=True)
  
  plot_map_png_from_file(
      out_dir + '/pathway_map',
      out_dir + '/pathway_map.png',
      '',
      parms['min_val'], parms['max_val'],
      'perturbing residue', 'responding residue', 0.01,
      is_diagonal=True)
  
  soup = showpdb.transformed_soup_from_pdb(
      out_dir + '/sim.pdb', 
      center, top, 800, 800, 
      frame_residues=frame_residues)
  make_imgs_from_map(
      out_dir+'/coupling_map', soup)  
  make_imgs_from_map(     
      out_dir + '/buried_coupling_map', soup)  
  make_imgs_from_map(     
      out_dir + '/tertiary_coupling_map', soup)  
  make_imgs_from_map(     
      out_dir + '/buried_tertiary_coupling_map', soup)
      
  attributes = {'name': name}

  attributes['top_part'] = \
      [['heatflow_map.png', 'pathway_map.png', 'coupling_map.png'],
       ['buried_coupling_map.png', 'tertiary_coupling_map.png',
        'buried_tertiary_coupling_map.png'],
       ['buried_coupling_map_profile.png',
        'tertiary_coupling_map_profile.png',
        'buried_tertiary_coupling_map_profile.png']]

  attributes['couplings'] = []

  raw_soup = None
  if os.path.isfile(out_dir + '/raw.pdb'):
    raw_soup = pdbstruct.Soup(out_dir + '/raw.pdb')
  coupling_map = util.read_array(out_dir+'/coupling_map') 
  pathway_map = util.read_array(out_dir+'/pathway_map') 

  for i in range(soup.n_residue()):
    partners = indices_above_threshold(coupling_map, i, 0.5)
    if partners:
      i_res = i+1
      print(">>>>> Residue", i_res, "couples to", \
            [k+1 for k in partners])
      make_rip_slice_png(
          '%s/%d-kin-slice.png' % (out_dir, i_res),
          out_dir + '/%s-kin.ave' % i_res, 
          [5], parms['max_val'], parms['min_val'], 0.04, i_res)
      pathway_soup = showpdb.transformed_soup_from_pdb(
          '%s/%d.pdb' % (out_dir, i_res),
          center, top, frame_residues=frame_residues)
      res = pathway_soup.residue(i)
      text = "%s-%s" % (res.type, showpdb.get_res_name(res))
      if raw_soup:
        text += " (%s)" % showpdb.get_res_name(raw_soup.residue(i))
      make_pathway_img(
          pathway_map, i, parms['min_val'], parms['max_val'],
          parms['multiplier'], 
          pathway_soup, 
          '%s/%d-pathway.png' % (out_dir, i_res))
      coupling_pngs = []
      for j in partners:
        png = '%d-coupling-%d.png' % (i_res, j+1)
        make_coupling_img(
            pathway_soup, i, j, '%s/%s' % (out_dir, png))
        coupling_pngs.append(png)
      del pathway_soup
      attributes['couplings'].append((text, i_res, coupling_pngs))

  plot.make_html(out_dir + '/index.html', template, attributes)



usage = """
pathway.py [-h] dir ps
"""


if __name__ == "__main__":
  import sys
  if len(sys.argv) == 1:
    print(usage)
    sys.exit(1)

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

  ref_pdb = sim_dir + '/../sim.pdb'
  out_dir = sim_dir + '/ca_rmsd.%dps' % i_ps
  name, parms = sim_dir.split('/')[-2:]
  out_dir = '%s/kin.%dps' % (sim_dir, i_ps)

  print(">>>>>>>>>>", out_dir)
  generate_datafiles(
      sim_dir, 'kin.ave', ref_pdb, out_dir, i_ps)            
  if is_html:
    div = " &nbsp; : &nbsp; "
    name = sim_dir + div + "slow RIP" 
    name += div + "T<sub>back</sub>=10K"
    name += div + "T<sub>rot</sub>=26K"
    generate_images_and_html(
        name, out_dir, center_res, top_res, [])

