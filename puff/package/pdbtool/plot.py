#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pylab
import matplotlib
import Image, ImageDraw, ImageFont
import util
import math
import os
import glob
import math
import numpy
import string
import re
import codecs

# array handling


def flip_x(data):
  result = data.copy()
  n_x, n_y = numpy.shape(result)
  for i in range(n_x):
    for j in range(n_y):
      k = n_x - 1 - i
      result[i, j] = data[k, j]
  return result


def truncate_array(data, max_val, min_val):
  result = data.copy()
  n_x, n_y = numpy.shape(data)
  for i in range(n_x):
    for j in range(n_y):
      result[i,j] = min(data[i,j], max_val)
      if result[i,j] < min_val:
        result[i,j] = 0.0
  return result
truncate_2d_array = truncate_array
  
  
def symmetric_array(array):
  result = array.copy()
  n_x, n_y = array.shape
  for i in range(n_x):
    for j in range(n_y):
      if i<=j:
        v = max(array[i,j], array[j,i])
        result[i,j] = v
        result[j,i] = v
  return result
symmetric_map = symmetric_array


def filter_array(array, filter_fn):
  n_x, n_y = array.shape
  for i in range(n_x):
    for j in range(n_y):
      if filter_fn(i, j):
        array[i,j] = 0.0
filter_map = filter_array


def filtered_array(array, filter_fn, save_fname=""):
  if save_fname and os.path.isfile(save_fname):
    return util.read_array(save_fname)
  result = array.copy()
  filter_array(result, filter_fn)
  if save_fname:
    util.write_array(save_fname, result)
  return result
filtered_map = filtered_array


def map_array(array, fn):
  n_x, n_y = array.shape
  result = numpy.zeros((n_x, n_y), float)
  for i in range(n_x):
    for j in range(n_y):
      result[i,j] = fn(array[i,j])
  return result


# some stat functions

def count_bins(vals, n_bin, val_min, val_max):
  "Returns x_centers and bins of vals"
  delta = (val_max - val_min) / float(n_bin)
  bins = [0 for i in range(n_bin)]   
  for val in vals:
    i = int((val-val_min) / delta)
    if i >= n_bin:
      i = n_bin - 1
    if i < 0:
      i = 0
    bins[i] += 1
  x_centers = [val_min + (i+0.5)*delta for i in range(n_bin)]
  return x_centers, bins


def calculate_entropy(bins):
  total = sum(bins)
  probs = [v/float(total) for v in bins]
  terms = [p*math.log(p) for p in probs if p > 0.0]  
  return -sum(terms)


def average(vals):
  return sum(vals)/float(len(vals))


def sorted_sets(*args):
  """
  Given a bunch of lists of equal lengths, will sort
  w.r.t. to the values of the first list and return
  the lists
  """
  zipped_sets = list(zip(*args))
  zipped_sets.sort()
  return [list(i) for i in zip(*zipped_sets)]
 

def extract_floats(txt):
  strs = re.findall(r'[+|-]?\d+\.?\d*', txt)
  return [float(s) for s in strs]


# curve-fitting
  
def interpolate(x, x1, y1, x2, y2):
  slope = (y2 - y1)/(x2 - x1)
  return y1 + (x - x1)*slope


def interpolated_fn(x_vals, y_vals):
  """
  Returns a function for x that interpolates 
  the y_vals given evenly-spaced x_vals
  """
  x_min, x_max = min(x_vals), max(x_vals)
  dx = x_vals[1] - x_vals[0]
  x_range = x_max - x_min + dx
  n_bin = len(y_vals)

  def fn(x):
    k = n_bin*(x-x_min)/x_range
    l, m = math.floor(k), math.ceil(k)
    if m == l: 
      m = l+1
    if l >= n_bin or m >= n_bin:
      l, m = n_bin-2, n_bin-1
    elif l < 0 or m < 0:
      l, m = 0, 1
    l, m = int(l), int(m)
    slope = (y_vals[l]-y_vals[m])/float(l-m)
    return y_vals[l] + (k-l)*slope

  return fn


def polyfit(x_vals, y_vals, n_order=15):
  soln = numpy.polyfit(x_vals, y_vals, n_order)
  return [numpy.polyval(soln, x) for x in x_vals]


def low_pass_filter(x_vals, y_vals, freq_cutoff=20):
  freqs = numpy.fft.rfft(y_vals)
  n_freq = len(freqs)
  for i in range(freq_cutoff, n_freq):
    freqs[i] = 0
  fft_y_vals = [y.real for y in numpy.fft.ifft(freqs)]
  # rescale to original y_vals
  diff_y = average(y_vals) - average(fft_y_vals)
  fft_y_vals += diff_y
  # rescale to width of original x_vals for fn
  x_min, x_max = min(x_vals), max(x_vals)
  x_delta = (x_max - x_min)/n_freq
  fft_x_vals = [i*x_delta + x_min for i in range(n_freq)]
  return fft_x_vals, fft_y_vals
  
  
def interpolated_array_fn(array, x_min, x_max):
  """
  Assuming square matrix, returns a fn(x,y) that maps
  x,y values to the array values
  """
  n = array.shape[0]
  x_delta = (x_max - x_min)/float(n)
  i_from_x = lambda x: (x - x_min)/x_delta
  
  def i_limits(x):
    "Returns scaled i and nearest integers l, u in array"
    i = i_from_x(x) 
    l, u = math.floor(i), math.ceil(i)
    if l < 0 or u < 0:
      l, u = 0, 1
    if l >= n or u >= n:
      l, u = n-2, n-1
    return i, l, u

  def fn(x1, x2):
    i1, l1, u1 = i_limits(x1) 
    i2, l2, u2 = i_limits(x2) 
    if i1 % 1 == 0 and i2 % 1 == 0:
      y = array[i1, i2]
    elif i1 % 1 == 0:
      y = interpolate(i2, l2, array[i1, l2], u2, array[i1, u2])
    elif i2 % 1 == 0:
      y = interpolate(i1, l1, array[l1, i2], u1, array[u1, i2])
    else:
      y = array[l1, l2]*(u1-i1)*(u2-i2) + \
          array[u1, u2]*(i1-l1)*(i2-l2) + \
          array[l1, u2]*(u1-i1)*(i2-l2) + \
          array[u1, l2]*(i1-l1)*(u2-i2)
    return y
    
  return fn
  

# interface with the Python Imaging Library

def im_to_width(im, new_w):
  w, h = im.size
  ratio = w / float(h)
  new_h = int(new_w / ratio)
  return im.resize((new_w, new_h), Image.ANTIALIAS)


def resize_png(png, width):
  im = Image.open(png)
  new_im = im_to_width(im, width)
  new_im.save(png)


def insert_png(big_im, png, width, coords):
  if width:
    part_im = im_to_width(Image.open(png), int(width))
  else:
    part_im = Image.open(png)
  big_im.paste(part_im, coords)
  del part_im


def construct_plot(
    out_png, parts, width, height, 
    label_size=120, label_fill=(0,0,0),
    bg_color="white"):
  im = Image.new("RGB", (width, height), bg_color)
  draw = ImageDraw.Draw(im)
  f = ImageFont.truetype(
      "/Library/Fonts/Arial.ttf", label_size)
  for i, (x, y, w, png, label_x, label_y) in enumerate(parts):
    if not os.path.isfile(png):
      raise IOError("can't find %s" % png)
    insert_png(im, png, w, (x, y))
    s = string.ascii_uppercase[i]
    draw.text((x + label_x, y + label_y), s, font=f, fill=label_fill)
  if out_png.endswith('.tiff'):
    im.save(out_png, dpi=(300,300))
  else:
    im.save(out_png)


def shrink_png(png, l=0.0, u=0.0, r=1.0, b=1.0):
  im = Image.open(png)
  (abs_l, abs_u, abs_r, abs_b) = im.getbbox()
  im.crop((int(abs_r*l), int(abs_b*u), 
           int(abs_r*r), int(abs_b*b))).save(png)


# interface with pylab/matplotlib

def locs_labels(n, first, major=50, minor=10):
  locs = []
  labels = []
  for i in range(n):
    val = i+first
    if i == 0 or (val % minor) == 0:
      locs.append(i)
      if i == 0 or (val % major) == 0:
        labels.append(val)
      else:
        labels.append("")
  return locs, labels  


def my_color_map():
  # color_map = pylab.get_cmap("OrRd")
  cdict = { \
      'red':   ((0.0, 1.0, 1.0),
                (1.0, 1.0, 1.0)),
      'green': ((0.0, 1.0, 1.0),
                (1.0, 0.0, 0.0)),
      'blue':  ((0.0, 1.0, 1.0),
                (1.0, 0.0, 0.0))}
  color_map = matplotlib.colors.LinearSegmentedColormap(
      'my_color_map', cdict, 256)
  return color_map
red_color_map = my_color_map()  
  
  
def pylab_array(
    raw_data, min_val, max_val, first_num=1, 
    major=50, minor=10, color_map=red_color_map,
    title='', xlabel=None, ylabel=None, gap_z=None, 
    is_legend=True, is_diagonal=False):

  pylab.rc('xtick', direction='in')
  pylab.rc('ytick', direction='in')
  data = flip_x(raw_data.transpose())
  data = truncate_array(data, max_val, min_val)
  pylab.imshow(
      data, cmap=color_map, 
      vmin=min_val, vmax=max_val, 
      interpolation='nearest')

  locs, labels = locs_labels(
      data.shape[0], first_num, major, minor)
  pylab.xticks(locs, labels)
  last = data.shape[0]-1
  new_locs = [last - l for l in locs]
  pylab.yticks(new_locs, labels)

  if is_legend:
    if gap_z:
      z = int(min_val/gap_z)*gap_z
      z_ticks = [z]
      while z<max_val:
        z_ticks.append(z)
        z += gap_z
      pylab.colorbar(ticks=z_ticks)
    else:
      pylab.colorbar()

  pylab.title(title)
  if xlabel: pylab.xlabel(xlabel)
  if ylabel: pylab.ylabel(ylabel)

  if is_diagonal:
    x_vals, y_vals = [], []
    n = raw_data.shape[0]
    size = int((350/float(n))**2)
    if size==0:
      size = 1
    for i in range(n):
      if sum(raw_data[i,:]) > 0:
        x_vals.append(i)
        y_vals.append(n-i-1)
    if x_vals:
      pylab.scatter(
          x_vals, y_vals, s=size, marker='o', 
          facecolor='g', edgecolor='none') 
      pylab.xlim([-0.5, n-0.5])
      pylab.ylim([n, 0.5])

  pylab.xlim([-0.5, raw_data.shape[0]-0.5])
  pylab.ylim([raw_data.shape[1]-0.5, -0.5])

pylab_2dmap = pylab_array


def intensity_map_png(
    raw_data, out_png, title='', 
    min_val=None, max_val=None, 
    xlabel=None, ylabel=None, gap_z=None, 
    first_num=1, is_legend=True,
    major=50, minor=10, color_map=red_color_map,
    is_diagonal=False):

  if not min_val: min_val = raw_data.min()
  if not max_val: max_val = raw_data.max()
  
  pylab.rc('xtick', direction='in')
  pylab.rc('ytick', direction='in')
  pylab.rc('axes', labelsize=25)
  pylab.rc('axes', titlesize=25)
  pylab.rc('xtick', labelsize=20)
  pylab.rc('ytick', labelsize=20)
  pylab.rc('legend', fontsize=20)

  pylab.clf()
  pylab.axes([0.08, 0.12, 0.8, 0.8])
  pylab_array(
    raw_data, min_val, max_val, first_num, 
    major, minor, color_map,
    title, xlabel, ylabel, gap_z, 
    is_legend, is_diagonal)
  pylab.savefig(out_png, dpi=300)

  im = Image.open(out_png)
  (l, u, r, b) = im.getbbox()
  new_im = im.crop((l, u, int(r*.9), b))
  new_im.save(out_png)


def pylab_bars(
    x_vals, y_vals, in_color='r', title=None,
    xlabel=None, ylabel=None):
  dx = x_vals[1]-x_vals[0]
  x_min, x_max = min(x_vals), max(x_vals)
  dx_half = 0.5*dx
  x_left_vals = [x-dx_half for x in x_vals]
  pylab.bar(
      x_left_vals, y_vals, dx, 
      edgecolor=in_color, facecolor=in_color)
  pylab.xlim([x_min-1.01*dx_half, x_max+1.01*dx_half])
  y_min = min(y_vals)
  y_max = max(y_vals)
  dy = (y_max - y_min)*0.001
  pylab.ylim([y_min-dy, y_max+dy])
  lines = pylab.getp(pylab.gca(), 'xticklines')
  toplines = [l for l in lines if pylab.getp(l, 'ydata') == (1,)]
  pylab.setp(toplines, visible=False)
  if title: pylab.title(title)
  if xlabel: pylab.xlabel(xlabel)
  if ylabel: pylab.ylabel(ylabel)
  

def pylab_bins(bins, first_num, y_max=None, in_color='r'):
  n_bin = len(bins)
  pylab_bars(range(n_bin), bins, in_color)
  pylab.xticks(*locs_labels(n_bin, first_num))
  if y_max is None:
    y_max = max(bins)
  if y_max is not None:
    pylab.ylim([0, 1.001*y_max])
    pylab.yticks([0, y_max])


def make_trace_png(
    out_png, x_vals, y_vals_list, xlabel=None, ylabel=None, 
    min_y=None, max_y=None, gap_y=None, legend_labels=None,
    xticks=None, xtick_labels=None, highlight_x=None):

  params = {'axes.labelsize': 15,
            'text.fontsize': 15,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            }
  pylab.rcParams.update(params)
  pylab.clf() 
  pylab.axes([.1, .1, .5, .2])

  for y_vals in y_vals_list:
    pylab.plot(x_vals, y_vals)

  # plot a minimium horizontal line
  if min_y:
    pylab.plot(x_vals, [min_y for j in range(len(x_vals))])

  pylab.xticks(xticks, xtick_labels)
  if gap_y:
    pylab.yticks([i*gap_y for i in range(1000)])

  if not max_y:
    max_y = max([max(y_vals) for y_vals in y_vals_list]) * 1.1
  pylab.ylim(0, max_y)

  pylab.xlim([min(x_vals), max(x_vals)])

  if highlight_x:
    y = max_y * 0.7 
    pylab.arrow(highlight_x, max_y*0.7, 0, -max_y*0.15, 
                edgecolor="red", facecolor="red", width=0.5)
    pylab.text(highlight_x, max_y*0.75, "RIP", color="red")

  if xlabel: pylab.xlabel(xlabel)
  if ylabel: pylab.ylabel(ylabel)

  if legend_labels:
    pylab.legend(legend_labels, loc=(.77, .7))
    legend = pylab.gca().get_legend()
    legend.draw_frame(False)

  pylab.savefig(out_png, dpi=300)
  shrink_png(out_png, 0, .65, .64, 1)


def pylab_save_and_show(png):
  pylab.savefig(png, dpi=72)
  os.system('open ' + png)
  

# some decorators to help process image generation

def skip_if_png_exists(fn):
  def new_fn(*args, **kwargs):
    for a in args:
      if isinstance(a, str):
        if a.endswith('png'):
          if os.path.isfile(a):
            return
    return fn(*args, **kwargs)
  return new_fn


def success_on_png(fn):
  def new_fn(*args, **kwargs):
    not_fnames = []
    for a in args:
      if isinstance(a, str):
        if not os.path.isfile(a):
          not_fnames.append(a)
    output = fn(*args, **kwargs)
    for f in not_fnames:
      if os.path.isfile(f):
        print "Made", os.path.basename(f)
    return output
  return new_fn
      

  
# Interface with MAKO to make html from template

def make_html(fname, template, attrs):
  if 'Template' not in globals():
    from mako.template import Template
  html = Template(template).render_unicode(attributes=attrs)
  codecs.open(fname, encoding='utf-8', mode='w').write(html)
  
  
  

if __name__ == "__main__":
  import sys, getopt

  opts, args = getopt.getopt(sys.argv[1:], "x:")
  if len(args) < 1:
    print "Usage plot.py [-x xval_col] data_fname [y_col1 y_col2...]"
    sys.exit(1)

  data_fname = args[0]
  if not os.path.isfile(data_fname):
    print "Can't find", data_fname
    sys.exit(1)
  
  data = util.read_array(data_fname)

  x_vals = range(data.shape[1])
  for opt, val in opts:
    if '-x' in opt:
      x_vals = data[int(val),:]      

  if len(args) > 1:
    for i in args[1:]:
      pylab.plot(x_vals, data[int(i),:])
  else:
    pylab.plot(x_vals, data[0,:])

  png = data_fname + '.png'
  pylab_save_and_show(png)

