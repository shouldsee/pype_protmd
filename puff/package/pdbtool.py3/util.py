import os
import stat
import copy
import tempfile
import re
import glob
import subprocess
import time
import numpy


def re_glob(dir_tag, reg_exp=""):
  fnames = glob.glob(dir_tag)
  return [f for f in fnames if re.search(reg_exp, f)]


def wait_for_stupid_nfs_system():
  txt = run_with_output('ls')
  while 'Stale NFS file handle' in txt:
    txt = run_with_output('ls')


def goto_dir(new_dir):
  save_dir = os.getcwd()
  wait_for_stupid_nfs_system()
  if not os.path.isdir(new_dir):
    os.makedirs(new_dir)
  os.chdir(new_dir)
  wait_for_stupid_nfs_system()


def insert_path(path, insert):
  if path.startswith('/'):
    return path
  else:
    return os.path.join(insert, path)


def temp_fname(suffix=''):
  fd, fname = tempfile.mkstemp(suffix, 'tmp-', '.')
  f = os.fdopen(fd, 'w')
  f.close()
  os.unlink(fname)
  return os.path.basename(fname)


def fname_variant(fname):
  root, ext = os.path.splitext(fname)
  i = 1
  new_fname = "%s-%d%s" % (root, i, ext)
  while os.path.isfile(new_fname):
    i += 1
    new_fname = "%s-%d%s" % (root, i, ext)
  return new_fname


def clean_fname(fname):
  try:
    os.remove(fname)
  except:
    pass


def run_with_output(cmd):
  p = subprocess.Popen(
      cmd, shell=True, stdout=subprocess.PIPE, 
      stderr=subprocess.PIPE)
  return p.stdout.read()


def run_with_output_file(cmd, out_fname=None):
  if not out_fname:
    txt = run_with_output(cmd)
    return
  sh_file = out_fname + '.sh'
  log_file = out_fname + '.log'
  cmd_log = cmd + ' >& ' + log_file
  open(sh_file, 'w').write(cmd_log)
  S_IRWXU = stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
  os.chmod(sh_file, S_IRWXU)
  stopwatch = Timer()
  os.system(cmd_log)
  stopwatch.stop()
  open(out_fname + '.time', 'w').write(stopwatch.str())
 

def get_floats_from_string(s):
  val_strs = re.findall(r'[-+]?([0-9]*\.[0-9]+|[0-9]+)', s)
  return [float(v) for v in val_strs]
  
  
def is_same_dict_in_file(parms, fname):
  try:
    saved_parms = eval(open(fname).read())
    result = saved_parms == parms
  except:
    result = False
  return result


def replace_dict(t, d, max_len = 0):
  """Replace $token in string t using a dictionary

  By default, converts floats to %.3f format:
    - t: string whose tokens will be replace_dict
    - d: dictionary of tokens (string keys) and their 
         replacements (values, may be of any type)
    - max_len: maximum length (in characters) of 
               replacement values.
  """
  s = copy.deepcopy(t)
  for (k, v) in d.items():
    token = "$%s" % k
    if max_len == 0:
      if type(v) is float:
        s = s.replace(token, "%.3f" % v)
      else:
        s = s.replace(token, str(v))
    else:
      if type(v) is float:
        s = s.replace(token, ("%.3f" % v).rjust(max_len))
      else:
        s = s.replace(token, str(v).rjust(max_len))
  return s


def write_dict(d, fname):
  """Outputs the dictionary to a file in repr format."""
  open(fname, "w").write(format_dict(d))


def read_dict(fname):
  return eval(open(fname, 'r').read())
  

def format_dict(d):
  """Makes repr() comaptible string of a dictionary"""
  s = "{ "
  n = len(d)
  items = sorted(d.items())
  max_len = max([len(repr(k)) for k in list(d.keys())])
  for i in range(n):
    if i > 0:
      s += "  "
    key, value = items[i]
    s += "%s : %s" % (repr(key).ljust(max_len), 
                      repr(value))
    if i < n-1:
      s += ", \n"
  s += " }"
  return s


def words_in_file(fname):
  result = []
  for line in open(fname).readlines():
    result.extend(line.split())
  return result


def read_parameters(fname):
  class DataHolder: pass
  f = open(fname, 'r')
  result = DataHolder()
  result.__dict__ = eval(f.read())
  f.close()
  return result


def elapsed_time_str(time):
  s = str(time) + ' '
  minute = time / 60.0
  if minute > 60:
    s += "%.f:%02.f:" % (minute / 60, minute % 60)
  elif minute >= 1:
    s += "%.f:" % minute
  sec = time % 60.0
  if sec < 0.01:
    s += "%07.4fs" % sec
  else:
    s += "%05.2fs" % sec
  return s
  
  
class Timer:
  def __init__(self):
    self._elapsed = 0;
    self._start = time.time()

  def start(self):
    self._start = time.time()
    self._elapsed = 0

  def stop(self):
    self._elapsed = time.time() - self._start

  def elapsed(self):
    if self._elapsed == 0:
      return time.time() - self._start
    else:
      return self._elapsed

  def str(self):
    elapsed_time = self.elapsed()
    return elapsed_time_str(elapsed_time)
    
  def __str__(self):
    return self.str()


def load_psyco_if_in_system():
  try:
    import psyco
    psyco.full()
    print("Loaded optional psyco JIT compiler")
  except:
    pass


def val_range(start, end, step):
  vals = []
  v = start
  while v <= end:
    vals.append(v)
    v += step
  return vals
 
 
def read_array(fname):
  def keep(l):
    if l.startswith("#"):
      return False
    if l.startswith("@"):
      return False
    if not len(l.strip()):
      return False
    return True
  lines = list(filter(keep, open(fname, "r").readlines()))
  n_y = len(lines)
  n_x = len(lines[0].split())
  data = numpy.zeros((n_x, n_y), float)
  for j, line in enumerate(lines):
    for i, word in enumerate(line.split()):
      data[i, j] = float(word)
  return data
    

def read_column(fname, i=0):
  words = [l.split()[i] for l in open(fname, "r").readlines()]
  return numpy.array([float(w) for w in words])


def write_column(fname, y_vals):
  s = '\n'.join([str(val) for val in y_vals])
  open(fname, 'w').write(s)


def write_array(fname, data):
  n_bin = data.shape[0]
  lines = []
  for j in range(n_bin):
    strs = [str(val) for val in data[:,j]]
    line = ' '.join(strs)
    lines.append(line)
  open(fname, 'w').write('\n'.join(lines))


def flip_x_array(data):
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


def filter_array(array, filter_fn):
  n_x, n_y = array.shape
  for i in range(n_x):
    for j in range(n_y):
      if filter_fn(i, j):
        array[i,j] = 0.0


def filtered_array(array, filter_fn, save_fname=""):
  if save_fname and os.path.isfile(save_fname):
    return util.read_array(save_fname)
  result = array.copy()
  filter_array(result, filter_fn)
  if save_fname:
    util.write_array(save_fname, result)
  return result


def map_array(array, fn):
  n_x, n_y = array.shape
  result = numpy.zeros((n_x, n_y), float)
  for i in range(n_x):
    for j in range(n_y):
      result[i,j] = fn(array[i,j])
  return result


def get_numbered_dirs(in_dir):
  return [int(d) for d in glob.glob(os.path.join(in_dir, '*')) 
          if d.isdigit() and os.path.isdir(d)]


def get_numbered_dir_data(data_fname, i_column):
  res_dirs = [d for d in glob.glob('*') if d.isdigit()]
  test_fname = os.path.join(res_dirs[0], data_fname)
  n_res = len(read_array(test_fname)[0,:])
  result = numpy.zeros((n_res, n_res), float)
  for i in range(n_res):
    try:
      fname = os.path.join(str(i+1), data_fname)
      data = read_array(fname)
      result[i,:] = data[i_column,:].copy()
    except:
      result[i,:] = numpy.zeros((n_res), float)
  return result


def count_n_vals_per_column(res_data, min_val):
  def block_array(data, min_val):
    result = data.copy()
    n_x, n_y = numpy.shape(data)
    for i in range(n_x):
      for j in range(n_y):
        if result[i,j] > min_val:
          result[i,j] = 1
        else:
          result[i,j] = 0
    return result
  data = block_array(res_data, min_val)
  n_bin = data.shape[0]
  for i in range(n_bin):
    data[i,i] = 0
  result = numpy.zeros(n_bin, int)
  for i in range(n_bin):
    result[i] = sum(data[i,:])
  return result

