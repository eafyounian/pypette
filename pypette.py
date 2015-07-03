from __future__ import print_function
import sys, re, datetime, pwd, os, subprocess, string, random, gc, itertools
from collections import defaultdict



###########
# LOGGING #
###########

def info(message):
	print(message, file=sys.stderr)
	
def error(message):
	print(message, file=sys.stderr)
	sys.exit(-1)



#########################
# FILESYSTEM MANAGEMENT #
#########################

def mkdir(path):
	try:
		os.makedirs(path)
	except OSError:
		pass
	
def zopen(path, mode='r'):
	gzip_executable = 'gzip'
	#gzip_executable = 'pigz'
	
	if path[0] == '~':
		path = os.path.expanduser(path)
	
	if path == '-':
		if mode == 'r':
			return sys.stdin
		elif mode == 'w':
			return sys.stdout
	
	if path.lower().endswith('.gz'):
		if mode == 'r':
			return subprocess.Popen('gunzip -c %s' % path,
				stdout=subprocess.PIPE, shell=True).stdout
		elif mode == 'w':
			return subprocess.Popen('%s -c > %s' % (gzip_executable, path),
				stdin=subprocess.PIPE, shell=True).stdin
	else:
		return open(path, mode)
	
def open_exclusive(path):
	try:
		fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL)
	except:
		return None
	return os.fdopen(fd, 'w')



######################
# PROCESS MANAGEMENT #
######################

def shell(command, **kwargs):
	try:
		subprocess.check_call(command, shell=True, executable='/bin/bash',
			**kwargs)
	except subprocess.CalledProcessError as e:
		print('Process returned with error %d.' % e.returncode)
		return False
	return True

def shell_stdin(command):
	return subprocess.Popen(command, stdin=subprocess.PIPE, shell=True,
		executable='/bin/bash').stdin

def shell_stdout(command):
	return subprocess.Popen(command, stdout=subprocess.PIPE, shell=True,
		executable='/bin/bash', bufsize=-1).stdout

def shell_stdinout(command):
	p = subprocess.Popen(command, stdin=subprocess.PIPE,
		stdout=subprocess.PIPE, shell=True, executable='/bin/bash',
		bufsize=-1)
	return (p.stdin, p.stdout)

def daemonize(close_outputs=False):
	pid = os.fork()
	if pid < 0: error('Fork error in daemonize().')
	if pid > 0: sys.exit(0)     # Parent has PID > 0
	os.setsid()                 # Stop listening to signals for parent process
	for fd in range(0, 3 if close_outputs else 1):
		try: os.close(fd)
		except: pass



###########################
# MISCELLANEOUS FUNCTIONS #
###########################

class Object(object):
	def __init__(self, **kwargs): self.__dict__.update(kwargs)
	def __repr__(self): return repr(self.__dict__)


class prioritydict(dict):
	def __init__(self, comparator):
		self.comparator = comparator
	
	def __setitem__(self, key, value):
		old = self.setdefault(key, value)
		if old and self.comparator(value) > self.comparator(old):
			dict.__setitem__(self, key, value)

def argsort(seq):
	return sorted(range(len(seq)), key = seq.__getitem__)

def natural_sorted(seq):
    alphanum_key = lambda key: [ int(c) if c.isdigit() else c.lower()
    	for c in re.split('([0-9]+)', key) ] 
    return sorted(seq, key = alphanum_key)

def flatten(seq_of_seq):
	return list(itertools.chain(*seq_of_seq))



####################
# DNA/RNA SEQUENCE #
####################

def read_fasta(fasta_path):
	entries = {}
	fasta = open(fasta_path)
	chr = None
	for line in fasta:
		if line[0] == '#': continue
		if line[0] == '>':
			if chr:
				entries[chr] = ''.join(sequence)
				gc.collect()
			chr = line[1:].rstrip()
			if chr in entries: error('Entry %s found twice in FASTA.' % chr)
			sequence = []
		else:
			sequence.append(line.strip())

	if chr:
		entries[chr] = ''.join(sequence)

	return entries

def read_flat_seq(flat_seq_dir):
	entries = {}
	files = [f for f in os.listdir(flat_seq_dir) if f.endswith('.seq')]
	for file in files:
		entries[file.replace('.seq', '')] = open(
			os.path.join(flat_seq_dir, file)).read()
	return entries

revcomp_translate = str.maketrans('acgtACGT', 'tgcaTGCA')

def revcomplement(seq):
	return seq.translate(revcomp_translate)[::-1]
	
def point_region_distance(x, region):
	return max([0, region[0] - x, x - region[1]])

def regions_from_bed(bed_path):
	regions = []
	bed = zopen(bed_path)
	for line in bed:
		if line[0] == '#': continue
		c = line.rstrip().split('\t')
		regions.append((c[0], c[5], int(c[1])+1, int(c[2])))
	bed.close()
	return regions

def regions_to_bed(regions, bed_path):
	regions_file = open(bed_path, 'w')
	for region in regions:
		regions_file.write('%s\t%d\t%d\n' % (region[0], region[1]-1, region[2]))
	regions_file.close()

class GenomicFeatures:
	def __init__(self, features):
		self.features = features
		
	def __iter__(self):
		return iter(self.features)
		
	@classmethod
	def from_bed(cls, bed_path):
		features = []
		bed = open(bed_path)
		for line in bed:
			if line[0] == '#': continue
			
			tokens = line[:-1].split('\t')
			if len(tokens) == 4: tokens += ('', '+')
				
			# The 'chr' prefix is trimmed from chromosome names.
			features.append((tokens[0].replace('chr', ''), tokens[5], 
				int(tokens[1])+1, int(tokens[2]), tokens[3]))
		bed.close()
		return cls(features)
	
	def near(self, chr, pos, max_distance):
		chr = chr.replace('chr', '')
		nearby = [(f[3], point_region_distance(pos, f[2]))
			for f in self.features if f[0] == chr]
		nearby = [f for f in nearby if f[1] <= max_distance]
		return sorted(nearby, key=lambda x: x[1])
	
	def on_boundary(self, chr, pos):
		chr = chr.replace('chr', '')
		return [f for f in self.features if f[0] == chr and
			(f[2] == pos or f[3] == pos)]
	
		
	

	
	

