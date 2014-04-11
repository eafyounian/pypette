#!/bin/env python

"""

A tool for running parallel jobs on the Merope cluster. Reads a list of items
from standard input and runs the specified command on each item. Variable $x
contains the item name in each command invocation. The command must be given
inside single tick marks. The output of each job is stored in a log file 
under ~/.jobs/.

Usage:
  parallel [options] <command>

Examples:
  echo *.bam | parallel 'samtools sort $x ${x%.bam} && samtools index $x'
  echo *_1.fq.gz | parallel -n8 -c8 -m10 \\
	'tophat2 tophat-indexes/hg19 $x ${x%_1.fq.gz}_2.fq.gz -o ${x%_1.fq.gz}'

Options:
  -h --help            Show this screen.
  -w --wait            Do not daemonize, wait for jobs to finish.
  -n --jobs=N          Maximum number of jobs to run in parallel [default: 4].
  -c --cpus=N          How many CPUs to allocate for each job [default: 1].
  -m --memory=N        How much memory (in gigabytes) to allocate per job
                       [default: 5].
  -t --time=N          Time limit in hours for the full analysis [default: 168].
  -J --job-name=NAME   Job name shown by SLURM [default: job].
  -P --partition=NAME  SLURM partition on which to run jobs [default: local].

Author: Matti Annala <matti.annala@tut.fi>

"""

from __future__ import print_function
import subprocess, sys, re, docopt, socket, os, datetime, time
from pypette import shell, daemonize, info, error, open_exclusive



def sanitize_path(path):
	path = path.replace('../', '')
	if path.startswith('./'): path = path[2:]
	return path.replace('/', '_')

def parallel(command, job_name, max_parallel, cpus, memory, partition,
	time_limit, wait_to_finish):
	
	if partition == 'local' and socket.gethostname() == 'merope.local':
		error('ERROR: Running jobs on the Merope gateway is not allowed.')
	
	# Allow splitting the command string onto multiple lines.
	command = command.replace('\n', '')
	
	if sys.stdin.isatty():
		# If the user did not provide any input, just run the command once.
		# The command must not contain $x.
		if '$x' in command or '${x' in command:
			error('Command requires targets but none provided.')
		
		tokens = ['']
		info('Running job "%s" on %s partition (with %d %s and '
			'%d GB of memory).' % (job_name, partition, cpus,
			'CPUs' if cpus != 1 else 'CPU', memory))
	else:
		# Parse whitespace-delimited target items from standard input.
		tokens = []
		for line in sys.stdin:
			tokens += line.split(' ')
		tokens = [t.replace('\n', '') for t in tokens]
		
		if not tokens: error('Command requires targets but none provided.')
	
		info('Distributing %d %s named "%s" on %s partition '
			'(with %d %s and %d GB of memory per job).' % (
			len(tokens), 'jobs' if len(tokens) != 1 else 'job',
			job_name, partition, cpus, 'CPUs' if cpus != 1 else 'CPU', memory))
	
	if max_parallel > len(tokens):
		max_parallel = len(tokens)
			
	script_dir = os.path.expanduser('~/.jobs/%s_%s' % (job_name, 
		datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')))
	os.makedirs(script_dir)
	scripts = ['%s/%s.queued' % (script_dir, sanitize_path(token))
		for token in tokens]

	if len(set(scripts)) < len(scripts):
		for s in scripts:
			if sum(s == x for x in scripts) > 1: print(s)
		error('Input list contains multiple instances of some items.')
		
	for token, script in zip(tokens, scripts):
		with open(script, 'w') as f:
			f.write('export x=%s; %s' % (token, command))
			
	if not wait_to_finish: daemonize()
	
	worker_cmd = ['parallel', 'nested', script_dir]
	if partition == 'local':
		workers = [subprocess.Popen(worker_cmd) for p in range(max_parallel)]
	else:
		# User wants to run the job steps on a SLURM cluster. Allocate resources
		# using salloc and then run the commands using srun.
		# Required memory is given in GB per job step. Convert to MB per CPU.
		mem_per_cpu = round(float(memory) / cpus * 1000)
		flags = '--partition=%s --job-name=%s --ntasks=1 ' \
			'--cpus-per-task=%d --mem-per-cpu=%d --time=%d' % \
			(partition, job_name, cpus, mem_per_cpu, 60 * time_limit)
		
		workers = [subprocess.Popen(['srun', '-Q'] + flags.split() + worker_cmd)
			for p in range(max_parallel)]
			
	while workers:
		workers = [w for w in workers if w.poll() == None]
		time.sleep(1)
		
			




def nested_parallel(script_dir):
	scripts = [os.path.join(script_dir, f) for f in os.listdir(script_dir)
		if f.endswith('.queued')]
	out_files = [re.sub('.queued$', '.out', f) for f in scripts]
		
	for script, out_file in zip(scripts, out_files):
		out = open_exclusive(out_file)
		if not out: continue
		
		with open(script) as f: command = f.read()
		os.remove(script)
		out.write('%s\n%s\n' % (command, '-'*80))
		out.flush()
		start_time = datetime.datetime.now()
		shell(command, stdout=out, stderr=out)
		end_time = datetime.datetime.now()
		out.write('%s\nJOB FINISHED. ELAPSED TIME WAS %s.\n' %
			('-'*80, end_time - start_time)) 
		out.close()
		
		






#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	# The "parallel nested" invocation is hidden from the user.
	if len(sys.argv) >= 3 and sys.argv[1] == 'nested':
		nested_parallel(script_dir=sys.argv[2])
		exit()
	
	args = docopt.docopt(__doc__)
	parallel(args['<command>'], job_name=args['--job-name'],
		max_parallel=int(args['--jobs']), cpus=int(args['--cpus']),
		memory=int(args['--memory']), partition=args['--partition'],
		time_limit=int(args['--time']), wait_to_finish=args['--wait'])

