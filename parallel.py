#!/bin/env python

"""

A tool for running parallel jobs on a SLURM cluster. Reads a list of items
from standard input and runs the specified command on each item. Variable $x
contains the item name in each command invocation. The command must be given
inside single tick marks. The output of each job is stored in a log file 
under ~/.jobs/.

Usage:
  parallel [options] <command>

Examples:
  echo *.bam | parallel 'samtools sort $x ${x/.bam/} && samtools index $x'
  echo *_1.fq.gz | parallel -n8 -c8 -m10 \\
	'tophat2 tophat-indexes/hg19 $x ${x/_1.fq/_2.fq} -o ${x/_1.fq.gz/}'

Options:
  -h --help            Show this screen.
  -n --workers=N       Number of workers to run in parallel [default: 1].
  -c --cpus=N          How many CPUs to allocate for each job [default: 1].
  -m --memory=N        How much memory (in gigabytes) to allocate per job
                       [default: 5].
  -t --time=N          Time limit in hours for the full analysis [default: 72].
  -J --job-name=NAME   Job name shown by SLURM [default: job].
  -P --partition=NAME  SLURM partition on which to run jobs [default: serial].

"""

from __future__ import print_function
import subprocess, sys, re, docopt, socket, os, datetime, time
from pypette import shell, info, error, open_exclusive, daemonize

sbatch_template = '''#!/bin/bash -l
#SBATCH -p %s
#SBATCH -J %s
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH -n 1
#SBATCH -c %d
#SBATCH --mem-per-cpu=%d
#SBATCH -t %d
parallel worker %s
'''

def sanitize_path(path):
	path = path.replace('../', '')
	if path.startswith('./'): path = path[2:]
	return path.replace('/', '_')

def parallel(command, job_name, max_workers, cpus, memory, partition,
	time_limit):
	
	if partition == 'local' and socket.gethostname() == 'merope.local':
		error('ERROR: Running jobs on the Merope gateway is not allowed.')

	# SLURM uses environment variables to know when "srun" is being run
	# from inside a running job. By removing these environment variables,
	# we force SLURM to create an independent new job.
	if 'SLURM_JOB_ID' in os.environ:
		for key in os.environ.keys():
			if key.startswith('SLURM_'): del os.environ[key]
		info('Hiding SLURM environment variables...')
	
	# Allow splitting the command string onto multiple lines.
	command = command.replace('\n', ' ')
	
	if sys.stdin.isatty():
		# If the user did not provide any input, just run the command once.
		# The command must not contain $x.
		if '$x' in command or '${x' in command:
			error('Command contains $x but no targets provided.')
		
		targets = ['']
		info('Running job "%s" on %s partition (with %d %s and '
			'%d GB of memory).' % (job_name, partition, cpus,
			'CPUs' if cpus != 1 else 'CPU', memory))
	else:
		# Parse whitespace-delimited target items from standard input.
		targets = []
		for line in sys.stdin:
			targets += line.split(' ')
		targets = [t.replace('\n', '') for t in targets]
		
		if not targets: error('Command requires targets but none provided.')
	
		info('Distributing %d %s named "%s" on %s partition '
			'(with %d %s and %d GB of memory per job).' % (
			len(targets), 'jobs' if len(targets) != 1 else 'job',
			job_name, partition, cpus, 'CPUs' if cpus != 1 else 'CPU', memory))
	
	if max_workers > len(targets): max_workers = len(targets)
			
	script_dir = os.path.expanduser('~/.jobs/%s_%s' % (job_name, 
		datetime.datetime.now().strftime('%Y-%m-%d_%H:%M:%S')))
	os.makedirs(script_dir)
	scripts = ['%s/%s.queued' % (script_dir, sanitize_path(target))
		for target in targets]

	if len(set(scripts)) < len(scripts):
		error('Target list contains multiple instances of the following targets:\n' + '\n'.join(s for s in set(scripts)
			if scripts.count(s) > 1))
		
	for target, script in zip(targets, scripts):
		with open(script, 'w') as f:
			f.write('export x=%s; %s' % (target, command))
		
	if partition == 'local':
		worker_cmd = ['parallel', 'worker', script_dir]
		workers = [subprocess.Popen(worker_cmd) for w in range(max_workers)]
		for w in workers: w.wait()
	else:
		# Run the job steps on a SLURM cluster using sbatch.
		# Required memory is given in GB per job step. Convert to MB per CPU.
		mem_per_cpu = round(float(memory) / cpus * 1000)
		sbatch_script = sbatch_template % (partition, job_name, cpus,
			mem_per_cpu, 60 * time_limit, script_dir)
		workers = [subprocess.Popen(['sbatch', '-Q'], stdin=subprocess.PIPE)
			for p in range(max_workers)]
		for w in workers:
			w.stdin.write(sbatch_script)
			w.stdin.close()
		for w in workers: w.wait()

	

		
			




def parallel_worker(script_dir):
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
	# Users do not need to know about the "parallel worker" invocation. 
	if len(sys.argv) >= 3 and sys.argv[1] == 'worker':
		parallel_worker(script_dir=sys.argv[2])
		exit()
	
	args = docopt.docopt(__doc__)
	parallel(args['<command>'], job_name=args['--job-name'],
		max_workers=int(args['--workers']), cpus=int(args['--cpus']),
		memory=int(args['--memory']), partition=args['--partition'],
		time_limit=int(args['--time']))

