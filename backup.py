#!/bin/env python

"""

A tool for performing incremental backups of high data volumes.

Usage:
  backup <rules_file>

Options:
  -h --help       Show this screen.

Author: Matti Annala <matti.annala@tut.fi>

"""

from __future__ import print_function
import subprocess, sys, re, docopt, os, getpass
from pypette import shell, shell_stdout, error, Object

def backup(rules_path):
	passwords = {}
	rules = []
	for line in open(rules_path):
		line = line.strip()
		if not line or line[0] == '#': continue
		tokens = line.strip().split()
		if len(tokens) != 2: error('Invalid rule: "%s"' % line)

		if not ':' in tokens[1]: error('Missing host: "%s"' % line)
		host, path = tokens[1].split(':')

		username = getpass.getuser()
		if '@' in host:
			username, host = host.split('@')

		if not os.path.isdir(tokens[0]):
			print('Directory %s does not exist. Ignoring rule...' % tokens[0])
			continue

		rule = Object()
		rule.src_dir = tokens[0]
		rule.dst_host = host
		rule.dst_dir = path
		rule.username = username
		rule.password = passwords[host] if host in passwords else \
			getpass.getpass('Password for %s: ' % host)
		passwords[host] = rule.password
		rules.append(rule)

	def lftp_mirror(rule, dry_run=False):

		cmds = open('.lftp_script', 'w')
		cmds.write('open -u %s,%s sftp://%s\n' % (
			rule.username, rule.password, rule.dst_host))
		cmds.write('mirror -Rae %s %s %s\n' % (
			'--dry-run' if dry_run else '-v', rule.src_dir, rule.dst_dir))
		cmds.close()

		if dry_run:
			userpass = rule.username + ':' + rule.password + '@'
			host = rule.dst_host

			out = shell_stdout('lftp -f .lftp_script')
			for line in out:
				if line.startswith('chmod'): continue
				if line.startswith('mkdir'): continue

				m = re.match('get -O sftp://(.+) file:/.+/(.+)', line)
				if m:
					dst = m.group(1)
					if dst.startswith(userpass): dst = dst[len(userpass):]
					if dst.startswith(host): dst = dst[len(host):]
					print('ADD %s/%s' % (dst, m.group(2)))
					continue

				m = re.match('get -e -O sftp://(.+) file:/.+/(.+)', line)
				if m:
					dst = m.group(1)
					if dst.startswith(userpass): dst = dst[len(userpass):]
					if dst.startswith(host): dst = dst[len(host):]
					print('UPDATE %s/%s' % (dst, m.group(2)))
					continue

				m = re.match('rm .*sftp://(.+)', line)
				if m:
					dst = m.group(1)
					if dst.startswith(userpass): dst = dst[len(userpass):]
					if dst.startswith(host): dst = dst[len(host):]
					print('DELETE %s' % dst)
					continue

				sys.stdout.write(line)
		else:
			shell('lftp -f .lftp_script')

		os.remove('.lftp_script')
	
	for rule in rules:
		lftp_mirror(rule, dry_run=True)
			
	if not raw_input('Proceed with backup? [y/N] ').lower() in ('y', 'yes'):
		error('Backup canceled.')
	
	for rule in rules:
		lftp_mirror(rule)






#######################
# COMMAND LINE PARSER #
#######################

if __name__ == '__main__':
	args = docopt.docopt(__doc__)
	backup(args['<rules_file>'])

