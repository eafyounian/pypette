
all:
	mkdir -p bin
	ln -sf ../breakfast.py bin/breakfast
	ln -sf ../fasta.py bin/fasta
	ln -sf ../gtf.py bin/gtf
	ln -sf ../sam.py bin/sam
	ln -sf ../variant.py bin/variant
	ln -sf ../coverage.py bin/coverage
	ln -sf ../parallel.py bin/parallel
	ln -sf ../expression.py bin/expression
	ln -sf ../backup.py bin/backup
	ln -sf ../cghub.py bin/cghub
	ln -sf ../swiss.py bin/swiss
	gcc -O3 -std=c99 -o compiled/spileup compiled/spileup.c

