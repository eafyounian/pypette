
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ALLELE_LEN 32
#define MAX_ALLELES 100
#define MAX_SAMPLES 1000

typedef struct {
	char seq[MAX_ALLELE_LEN];
	int high_count;			// High quality reads
	int low_count;			// Low quality reads
} Allele;

typedef struct {
	Allele alleles[MAX_ALLELES];
	int total;
} AlleleList;

int min_mapq = 0;

void count_allele(char* allele, int len, char quality, AlleleList* list) {
	quality = quality - 33;

	// Zero terminate temporarily.
	char ochar = allele[len]; allele[len] = 0;
	
	for (int i = 0; i < list->total; i++) {
		if (strcmp(allele, list->alleles[i].seq) == 0) {
			if (quality >= min_mapq)
				list->alleles[i].high_count += 1;
			else
				list->alleles[i].low_count += 1;
			goto cleanup;
		}
	}
	
	// Did not find our allele. Add a new one.
	list->alleles[list->total].high_count = (quality >= min_mapq);
	list->alleles[list->total].low_count = (quality < min_mapq);
	strncpy(list->alleles[list->total].seq, allele, MAX_ALLELE_LEN); 
	list->total += 1;
	
cleanup:
	allele[len] = ochar;   // Undo zero termination
}

void print_alleles(AlleleList* list) {
	if (list->total == 0) return;
	int i;
	for (i = 0; i < list->total - 1; i++) {
		printf("%s %d %d ", list->alleles[i].seq, list->alleles[i].high_count,
			list->alleles[i].low_count);
	}
	printf("%s %d %d", list->alleles[i].seq, list->alleles[i].high_count,
		list->alleles[i].low_count);
}

int parse_pileup(char* bases, char* quality, AlleleList* alleles) {
	alleles->total = 0;
	
	//printf("Bases: %s\n", bases);
	//printf("Quality: %s\n", quality);
	
	// First we need to convert everything into lowercase.
	for (int i = 0; bases[i]; i++) {
		if (bases[i] == ',') {
			bases[i] = '.';
		} else {
			bases[i] = tolower(bases[i]);
		}
	}
	
	int j = 0;
	for (int i = 0; bases[i]; i++) {
		if (bases[i] == '^') {
			i++;
		} else if (bases[i] == '*' || bases[i] == '$' || bases[i] == '<' ||
			bases[i] == '>') {
			// Ignore.
		} else if (bases[i+1] == '+' || bases[i+1] == '-') {
			char* allele = 0;  // First char after number stored here.
			int len = strtol(bases + i + 2, &allele, 10);
			allele[-1] = bases[i+1]; allele[-2] = bases[i];
			allele -= 2; len += 2;  // Include pre-base and sign
			count_allele(allele, len, quality[j++], alleles);
			i += len;
		} else {
			count_allele(bases + i, 1, quality[j++], alleles);
		}
	}
}

int main(int argc, char** argv) {
	if (argc != 3) {
		printf("Usage: spileup <min_alt_reads> <min_mapq>\n");
		return -1;
	}

	// Parse command line arguments.
	int min_alt_reads = strtol(argv[1], NULL, 10);
	min_mapq = strtol(argv[2], NULL, 10);
	
	int bufsize = 1000 * 1000;
	char line[bufsize];
	
	int tabs[2*MAX_SAMPLES];    // Tab character locations
	
	AlleleList sample_alleles[MAX_SAMPLES];
	
	while (fgets(line, bufsize, stdin)) {
		// Find all tab characters.
		int ntabs = 0;
		for (int i = 0; line[i]; i++) {
			if (line[i] == '\t') {
				tabs[ntabs++] = i;
				if (ntabs >= 3) line[i] = '\0';
			}
		}
		if (ntabs < 3) continue;

		int S = 0, i = 2;
		while (i < ntabs) {
			if (line[tabs[i]+1] == '0') { i += 3; }
			else {
				parse_pileup(line + tabs[i+1] + 1, line + tabs[i+3] + 1,
					&sample_alleles[S]);
				i += 4;
			}
			S += 1;
		}
		
		// If all samples only show reference alleles, don't print anything.
		int best_alt_reads = 0;
		for (int s = 0; s < S; s++) {
			AlleleList* alleles = &sample_alleles[s];
			int alt_alleles = 0;
			for (int a = 0; a < alleles->total; a++) {
				if (alleles->alleles[a].seq[0] == '.') continue;
				if (alleles->alleles[a].high_count > best_alt_reads)
					best_alt_reads = alleles->alleles[a].high_count;
			}
		}
		
		if (best_alt_reads >= min_alt_reads) {
			fwrite(line, 1, tabs[2], stdout);  // First three columns.
			for (int s = 0; s < S; s++) {
				printf("\t");
				print_alleles(&sample_alleles[s]);
			}
			printf("\n");
			fflush(stdout);
		}

	}
	return 0;
}


