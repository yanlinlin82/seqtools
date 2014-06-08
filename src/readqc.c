#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "fastq.h"
#include "readqc.h"

#define LOW_BASE 0.5
#define LOW_QUAL 5
#define MAX_CHECK_PAIRS 1000000
#define MIN_CHECK_PAIRS 10000

static float low_base = LOW_BASE;
static int low_qual = LOW_QUAL;
static long long max_check_pairs = MAX_CHECK_PAIRS;

static inline bool is_read_low_qual(int low_qual_base_count, int base_count)
{
	if (low_base >= 1) {
		return (low_qual_base_count > low_base);
	} else {
		return (low_qual_base_count > (int)(base_count * low_base));
	}
}

static int qc_check(const char *file1, const char *file2, int base)
{
	struct fastq *fq;
	struct seq_pair pair;
	long long read_count = 0, low_qual_read_count = 0;
	int i, base_count, low_qual_base_count;

	if (!(fq = fastq_open(file1, file2, base))) {
		return 1;
	}
	while (read_count < max_check_pairs && fastq_read(fq, &pair)) {
		++read_count;

		base_count = 0;
		low_qual_base_count = 0;
		for (i = 0; pair.r1.qual[i]; ++i) {
			++base_count;
			low_qual_base_count += ((pair.r1.qual[i] - fq->base) <= low_qual);
		}
		if (is_read_low_qual(low_qual_base_count, base_count)) {
			++low_qual_read_count;
			continue;
		}

		base_count = 0;
		low_qual_base_count = 0;
		for (i = 0; pair.r2.qual[i]; ++i) {
			++base_count;
			low_qual_base_count += ((pair.r2.qual[i] - fq->base) <= low_qual);
		}
		if (is_read_low_qual(low_qual_base_count, base_count)) {
			++low_qual_read_count;
			continue;
		}
	}
	if (fastq_has_error(fq)) {
		fastq_close(fq);
		return 1;
	}
	fprintf(stdout, "%lld\t%lld\n", low_qual_read_count, read_count);
	fastq_close(fq);
	return 0;
}

static void print_usage(const char *cmd)
{
	fprintf(stderr, "\n"
			"Usage: seqtools %s [options] <1.fq> <2.fq>\n"
			"\n"
			"Options:\n"
			"   -v         Show verbose messages\n"
			"   -q <INT>   Quality base value, 33 or 64, default to auto-detect\n"
			"   -L <INT>   Maximum low base quality value, default: %d\n"
			"   -n <NUM>   Minimum bad bases for a low quality read, default: %.1f\n"
			"   -N <INT>   Maximum pairs to check, default: %lld\n"
			"\n", cmd, low_qual, low_base, max_check_pairs);
}

int readqc_main(int argc, char * const *argv)
{
	int c, base = 0;

	while ((c = getopt(argc, argv, "hvq:L:n:N:")) != -1) {
		switch (c) {
		case 'h':
			print_usage(argv[0]);
			return 1;
		case 'v':
			++verbose;
			break;
		case 'q':
			base = atoi(optarg);
			if (base != 33 && base != 64) {
				fprintf(stderr, "Qualtiy base value (-q) should be 33 or 64!\n");
				return 1;
			}
			break;
		case 'L':
			low_qual = atoi(optarg);
			break;
		case 'n':
			low_base = atof(optarg);
			break;
		case 'N':
			max_check_pairs = strtoll(optarg, NULL, 10);
			if (max_check_pairs < MIN_CHECK_PAIRS) {
				fprintf(stderr, "Maximum pairs to check (-N) should be at least %d!\n", MIN_CHECK_PAIRS);
				return 1;
			}
			break;
		default:
			fprintf(stderr, "Unknown option '-%c'!\n", c);
			return 1;
		}
	}
	if (optind + 2 != argc) {
		print_usage(argv[0]);
		return 1;
	}

	return qc_check(argv[optind], argv[optind + 1], base);
}
