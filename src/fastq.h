#ifndef __FASTQ_H__
#define __FASTQ_H__

#include <stdbool.h>
#include <zlib.h>

struct seq {
	char *name;
	char *base;
	char *qual;
};

struct seq_pair {
	struct seq r1;
	struct seq r2;
};

struct cache {
	size_t capacity;
	size_t size;
	struct seq *data;
	size_t skip;
};

struct fastq {
	gzFile fp1;
	gzFile fp2;
	int base;
	const char *file1;
	const char *file2;
	long long lineno;
	struct cache cache;
	bool error;
};

struct fastq *fastq_open(const char *file1, const char *file2, int base);
bool fastq_has_error(struct fastq *fq);
bool fastq_read(struct fastq *fq, struct seq_pair *pair);
void fastq_close(struct fastq *fq);

#endif /* __FASTQ_H__ */
