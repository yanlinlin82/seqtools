#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "common.h"
#include "fastq.h"

#define MAX_CHAR_IN_LINE 1024

#define QBASE33_MIN 33
#define QBASE33_MAX 126
#define QBASE64_MIN 59
#define QBASE64_MAX 126

static int check_qual_base(const char *buf)
{
	int i, min_qual, max_qual;

	for (i = 1, min_qual = max_qual = (unsigned char)buf[0]; buf[i] && buf[i] != '\n'; ++i) {
		if (min_qual > (unsigned char)buf[i]) {
			min_qual = (unsigned char)buf[i];
		}
		if (max_qual < (unsigned char)buf[i]) {
			max_qual = (unsigned char)buf[i];
		}
	}
	if (min_qual < QBASE33_MIN) {
		fprintf(stderr, "Error: Invalid base quality value: %d (in ascii)\n", min_qual);
		return -1;
	}
	if (max_qual > QBASE64_MAX) {
		fprintf(stderr, "Error: Invalid base quality value: %d (in ascii)\n", max_qual);
		return -1;
	}
	if (min_qual < QBASE64_MIN && max_qual > QBASE33_MAX) {
		fprintf(stderr, "Error: Invalid base quality range: %d - %d (in ascii)\n", min_qual, max_qual);
		return -1;
	}
	if (max_qual > QBASE33_MAX) {
		return 64; /* 64-based: Solexa, Illumina 1.3+, Illumina 1.5+ */
	} else if (min_qual < QBASE64_MIN) {
		return 33; /* 33-based: Sanger, Illumina 1.8+ */
	} else {
		return 0; /* Unsured */
	}
}

static bool open_and_check_fastq(gzFile *fp, const char *file)
{
	int c;

	if (!(*fp = gzopen(file, "rb"))) {
		fprintf(stderr, "Error: Can not open read file '%s'!\n", file);
		return false;
	}
	if ((c = gzgetc(*fp)) != '@') {
		fprintf(stderr, "Error: First read file '%s' is not in FASTQ format!\n", file);
		gzclose(*fp);
		return false;
	}
	gzungetc(c, *fp);
	return true;
}

#define MIN_INCREMENT 16
#define MAX_INCREMENT 256

static bool cache_reserve(struct cache *cache, size_t size)
{
	if (cache->capacity < size) {
		struct seq *p;
		size_t capacity = cache->capacity;
		while (capacity < size) {
			if (capacity < MIN_INCREMENT) {
				capacity += MIN_INCREMENT;
			} else if (capacity > MAX_INCREMENT) {
				capacity += MAX_INCREMENT;
			} else {
				capacity += capacity;
			}
		}
		p = malloc(sizeof(struct seq) * capacity);
		if (!p) {
			fprintf(stderr, "Error: Allocate memory for FASTQ cache failed!\n");
			return false;
		}
		if (cache->size > 0) {
			memcpy(p, cache->data, sizeof(struct seq) * cache->size);
		}
		free(cache->data);
		cache->data = p;
		cache->capacity = capacity;
	}
	return true;
}

static bool cache_append(struct cache *cache, const struct seq *seq)
{
	if (!cache_reserve(cache, cache->size + 1)) {
		return false;
	}
	memcpy(cache->data + cache->size, seq, sizeof(struct seq));
	++cache->size;
	return true;
}

static void cache_free(struct cache *cache)
{
	free(cache->data);
	cache->data = NULL;
	cache->size = 0;
	cache->capacity = 0;
}

bool load_one_read(gzFile fp, const char *file, long long line, struct seq *seq, int *base, bool *error)
{
	char buf[MAX_CHAR_IN_LINE];
	int i;
	size_t len;

	for (i = 0; i < 4; ++i) {
		if (!gzgets(fp, buf, sizeof(buf))) {
			return false;
		}
		len = strlen(buf);
		if (len > 0 && buf[len - 1] == '\n') {
			buf[--len] = '\0';
		}
		if (i == 0) {
			if (buf[0] != '@') {
				fprintf(stderr, "Error: '@' is expected in line %lld of file '%s'!\n",
						line + 1, file);
				*error = true;
				return false;
			}
			if (len <= 1) {
				fprintf(stderr, "Error: Unexpected empty read name in line %lld of file '%s'!\n",
						line + 1, file);
				*error = true;
				return false;
			}
			seq->name = strdup(buf);
		} else if (i == 1) {
			if (len <= 0) {
				fprintf(stderr, "Error: Unexpected empty sequence in line %lld of file '%s'!\n",
						line + 1, file);
				*error = true;
				return false;
			}
			seq->base = strdup(buf);
		} else if (i == 2) {
			if (buf[0] != '+') {
				fprintf(stderr, "Error: '+' is expected in line %lld of file '%s'!\n",
						line + 3, file);
				*error = true;
				return false;
			}
		} else {
			seq->qual = strdup(buf);
			if (strlen(seq->base) != strlen(seq->qual)) {
			}
			if (*base == 0) {
				*base = check_qual_base(seq->qual);
				if (*base < 0) {
					*error = true;
					return false;
				}
				if (verbose) {
					fprintf(stderr, "Qualit base value detected: %d\n", *base);
				}
			} else {
				int new_base = check_qual_base(seq->qual);
				if (new_base < 0) {
					*error = true;
					return false;
				} else if (new_base > 0 && *base != new_base) {
					fprintf(stderr, "Error: Unexpected base quality!\n");
					*error = true;
					return false;
				}
			}
		}
	}
	return true;
}

struct fastq *fastq_open(const char *file1, const char *file2, int base)
{
	struct fastq *fq = malloc(sizeof(struct fastq));
	if (!fq) {
		fprintf(stderr, "Error: Allocate memory for FASTQ loader failed!\n");
		return NULL;
	}
	memset(fq, 0, sizeof(struct fastq));

	if (!open_and_check_fastq(&fq->fp1, file1)) {
		free(fq);
		return NULL;
	}
	if (!open_and_check_fastq(&fq->fp2, file2)) {
		gzclose(fq->fp1);
		free(fq);
		return NULL;
	}

	fq->file1 = file1;
	fq->file2 = file2;
	fq->base = base;

	while (fq->base == 0) {
		struct seq r;
		if (!load_one_read(fq->fp1, file1, fq->lineno, &r, &fq->base, &fq->error)) {
			goto failed;
		}
		cache_append(&fq->cache, &r);

		if (!load_one_read(fq->fp2, file2, fq->lineno, &r, &fq->base, &fq->error)) {
			goto failed;
		}
		cache_append(&fq->cache, &r);
	}
	return fq;
failed:
	cache_free(&fq->cache);
	gzclose(fq->fp2);
	gzclose(fq->fp1);
	free(fq);
	return NULL;
}

bool fastq_has_error(struct fastq *fq)
{
	return fq->error;
}

bool fastq_read(struct fastq *fq, struct seq_pair *pair)
{
	if (fq->cache.skip < fq->cache.size) {
		assert(fq->cache.skip + 1 < fq->cache.size);
		memcpy(&pair->r1, fq->cache.data + fq->cache.skip, sizeof(struct seq));
		memcpy(&pair->r2, fq->cache.data + fq->cache.skip + 1, sizeof(struct seq));
		fq->cache.skip += 2;
	} else {
		if (!load_one_read(fq->fp1, fq->file1, fq->lineno, &pair->r1, &fq->base, &fq->error)) {
			return false;
		}
		if (!load_one_read(fq->fp2, fq->file2, fq->lineno, &pair->r2, &fq->base, &fq->error)) {
			fq->error = true;
			return false;
		}
	}
	return true;
}

void fastq_close(struct fastq *fq)
{
	if (fq) {
		cache_free(&fq->cache);
		gzclose(fq->fp2);
		gzclose(fq->fp1);
		free(fq);
	}
}
