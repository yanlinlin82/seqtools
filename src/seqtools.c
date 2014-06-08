#include <stdio.h>
#include <string.h>
#include "version.h"
#include "common.h"
#include "readqc.h"

int verbose = 0;

struct command_item {
	const char *cmd;
	int (*func)(int argc, char * const *argv);
	const char *desc;
};

static struct command_item commands[] = {
	{ "readqc", readqc_main, "Check quality of sequencing reads" },
	{ },
};

static void print_usage(void)
{
	struct command_item *p;

	fprintf(stderr, "\n"
			"Program: seqtools (Tools for manipulating high throughput sequencing data)\n"
			"Version: "VERSION"\n"
			"Author : Linlin Yan (yanll<at>mail.cbi.pku.edu.cn)\n"
			"Copyright: 2014, Centre for Bioinformatics, Peking University, China\n"
			"\n"
			"Usage: seqtools <command> [options]\n"
			"\n"
			"Commands:\n");

	for (p = commands; p->cmd; ++p) {
		fprintf(stderr, "   %s    %s\n", p->cmd, p->desc);
	}
	fprintf(stderr, "\n");
}

int main(int argc, char * const *argv)
{
	struct command_item *p;

	if (argc < 2) {
		print_usage();
		return 1;
	}

	for (p = commands; p->cmd; ++p) {
		if (strcmp(p->cmd, argv[1]) == 0) {
			return (*p->func)(argc - 1, argv + 1);
		}
	}
	fprintf(stderr, "[%s] Unknown command '%s'!\n", __FUNCTION__, argv[1]);
	return 1;
}
