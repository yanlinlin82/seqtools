CC      := gcc
CFLAGS  := -O2 -Wall

ifneq ("${MAKECMDGOALS}", "clean")
VERSION := $(shell ./src/update_version.sh ./src/version.h)
endif

TARGET  := seqtools ${VERSION}
MODULES := seqtools readqc

all: ${TARGET}

clean:
	@rm -vrf ${TARGET} tmp/ src/version.h

seqtools: ${MODULES:%=tmp/%.o}
	${CC} ${CFLAGS} -o $@ $^

tmp/%.o: src/%.c
	${CC} ${CFLAGS} -c -o $@ $<

ifneq ("${MAKECMDGOALS}", "clean")
sinclude ${MODULES:%=tmp/%.d}
tmp/%.d: src/%.c
	@echo "Parsing dependency for '$<'"
	@[ -d ${@D} ] || mkdir -p ${@D}
	@${CC} -MM $< -MT ${@:%.d=%.o} | sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' > $@
endif
