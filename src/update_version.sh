#!/bin/bash

if [ -z "$1" ]; then
	echo "Usage: $0 <version.h>"
	exit 1
fi

GIT_VER=$(git rev-parse --short HEAD)

if [ -e "$1" ]; then
	if (cat $1.in | sed "s,GIT_VER,${GIT_VER},g" | cmp -s - $1); then
		exit 0
	fi
fi

cat $1.in | sed "s,GIT_VER,${GIT_VER},g" > $1
