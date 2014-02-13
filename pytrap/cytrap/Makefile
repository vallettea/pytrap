CC=gcc
CFLAGS=-Wall -Wextra -O3 -DHAVE_INLINE

LDFLAGS+=$(shell pkg-config --libs gsl) -s

all: solveode

solveode: solveode.o conformal_map.o