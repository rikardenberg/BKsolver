# Position of sources and libraries
VPATH= ../src 
GSL_HOME = /usr/local
GSLLIB   = -L$(GSL_HOME)/lib -lgsl -lgslcblas -lm
GSLINC   = -I$(GSL_HOME)/include


# Files to be compiled
OBJECTS  = bkeqn.o bksolve.o bksetup.o fft4g.o parsecfg.o
SATSCALE = qs.o


# Compiler and choice of compilation flags
CC = gcc
CFLAGS = -O3 -g -Wall -W -Wshadow

OPTFLAGS = -O3 -g -Wall -W -Wshadow \
-march=i686 -malign-double -fshort-enums -fomit-frame-pointer

BIGFLAGS = -ansi -Wall -W \
-Wmissing-prototypes -Wstrict-prototypes \
-Wmissing-declarations -Wredundant-decls \
-Wconversion -Wshadow \
-Wpointer-arith -Wcast-qual -Wcast-align \
-Wnested-externs -Wfloat-equal \
-Wpointer-arith \
-fshort-enums -fno-common -Dinline= -g -O3

FLAGG = $(CFLAGS)



all: bk

bk: $(OBJECTS) 
	$(CC) $(FLAGG) -o bk $(OBJECTS) $(GSLLIB)
	
qs: $(SATSCALE) 
	$(CC) $(FLAGG) -o qs $(SATSCALE) $(GSLLIB) 
	
.c.o:
	$(CC) -c $(FLAGG) $(GSLINC) $<

clean:
	rm -f *.BAK *.bak core *~ *# *.o

