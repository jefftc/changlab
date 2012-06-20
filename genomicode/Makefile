# This is a semi-generic makefile for compiling Python shared
# libraries.
# Created 5/6/2002, jchang

# This should be SOLARIS, IRIX, OSX, OSX4, or LINUX.
#PLATFORM=LINUX
PLATFORM=OSX

# Set this to point to the include files for your current installation
# of python.
#PYINCLUDE=/usr/local/include/python2.7
#PYINCLUDE=/home/jefftc/include/python2.5
#PYINCLUDE=/opt/local/include/python2.4
#PYINCLUDE=/opt/local/Library/Frameworks/Python.framework/Versions/2.4/include/python2.4
PYINCLUDE=/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7

# Use this to pass any other flags to the C compiler.  This should be
# at least -I$(PYINCLUDE).
MORECFLAGS=-I$(PYINCLUDE)

# Use this to pass any other flags to the linker.
MORELDFLAGS=

# Specify the files to make.
all: ciolib.so cMatrix.so cjmath.so cgeolib.so ccoherencelib.so cGeneSet.so
.PHONY: all

cjmath.so: cjmathmodule.o mathlib.o

ccoherencelib.so: ccoherencelibmodule.o mathlib.o

mathlib.o: mathlib.c
	${CC} -c -Wall -fPIC -march=nocona -O2 ${MORECFLAGS} $^ -o $@

#cMatrixmodule.o: cMatrixmodule.c
#	${CC} -c -Wall -fPIC -DHASH_FUNCTION=HASH_OAT -O2 ${MORECFLAGS} $^ -o $@
#	${CC} -c -Wall -fPIC -DHASH_BLOOM=24 -O2 ${MORECFLAGS} $^ -o $@


# You should not have to change things past here, unless you want to
# tweak flags for your platform.

# Set platform-specific values for these variables.
SOLARIS_CC=gcc
SOLARIS_LD=ld
SOLARIS_CFLAGS=-O2 -Wall
SOLARIS_LDFLAGS=-G

IRIX_CC=gcc
IRIX_LD=ld
IRIX_CFLAGS=-O2 -Wall -g -shared
IRIX_LDFLAGS=-shared -all

# Works for OS X 10.4 and below.
OSX4_CC=cc
OSX4_LD=cc
OSX4_CFLAGS=-O2 -Wall -g -traditional-cpp
OSX4_LDFLAGS=-bundle -undefined suppress -flat_namespace

OSX_CC=cc
OSX_LD=cc
OSX_CFLAGS=-O3 -Wall -g -dynamic
OSX_LDFLAGS=-bundle -undefined dynamic_lookup

LINUX_CC=gcc
LINUX_LD=ld
LINUX_CFLAGS=-O2 -Wall -fPIC
LINUX_LDFLAGS=-G

CC=$($(PLATFORM)_CC)
LD=$($(PLATFORM)_LD)
CFLAGS=$($(PLATFORM)_CFLAGS)
LDFLAGS=$($(PLATFORM)_LDFLAGS)

# Export these variables to sub-makes.
export PYINCLUDE
export MORECFLAGS
export MORELDFLAGS
export CC
export LD
export CFLAGS
export LDFLAGS

# Implicit in make, except for the MORECFLAGS thing.
%.o: %.c
	${CC} -c ${CFLAGS} ${MORECFLAGS} $^ -o $@

%.so: %module.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $^ -o $@

%.so: %.o
	${LD} ${LDFLAGS} ${MORELDFLAGS} $^ -o $@

.PHONY : clean
clean:
	@find . -name '*.py[oc]' -print -exec rm -f {} \;
	@find . -name '*~' -print -exec rm -f {} \;
	@find . -name '*.so' -print -exec rm -f {} \;
	@find . -name '*.o' -print -exec rm -f {} \;
	@rm -f *_wrap.c
#	$(MAKE) clean -C abbreviation
