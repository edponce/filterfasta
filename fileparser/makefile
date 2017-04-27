# makefile for filterfasta program.
# CC - compiler to use
CC=mpicc

# LDFLAGS - linker options
LDFLAGS=

# CFLAGS - compiler options
# -Wall = turn on most warnings
# -Wextra = warn about type limits
# -g = compiles with debug info
# -pg = create gmon.out for gprof
CFLAGS=-Wall -Wextra -g -O3 -DBCAST_INFILES -DBCAST_OUTFILES

# INCLUDES - define directories containing header files in addition to /usr/include
# Example: INCLUDES=-I/dir1 -I/dir2
INCLUDES=

# LFLAGS - define library paths in addition to /usr/lib
# Example: LFLAGS=-L/dir1 -L/dir2
#LFLAGS=-L/sw/beacon/fpmpi/1.1/centos6.2_icc2013.1.117/lib
LFLAGS=

# LIBS - define libraries to link into executable
# -lm = math library
#LIBS=-lm -lfpmpi -lirc
LIBS=-lm

# SOURCES - source files to compile
SOURCES=src/mpifilterfasta.c

# OBJECTS - object files to link
# This uses a suffix replacement rule within a macro:
#  $(name:string1=string2)
#	For each word in 'name' replace 'string1' with 'string2'
OBJECTS=$(SOURCES:.c=.o)

# TARGETDIR - target directory to place executable file
TARGETDIR=bin/

# EXE - executable file
EXE=filterfasta

all: $(EXE)
	@echo "Build complete: $(TARGETDIR)$(EXE)"

# $@ - name of target rule (in this case, $(EXE))
$(EXE): $(TARGETDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) $(INCLUDES) $(LFLAGS) $(SOURCES) $(LIBS) -o $(TARGETDIR)$@
	@echo "Compilation and linking complete"

# This is a suffix replacement rule for building .o's from .c's
# It uses automatic variables
#	$< - name of the first prerequisite of the rule (.c file)
#	$@ - name of target rule (.o file)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
	@echo "Compilation complete"

$(TARGETDIR):
	@mkdir $(TARGETDIR)

clean:
	@rm $(OBJECTS) $(TARGETDIR)$(EXE)

rebuild: clean all
