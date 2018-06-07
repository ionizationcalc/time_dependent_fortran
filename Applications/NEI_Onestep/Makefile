# The compiler
FC = pgfortran
# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
FCFLAGS = -O2
# flags forall (e.g. look for system .mod files, required in gfortran)
#FCFLAGS += -I/usr/include

# libraries needed for linking, unused in the examples
#LDFLAGS = -li_need_this_lib

# List of executables to be built within the package
PROGRAMS = main

# "make" builds all
all: $(PROGRAMS)

objs = sub_read_eigen_matrix.o sub_read_prameter.o sub_solve_ioni_onestep.o func_equilibrium_eigen.o func_solveionization_eigen.o

main: main.o mod_eigen_matrix.o ${objs}
main.o: mod_eigen_matrix.o ${objs}
sub_read_eigen_matrix.o: mod_eigen_matrix.o
sub_read_prameter.o: mod_eigen_matrix.o
sub_solve_ioni_onestep.o: mod_eigen_matrix.o
func_equilibrium_eigen.o: mod_eigen_matrix.o
func_solveionization_eigen.o: mod_eigen_matrix.o


# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================
# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)