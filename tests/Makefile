#NAME = enomo
FC = gfortran
FFLAGS = -O2 -fconvert=big-endian -I../
CC = gfortran
LDFLAGS = -L.. -L${HOME}/local/lib
LIBNAME = enomo
LDLIBS = -l$(LIBNAME) -lspherepack
PROGS = test_alf test_alff test_alfq test_alfx test_alfxq test_glatwgt test_glatwgtq test_sum test_xreal test_shtrans test_gaqd

all : $(PROGS)

clean :
	rm -f $(PROGS) *.dat *.txt

%.o : %.f90
	${FC} ${FFLAGS} $< -c

$(PROGS) : ../lib$(LIBNAME).a
