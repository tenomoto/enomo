#NAME = enomo
FC = gfortran-mp-5
FFLAGS = -O2 -fconvert=big-endian -I../
CC = gfortran-mp-5
LDFLAGS = -L..
LIBNAME = enomo_gfortran
LDLIBS = -l$(LIBNAME)
PROGS = test_alf test_alff test_alfq test_alfx test_alfxq test_glatwgt test_glatwgtq test_sum test_xreal

all : $(PROGS)

clean :
	rm -f $(PROGS) *.dat *.txt

%.o : %.f90
	${FC} ${FFLAGS} $< -c

$(PROGS) : ../lib$(LIBNAME).a
