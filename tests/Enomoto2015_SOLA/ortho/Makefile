FC = gfortran-mp-5
FFLAGS = -O3 -I${HOME}/local/include/enomo_gfortran
LDFLAGS = -L${HOME}/local/lib
LIBS = -lenomo_gfortran
TARGET = calc ortho

all : $(TARGET)

calc : calc.o
	$(FC) $(LDFLAGS) $^ $(LIBS) -o $@

ortho : ortho.o
	$(FC) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(TARGET) *.o

%.o : %.f90
	$(FC) $(FFLAGS) $^ -c
