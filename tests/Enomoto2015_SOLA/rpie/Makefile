FC = mpif90
FFLAGS = -O2 -I${HOME}/local/include/enomo_gfortran
LDFLAGS = -L${HOME}/local/lib
LIBS = -lenomo_gfortran
TARGET = rpie2

$(TARGET) : $(TARGET).o
	$(FC) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(TARGET) *.o

%.o : %.f90
	$(FC) $(FFLAGS) $^ -c
