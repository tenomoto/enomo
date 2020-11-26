FC = gfortran
FFLAGS = -O3 -I../../..
LDFLAGS = -L../../..
LIBS = -lenomo
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
