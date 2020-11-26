FC = mpif90
FFLAGS = -O2 -I../../..
LDFLAGS = -L../../..
LIBS = -lenomo
TARGET = rpie2

$(TARGET) : $(TARGET).o
	$(FC) $(LDFLAGS) $^ $(LIBS) -o $@

clean :
	rm -f $(TARGET) *.o

%.o : %.f90
	$(FC) $(FFLAGS) $^ -c
