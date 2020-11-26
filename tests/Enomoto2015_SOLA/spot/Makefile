FC = gfortran
FFLAGS = -O2 -fconvert=big-endian -I../../..
NAME = enomo
LDFLAGS = -L../../.. -l$(NAME)
TARGET = spot

all : $(TARGET)

clean :
	rm -f $(TARGET) *.o

% : %.f90
	${FC} ${FFLAGS} ${LDFLAGS} $< -o $@

