PREFIX = ${HOME}/local
FC = mpif90
FFLAGS = -O2 -fconvert=big-endian -I$(HOME)/local/include/$(NAME)
NAME = enomo_gfortran
LDFLAGS = -L$(HOME)/local/lib -l$(NAME)
TARGET = calc_checksum

all : $(TARGET)

clean :
	rm -f $(TARGET) *.o

% : %.f90
	${FC} ${FFLAGS} ${LDFLAGS} $< -o $@

