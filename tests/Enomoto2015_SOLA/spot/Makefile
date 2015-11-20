PREFIX = ${HOME}/local
FC = gfortran-mp-5
FFLAGS = -O2 -fconvert=big-endian -I$(HOME)/local/include/$(NAME)
NAME = enomo_gfortran
LDFLAGS = -L$(HOME)/local/lib -l$(NAME)
TARGET = spot

all : $(TARGET)

clean :
	rm -f $(TARGET) *.o

% : %.f90
	${FC} ${FFLAGS} ${LDFLAGS} $< -o $@

