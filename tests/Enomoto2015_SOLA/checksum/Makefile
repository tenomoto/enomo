FC = mpif90
FFLAGS = -O2 -fconvert=big-endian -I../../..
NAME = enomo
LDFLAGS = -L../../.. -l$(NAME)
TARGET = calc_checksum

all : $(TARGET)

clean :
	rm -f $(TARGET) *.o

% : %.f90
	${FC} ${FFLAGS} ${LDFLAGS} $< -o $@

