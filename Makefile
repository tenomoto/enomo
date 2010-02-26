PREFIX = ${HOME}/local
FC = g95
FFLAGS = -O2 -fendian=big -I/opt/local/include -framework accelerate
AR = ar
ARFLAGS = -cru
TARGET = libenomo.a
SRC = kind_module.f90 math_module.f90 earth_module.f90 air_module.f90 water_module.f90 \
      sphere_module.f90 glatwgt_module.f90 interpolate_module.f90 slp_module.f90       \
      confmap_module.f90 regrid_module.f90 upstream_module.f90 fft_module.f90          \
      vectrans_module.f90 calendar_module.f90 moist_module.f90 stability_module.f90    \
      string_module.f90 grads_module.f90 besttrack_module.f90 \
      time_module.f90 udunits_module.f90 search_module.f90 svd_module.f90 matrix_module.f90
OBJ = ${SRC:%.f90=%.o}
MOD = ${SRC:%.f90=%.mod}

${TARGET} : ${OBJ}
	${AR} ${ARFLAGS} $@ $^ 

install :
	install -m 644 ${TARGET} ${PREFIX}/lib/
	install -d ${PREFIX}/include/enomo
	install -m 644 ${MOD} ${PREFIX}/include/enomo

uninstall :
	rm -f ${PREFIX}/lib/${TARGET}
	rm -f ${PREFIX}/include/enomo/*.mod

clean :
	rm -f ${TARGET} ${OBJ} ${MOD}

%.o : %.f90
	${FC} ${FFLAGS} $< -c

math_module.o : kind_module.o
earth_module.o : kind_module.o math_module.o
air_module.o : kind_module.o math_module.o
water_module.o : kind_module.o air_module.o
slp_module.o : kind_module.o earth_module.o air_module.o
sphere_module.o : kind_module.o math_module.o glatwgt_module.o
glatwgt_module.o : kind_module.o math_module.o
interpolate_module.o : kind_module.o math_module.o
confmap_module.o : kind_module.o math_module.o sphere_module.o
regrid_module.o : kind_module.o math_module.o search_module.o
upstream_module.o : kind_module.o math_module.o earth_module.o sphere_module.o regrid_module.o
fft_module.o : kind_module.o
vectrans_module.o : kind_module.o
calendar_module.o : kind_module.o string_module.o
moist_module.o : kind_module.o math_module.o air_module.o water_module.o
stability_module.o : kind_module.o math_module.o air_module.o water_module.o
io_module.o : kind_module.o
string_module.o : kind_module.o
grads_module.o : kind_module.o string_module.o
besttrack_module.o : kind_module.o time_module.o string_module.o
time_module.o : kind_module.o
udunits_module.o : kind_module.o
#vis5d_module.o : kind_module.o time_module.o calendar_module.o
search_module.o : kind_module.o
svd_module : kind_module.o
matrix_module : kind_module.o
