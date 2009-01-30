PREFIX = ${HOME}/local
FC = g95
FFLAGS = -O2 -fendian=big
AR = ar
ARFLAGS = -cru
TARGET = libenomo.a
SRC = type_module.f90 math_module.f90 earth_module.f90 air_module.f90 water_module.f90 \
      sphere_module.f90 glatwgt_module.f90 interpolate_module.f90 slp_module.f90       \
      confmap_module.f90 regrid_module.f90 upstream_module.f90
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

math_module.o : type_module.o
earth_module.o : type_module.o math_module.o
air_module.o : type_module.o math_module.o
water_module.o : type_module.o air_module.o math_module.o
slp_module.o : type_module.o earth_module.o air_module.o
sphere_module.o : type_module.o math_module.o glatwgt_module.o
glatwgt_module.o : type_module.o math_module.o
interpolate_module.o : type_module.o math_module.o
confmap_module.o : type_module.o math_module.o sphere_module.o
regrid_module.o : type_module.o math_module.o
upstream_module.o : type_module.o math_module.o earth_module.o sphere_module.o regrid_module.o
