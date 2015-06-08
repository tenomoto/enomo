PREFIX = ${HOME}/local
FC = g95
#FC = gfortran-mp-4.5
NAME = enomo
FFLAGS = -O2 -fendian=big -I/opt/local/include
#FFLAGS = -O2 -fconvert=big-endian -I/opt/local/include
AR = ar
ARFLAGS = -cru
TARGET = lib${NAME}.a
SRC = \
air_module.f90 \
alf_module.f90 \
alff_module.f90 \
alfx_module.f90 \
besttrack_module.f90 \
calendar_module.f90 \
confmap_module.f90 \
earth_module.f90 \
fft_module.f90 \
glatwgt_module.f90 \
grads_module.f90 \
grmsm_sig_module.f90 \
grmsm_sfc_module.f90 \
integer_module.f90 \
interpolate_module.f90 \
kind_module.f90 \
math_module.f90 \
matrix_module.f90 \
moist_module.f90 \
regrid_module.f90 \
search_module.f90 \
semiimplicit_module.f90 \
shregrid_module.f90 \
shtrans_module.f90 \
sigma_module.f90 \
sigmap_module.f90 \
slp_module.f90 \
sphere_module.f90 \
stability_module.f90 \
string_module.f90 \
svd_module.f90 \
time_module.f90 \
udunits_module.f90 \
upstream_module.f90 \
vectrans_module.f90 \
vectransf_module.f90 \
water_module.f90 \
xreal_module.f90 \
machine_module.f90 \
sum_module.f90
#monotonecubic_module.f90 \
#vis5d_module.f90
OBJ = ${SRC:%.f90=%.o}
MOD = ${SRC:%.f90=%.mod}

${TARGET} : ${OBJ}
	${AR} ${ARFLAGS} $@ $^ 

install :
	install -m 644 ${TARGET} ${PREFIX}/lib/
	install -d ${PREFIX}/include/${NAME}
	install -m 644 ${MOD} ${PREFIX}/include/${NAME}

uninstall :
	rm -f ${PREFIX}/lib/${TARGET}
	rm -f ${PREFIX}/include/${NAME}/*.mod

test : test.o ${TARGET}
	${FC}  $< ./${TARGET} -o test

clean :
	rm -f ${TARGET} ${OBJ} ${MOD}

%.o : %.f90
	${FC} ${FFLAGS} $< -c

air_module.o : kind_module.o math_module.o
alf_module.o : kind_module.o math_module.o glatwgt_module.o
alff_module.o : kind_module.o math_module.o integer_module.o glatwgt_module.o alf_module.o
alfx_module.o : kind_module.o math_module.o xreal_module.o glatwgt_module.o alf_module.o
besttrack_module.o : kind_module.o time_module.o string_module.o
calendar_module.o : kind_module.o string_module.o
confmap_module.o : kind_module.o math_module.o sphere_module.o
earth_module.o : kind_module.o math_module.o
fft_module.o : kind_module.o
glatwgt_module.o : kind_module.o math_module.o machine_module.o
grads_module.o : kind_module.o string_module.o
grmsm_sig_module.o : kind_module.o sigma_module.o
grmsm_sfc_module.o : kind_module.o
integer_module.o : kind_module.o
interpolate_module.o : kind_module.o math_module.o
io_module.o : kind_module.o
math_module.o : kind_module.o
matrix_module.o : kind_module.o
moist_module.o : kind_module.o math_module.o air_module.o water_module.o
regrid_module.o : kind_module.o math_module.o search_module.o
search_module.o : kind_module.o
semiimplicit_module.o : kind_module.o
shregrid_module.o : kind_module.o
shtrans_module.o : kind_module.o
sigma_module.o : kind_module.o
sigmap_module.o : kind_module.o
slp_module.o : kind_module.o earth_module.o air_module.o
sphere_module.o : kind_module.o math_module.o glatwgt_module.o
stability_module.o : kind_module.o math_module.o air_module.o water_module.o
string_module.o : kind_module.o
svd_module.o : kind_module.o matrix_module.o
time_module.o : kind_module.o calendar_module.o
udunits_module.o : kind_module.o
upstream_module.o : kind_module.o math_module.o earth_module.o sphere_module.o regrid_module.o
vectrans_module.o : kind_module.o shtrans_module.o
vectransf_module.o : kind_module.o
#vis5d_module.o : kind_module.o time_module.o calendar_module.o
water_module.o : kind_module.o air_module.o
xreal_module.o : kind_module.o
machine_module.o : kind_module.o
sum_module.o : kind_module.o
