program spot
  use kind_module, only: i4b, dp, qp
  use math_module, only: pi => math_pi, piq => math_piq
  use alfq_module, only: alfq_init, alfq_calc_m, alfq_clean
  use alfx_module, only: alfx_init, alfx_calc_m, alfx_clean
  use alff_module, only: alff_init, alff_calc_m, alff_clean
  use alfxq_module, only: alfxq_init, alfxq_calc_m, alfxq_clean
  implicit none

  integer(kind=i4b), parameter :: ntrunc = 21600, nstart = 21598, &
    mstart = 10798, mend = 10800
!  integer(kind=i4b), parameter :: ntrunc = 39, nstart = 0, nend = 4, &
!    mstart = 0, mend = 4
  integer(kind=i4b) :: n, m
  real(kind=dp), dimension(1) :: lat
  real(kind=dp), dimension(0:ntrunc, 1) :: alfxm, alffm
  real(kind=qp), dimension(1) :: latq
  real(kind=qp), dimension(0:ntrunc, 1) :: alfqm, alfxqm

  lat(1) = pi/6.0_dp
  latq(1) = piq/6.0_qp
  call alfx_init(ntrunc)
  call alff_init(ntrunc)
  call alfq_init(ntrunc)
  call alfxq_init(ntrunc)
  do m = mstart, mend
    call alfx_calc_m(m, lat, alfxm)
    call alff_calc_m(m, lat, alffm)
    call alfq_calc_m(m, latq, alfqm)
    call alfxq_calc_m(m, latq, alfxqm)
    do n = nstart, ntrunc
!    do n = m, nend
      print *, n, m, alfxm(n,1)*2, alffm(n,1)*2, alfqm(n,1)*2, alfxqm(n,1)*2
    end do
  end do
  call alfx_clean()
  call alff_clean()
  call alfq_clean()
  call alfxq_clean()

end program spot
