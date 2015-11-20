program ortho
  use kind_module, only: i4b, dp, qp
  use math_module, only: rad2deg=>math_rad2deg
  use glatwgt_module, only: glatwgt_calc
  use glatwgtq_module, only: glatwgtq_calc
  implicit none
  
  integer :: iargc

  integer(kind=i4b), parameter :: unx = 41, unf = 42, unq = 43
  character(len=*), parameter :: datadir = "/Volumes/Pegasus/alfx/run/"

  integer(kind=i4b) :: ntrunc, nlat, nlath, i, j, m, n, nn
  character(len=5) :: ntrunc_str
  real(kind=dp) :: lat(1), smaxx, smaxf, t, y
  real(kind=dp), dimension(:), allocatable :: glat, gwgt
  real(kind=dp), dimension(:, :), allocatable :: s, c
  real(kind=dp), dimension(:,:,:), allocatable :: pnm
  real(kind=qp) :: latq(1), smaxq, tq, yq
  real(kind=qp), dimension(:), allocatable :: glatq, gwgtq
  real(kind=qp), dimension(:,:), allocatable :: sq, cq
  real(kind=qp), dimension(:,:,:), allocatable :: pnmq

  if (iargc() < 1) then
    print *, "Usage :: ortho ntrunc"
    stop
  end if

  call getarg(1,ntrunc_str)
  read(unit=ntrunc_str,fmt=*) ntrunc
  nlat = (ntrunc+1)*3/2
  nlath = nlat/2
!  print *, "ntrunc=", ntrunc, " nlat=", nlat, " jmax=", jmax

  allocate(pnm(0:ntrunc,0:ntrunc,1), pnmq(0:ntrunc,0:ntrunc,1))

  allocate(glat(nlat), gwgt(nlat))
  call glatwgt_calc(glat,gwgt)
  allocate(glatq(nlat), gwgtq(nlat))
  call glatwgtq_calc(glatq,gwgtq)

  open(unit=unx, file=datadir//"alfx_T"//trim(ntrunc_str)//".dat", &
    access="direct", recl=size(pnm)*8, &
    status="old",action="read")
  open(unit=unf, file=datadir//"alff_T"//trim(ntrunc_str)//".dat", &
    access="direct", recl=size(pnm)*8, &
    status="old",action="read")
!  open(unit=unq, file=datadir//"alfq_T"//trim(ntrunc_str)//".dat", &
  open(unit=unq, file=datadir//"alfxq_T"//trim(ntrunc_str)//".dat", &
    access="direct", recl=size(pnm)*16, &
    status="old",action="read")

  allocate(s(0:ntrunc-1,0:ntrunc-1), sq(0:ntrunc-1,0:ntrunc-1), &
           c(0:ntrunc-1,0:ntrunc-1), cq(0:ntrunc-1,0:ntrunc-1))
  do m = 0, ntrunc-1
! alfx
    s(:,:) = 0.0_dp
    c(:,:) = 0.0_dp
    do j = 1, nlath
      read(unit=unx, rec=j) pnm
      do n = m, ntrunc-1
        do nn = n+1, ntrunc-1
          y = pnm(n,m,1) * pnm(nn,m,1) * (1.0_dp + (-1.0_dp)**(n+nn)) * gwgt(j) - c(n,nn)
          t = s(n,nn) + y
          c(n,nn) = (t - s(n,nn)) - y
          s(n,nn) = t
        end do
      end do
    end do
    smaxx = maxval(abs(s))
! alff
    s(:,:) = 0.0_dp
    c(:,:) = 0.0_dp
    do j = 1, nlath
      read(unit=unf, rec=j) pnm
      do n = m, ntrunc-1
        do nn = n+1, ntrunc-1
          y = pnm(n,m,1) * pnm(nn,m,1) * (1.0_dp + (-1.0_dp)**(n+nn)) * gwgt(j) - c(n,nn)
          t = s(n,nn) + y
          c(n,nn) = (t - s(n,nn)) - y
          s(n,nn) = t
        end do
      end do
    end do
    smaxf = maxval(abs(s))
! alff
    sq(:,:) = 0.0_dp
    cq(:,:) = 0.0_dp
    do j = 1, nlath
      read(unit=unq, rec=j) pnmq
      do n = m, ntrunc-1
        do nn = n+1, ntrunc-1
          yq = pnmq(n,m,1) * pnmq(nn,m,1) * (1.0_qp + (-1.0_qp)**(n+nn)) * gwgtq(j) - cq(n,nn)
          tq = sq(n,nn) + yq
          cq(n,nn) = (tq - sq(n,nn)) - yq
          sq(n,nn) = tq
        end do
      end do
    end do
    smaxq = maxval(abs(sq))
    print *, m, smaxx, smaxf, smaxq
  end do

  close(unit=unx)
  close(unit=unf)
  close(unit=unq)

  deallocate(glat, gwgt, glatq, gwgtq)
  deallocate(s, c, sq, cq, pnm, pnmq)

end program ortho
