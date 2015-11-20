program rpie2
  use kind_module, only: i4b, dp, qp
  use math_module, only: rad2deg=>math_rad2deg
  use glatwgt_module, only: glatwgt_calc
  use glatwgtq_module, only: glatwgtq_calc
  use alfx_module, only: alfx_init, alfx_clean, alfx_calc_inline
!  use alfx_module, only: alfx_init, alfx_clean, alfx_calc
  use alff_module, only: alff_init, alff_clean, alff_calc
!  use alfq_module, only: alfq_init, alfq_clean, alfq_calc
  use alfxq_module, only: alfxq_init, alfxq_clean, alfxq_calc_inline
!  use alfxq_module, only: alfxq_init, alfxq_clean, alfxq_calc
  use mpi
  implicit none
  
  integer :: iargc

  integer(kind=i4b), parameter :: oun = 41, ne = 5
  character(len=*), parameter :: ofile = "rpie2.txt"

  integer(kind=i4b) :: ntrunc, nlat, nlath, i, j, jj, k, jg, jmax
  character(len=5) :: ntrunc_str
  real(kind=dp) :: sx, sf, lat(1)
  real(kind=dp), dimension(:), allocatable :: &
    glat, gwgt, sbuf, rbuf
  real(kind=dp), dimension(:, :), allocatable :: e
  real(kind=dp), dimension(:,:,:), allocatable :: pnm
  real(kind=qp) :: tx, tf, tq, sq, latq(1)
  real(kind=qp), dimension(:), allocatable :: glatq, gwgtq
  real(kind=qp), dimension(:,:,:), allocatable :: pnmq

  integer(kind=i4b) :: ierr, nproc, myrank

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

  if (iargc() < 1) then
    print *, "Usage :: rpie2 ntrunc"
    stop
  end if

  call getarg(1,ntrunc_str)
  read(unit=ntrunc_str,fmt=*) ntrunc
  nlat = (ntrunc+1)*3/2
  nlath = nlat/2
  jmax = nlath/nproc
!  print *, "ntrunc=", ntrunc, " nlat=", nlat, " jmax=", jmax

  if (myrank == 0) then
    open(unit=oun,file=ofile,status="replace",action="write")
    allocate(rbuf(nlath*ne), e(ne, nlath))
  end if

  allocate(glat(nlat), gwgt(nlat))
print *, "glatwgt_calc"
  call glatwgt_calc(glat,gwgt)
  allocate(glatq(nlat), gwgtq(nlat))
print *, "glatwgtq_calc"
  call glatwgtq_calc(glatq,gwgtq)
  deallocate(gwgt, gwgtq)

  allocate(pnm(0:ntrunc,0:ntrunc,1))
print *, "alfx_init"
  call alfx_init(ntrunc)
print *, "alff_init"
  call alff_init(ntrunc)
  allocate(pnmq(0:ntrunc,0:ntrunc,1))
!print *, "alfq_init"
!  call alfq_init(ntrunc)
print *, "alfxq_init"
  call alfxq_init(ntrunc)

  allocate(sbuf(jmax * ne))
  
  do j=1, jmax
print *, j, "/", jmax, j*100.0_dp/jmax
    jg = jmax * myrank + j
    lat(1) = glat(jg)
    latq(1) = glatq(jg)
!    call alfq_calc(latq,pnmq)
    call alfxq_calc_inline(latq,pnmq)
!    call alfxq_calc(latq,pnmq)
    sq = nmsumq2(pnmq)
    tq = nmsumq(pnmq)
    call alfx_calc_inline(lat,pnm)
!    call alfx_calc(lat,pnm)
    tx = dnmsumq(pnm,pnmq)
    sx = nmsum2(pnm)
    call alff_calc(lat,pnm)
    tf = dnmsumq(pnm,pnmq)
    sf = nmsum2(pnm)
    jj = (j - 1) * ne 
    sbuf(jj + 1) = real(tx/tq, kind=dp)
    sbuf(jj + 2) = real(tf/tq, kind=dp)
    sbuf(jj + 3) = abs(sx / (ntrunc + 1.0_dp)**2 - 1.0_dp)
    sbuf(jj + 4) = abs(sf / (ntrunc + 1.0_dp)**2 - 1.0_dp)
    sbuf(jj + 5) = real(abs(sq / (ntrunc + 1.0_qp)**2 - 1.0_qp), kind=dp)
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Gather(sbuf, jmax * ne, MPI_REAL8, rbuf(jmax * ne * myrank + 1), jmax * ne, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  deallocate(pnm, pnmq, glatq, sbuf)
  if (myrank == 0) then
    do i = 0, nproc - 1
      do j = 1, jmax
        jg = jmax * i + j
        e(:, jg) = rbuf((jg - 1) * ne + 1:jg * ne)
      end do
    end do
    deallocate(rbuf)
    do j = 1, nlath 
      write(unit=oun,fmt=*), j, glat(j)*rad2deg, (e(i,j), i = 1, ne)
    end do
    deallocate(e)
  end if

  deallocate(glat)
  call alfx_clean()
  call alff_clean()
!  call alfq_clean()
  call alfxq_clean()

  close(unit=oun)

  call MPI_Finalize(ierr)

contains

  function nmsumq(pnmq) result (s)
    real(kind=qp), dimension(0:,0:, :) :: pnmq
    integer(kind=i4b) :: m, n, ntrunc
    real(kind=qp) :: s, t, c, y

    ntrunc = size(pnmq, 1) - 1
    s = 0.0_dp
    c = 0.0_dp
    do m = ntrunc, 0, -1
      do n = m, ntrunc
        y = abs((pnmq(n, m, 1))) - c ! apply correction
        t = s + y            ! add the corrected input
        c = (t - s) - y      ! calculate correction
        s = t
      end do
    end do

  end function nmsumq

  function dnmsumq(pnm,pnmq) result (s)
    real(kind=dp), dimension(0:,0:, :) :: pnm
    real(kind=qp), dimension(0:,0:, :) :: pnmq
    integer(kind=i4b) :: m, n, ntrunc
    real(kind=qp) :: s, t, c, y

    ntrunc = size(pnm, 1) - 1
    s = 0.0_dp
    c = 0.0_dp
    do m = ntrunc, 0, -1
      do n = m, ntrunc
        y = abs((pnm(n, m, 1) - pnmq(n, m, 1))) - c ! apply correction
        t = s + y            ! add the corrected input
        c = (t - s) - y      ! calculate correction
        s = t
      end do
    end do

  end function dnmsumq

  function nmsum2(pnm) result (s)
    real(kind=dp), dimension(0:,0:, :) :: pnm
    integer(kind=i4b) :: m, n, ntrunc
    real(kind=dp) :: s, t, u, c, y

    ntrunc = size(pnm, 1) - 1
    s = 0.0_dp
    c = 0.0_dp
    do m = ntrunc, 1, -1
      u = 0.0_dp
      do n = m, ntrunc
        y = pnm(n, m, 1)**2 - c ! apply correction
        t = u + y               ! add the corrected input
        c = (t - u) - y         ! calculate correction
        u = t
      end do
      s = s + 4.0_dp * u
    end do
    c = 2.0_dp * c
    u = 0.0_dp
    do n = 0, ntrunc
      y = pnm(n, 0, 1)**2 - c ! apply correction
      t = u + y               ! add the corrected input
      c = (t - u) - y         ! calculate correction
      u = t
    end do
    s = s + 2.0_dp * u

  end function nmsum2

  function nmsumq2(pnm) result (s)
    real(kind=qp), dimension(0:,0:, :) :: pnm
    integer(kind=i4b) :: m, n, ntrunc
    real(kind=qp) :: s, t, u, c, y

    ntrunc = size(pnm, 1) - 1
    s = 0.0_qp
    c = 0.0_qp
    do m = ntrunc, 1, -1
      u = 0.0_qp
      do n = m, ntrunc
        y = pnm(n, m, 1)**2 - c ! apply correction
        t = u + y               ! add the corrected input
        c = (t - u) - y         ! calculate correction
        u = t
      end do
      s = s + 4.0_qp * u
    end do
    c = 2.0_qp * c
    u = 0.0_qp
    do n = 0, ntrunc
      y = pnm(n, 0, 1)**2 - c ! apply correction
      t = u + y               ! add the corrected input
      c = (t - u) - y         ! calculate correction
      u = t
    end do
    s = s + 2.0_qp * u

  end function nmsumq2

end program rpie2
