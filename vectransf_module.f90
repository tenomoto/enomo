module vectransf_module
  use kind_module, only: dp, i4b
  use math_module, only: pi=>math_pi
  private

  real(kind=dp), private :: kfil
  real(kind=dp), dimension(:), allocatable, private :: wvts, wvhaes, wshses, work
  real(kind=dp), dimension(:,:,:), allocatable, private :: v, w, br, bi, cr, ci
  integer(kind=i4b), private :: &
    nlon, nlat, nt, ntrunc, lvhaes, lwvts, lshses, mdab, ndab, idvw, jdvw, ierror, lwork

  public vectransf_init, vectransf_clean, vectransf_dlat, vectransf_uv2vrdv

contains

  subroutine vectransf_init(nx,ny,n)
    implicit none

    integer(kind=i4b), intent(in) :: nx, ny, n

    real(kind=dp), dimension(:), allocatable :: dwork
    integer(kind=i4b) :: l1, l2, ldwork

    nlon = nx
    nlat = ny
    nt = n

    l1 = min0(nlat,nlon/2) 
    l2 = nlat/2
    lvhaes = l1*l2*(nlat+nlat-l1+1)+nlon+15
    ldwork = (3*nlat*(nlat+3)+2)/2
    lwork = 3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+5*l2*nlat
    allocate(wvhaes(lvhaes),dwork(ldwork),work(lwork))
    call vhaesi(nlat,nlon,wvhaes,lvhaes,work,lwork,dwork,ldwork,ierror)
    if (ierror/=0) then
      print *, "*** error in vhaesi: ", ierror
      stop
!    else
!      print *, "nlat,nlon:", nlat, nlon
!      print *, "wvhaes:", minval(wvhaes), maxval(wvhaes), size(wvhaes), lvhaes
!      print *, "dwork:", size(dwork)
    end if
    deallocate(dwork)

    lwvts  = l1*l2*(nlat+nlat-l1+1)+nlon+15
    ldwork = 2*(nlat+1)
    allocate(wvts(lwvts),dwork(ldwork))
    call vtsesi(nlat,nlon,wvts,lwvts,work,lwork,dwork,ldwork,ierror)
    if (ierror/=0) then
      print *, "*** error in vtsesi: ", ierror
      stop
!    else
!      print *, "wvts:", minval(wvts), maxval(wvts)
    end if
    deallocate(dwork,work)

    lshses = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
    lwork = 5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
    ldwork = nlat+1
    allocate(wshses(lshses),work(lwork),dwork(ldwork))
    call shsesi(nlat,nlon,wshses,lshses,work,lwork,dwork,ldwork,ierror)
    if (ierror/=0) then
      print *, "*** error in shsesi: ", ierror
      stop
    end if

    mdab = min0(nlat,nlon/2)
    ndab = nlat
    allocate(br(mdab,ndab,nt),bi(mdab,ndab,nt),cr(mdab,ndab,nt),ci(mdab,ndab,nt))

    idvw = nlat
    jdvw = nlon
    allocate(v(idvw,jdvw,nt), w(idvw,jdvw,nt))

    deallocate(dwork,work)
    lwork = max(3*nlat*(nlat+1)+2,(2*nt+1)*nlat*nlon)
    lwork = max(lwork,(2*nt+1)*nlat*nlon)
    lwork = max(lwork,nlat*((nt+1)*nlon+2*nt*l1+1))
    allocate(work(lwork))

! Hoskins filter
    ntrunc = nlon/3
    if (mod(nlon,3)/=0) then
      ntrunc = ntrunc - 1
    end if
    kfil = -log(0.1_dp)/(ntrunc*(ntrunc+1))**2

  end subroutine vectransf_init

  subroutine vectransf_clean
    implicit none

    deallocate(wvhaes,wvts,wshses,work,br,bi,cr,ci,v,w)

  end subroutine vectransf_clean

  subroutine vectransf_dlat(ui,vi,uy,vy,np2sp,lfil)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: ui, vi
    real(kind=dp), dimension(:,:,:), intent(inout) :: uy, vy
    logical, optional, intent(in) :: np2sp,lfil

    integer(kind=i4b), parameter :: ityp = 0

    real(kind=dp) :: fil
    integer(kind=i4b) :: i, j, jr, k, n, m

! spherepack uses colatitudes (colat = pi/2 - lat)
!   - colat = pi/2 - lat, colat = 0 at NP
!   - data is stored from NP->SP
!   - meridional wind is reveresed -v
!   - merional derivative is reversed -d/dtheta
    if (present(np2sp).and.np2sp) then
      do k=1, nt
        do j=1, nlat
          v(j,1:nlon,k) = -vi(1:nlon,j,k)
          w(j,1:nlon,k) =  ui(1:nlon,j,k)
        end do
      end do
    else
      do k=1, nt
        do j=1, nlat
          jr = nlat - j +1
          v(jr,1:nlon,k) = -vi(1:nlon,j,k)
          w(jr,1:nlon,k) =  ui(1:nlon,j,k)
        end do
      end do
    end if
    call vhaes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvhaes,lvhaes,work,lwork,ierror)
    if (ierror/=0) then
      print *, "*** error in vhaes: ", ierror
      stop
    end if

! triangular trunction
    br(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp
    bi(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp
    cr(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp
    ci(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp

    if (present(lfil).and.lfil) then
      do n=0, ntrunc
        fil = exp(-kfil*(n*n*(n+1)*(n+1)))
        do m=0, n
          br(m+1,n+1,:) = fil*br(m+1,n+1,:)
          bi(m+1,n+1,:) = fil*bi(m+1,n+1,:)
          cr(m+1,n+1,:) = fil*cr(m+1,n+1,:)
          ci(m+1,n+1,:) = fil*ci(m+1,n+1,:)
        end do
      end do
    end if

!! debug
!    print *, "v:", minval(v), maxval(v)
!    print *, "w:", minval(w), maxval(w)
!    print *, "wvhaes:", minval(wvhaes), maxval(wvhaes)
!    print *, "br:", minval(br), maxval(br)
!    print *, "bi:", minval(bi), maxval(bi)
!    print *, "cr:", minval(cr), maxval(cr)
!    print *, "ci:", minval(ci), maxval(ci)
!stop

    call vtses(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvts,lwvts,work,lwork,ierror)
    if (ierror==0) then
      if (present(np2sp).and.np2sp) then
        do k=1, nt
          do j=1, nlat
            vy(1:nlon,j,k) =  v(j,1:nlon,k)
            uy(1:nlon,j,k) = -w(j,1:nlon,k)
          end do
        end do
      else
        do k=1, nt
          do j=1, nlat
            jr = nlat - j + 1
            vy(1:nlon,jr,k) =  v(j,1:nlon,k)
            uy(1:nlon,jr,k) = -w(j,1:nlon,k)
          end do
        end do
      end if
    else
      print *, "*** error in vtses: ", ierror
      stop
    end if

  end subroutine vectransf_dlat

  subroutine vectransf_uv2vrdv(ui,vi,vr,dv,np2sp,lfil)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: ui, vi
    real(kind=dp), dimension(:,:,:), intent(inout) :: vr, dv
    logical, optional, intent(in) :: np2sp, lfil

    integer(kind=i4b), parameter :: ityp = 0

    integer(kind=i4b) :: i, j, jr, k, n, m
    real(kind=dp) :: fil

! spherepack uses colatitudes (colat = pi/2 - lat)
!   - colat = pi/2 - lat, colat = 0 at NP
!   - data is stored from NP->SP
!   - meridional wind is reveresed -v
!   - merional derivative is reversed -d/dtheta
    if (present(np2sp).and.np2sp) then
      do k=1, nt
        do j=1, nlat
          v(j,1:nlon,k) = -vi(1:nlon,j,k)
          w(j,1:nlon,k) =  ui(1:nlon,j,k)
        end do
      end do
    else
      do k=1, nt
        do j=1, nlat
          jr = nlat - j +1
          v(jr,1:nlon,k) = -vi(1:nlon,j,k)
          w(jr,1:nlon,k) =  ui(1:nlon,j,k)
        end do
      end do
    end if
    call vhaes(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvhaes,lvhaes,work,lwork,ierror)
    if (ierror/=0) then
      print *, "*** error in vhaes: ", ierror
      stop
    end if

! triangular trunction
    br(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp
    bi(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp
    cr(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp
    ci(ntrunc+2:mdab,ntrunc+2:ndab,:) = 0.0_dp

    if (present(lfil).and.lfil) then
      do n=1, ntrunc
!        fil = n*pi/ntrunc
!        fil = sin(fil)/(fil)
        fil = exp(-kfil*(n*n*(n+1)*(n+1)))
        do m=1, n
          br(m+1,n+1,:) = fil*br(m+1,n+1,:)
          bi(m+1,n+1,:) = fil*bi(m+1,n+1,:)
          cr(m+1,n+1,:) = fil*cr(m+1,n+1,:)
          ci(m+1,n+1,:) = fil*ci(m+1,n+1,:)
        end do
      end do
    end if

!! debug
!    print *, "v:", minval(v), maxval(v)
!    print *, "wvhaes:", minval(wvhaes), maxval(wvhaes)
!    print *, "br:", minval(br), maxval(br)
!    print *, "bi:", minval(bi), maxval(bi)
!    print *, "cr:", minval(cr), maxval(cr)
!    print *, "ci:", minval(ci), maxval(ci)
!stop

    call vrtes(nlat,nlon,ityp,nt,v,idvw,jdvw,cr,ci, &
               mdab,ndab,wshses,lshses,work,lwork,ierror)
    if (ierror==0) then
      if (present(np2sp).and.np2sp) then
        do k=1, nt
          do j=1, nlat
            vr(1:nlon,j,k) =  v(j,1:nlon,k)
          end do
        end do
      else
        do k=1, nt
          do j=1, nlat
            jr = nlat - j + 1
            vr(1:nlon,jr,k) =  v(j,1:nlon,k)
          end do
        end do
      end if
    else
      print *, "*** error in vrtes : ", ierror
      stop
    end if
    call dives(nlat,nlon,ityp,nt,w,idvw,jdvw,br,bi, &
               mdab,ndab,wshses,lshses,work,lwork,ierror)
    if (ierror==0) then
      if (present(np2sp).and.np2sp) then
        do k=1, nt
          do j=1, nlat
            dv(1:nlon,j,k) =  w(j,1:nlon,k)
          end do
        end do
      else
        do k=1, nt
          do j=1, nlat
            jr = nlat - j + 1
            dv(1:nlon,jr,k) =  w(j,1:nlon,k)
          end do
        end do
      end if
    else
      print *, "*** error in dives: ", ierror
      stop
    end if

  end subroutine vectransf_uv2vrdv

end module vectransf_module
