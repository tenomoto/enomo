module vectrans_module
  use kind_module, only: dp, i4b
  private

  real(kind=dp), dimension(:), allocatable, private :: wvts, wvhags, work, workt
  real(kind=dp), dimension(:,:,:), allocatable, private :: v, w, br, bi, cr, ci
  integer(kind=i4b), private :: &
    nlon, nlat, nt, lvhags, lwvts, mdab, ndab, idvw, jdvw, ierror, lwork, lworkt

  public vectrans_init, vectrans_clean, vectrans_dlat

contains

  subroutine vectrans_init(nx,ny,n)
    implicit none

    integer(kind=i4b), intent(in) :: nx, ny, n

    real(kind=dp), dimension(:), allocatable :: dwork, dworkt
    integer(kind=i4b) :: l1, l2, ldwork, ldworkt

    nlon = nx
    nlat = ny
    nt = n

    l1 = min0(nlat,nlon/2) 
    l2 = nlat/2
    lvhags = l1*l2*(nlat+nlat-l1+1)+nlon+15
    ldwork = (3*nlat*(nlat+3)+2)/2
    allocate(wvhags(lvhags),dwork(ldwork))
    call vhagsi(nlat,nlon,wvhags,lvhags,dwork,ldwork,ierror)
    if (ierror/=0) then
      print *, "*** error in vhagsi: ", ierror
      stop
!    else
!      print *, "nlat,nlon:", nlat, nlon
!      print *, "wvhags:", minval(wvhags), maxval(wvhags), size(wvhags), lvhags
!      print *, "dwork:", size(dwork)
    end if

    lwvts  = l1*l2*(nlat+nlat-l1+1)+nlon+15
    ldworkt = 3*nlat+2
    lworkt = 3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+(5*l2+2)*nlat
    allocate(wvts(lwvts),workt(lworkt),dworkt(ldworkt))
    call vtsgsi(nlat,nlon,wvts,lwvts,workt,lworkt,dworkt,ldworkt,ierror)
    if (ierror/=0) then
      print *, "*** error in vtsgsi: ", ierror
      stop
!    else
!      print *, "wvts:", minval(wvts), maxval(wvts)
    end if

    mdab = min0(nlat,nlon/2)
    ndab = nlat
    allocate(br(mdab,ndab,nt),bi(mdab,ndab,nt),cr(mdab,ndab,nt),ci(mdab,ndab,nt))

    idvw = nlat
    jdvw = nlon
    allocate(v(idvw,jdvw,nt), w(idvw,jdvw,nt))

    deallocate(dwork,dworkt,workt)
    lwork = max(3*nlat*(nlat+1)+2,(2*nt+1)*nlat*nlon)
    lworkt = (2*nt+1)*nlat*nlon
    allocate(work(lwork),workt(lworkt))

  end subroutine vectrans_init

  subroutine vectrans_clean
    implicit none

    deallocate(wvhags,wvts,work,br,bi,cr,ci,v,w)

  end subroutine vectrans_clean

  subroutine vectrans_dlat(ui,vi,uy,vy,np2sp)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: ui, vi
    real(kind=dp), dimension(:,:,:), intent(inout) :: uy, vy
    logical, optional, intent(in) :: np2sp

    integer(kind=i4b), parameter :: ityp = 0

    integer(kind=i4b) :: i, j, jr, k

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
    call vhags(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvhags,lvhags,work,lwork,ierror)
    if (ierror/=0) then
      print *, "*** error in vhags: ", ierror
      stop
    end if

!! debug
!    print *, "v:", minval(v), maxval(v)
!    print *, "w:", minval(w), maxval(w)
!    print *, "wvhags:", minval(wvhags), maxval(wvhags)
!    print *, "br:", minval(br), maxval(br)
!    print *, "bi:", minval(bi), maxval(bi)
!    print *, "cr:", minval(cr), maxval(cr)
!    print *, "ci:", minval(ci), maxval(ci)
!stop

    call vtsgs(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvts,lwvts,workt,lworkt,ierror)
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
      print *, "*** error in vtsgs: ", ierror
      stop
    end if

  end subroutine vectrans_dlat

end module vectrans_module
