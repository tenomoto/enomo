module vectrans_module
  use shtrans_module, only: shtrans_verbose, shtrans_isinit, shtrans_init, & 
    shtrans_transf, shtrans_transb, shtrans_hoskinsfilter, shtrans_truncate
  use kind_module, only: dp, i4b
  private

  logical, public :: vectrans_verbose = .false.

  integer(kind=i4b), parameter, private :: ityp = 0, isym = 0
  real(kind=dp), dimension(:), allocatable, private :: wvts, wvhags, wvhsgs, work
  real(kind=dp), dimension(:,:,:), allocatable, private :: v, w, br, bi, cr, ci, a, b, g
  integer(kind=i4b), private :: &
    nlon, nlat, ntrunc, mdab, ndab, idvw, jdvw
  logical, private :: lyrev

  public :: vectrans_init, vectrans_clean, vectrans_dlat, vectrans_uv2vrdv, vectrans_gradient
  private :: checkerror, checkval

contains

  subroutine vectrans_init(nx,ny,nt,np2sp)
    implicit none

    integer(kind=i4b), intent(in) :: nx, ny, nt
    logical, optional, intent(in) :: np2sp

    real(kind=dp), dimension(:), allocatable :: dwork
    integer(kind=i4b) :: l1, l2, ldwork, lvhags, lvhsgs, lwvts, lwork, ierror

    nlon = nx
    nlat = ny
    ntrunc = ceiling(nlon/3.0)-1
    lyrev = .true.
    if (present(np2sp).and.np2sp) then
      lyrev = .false.
    end if
    if (vectrans_verbose) then
      print *, "nlat=", nlat, " nlon=", nlon, " ntrunc=", ntrunc, &
        " nt=", nt, " lyrev=", lyrev
    end if

! spherepack uses colatitudes (colat = pi/2 - lat)
!   - colat = pi/2 - lat, colat = 0 at NP
!   - data is stored from NP->SP
!   - meridional wind is reveresed -v
!   - merional derivative is reversed -d/dtheta

! init analysis
    l1 = min0(nlat,nlon/2) 
    l2 = nlat/2
    lvhags = l1*l2*(nlat+nlat-l1+1)+nlon+15
    ldwork = (3*nlat*(nlat+3)+2)/2
    allocate(wvhags(lvhags),dwork(ldwork))
    call vhagsi(nlat,nlon,wvhags,lvhags,dwork,ldwork,ierror)
    call checkerror(ierror,"vhagsi")

! init synthesis
    lvhsgs = l1*l2*(nlat+nlat-l1+1)+nlon+15+2*nlat
    allocate(wvhsgs(lvhsgs))
    call vhsgsi(nlat,nlon,wvhsgs,lvhsgs,dwork,ldwork,ierror)
    call checkerror(ierror,"vhsgsi")
    deallocate(dwork)

! init colatitude derivatives
    lwvts  = l1*l2*(nlat+nlat-l1+1)+nlon+15
    ldwork = 3*nlat+2
    lwork = 3*(max0(l1-2,0)*(nlat+nlat-l1-1))/2+(5*l2+2)*nlat
    allocate(wvts(lwvts),work(lwork),dwork(ldwork))
    call vtsgsi(nlat,nlon,wvts,lwvts,work,lwork,dwork,ldwork,ierror)
    call checkerror(ierror,"vtsgsi")
    deallocate(dwork,work)

    mdab = min0(nlat,nlon/2)
    ndab = nlat
    allocate(br(mdab,ndab,nt),bi(mdab,ndab,nt),cr(mdab,ndab,nt),ci(mdab,ndab,nt))
    allocate(g(nlat,nlon,nt),a(mdab,ndab,nt),b(mdab,ndab,nt))

    idvw = nlat
    jdvw = nlon
    allocate(v(idvw,jdvw,nt), w(idvw,jdvw,nt))

! vhags, vtsgs=vhsgs
    lwork = max(3*nlat*(nlat+1)+2,(2*nt+1)*nlat*nlon)
! vrtgs
    lwork = max(lwork,nlat*((nt+1)*nlon+2*nt*l1+1))
! divgs
    l1 = min0(nlat,(nlon+2)/2)
    lwork = max(lwork,nlat*((nt+1)*nlon+2*nt*l1+1))
    allocate(work(lwork))

  end subroutine vectrans_init

  subroutine vectrans_clean
    use shtrans_module, only: shtrans_clean
    implicit none

    deallocate(wvhags,wvhsgs,wvts,work,br,bi,cr,ci,v,w,g,a,b)
    if (shtrans_isinit) then
      call shtrans_clean()
    end if

  end subroutine vectrans_clean

  subroutine vectrans_dlat(ui,vi,uy,vy,lfil)
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: ui, vi
    real(kind=dp), dimension(:,:,:), intent(inout) :: uy, vy
    logical, optional, intent(in) :: lfil

    real(kind=dp) :: fil
    integer(kind=i4b) :: i, j, jr, k, n, m, nt, ierror

    nt = size(ui,3)
    call shtrans_transf(vi,v,lyrev)
    call shtrans_transf(ui,w,lyrev)
    v(:,:,:) = -v(:,:,:)
    call vhags(nlat,nlon,isym,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvhags,size(wvhags),work,size(work),ierror)
    call checkerror(ierror,"vhags")

    call shtrans_truncate(br,ntrunc)
    call shtrans_truncate(bi,ntrunc)
    call shtrans_truncate(cr,ntrunc)
    call shtrans_truncate(ci,ntrunc)

    if (present(lfil).and.lfil) then
      call shtrans_hoskinsfilter(br,ntrunc)
      call shtrans_hoskinsfilter(bi,ntrunc)
      call shtrans_hoskinsfilter(cr,ntrunc)
      call shtrans_hoskinsfilter(ci,ntrunc)
    end if

    call vtsgs(nlat,nlon,isym,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvts,size(wvts),work,size(work),ierror)
    call checkerror(ierror,"vtsgs")
    call shtrans_transb(w,uy,lyrev)
    call shtrans_transb(v,vy,lyrev)

    if (vectrans_verbose) then
      call checkval(ui,"ui")
      call checkval(vi,"vi")
      call checkval(br,"br")
      call checkval(bi,"bi")
      call checkval(cr,"cr")
      call checkval(ci,"ci")
      call checkval(uy,"uy")
      call checkval(vy,"vy")
    end if

  end subroutine vectrans_dlat

  subroutine vectrans_uv2vrdv(ui,vi,vr,dv,lfil)
    use shtrans_module, only: wshsgs=>shtrans_wshsgs
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: ui, vi
    real(kind=dp), dimension(:,:,:), intent(inout) :: vr, dv
    logical, optional, intent(in) :: lfil

    real(kind=dp) :: fil
    integer(kind=i4b) :: i, j, jr, k, n, m, nx, ny, nt, ierror
    
    nx = size(ui,1)
    ny = size(ui,2)
    nt = size(ui,3)
    if (.not.shtrans_isinit) then
      shtrans_verbose = vectrans_verbose
      call shtrans_init(nx,ny,nt,.not.lyrev)
    end if

    call shtrans_transf(ui,w,lyrev)
    call shtrans_transf(vi,v,lyrev)
    v(:,:,:) = -v(:,:,:)
    call vhags(nlat,nlon,ityp,nt,v,w,idvw,jdvw,br,bi,cr,ci, &
               mdab,ndab,wvhags,size(wvhags),work,size(work),ierror)
    call checkerror(ierror,"vhags")

    call shtrans_truncate(br,ntrunc)
    call shtrans_truncate(bi,ntrunc)
    call shtrans_truncate(cr,ntrunc)
    call shtrans_truncate(ci,ntrunc)

    if (present(lfil).and.lfil) then
      call shtrans_hoskinsfilter(br,ntrunc)
      call shtrans_hoskinsfilter(bi,ntrunc)
      call shtrans_hoskinsfilter(cr,ntrunc)
      call shtrans_hoskinsfilter(ci,ntrunc)
    end if

    call vrtgs(nlat,nlon,isym,nt,v,idvw,jdvw,cr,ci, &
               mdab,ndab,wshsgs,size(wshsgs),work,size(work),ierror)
    call checkerror(ierror,"vrtgs")
    call shtrans_transb(v,vr,lyrev)
    call divgs(nlat,nlon,isym,nt,w,idvw,jdvw,br,bi, &
               mdab,ndab,wshsgs,size(wshsgs),work,size(work),ierror)
    call checkerror(ierror,"divgs")
    call shtrans_transb(w,dv,lyrev)

    if (vectrans_verbose) then
      call checkval(ui,"ui")
      call checkval(vi,"vi")
      call checkval(br,"br")
      call checkval(bi,"bi")
      call checkval(cr,"cr")
      call checkval(ci,"ci")
      call checkval(vr,"vr")
      call checkval(dv,"dv")
    end if

  end subroutine vectrans_uv2vrdv

  subroutine vectrans_gradient(gin,gx,gy,lfil)
    use shtrans_module, only: wshags=>shtrans_wshags
    implicit none

    real(kind=dp), dimension(:,:,:), intent(in) :: gin
    real(kind=dp), dimension(:,:,:), intent(inout) :: gx, gy
    logical, optional, intent(in) :: lfil

    integer(kind=i4b) :: nx, ny, nt, ierror

    nx = size(gin,1)
    ny = size(gin,2)
    nt = size(gin,3)
    if (.not.shtrans_isinit) then
      shtrans_verbose = vectrans_verbose
      call shtrans_init(nx,ny,nt,.not.lyrev)
    end if
    call shtrans_transf(gin,g,lyrev)
    call shags(nlat,nlon,isym,nt,g,idvw,jdvw,a,b,mdab,ndab, &
               wshags,size(wshags),work,size(work),ierror)
    call checkerror(ierror,"shags")

    call shtrans_truncate(a,ntrunc)
    call shtrans_truncate(b,ntrunc)
    if (present(lfil).and.lfil) then
      call shtrans_hoskinsfilter(a,ntrunc)
      call shtrans_hoskinsfilter(b,ntrunc)
    end if
    call gradgs(nlat,nlon,isym,nt,v,w,idvw,jdvw,a,b,mdab,ndab, &
                wvhsgs,size(wvhsgs),work,size(work),ierror)
    call checkerror(ierror,"gradgs")

    v(:,:,:) = -v(:,:,:)
    call shtrans_transb(v,gy,lyrev)
    call shtrans_transb(w,gx,lyrev)

    if (vectrans_verbose) then
      call checkval(gin,"gin")
      call checkval(a,"a")
      call checkval(b,"b")
      call checkval(gx,"gx")
      call checkval(gy,"gy")
    end if

  end subroutine vectrans_gradient

  subroutine checkerror(ierror,name)
    implicit none

    integer(kind=i4b), intent(in) :: ierror
    character(len=*), intent(in) :: name

    if (ierror/=0) then
      print *, "*** error in ",name,": ", ierror
      stop
    end if

  end subroutine checkerror

  subroutine checkval(x,name)
    implicit none

    real(kind=dp), dimension(:,:,:) :: x
    character(len=*), intent(in) :: name

    print *, name, " min=", minval(x), " at ", minloc(x), &
             name, " max=", maxval(x), " at ", maxloc(x)

  end subroutine checkval

end module vectrans_module
