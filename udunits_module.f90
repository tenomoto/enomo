module udunits_module
  use kind_module, only : i4b, i8b, sp, dp
  use time_module, only : time_type
  implicit none
  private

  integer(kind=i8b), private :: unitptr
  logical, private :: linit = .false. 

  public :: udunits_init, udunits_clean, udunits_hours, udunits_dayofyear, udunits_utc

contains

  subroutine udunits_init()
    implicit none

    integer(kind=i4b) :: utopen

    if (utopen("").eq.0) then
      linit = .true.
    else
      linit = .false.
      print *, "### error in udunits_init."
      print *, "Try setting UDUNITS_PATH"
    end if

  end subroutine udunits_init

  subroutine  udunits_utc(h, time)
    implicit none

    real(kind=dp), intent(in) :: h
    type(time_type), intent(inout) :: time

    character*80, parameter :: unitstr = "hours since 1-1-1 00:00:0.0"

    integer(kind=i4b) :: utdec, utcaltime 
    integer(kind=i8b) :: utmake

    if (.not.linit) then
      call udunits_init()
    end if

    unitptr = utmake()
    if (utdec(unitstr, unitptr)==0) then
      if (utcaltime(h, unitptr, time % year, time % month, time % day, &
        time % hour, time % minute, time % second)/=0) then
        print *, "### error in utcaltime()"
        print *, "setting time = 0"
        time = time_type(0,0,0,0,0,0.0)
      end if
    else
      print *, "### error in utdec()"
      print *, "setting time = 0"
      time = time_type(0,0,0,0,0,0.0)
    end if

  end subroutine udunits_utc

  subroutine  udunits_hours(time, h)
    implicit none

    type(time_type), intent(in) :: time
    real(kind=dp), intent(out) :: h

    character*80, parameter :: unitstr = "hours since 1-1-1 00:00:0.0"

    integer(kind=i4b) :: utdec, uticaltime 
    integer(kind=i8b) :: utmake

    if (.not.linit) then
      call udunits_init()
    end if

    unitptr = utmake()
    if (utdec(unitstr, unitptr)==0) then
      if (uticaltime(time % year, time % month, time % day, &
        time % hour, time % minute, time % second, unitptr, h)/=0) then
        print *, "### error in uticaltime()"
        print *, "setting h = 0"
        h = 0.0d0
      end if
    else
      print *, "### error in utdec()"
      print *, "setting h = 0"
      h = 0.0d0
    end if

  end subroutine udunits_hours

  subroutine  udunits_dayofyear(time, day)
    implicit none

    type(time_type), intent(in) :: time
    real(kind=dp), intent(out) :: day

    real(kind=sp) :: sc
    integer(kind=i4b) :: utdec, uticaltime 
    integer(kind=i8b) :: utmake

    character(len=*), parameter :: &
      unitprefix = "days since ", unitpostfix="-1-1 00:00:0.0"

    character(len=4) :: ystr
    character(len=80) :: unitstr

    write(ystr, "(i4)") time % year
    unitstr = unitprefix//ystr//unitpostfix
!   print *, "ystr=",ystr, " time % year=",time % year, " unitstr=",unitstr
  
    if (.not.linit) then
      call udunits_init()
    end if

    unitptr = utmake()
    if (utdec(unitstr, unitptr)==0) then
      if (uticaltime(time % year, time % month, time % day, &
        time % hour, time % minute, time % second, unitptr, day)/=0) then
        print *, "### error in uticaltime()"
        print *, "setting day = 0"
        day = 0.0d0
      end if
    else
      print *, "### error in utdec()"
      print *, "setting day = 0"
      day = 0.0d0
    end if

  end subroutine udunits_dayofyear

  subroutine udunits_clean()
    implicit none

    if (linit) then
      call utfree(unitptr)
      call utcls()
    end if

  end subroutine udunits_clean

end module udunits_module

