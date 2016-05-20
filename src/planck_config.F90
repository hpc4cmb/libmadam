!> Module containing non-standard, platform-dependent routines and
!> type definitions used by all LevelS modules.
!> The Fortran 2003 standard was emulated where possible
!> \author Martin Reinecke
module planck_config
#if defined (PLANCK_NAGF95)
  use f90_unix, only:exit
#endif
implicit none
public

  !> 1-byte signed integer kind
  integer, parameter, public :: i1b = selected_int_kind(2)
  !> 2-byte signed integer kind
  integer, parameter, public :: i2b = selected_int_kind(4)
  !> 4-byte signed integer kind
  integer, parameter, public :: i4b = selected_int_kind(9)
  !> 8-byte signed integer kind
  integer, parameter, public :: i8b = selected_int_kind(18)
  !> single precision real kind
  integer, parameter, public :: sp  = kind(1.0e0)
  !> double precision real kind
  integer, parameter, public :: dp  = kind(1.0d0)

  ! Numerical Constants (double precision)
  real(dp), parameter, public :: &
    pi        = 3.141592653589793238462643383279502884197_dp, &
    twopi     = pi*2._dp, &
    fourpi    = pi*4._dp, &
    halfpi    = pi/2._dp, &
    sqrt2     = 1.41421356237309504880168872420969807856967_dp, &
    twothird  = 2._dp/3._dp, &
    sec_58_70 = 378691200._dp ! seconds between Jan 1, 1958 and Jan 1, 1970

  ! useful parameters
  integer, parameter, public :: &
    filenamelen = 1024

contains

#if defined (__INTEL_COMPILER)
#if (__INTEL_COMPILER<810)
ERROR_COMPILER_TOO_OLD
#elif (__INTEL_COMPILER<900)
subroutine get_command_argument(number,value)
  use ifposix
  integer, intent(in) :: number
  character (len=*), intent(out) :: value
  integer l, err

  call pxfgetarg(number,value,l,err)
end subroutine
#endif
#endif

#if defined (PLANCK_NAGF95)
function command_argument_count ()
  use f90_unix
  integer command_argument_count

  command_argument_count=iargc()
end function

subroutine get_command_argument(number,value)
  use f90_unix
  integer, intent(in) :: number
  character (len=*), intent(out) :: value

  call getarg(number,value)
end subroutine
#endif

#if (defined(pgiFortran) || defined(__SUNPRO_F90) || defined(__IBMC__) || defined(PLANCK_F95F))
function command_argument_count ()
  integer command_argument_count
  interface
  function iargc()
    integer iargc
  end function
  end interface
  command_argument_count=iargc()
end function

subroutine get_command_argument(number,value)
  integer, intent(in) :: number
  character (len=*), intent(out) :: value

  call getarg(number,value)
end subroutine
#endif

#if defined (__IBMC__)
subroutine exit(code)
  integer, intent(in) :: code
  call exit_(code)
end subroutine
#endif

!> prints the error message in \a msg and tries to exit with a return status
!> of \a code. Note that this is not possible on all platforms.
subroutine exit_with_status (code, msg)
  integer, intent(in) :: code
  character (len=*), intent(in), optional :: msg

  if (present(msg)) print *,trim(msg)
  print *,'program exits with exit code ', code

  call exit (code)
end subroutine exit_with_status

end module planck_config
