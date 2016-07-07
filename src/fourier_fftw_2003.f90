MODULE fourier

  use, intrinsic :: iso_c_binding

  implicit none

  private

  include "fftw3.f03"

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: nof, info=0
  real(C_DOUBLE), pointer :: in(:) => NULL()
  complex(C_DOUBLE_COMPLEX), pointer :: out(:) => NULL()
  complex(C_DOUBLE_COMPLEX), allocatable :: fyy_backup(:)

  type(C_PTR) :: plan = C_NULL_PTR, plan_inv = C_NULL_PTR
  integer ::  ierr

  ! integer               :: fftw_planning_strategy = FFTW_MEASURE
  integer               :: fftw_planning_strategy = FFTW_ESTIMATE

  public init_fourier, close_fourier, dfft, dfftinv

CONTAINS

  SUBROUTINE init_fourier(nof_in)

    integer :: nof_in
    type(C_PTR) :: pin, pout

    if (nof_in == 0) return

    nof = nof_in

    allocate(fyy_backup(nof/2+1), stat=ierr)
    if (ierr /= 0) stop 'No room for fyy_backup'
    
    pin = fftw_alloc_real( int(nof, C_SIZE_T) )
    call c_f_pointer( pin, in, [nof] )

    pout = fftw_alloc_complex( int(nof/2+1, C_SIZE_T) )
    call c_f_pointer( pout, out, [nof/2+1] )

    if ( .not. associated(in) .or. .not. associated(out) ) then
       stop "Unable to allocate working space for Fourier transform"
    end if

    ierr = fftw_import_system_wisdom()

    if (info /= 0) then
       if (ierr /= 0) then
          write (*,*) 'FFTW: System wide wisdom loaded'
       else
          write (*,*) 'FFTW: Unable to load system wisdom'
       end if
    end if

    if (nof.le.65536) fftw_planning_strategy = FFTW_MEASURE

    plan = fftw_plan_dft_r2c_1d(nof, in, out, ior(fftw_planning_strategy,FFTW_PRESERVE_INPUT))
    plan_inv = fftw_plan_dft_c2r_1d(nof, out, in, ior(fftw_planning_strategy,FFTW_DESTROY_INPUT))

    call fftw_free(pin)
    call fftw_free(pout)
    in => NULL()
    out => NULL()

  END SUBROUTINE init_fourier

  !-----------------------------------------------------------------------

  SUBROUTINE close_fourier

    if (c_associated(plan)) then
       call fftw_destroy_plan( plan )
       plan = C_NULL_PTR
    end if
    if (c_associated(plan_inv)) then
       call fftw_destroy_plan( plan_inv )
       plan_inv = C_NULL_PTR
    end if

    if (allocated(fyy_backup)) deallocate( fyy_backup )

  END SUBROUTINE close_fourier

  !-----------------------------------------------------------------------

  SUBROUTINE dfft(fyy, yy)
    ! Fourier transform of a real vector.

    complex(C_DOUBLE_COMPLEX), intent(out) :: fyy(nof/2+1)
    real(C_DOUBLE) :: yy(nof)

    call fftw_execute_dft_r2c(plan, yy, fyy)

  END SUBROUTINE dfft

  !-----------------------------------------------------------------------

  SUBROUTINE dfftinv(yy, fyy)
    ! Inverse Fourier transform complex -> real.

    real(C_DOUBLE), intent(out) :: yy(nof)
    complex(C_DOUBLE_COMPLEX) :: fyy(nof/2+1)

    ! The plan destroys the input

    fyy_backup = fyy

    call fftw_execute_dft_c2r(plan_inv, fyy, yy)
    
    fyy = fyy_backup

    yy = yy/nof

  END SUBROUTINE dfftinv

END MODULE fourier
