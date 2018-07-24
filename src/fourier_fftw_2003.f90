MODULE fourier

  use, intrinsic :: iso_c_binding

  implicit none

  private

  include "fftw3.f03"

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)

  integer :: nof = 0, info = 0
  real(dp) :: nofinv

  type(C_PTR) :: plan = C_NULL_PTR, plan_inv = C_NULL_PTR
  integer :: ierr

  public init_fourier, close_fourier, dfft, dfftinv

CONTAINS

  SUBROUTINE init_fourier(nof_in)

    integer, intent(in) :: nof_in
    integer :: fftw_planning_strategy
    real(C_DOUBLE), allocatable :: in(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: out(:)

    if (nof_in == 0 .or. nof_in == nof) return

    if (nof > 0) call close_fourier()

    nof = nof_in
    nofinv = 1 / dble(nof)

    allocate(in(nof), out(nof/2+1), stat=ierr)
    if (ierr /= 0) then
       stop 'Unable to allocate working space for Fourier transform'
    end if

    ierr = fftw_import_system_wisdom()

    if (info /= 0) then
       if (ierr /= 0) then
          write (*,*) 'FFTW: System wide wisdom loaded'
       else
          write (*,*) 'FFTW: Unable to load system wisdom'
       end if
    end if

    if (nof <= 65536) then
       fftw_planning_strategy = FFTW_MEASURE
    else
       fftw_planning_strategy = FFTW_ESTIMATE
    end if

    fftw_planning_strategy = ior(fftw_planning_strategy, FFTW_PRESERVE_INPUT)
    fftw_planning_strategy = ior(fftw_planning_strategy, FFTW_UNALIGNED)

    plan = fftw_plan_dft_r2c_1d(nof, in, out, fftw_planning_strategy)
    plan_inv = fftw_plan_dft_c2r_1d(nof, out, in, fftw_planning_strategy)

    deallocate(in, out)

  END SUBROUTINE init_fourier

  !-----------------------------------------------------------------------

  SUBROUTINE close_fourier

    ! Destroying the FFTW plans was disabled for a while, after
    ! observing inexplicable segfaults in some environments.

    if (c_associated(plan)) then
       call fftw_destroy_plan(plan)
       plan = C_NULL_PTR
    end if

    if (c_associated(plan_inv)) then
       call fftw_destroy_plan(plan_inv)
       plan_inv = C_NULL_PTR
    end if

    nof = 0

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

    call fftw_execute_dft_c2r(plan_inv, fyy, yy)

    yy = yy * nofinv

  END SUBROUTINE dfftinv

END MODULE fourier
