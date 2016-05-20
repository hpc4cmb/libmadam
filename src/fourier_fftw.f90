MODULE fourier

  implicit none
  private

  include "fftw3.f"

  integer, parameter :: sp = kind(1.0)
  integer, parameter :: dp = kind(1.0d0)

  integer       :: nof, info=0
  real(dp), pointer     :: in(:)
  complex(dp), pointer  :: out(:)

  integer      :: plan(8), plan_inv(8), ierr

  ! integer               :: fftw_planning_strategy = FFTW_MEASURE
  integer               :: fftw_planning_strategy = FFTW_ESTIMATE

  public init_fourier, close_fourier, dfft, dfftinv

CONTAINS

  SUBROUTINE init_fourier(nof_in)

    integer :: nof_in

    nof = nof_in

    allocate(in(nof), out(nof/2+1), stat = ierr)

    if (ierr /= 0) then
       stop "Unable to allocate working space for Fourier transform"
    end if

    if (nof==0) return

    call dfftw_import_system_wisdom(ierr)

    if (info /= 0) then
       if (ierr /= 0) then
          write (*,*) 'FFTW: System wide wisdom loaded'
       else
          write (*,*) 'FFTW: Unable to load system wisdom'
       end if
    end if

    if (nof.le.65536) fftw_planning_strategy = FFTW_MEASURE

    call dfftw_plan_dft_r2c_1d(plan, nof_in, in, out, fftw_planning_strategy)
    call dfftw_plan_dft_c2r_1d(plan_inv, nof_in, out, in, fftw_planning_strategy)

  END SUBROUTINE init_fourier

  !-----------------------------------------------------------------------

  SUBROUTINE close_fourier

    deallocate(in, out)

  END SUBROUTINE close_fourier

  !-----------------------------------------------------------------------

  SUBROUTINE dfft(fyy,yy)
    ! Fourier transform of a real vector.

    complex(dp),intent(out) :: fyy(nof/2+1)
    real(dp),   intent(in)  :: yy(nof)
    integer                 :: i

    in  = yy

    call dfftw_execute(plan)

    fyy = out

  END SUBROUTINE dfft

  !-----------------------------------------------------------------------

  SUBROUTINE dfftinv(yy,fyy)
    ! Inverse Fourier transform complex -> real.

    real(dp),   intent(out) :: yy(nof)
    complex(dp),intent(in)  :: fyy(nof/2+1)
    integer                 :: i

    out = fyy

    call dfftw_execute(plan_inv)

    yy = in/nof

  END SUBROUTINE dfftinv

END MODULE fourier
