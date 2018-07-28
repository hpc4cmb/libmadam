! Small utility module for covariance matrices
! 2007-09-07 RK

MODULE covmat_util

  USE planck_config, ONLY : i4b, i8b, sp, dp
  !USE healpix_types

  use mpi_wrappers, only : mpi_wtime

  IMPLICIT NONE

  REAL(dp), private :: time_start, time_end, tic_times(1000)
  INTEGER, private :: ID = 0


CONTAINS



  FUNCTION wtime()
    REAL(dp) :: wtime
    !INTEGER :: count, rate
    !call cpu_time(wtime)
    !call system_clock(count, rate)
    !wtime = dble(count)/rate

    wtime =  mpi_wtime()
  END FUNCTION wtime



  SUBROUTINE tic(index)
    INTEGER, OPTIONAL :: index
    
    IF (PRESENT(index)) THEN
       tic_times(index) = wtime()
    ELSE
       time_start = wtime()
    END IF

  END SUBROUTINE tic



  SUBROUTINE toc(label, index, threshold)
    CHARACTER(len=*), OPTIONAL :: label
    INTEGER, OPTIONAL :: index
    REAL(dp), OPTIONAL :: threshold
    CHARACTER(len=512) :: middle_string
    REAL(dp) :: elapsed, time_now
    REAL(dp), SAVE :: reporting_threshold=1e-3

    IF (PRESENT(threshold)) reporting_threshold = threshold

    time_now = wtime()

    IF (PRESENT(index)) THEN
       elapsed = time_now - tic_times(index)
    ELSE
       elapsed = time_now - time_start
    END IF

    IF (elapsed < reporting_threshold) RETURN

    IF (PRESENT(label)) THEN
       middle_string = ' : '//TRIM(label)//' completed in '
    ELSE
       middle_string = ' : elapsed time : '
    END IF

    WRITE (*,'(i3,a,f8.3,a)') &
         ID, TRIM(middle_string), elapsed, ' s'

  END SUBROUTINE toc


END MODULE covmat_util
