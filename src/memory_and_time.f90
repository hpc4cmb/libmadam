MODULE memory_and_time

  use commonparam
  use mpi_wrappers

  implicit none
  private

  real(sp), public :: memory_total = 0.0
  real(sp) :: time1 =0.0, time2 =0.0
  real(sp), public :: time_cum =0.0
  real(sp), public :: time2_sum =0.0

  public write_memory, write_time, check_stat

CONTAINS

  !------------------------------------------------------------------------------


  SUBROUTINE check_stat(allocstat, name)

    integer :: allocstat
    character(len=*) :: name

    if (allocstat.ne.0) then
       write(*,*) 'ERROR: out of memory allocating ',name
       call exit_with_status(1)
    endif

  END SUBROUTINE check_stat


  !------------------------------------------------------------------------------


  SUBROUTINE write_memory(text,memory_in)

    character(len=*),intent(in)          :: text
    real,            intent(in),optional :: memory_in
    real                                 :: memory_sum

    if (text=='Total') then
       if (ID==0) then
          write(*, '(x,a,t26,f9.1," MB")') text, memory_total
          write(*, '(x,a,t26,f9.1," GB")') ' ',  memory_total/1024.
       endif
    else
       memory_sum = memory_in/1024./1024.
       call sum_mpi(memory_sum)

       memory_total = memory_total+memory_sum

       if (ID==0.and.memory_sum.gt.0.05)  &
            write(*, '(x,a,t26,f9.1)') text, memory_sum
    endif

  END SUBROUTINE write_memory


  !------------------------------------------------------------------------------


  SUBROUTINE write_time(text,time)

    character(len=*),intent(in)          :: text
    real,            intent(in),optional :: time

    real :: time_min, time_max ! -RK

    if (text=='- Other') then
       time2 = time1-time2_sum

       if (ID==0.and.time2.gt.0.1)    &
            write(*,'(3x,a,t32,f7.1)') text, time2

    elseif (text(1:1)=='-') then
       time2 = time

       time_min = time2; time_max = time2 ! -RK
       call min_mpi(time_min); call max_mpi(time_max) ! -RK

       call sum_mpi(time2)

       time2 = time2/ntasks
       time2_sum = time2_sum+time2

       !if (ID==0.and.time2.gt.0.05) write(*,'(3x,a,t32,f8.1)') text, time2
       if (ID==0.and.time2.gt.0.05) write(*,'(3x,a,t32,3f8.1)') & ! -RK
            text, time2, time_min, time_max ! -RK

    elseif (text=='Total') then
       time1 = time

       call sum_mpi(time1)
       time1 = time1/ntasks

       if (ID==0) write(*,'(x,a,t28,f8.1," s  (",f8.2," CPU hours)")')  &
            text, time1, time1*ntasks*nthreads_max/3600.

    else
       time1 = time

       time_min = time1; time_max = time1 ! -RK
       call min_mpi(time_min); call max_mpi(time_max) ! -RK

       time_cum = time_cum+time1

       call sum_mpi(time1)
       time1 = time1/ntasks
       time2_sum = 0.0

       !if (ID==0.and.time1.gt.0.05) write (*,'(x,a,t28,f8.1)') text, time1 -RK
       if (ID==0.and.time1.gt.0.05) write (*,'(x,a,t28,3f8.1)') &
            text, time1, time_min, time_max

    endif

  END SUBROUTINE write_time


  !-----------------------------------------------------------------------------

END MODULE memory_and_time
