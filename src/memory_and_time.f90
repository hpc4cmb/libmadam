MODULE memory_and_time

  use commonparam
  use mpi_wrappers

  implicit none
  private

  real(dp), public :: memory_total = 0
  real(dp) :: time1 = 0, time2 = 0
  real(dp), public :: time_cum = 0
  real(dp), public :: time2_sum = 0

  public write_memory, write_time, check_stat

CONTAINS

  !------------------------------------------------------------------------------


  SUBROUTINE check_stat(allocstat, name)

    integer :: allocstat
    character(len=*) :: name

    if (allocstat /= 0) then
       write(*,*) 'ERROR: out of memory allocating ',name
       call exit_with_status(1)
    endif

  END SUBROUTINE check_stat


  !------------------------------------------------------------------------------


  function mem_to_string(mem) result(memstring)
    ! Return a string representation of memory value
    real(dp) :: mem  ! memory in MB
    character(len=12) :: memstring
    character(len=2) :: unit
    
    if (mem > 2 ** 30) then
       mem = mem / 2 ** 30
       unit = 'PB'
    else if (mem > 2 ** 20) then
       mem = mem / 2 ** 20
       unit = 'TB'
    else if (mem > 2 ** 10) then
       mem = mem / 2 ** 10
       unit = 'GB'
    else if (mem > 1) then
       unit = 'MB'
    else if (mem > 1 / 2 ** 10) then
       mem = mem * 2 ** 10
       unit = 'kB'
    else
       mem = mem ** 2 ** 20
       unit = 'B '
    end if
    write(memstring, '(f9.2,x,a2)') mem, unit
  end function mem_to_string


  SUBROUTINE write_memory(text, memory_in)
    character(len=*), intent(in) :: text  ! prefix the memory string
    real(dp), intent(in), optional :: memory_in  ! memory in bytes
    real(dp) :: memory, memory_sum, memory_min, memory_max

    if (text == 'Total') then
       memory = memory_total
    else
       memory = memory_in / 2 ** 20
       memory_total = memory_total + memory
    end if

    memory_min = memory
    memory_max = memory
    memory_sum = memory
    call min_mpi(memory_min)
    call max_mpi(memory_max)
    call sum_mpi(memory_sum)
    if (ID == 0 .and. info > 0) then
       write(*, '(x,a,t28,"min = ",a,"   max =",a,"   total =",a)') text, &
            mem_to_string(memory_min), mem_to_string(memory_max), &
            mem_to_string(memory_sum)
    end if

  END SUBROUTINE write_memory


  !------------------------------------------------------------------------------


  SUBROUTINE write_time(text, time)

    character(len=*), intent(in) :: text
    real(dp), intent(in), optional :: time

    real(dp) :: time_min, time_max

    if (text(1:1) == '-') then
       if (text == '- Other') then
          ! Difference between top level timer and the sum of the sub timers
          time2 = time1 - time2_sum
       else
          ! sub timer
          time2 = time
       end if
       time_min = time2
       time_max = time2
       call min_mpi(time_min)
       call max_mpi(time_max)
       call sum_mpi(time2)
       time2 = time2 / ntasks
       time2_sum = time2_sum + time2

       if (ID == 0 .and. time2 > 0.05) then
          write (*,'(3x,a,t34,"mean =",f8.1,"  min =",f8.1,"  max =",f8.1)') &
               text, time2, time_min, time_max
       end if
    else if (text == 'Total') then
       ! zero level timer
       time1 = time

       call sum_mpi(time1)
       time1 = time1 / ntasks

       if (ID == 0) then
          write(*,'(x,a,t28,f8.1," s  (",f12.2," CPU hours)")') &
               text, time1, time1 * ntasks * nthreads_max / 3600.
       end if
    else
       ! New top level timer, possibly followed by sub timers
       time1 = time
       time_min = time1
       time_max = time1
       call min_mpi(time_min)
       call max_mpi(time_max)

       time_cum = time_cum + time1

       call sum_mpi(time1)
       time1 = time1 / ntasks
       time2_sum = 0

       if (ID == 0 .and. time1 > 0.05) then
          write (*,'(x,a,t34,"mean =",f8.1,"  min =",f8.1,"  max =",f8.1)') &
               text, time1, time_min, time_max
       end if
    end if

  END SUBROUTINE write_time


  !-----------------------------------------------------------------------------

END MODULE memory_and_time
