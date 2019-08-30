MODULE tod_storage
  !
  ! Allocate storage space for the TOD
  !
  use iso_c_binding

  use commonparam
  use mpi_wrappers
  use memory_and_time, only : check_stat, write_memory

  implicit none
  private

  real(c_double), pointer, public :: tod(:, :)
  real(c_double), pointer, public :: sampletime(:)

  real(dp), save, public :: memory_tod = 0

  character(len=40), parameter :: mstr='(x,a,t32,f9.1," MB")', &
       mstr3='(x,a,t32,3(f10.1," MB"))'

  public allocate_tod, free_tod

CONTAINS

  !---------------------------------------------------------------------------


  SUBROUTINE allocate_tod

    integer  :: allocstat
    real(sp) :: memsum, mem_min, mem_max

    ! Stored signal
    memory_tod = nosamples_proc * nodetectors * 8.

    ! Time stamps
    memory_tod = memory_tod + nosamples_proc * 8.

    allocate(surveyflags(nosamples_proc), stat=allocstat)
    call check_stat(allocstat, 'surveyflags')

    memory_tod = memory_tod + nosamples_proc
    surveyflags = .true.

    call write_memory("TOD memory", memory_tod)

  END SUBROUTINE allocate_tod


  !---------------------------------------------------------------------------


  SUBROUTINE free_tod()
    
    if (allocated(surveyflags)) deallocate(surveyflags)

  END SUBROUTINE free_tod


  !---------------------------------------------------------------------------

END MODULE tod_storage
