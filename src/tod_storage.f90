MODULE tod_storage
  !
  ! Allocate storage space for the TOD
  !
  use iso_c_binding

  use commonparam
  use mpi_wrappers
  use memory_and_time, only : check_stat

  implicit none
  private

  real(c_double), pointer, public :: tod(:, :)

  real(sp), save, public :: memory_tod = 0

  real(c_double), pointer, public :: sampletime(:)

  character(len=40), parameter :: mstr='(x,a,t32,f9.1," MB")', &
       mstr3='(x,a,t32,3(f10.1," MB"))'

  public allocate_tod, free_tod

CONTAINS

  !---------------------------------------------------------------------------


  SUBROUTINE allocate_tod

    integer  :: allocstat
    real(sp) :: memsum, mem_min, mem_max

    ! Stored signal
    memory_tod = nosamples_proc*nodetectors*8.

    ! Time stamps
    memory_tod = memory_tod + nosamples_proc*8.

    allocate(surveyflags(nosamples_proc), stat=allocstat)
    call check_stat(allocstat, 'surveyflags')

    memory_tod = memory_tod + nosamples_proc
    surveyflags = .true.

    memsum = memory_tod/1024./1024.

    mem_min = memsum; mem_max = memsum
    call min_mpi(mem_min); call max_mpi(mem_max)
    call sum_mpi(memsum)

    if (id == 0 .and. info >= 1) &
         write(*,mstr3) 'Allocated memory for TOD:', memsum, mem_min, mem_max

  END SUBROUTINE allocate_tod


  !---------------------------------------------------------------------------


  SUBROUTINE free_tod()
    
    if (allocated(surveyflags)) deallocate(surveyflags)

  END SUBROUTINE free_tod


  !---------------------------------------------------------------------------

END MODULE tod_storage
