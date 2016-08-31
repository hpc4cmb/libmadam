MODULE  tod_storage
  !
  ! Allocate storage space for the TOD
  !
  use iso_c_binding

  use commonparam
  use mpi_wrappers
  use memory_and_time, only : check_stat

  implicit none
  private

  real(c_double), pointer, public :: tod_stored(:,:)
  real(c_double), pointer, public :: tod(:)

  real(sp),save,public :: memory_tod = 0.0

  real(c_double), pointer, public :: sampletime(:)
  integer, allocatable, target, public :: qualityFlag(:,:)

  real(dp), pointer, public :: time(:)

  character(len=40), parameter :: mstr='(x,a,t32,f9.1," MB")',  mstr3='(x,a,t32,3(f10.1," MB"))'

  public allocate_tod, free_tod, free_sampletime

CONTAINS

  !---------------------------------------------------------------------------


  SUBROUTINE allocate_tod

    integer  :: allocstat
    real(sp) :: memsum, mem_min, mem_max

    !allocate(tod_stored(nosamples_proc,nodetectors),stat=allocstat)
    !call check_stat(allocstat)
    !tod_stored = 0.0
    memory_tod = nosamples_proc*nodetectors*8.

    !if ( .not. allocated( sampletime ) ) then
    !   allocate(sampletime(nosamples_proc), stat=allocstat)
    !   call check_stat(allocstat)
    !   sampletime = 0.0
    !end if
    memory_tod = memory_tod + nosamples_proc*8.

    allocate(surveyflags(nosamples_proc), stat=allocstat)
    call check_stat(allocstat, 'surveyflags')

    memory_tod = memory_tod + nosamples_proc
    surveyflags = .true.

    tod => NULL()
    time => NULL()

    memsum = memory_tod/1024./1024.

    mem_min = memsum; mem_max = memsum ! -RK
    call min_mpi(mem_min); call max_mpi(mem_max) ! -RK

    call sum_mpi(memsum)

    if (id == 0 .and. info >= 1) &
         write(*,mstr3) 'Allocated memory for TOD:',memsum, mem_min, mem_max

  END SUBROUTINE allocate_tod


  !---------------------------------------------------------------------------


  SUBROUTINE free_tod()
    
    !if (allocated(tod_stored)) deallocate(tod_stored)
    if (allocated(surveyflags)) deallocate(surveyflags)

    tod =>NULL()
    time =>NULL()

  END SUBROUTINE free_tod


  !---------------------------------------------------------------------------


  SUBROUTINE free_sampletime()

    !if (allocated(sampletime)) deallocate(sampletime)

  END SUBROUTINE free_sampletime


  !---------------------------------------------------------------------------

END MODULE tod_storage
