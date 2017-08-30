MODULE pointing
  !
  ! Routines for storing and handling of pointing data
  use iso_c_binding

  use commonparam
  use mpi_wrappers
  use memory_and_time, only : check_stat

  implicit none
  private

  logical, allocatable, public :: kpolarized(:)

  integer(c_long), pointer, public :: pixels(:, :)
  real(c_double), pointer, public :: weights(:, :, :)

  ! Next two are used for sub ring maps

  integer(i2b), allocatable, target, public :: subchunk(:)
  integer(i2b), pointer, public :: subchunkpp(:)

  integer :: buffersize = 0

  real(sp),public :: memory_pointing = 0.0

  logical, allocatable, public :: ksubmap(:)
  integer, allocatable, public :: subtable1(:), subtable2(:)

  integer, public :: dummy_pixel = -1

  public init_pointing, close_pointing, reduce_pixels, restore_pixels

CONTAINS

  !--------------------------------------------------------------------------


  SUBROUTINE init_pointing
    !
    !Initialize the pointing module and allocate memory for pointing data.
    !
    integer  :: idet, k, allocstat
    real(dp) :: psi
    real(sp) :: memory, mem_min, mem_max

    if (id == 0 .and. info > 3) write (*,'(a)') ' Initializing pointing'

    allocate(kpolarized(nodetectors), stat=allocstat)
    call check_stat(allocstat, 'kpolarized')

    if (temperature_only) detectors%kpolar = .false.

    allocate(subchunk(nosamples_proc), stat=allocstat)
    call check_stat(allocstat, 'subchunk')
    subchunk = 0 ! initialize

    memory_pointing = nosamples_proc*nodetectors*4.
    memory_pointing = memory_pointing + nosamples_proc*nodetectors*24.

    allocate(ksubmap(0:nosubmaps_tot))
    allocate(subtable1(0:nosubmaps_tot))
    allocate(subtable2(0:nosubmaps_tot))
    ksubmap  = .true.
    subtable1 = 0
    subtable2 = 0

    dummy_pixel = 12*nside_max**2

    memory = memory_pointing / 2d0**20

    mem_min = memory; mem_max = memory
    call min_mpi(mem_min); call max_mpi(mem_max)

    call sum_mpi(memory)

    if (id == 0 .and. info > 0) write(*,'(a,t32,3(f12.1," MB"))')   &
         ' Allocated memory for pointing:', memory, mem_min, mem_max

  END SUBROUTINE init_pointing


  !------------------------------------------------------------------------


  SUBROUTINE close_pointing

    integer :: k

    if (allocated(kpolarized)) deallocate(kpolarized)

    if (allocated(ksubmap)) deallocate(ksubmap)
    if (allocated(subtable1)) deallocate(subtable1)
    if (allocated(subtable2)) deallocate(subtable2)

    if (allocated(subchunk)) deallocate(subchunk)

    ! Free various arrays not directly associated with pointing

    if (allocated(basis_functions)) then
       do k = 1, noba_short
          if (.not. basis_functions(k)%copy) deallocate(basis_functions(k)%arr)
       end do
       deallocate(basis_functions)
    end if

    if (allocated(id_submap)) deallocate(id_submap)

    if (allocated(detectors)) deallocate(detectors)

    buffersize = 0

  END SUBROUTINE close_pointing


  !------------------------------------------------------------------------


  SUBROUTINE reduce_pixels
    ! Reduce pixel numbers so that they point to locmap
    !
    integer :: i, idet, k, ip, ierr

    if (info > 4) write(*,idf) id, 'Reduce pixel numbers...'

    ksubmap = .false.

    do idet = 1, nodetectors
       do i = 1, nosamples_proc
          if (isubchunk /= 0 .and. subchunkpp(i) /= isubchunk) cycle
          ip = pixels(i, idet) / nosubpix_max
          ksubmap(ip) = .true.
       end do
    end do
    ksubmap(nosubmaps_tot) = .true.

    if (allreduce) then
       ! Flag all hit submaps on every process
       call mpi_allreduce(MPI_IN_PLACE, ksubmap, nosubmaps_tot, MPI_LOGICAL, &
            MPI_LOR, comm, ierr)
       if (ierr /= 0) call abort_mpi('Reducing ksubmaps failed.')
    end if

    subtable1 = -1
    subtable2 = -1
    k = -1
    do i = 0, nosubmaps_tot
       if (ksubmap(i)) then
          k = k+1
          subtable1(i) = i-k
          subtable2(k) = i-k
       end if
    end do

    do idet = 1, nodetectors
       do i = 1, nosamples_proc
          if (isubchunk /= 0 .and. subchunkpp(i) /= isubchunk) cycle
          ip = pixels(i, idet)/nosubpix_max
          pixels(i, idet) = pixels(i, idet) - subtable1(ip)*nosubpix_max
       end do
    end do

    dummy_pixel = (nosubmaps_tot-subtable1(nosubmaps_tot)) * nosubpix_max

    if (info > 4) write(*,idf) id, 'Done'

  END SUBROUTINE reduce_pixels


  !---------------------------------------------------------------------------

  SUBROUTINE restore_pixels
    ! restore original pixel numbers
    !
    integer :: i, idet, ip

    if (info.ge.5) write(*,idf) id, 'Restore pixel numbers (a)...'

    do idet = 1, nodetectors
       do i = 1, nosamples_proc
          if (isubchunk /= 0 .and. subchunkpp(i) /= isubchunk) cycle
          ip = pixels(i, idet) / nosubpix_max
          pixels(i, idet) = pixels(i, idet) + subtable2(ip)*nosubpix_max
       end do
    end do

    ksubmap = .true.
    subtable1 = 0
    subtable2 = 0
    dummy_pixel = 12*nside_max**2

    if (info > 4) write(*,idf) id, 'Done'

  END SUBROUTINE restore_pixels


  !---------------------------------------------------------------------------

END MODULE pointing
