MODULE pointing
  !
  ! Routines for storing and handling of pointing data
  use iso_c_binding

  use commonparam
  use mpi_wrappers

  implicit none
  private

  real(dp),allocatable,public :: cospsi(:),   &
       sinpsi(:)

  logical, allocatable, public :: kpolarized(:)

  integer, allocatable, public :: detbits(:)

  integer(c_long), pointer, public :: pixels(:,:)
  real(c_double), pointer, public :: weights(:,:,:)


  ! Next two are used for sub ring maps

  integer(i2b), allocatable,target,public :: subchunk(:)
  integer(i2b), pointer, public :: subchunkpp(:)

  integer :: buffersize = 0

  integer, allocatable,public :: ipix(:)
  logical, allocatable,public :: kpix(:)
  real(dp),allocatable,public :: pweight(:,:)

  real(sp),public :: memory_pointing = 0.0

  logical,allocatable,public :: ksubmap(:)
  integer,allocatable, public :: subtable1(:), subtable2(:)

  integer, public :: dummy_pixel = -1

  public init_pointing,      &
       close_pointing,       &
       allocate_pixbuffer

  public reduce_pixels_buff, &
       restore_pixels_buff,  &
       reduce_pixels_a,      &
       restore_pixels_a

CONTAINS

  !--------------------------------------------------------------------------


  SUBROUTINE init_pointing
    !
    !Initialize the pointing module and allocate memory for pointing data.
    !
    integer  :: idet, k, allocstat, detbit
    real(dp) :: psi
    real(sp) :: memory, mem_min, mem_max

    if (id == 0 .and. info > 3) write (*,'(a)') ' Initializing pointing'

    allocate(cospsi(nodetectors), sinpsi(nodetectors), &
         kpolarized(nodetectors), detbits(nodetectors), stat=allocstat)
    call check_stat(allocstat)

    if (temperature_only) detectors%kpolar = .false.

    detbit = 1
    do idet = 1,nodetectors
       detbits(idet) = detbit
       detbit = 2*detbit
    enddo

    allocate(subchunk(nosamples_proc), stat=allocstat)
    call check_stat(allocstat)
    subchunk = 0 ! initialize

    !allocate(pixels(nosamples_proc, nodetectors), stat=allocstat)
    !call check_stat(allocstat)
    memory_pointing = nosamples_proc*nodetectors*4.
    !pixels = 0

    !if (nmap /= 1) then
    !   allocate(weights(nmap, nosamples_proc, nodetectors), stat=allocstat)
    !   call check_stat(allocstat)
    memory_pointing = memory_pointing + nosamples_proc*nodetectors*24.
    !   weights = 0
    !endif

    call allocate_pixbuffer(1024)

    allocate(ksubmap(0:nosubmaps_tot))
    allocate(subtable1(0:nosubmaps_tot))
    allocate(subtable2(0:nosubmaps_tot))
    ksubmap  = .true.
    subtable1 = 0
    subtable2 = 0

    dummy_pixel = 12*nside_max**2

    memory = memory_pointing/1024./1024.

    mem_min = memory; mem_max = memory ! -RK
    call min_mpi(mem_min); call max_mpi(mem_max) ! -RK

    call sum_mpi(memory)

    if (id == 0 .and. info > 0) write(*,'(a,t32,3(f12.1," MB"))')   &
         ' Allocated memory for pointing:', memory, mem_min, mem_max

  END SUBROUTINE init_pointing


  !------------------------------------------------------------------------


  SUBROUTINE close_pointing

    integer :: k

    if (allocated(cospsi)) deallocate(cospsi, sinpsi, kpolarized, detbits)

    if (allocated(ipix))      deallocate(ipix)
    if (allocated(kpix))      deallocate(kpix)
    if (allocated(pweight)) deallocate(pweight)

    if (allocated(ksubmap))   deallocate(ksubmap)
    if (allocated(subtable1)) deallocate(subtable1)
    if (allocated(subtable2)) deallocate(subtable2)

    if (allocated(subchunk)) deallocate(subchunk)

    ! Free various arrays not directly associated with pointing

    if (allocated(basis_functions)) then
       do k = 1, noba_short
          if ( .not. basis_functions(k)%copy ) deallocate(basis_functions(k)%arr)
       end do
       deallocate(basis_functions)
    end if

    if (allocated(id_submap)) deallocate(id_submap)

    if (allocated(detectors)) deallocate(detectors)

    buffersize = 0

  END SUBROUTINE close_pointing


  !------------------------------------------------------------------------


  SUBROUTINE allocate_pixbuffer(n)

    integer, intent(in) :: n

    if (n.le.buffersize) return

    if (allocated(ipix))    deallocate(ipix)
    if (allocated(kpix))    deallocate(kpix)
    if (allocated(pweight)) deallocate(pweight)

    buffersize = 1
    do
       if (buffersize.ge.n) exit
       buffersize = 2*buffersize
    enddo

    allocate(ipix(buffersize))
    allocate(kpix(buffersize))
    allocate(pweight(nmap,buffersize))

    ipix = -1
    kpix = .true.
    pweight = 0.0
    pweight(1,:) = 1.d0

  END SUBROUTINE allocate_pixbuffer


  !--------------------------------------------------------------------------


  SUBROUTINE reduce_pixels_buff
    ! Reduce pixel numbers so that they point to locmap
    !
    integer :: i, k, ip, idet

    if (info > 4) write(*,idf) id,'Reduce pixel numbers...'

    ksubmap = .false.

    do idet = 1,nodetectors
       if ( .not. detflags(idet) ) cycle
       do i = 1,nosamples_proc
          if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
          if ( .not. surveyflags(i) ) cycle
          ip = pixels(i,idet) / nosubpix_max
          ksubmap(ip) = .true.
       end do
    end do
    ksubmap(nosubmaps_tot) = .true.  ! always allocate the dummy pixel

    subtable1 = -1
    subtable2 = -1
    k = -1
    do i = 0,nosubmaps_tot
       if (ksubmap(i)) then
          k = k + 1
          subtable1(i) = i - k
          subtable2(k) = i - k
       endif
    end do

    do idet = 1,nodetectors
       if ( .not. detflags(idet) ) cycle
       do i = 1,nosamples_proc
          if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
          if ( .not. surveyflags(i) ) cycle
          ip = pixels(i,idet) / nosubpix_max
          pixels(i,idet) = pixels(i,idet) - subtable1(ip)*nosubpix_max
       end do
    end do

    dummy_pixel = (nosubmaps_tot - subtable1(nosubmaps_tot)) * nosubpix_max

    if (info > 4) write(*,idf) id,'Done'

  END SUBROUTINE reduce_pixels_buff


  !-------------------------------------------------------------------------


  SUBROUTINE restore_pixels_buff
    ! restore original pixel numbers
    !
    integer :: i, ip, idet

    if (info > 4) write(*,idf) id,'Restore pixel numbers...'

    do idet = 1,nodetectors
       if ( .not. detflags(idet) ) cycle
       do i = 1,nosamples_proc
          if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
          if ( .not. surveyflags(i) ) cycle
          ip = pixels(i,idet) / nosubpix_max
          pixels(i,idet) = pixels(i,idet) + subtable2(ip)*nosubpix_max            
       end do
    end do

    ksubmap = .true.
    subtable1 = 0
    subtable2 = 0
    dummy_pixel = 12*nside_max**2

    if (info > 4) write(*, idf) ID,'Done'

  END SUBROUTINE restore_pixels_buff


  !------------------------------------------------------------------------


  SUBROUTINE reduce_pixels_a
    ! Reduce pixel numbers so that they point to locmap
    !
    integer :: i, idet, k, ip

    if (info.ge.5) write(*,idf) ID,'Reduce pixel numbers...'

    ksubmap = .false.

    do idet = 1, nodetectors
       do i = 1,nosamples_proc
          if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
          ip = pixels(i,idet) / nosubpix_max
          ksubmap(ip) = .true.
       enddo
    enddo
    ksubmap(nosubmaps_tot) = .true.

    subtable1 = -1
    subtable2 = -1
    k = -1
    do i = 0,nosubmaps_tot
       if (ksubmap(i)) then
          k = k+1
          subtable1(i) = i-k
          subtable2(k) = i-k
       endif
    enddo

    do idet = 1, nodetectors
       do i = 1,nosamples_proc
          if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
          ip = pixels(i, idet)/nosubpix_max
          pixels(i, idet) = pixels(i,idet) -subtable1(ip)*nosubpix_max
       enddo
    enddo

    dummy_pixel = (nosubmaps_tot-subtable1(nosubmaps_tot))*nosubpix_max

    if (info > 4) write(*,idf) ID,'Done'

  END SUBROUTINE reduce_pixels_a


  !---------------------------------------------------------------------------

  SUBROUTINE restore_pixels_a
    ! restore original pixel numbers
    !
    integer :: i, idet, ip

    if (info.ge.5) write(*,idf) ID,'Restore pixel numbers (a)...'

    do idet = 1, nodetectors
       do i = 1,nosamples_proc
          if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
          ip = pixels(i, idet) / nosubpix_max
          pixels(i, idet) = pixels(i,idet) + subtable2(ip)*nosubpix_max
       enddo
    enddo
       
    ksubmap = .true.
    subtable1 = 0
    subtable2 = 0
    dummy_pixel = 12*nside_max**2

    if (info > 4) write(*,idf) ID,'Done'

  END SUBROUTINE restore_pixels_a


  !--------------------------------------------------------------------------


  SUBROUTINE check_stat(allocstat)

    integer :: allocstat

    if (allocstat.ne.0) then
       write(*,*) 'ERROR: out of memory.'
       call exit_with_status(1)
    endif

  END SUBROUTINE check_stat


  !---------------------------------------------------------------------------

END MODULE pointing
