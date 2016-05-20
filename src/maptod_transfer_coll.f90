MODULE maptod_transfer
  ! E. Keihänen and M. Reinecke
  ! This modules uses the collective mpi_bcast and mpi_reduce routines,
  !  instead of process-to-process send. May be faster on some platforms.

  use commonparam
  use mpi_wrappers
  implicit none
  private

  real(dp),allocatable,public :: locmap(:,:)
  real(dp),allocatable,public :: loccc(:,:)

  integer, allocatable,public :: locmask(:)
  integer, allocatable,public :: lochits(:)

  ! Flag ksubmap_table(i,k) tells if the TOD handled by process k
  !  contributes to submap i
  logical,allocatable :: ksubmap_table(:,:)
  logical,allocatable :: emptymap(:)
  integer, save, public :: nsize_locmap = -1

  public  update_maptod_transfer,  &
       collect_map,             &
       collect_cc,              &
       collect_hits,            &
       scatter_map,             &
       scatter_mask, &
       free_locmaps ! -RK

  real(sp),save,public :: memory_locmap = 0.0

CONTAINS

  !-------------------------------------------------------------------------------


  SUBROUTINE update_maptod_transfer(ksubmap)
    !
    ! Upadte parameters nolocmaps, nolocpix, ksubmap_table
    !   which control the MPI communication.
    ! Ksubmap tells which submaps are hit by the TOD handled by this process.
    !
    logical,intent(in) :: ksubmap(0:nosubmaps_tot-1)
    logical            :: kbuffer(0:nosubmaps_tot-1)
    integer            :: id_send, i

    nolocmaps = count(ksubmap)

    nolocpix = nolocmaps*nosubpix_max

    if (.not.allocated(ksubmap_table))   &
         allocate(ksubmap_table(0:nosubmaps_tot-1,0:ntasks-1))

    if (.not.allocated(emptymap)) allocate(emptymap(0:nosubmaps_tot-1))

    do id_send = 0,ntasks-1
       if (ID==id_send) kbuffer=ksubmap

       call broadcast_mpi(kbuffer,nosubmaps_tot,id_send)

       ksubmap_table(:,id_send) = kbuffer
    enddo

    do i = 0,nosubmaps_tot-1
       emptymap(i) = (all(.not.ksubmap_table(i,:)))
    enddo

    if (nsize_locmap.lt.nolocpix) then
       if (nsize_locmap.ge.0) deallocate(locmap,loccc,locmask,lochits)

       nsize_locmap = nolocpix+100
       allocate(locmap(nostokes,0:nsize_locmap))
       allocate(loccc(ncc,0:nsize_locmap))
       allocate(locmask(0:nsize_locmap))
       allocate(lochits(0:nsize_locmap))

       memory_locmap = (nsize_locmap+1)*(nostokes*8.+ncc*8.+8)
    endif

    locmap = 0.0
    loccc  = 0.0
    locmask = 0
    lochits = 0

  END SUBROUTINE update_maptod_transfer


  !------------------------------------------------------------------------------


  SUBROUTINE free_locmaps

    if (allocated(locmap))  deallocate(locmap)
    if (allocated(loccc))   deallocate(loccc)
    if (allocated(locmask)) deallocate(locmask)
    if (allocated(lochits)) deallocate(lochits)
    if (allocated(ksubmap_table)) deallocate(ksubmap_table)
    if (allocated(emptymap)) deallocate(emptymap)
    nsize_locmap = -1

  END SUBROUTINE free_locmaps


  !------------------------------------------------------------------------------


  SUBROUTINE collect_map(map,nosubpix)
    !
    ! locmap->map
    integer, intent(in)    :: nosubpix
    real(dp),intent(inout) :: map(nostokes,nosubpix,nosubmaps)
    integer                :: i, j, k, m, ndegrade, lidx, ism, mapidx, is
    real(dp),allocatable   :: buffer(:,:,:), rbuffer(:,:,:)

    allocate(buffer(nostokes,nosubpix,nosubmaps_max))
    allocate(rbuffer(nostokes,nosubpix,nosubmaps_max))
    !Use same indexing with map to avoid confusion

    ndegrade = nosubpix_max/nosubpix

    do i = 0,ntasks-1

       ! build a local representation of process #i's "map", then reduce
       ! Do not distribute empty maps
       buffer = 0
       rbuffer = 0

       lidx = -1
       is = 0
       do mapidx = 0,nosubmaps_tot-1

          if (emptymap(mapidx)) cycle

          if (ksubmap_table(mapidx,id)) lidx=lidx+1

          if (id_submap(mapidx)==i) then
             is = is+1

             if (ksubmap_table(mapidx,id)) then
                m = lidx*nosubpix_max

                do k = 1,nosubpix
                   do j = 1,ndegrade
                      buffer(:,k,is) = buffer(:,k,is) +locmap(:,m)
                      m = m+1
                   enddo
                enddo
             endif

          endif
       enddo

       call collect_mpi_dp(buffer,rbuffer,nostokes*nosubpix*is,i)

       if (ID.ne.i) cycle

       ism = 0
       is = 0
       do mapidx = 0,nosubmaps_tot-1

          if (id_submap(mapidx)==i) then
             ism = ism+1

             if (.not.emptymap(mapidx)) then
                is = is+1
                map(:,:,ism) = map(:,:,ism) +rbuffer(:,:,is)
             endif
          endif
       enddo

    enddo

    deallocate(buffer,rbuffer)

  END SUBROUTINE collect_map


  !------------------------------------------------------------------------------


  SUBROUTINE collect_cc(cc,nosubpix)
    !
    ! locmap->map
    integer, intent(in)    :: nosubpix
    real(dp),intent(inout) :: cc(ncc,nosubpix,nosubmaps)
    integer                :: i, j, k, m, ndegrade, lidx, ism, mapidx, is
    real(dp),allocatable   :: buffer(:,:,:), rbuffer(:,:,:)

    allocate(buffer(ncc,nosubpix,nosubmaps_max))
    allocate(rbuffer(ncc,nosubpix,nosubmaps_max))
    !Use same indexing with map to avoid confusion

    ndegrade = nosubpix_max/nosubpix

    do i = 0,ntasks-1

       ! build a local representation of process #i's "map", then reduce
       ! Do not distribute empty maps
       buffer = 0
       rbuffer = 0

       lidx = -1
       is = 0
       do mapidx = 0,nosubmaps_tot-1

          if (emptymap(mapidx)) cycle

          if (ksubmap_table(mapidx,id)) lidx=lidx+1

          if (id_submap(mapidx)==i) then
             is = is+1

             if (ksubmap_table(mapidx,id)) then
                m = lidx*nosubpix_max

                do k = 1,nosubpix
                   do j = 1,ndegrade
                      buffer(:,k,is) = buffer(:,k,is) +loccc(:,m)
                      m = m+1
                   enddo
                enddo
             endif

          endif
       enddo

       call collect_mpi_dp(buffer,rbuffer,ncc*nosubpix*is,i)

       if (ID.ne.i) cycle

       ism = 0
       is = 0
       do mapidx = 0,nosubmaps_tot-1

          if (id_submap(mapidx)==i) then
             ism = ism+1

             if (.not.emptymap(mapidx)) then
                is = is+1
                cc(:,:,ism) = cc(:,:,ism) +rbuffer(:,:,is)
             endif
          endif
       enddo

    enddo

    deallocate(buffer,rbuffer)

  END SUBROUTINE collect_cc


  !------------------------------------------------------------------------------


  SUBROUTINE collect_hits(hits,nosubpix)
    !
    ! locmap->map
    integer, intent(in)    :: nosubpix
    integer, intent(inout) :: hits(nosubpix,nosubmaps)
    integer                :: i, j, k, m, ndegrade, lidx, ism, mapidx, is
    integer,allocatable    :: buffer(:,:), rbuffer(:,:)

    allocate(buffer(nosubpix,nosubmaps_max))
    allocate(rbuffer(nosubpix,nosubmaps_max))
    !Use same indexing with map to avoid confusion

    ndegrade = nosubpix_max/nosubpix

    do i = 0,ntasks-1

       ! build a local representation of process #i's "map", then reduce
       ! Do not distribute empty maps
       buffer = 0
       rbuffer = 0

       lidx = -1
       is = 0
       do mapidx = 0,nosubmaps_tot-1

          if (emptymap(mapidx)) cycle

          if (ksubmap_table(mapidx,id)) lidx=lidx+1

          if (id_submap(mapidx)==i) then
             is = is+1

             if (ksubmap_table(mapidx,id)) then
                m = lidx*nosubpix_max

                do k = 1,nosubpix
                   do j = 1,ndegrade
                      buffer(k,is) = buffer(k,is) +lochits(m)
                      m = m+1
                   enddo
                enddo
             endif

          endif
       enddo

       call collect_mpi_int(buffer,rbuffer,nosubpix*is,i)

       if (ID.ne.i) cycle

       ism = 0
       is = 0
       do mapidx = 0,nosubmaps_tot-1

          if (id_submap(mapidx)==i) then
             ism = ism+1

             if (.not.emptymap(mapidx)) then
                is = is+1
                hits(:,ism) = hits(:,ism) +rbuffer(:,is)
             endif
          endif
       enddo

    enddo

    deallocate(buffer,rbuffer)

  END SUBROUTINE collect_hits


  !------------------------------------------------------------------------------


  SUBROUTINE scatter_map(map,nosubpix)

    integer, intent(in)  :: nosubpix
    real(dp),intent(in)  :: map(nostokes,nosubpix,nosubmaps)
    integer              :: i, j, k, m, ndegrade, ism, lidx, mapidx, is
    real(dp),allocatable :: buffer(:,:,:)

    ndegrade = nosubpix_max/nosubpix
    locmap = 0.0

    allocate(buffer(nostokes,nosubpix,nosubmaps_max))
    buffer = 0.0

    do i = 0,ntasks-1
       ! broadcast all submaps of process #i
       ! skip empty maps
       ism = 0
       is = 0
       do mapidx = 0,nosubmaps_tot-1
          if (id_submap(mapidx)==i) then
             ism = ism+1

             if (.not.emptymap(mapidx)) then
                is = is+1
                if (ID==i) buffer(:,:,is)=map(:,:,ism)
             endif
          endif
       enddo

       call broadcast_mpi_vec_dp(buffer,nostokes*nosubpix*is,i)

       is = 0
       lidx = -1
       do mapidx = 0,nosubmaps_tot-1

          if (ksubmap_table(mapidx,id)) lidx=lidx+1

          if (id_submap(mapidx)==i.and..not.emptymap(mapidx)) then
             is = is+1

             if (ksubmap_table(mapidx,id)) then
                m = lidx*nosubpix_max

                do k = 1,nosubpix
                   do j = 1,ndegrade
                      locmap(:,m) = buffer(:,k,is)
                      m = m+1
                   enddo
                enddo
             endif

          endif

       end do
    end do

    deallocate(buffer)

  END SUBROUTINE scatter_map


  !------------------------------------------------------------------------------


  SUBROUTINE scatter_mask(mask,nosubpix)

    integer, intent(in)  :: nosubpix
    integer, intent(in)  :: mask(nosubpix,nosubmaps)
    integer              :: i, j, k, m, ndegrade, ism, lidx, mapidx, is
    integer, allocatable :: buffer(:,:)

    ndegrade = nosubpix_max/nosubpix
    locmask = 0.0

    allocate(buffer(nosubpix,nosubmaps_max))
    buffer = 0.0

    do i = 0,ntasks-1
       ! broadcast all submaps of process #i
       ! skip empty maps
       ism = 0
       is = 0
       do mapidx = 0,nosubmaps_tot-1
          if (id_submap(mapidx)==i) then
             ism = ism+1

             if (.not.emptymap(mapidx)) then
                is = is+1
                if (ID==i) buffer(:,is)=mask(:,ism)
             endif
          endif
       enddo

       call broadcast_mpi_vec_int(buffer,nosubpix*is,i)

       is = 0
       lidx = -1
       do mapidx = 0,nosubmaps_tot-1

          if (ksubmap_table(mapidx,id)) lidx=lidx+1

          if (id_submap(mapidx)==i.and..not.emptymap(mapidx)) then
             is = is+1

             if (ksubmap_table(mapidx,id)) then
                m = lidx*nosubpix_max

                do k = 1,nosubpix
                   do j = 1,ndegrade
                      locmask(m) = buffer(k,is)
                      m = m+1
                   enddo
                enddo
             endif

          endif

       end do
    end do

    deallocate(buffer)

  END SUBROUTINE scatter_mask


  !-------------------------------------------------------------------------------

END MODULE maptod_transfer
