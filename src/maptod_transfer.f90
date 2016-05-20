MODULE maptod_transfer

  use commonparam
  use mpi_wrappers
  implicit none
  private

  real(dp), allocatable, public :: locmap(:,:)
  real(dp), allocatable, public :: loccc(:,:,:)

  real(dp), allocatable, target, public :: submaps_send_map(:,:,:), submaps_recv_map(:,:,:)
  real(dp), pointer, public :: submaps_send_cross(:,:,:), submaps_recv_cross(:,:,:)
  integer, allocatable, target, public :: submaps_send_int_map(:,:), submaps_recv_int_map(:,:)
  integer, pointer, public :: submaps_send_int_cross(:,:), submaps_recv_int_cross(:,:)
  real(dp), allocatable, target, public :: subcc_send_map(:,:,:,:), subcc_recv_map(:,:,:,:)
  real(dp), pointer, public :: subcc_send_cross(:,:,:,:), subcc_recv_cross(:,:,:,:)
  integer, allocatable, public :: submaps_send_ind(:), submaps_recv_ind(:)
  integer, allocatable, public :: sendcounts(:), sendoffs(:), recvcounts(:), recvoffs(:)
  integer :: nsend_submap, nrecv_submap

  integer, allocatable,public :: locmask(:)
  integer, allocatable,public :: lochits(:)

  ! Flag ksubmap_table(i,k) tells if the TOD handled by process k
  !  contributes to submap i
  logical,allocatable :: ksubmap_table(:,:)
  integer,save, public :: nsize_locmap = -1

  public  update_maptod_transfer,  &
       collect_map,             &
       collect_cc,              &
       collect_hits,            &
       scatter_map,             &
       scatter_mask, &
       free_locmaps, initialize_alltoallv, assign_submaps ! -RK

  real(sp),save,public :: memory_locmap = 0.0, memory_all2all = 0.0

CONTAINS

  !---------------------------------------------------------------------------


  SUBROUTINE update_maptod_transfer(ksubmap)
    !
    ! Update parameters nolocmaps, nolocpix, ksubmap_table
    !   which control the MPI communication.
    ! Ksubmap tells which submaps are hit by the TOD handled by this process.
    !
    logical,intent(in) :: ksubmap(0:nosubmaps_tot-1)
    logical :: kbuffer(0:nosubmaps_tot-1)
    integer :: id_send
    integer :: ierr

    nolocmaps = count(ksubmap)

    nolocpix = nolocmaps*nosubpix_max

    if (.not. allocated(ksubmap_table))   &
         allocate(ksubmap_table(0:nosubmaps_tot-1,0:ntasks-1))

    do id_send = 0,ntasks-1
       if (id == id_send) kbuffer = ksubmap

       call broadcast_mpi(kbuffer, nosubmaps_tot, id_send)

       ksubmap_table(:,id_send) = kbuffer
    end do

    if (nsize_locmap < nolocpix) then
       if (nsize_locmap >= 0) deallocate(locmap, loccc, locmask, lochits)

       nsize_locmap = nolocpix !+ 100 This margin was not used anywhere. Disabled on 2014-05-09 -RK
       allocate(locmap(nmap,0:nsize_locmap), loccc(nmap,nmap,0:nsize_locmap), &
            locmask(0:nsize_locmap), lochits(0:nsize_locmap), stat=ierr)
       if (ierr /= 0) call abort_mpi('Failed to allocate locmap')

       memory_locmap = (nsize_locmap+1)*(nmap*8.+nmap**2*8.+8)
    endif

    locmap = 0.0
    loccc  = 0.0
    locmask = 0
    lochits = 0

  END SUBROUTINE update_maptod_transfer


  !---------------------------------------------------------------------------


  SUBROUTINE initialize_alltoallv()
    !
    ! set up buffers and index tables to perform map scatter/gather using alltoallv
    ! rather than point-to-point communication
    !
    ! Should be called after last call to update_maptod_transfer, before any calls
    ! to collect_map/scatter_map
    ! User may want to update id_submap table calling assign submaps before this routine.
    !
    integer :: ierr
    integer :: nsend, nrecv, itask, offset, ioffset, i

    if ( .not. concatenate_messages ) return

    ! allocate memory for collective alltoallv operations during CG iteration

    nsend_submap = count( ksubmap_table(:,id) )
    nrecv_submap = 0
    do i = 0, nosubmaps_tot-1
       if ( id_submap(i) == id ) nrecv_submap = nrecv_submap + count( ksubmap_table(i,:) )
    end do
    if ( allocated(submaps_send_map) ) then
       deallocate( submaps_send_map, submaps_recv_map )
       if ( use_inmask .or. do_hits ) then
          deallocate( submaps_send_int_map, submaps_recv_int_map )
          if ( nosubpix_cross /= nosubpix_map ) then
             deallocate( submaps_send_int_cross, submaps_recv_int_cross )
          end if
       end if
       deallocate( subcc_send_map, subcc_recv_map )
       if ( nosubpix_cross /= nosubpix_map ) then
          deallocate( submaps_send_cross, submaps_recv_cross )
          deallocate( subcc_send_cross, subcc_recv_cross )
       end if
       deallocate( submaps_send_ind, submaps_recv_ind )
       deallocate( sendcounts, sendoffs, recvcounts, recvoffs )
    end if

    allocate( submaps_send_ind(nsend_submap), submaps_recv_ind(nrecv_submap), stat=ierr )
    if ( ierr /= 0 ) stop 'No room to catenate messages'
    memory_all2all = nmap*nosubpix_cross*(nsend_submap + nrecv_submap)*4.

    allocate( submaps_send_map(nmap,nosubpix_map,nsend_submap), &
         submaps_recv_map(nmap,nosubpix_map,nrecv_submap), stat=ierr )
    if ( ierr /= 0 ) stop 'No room to catenate messages'
    memory_all2all = memory_all2all + nmap*nosubpix_map*(nsend_submap + nrecv_submap)*8.

    if ( use_inmask .or. do_hits ) then
       allocate( submaps_send_int_map(nosubpix_map,nsend_submap), &
            submaps_recv_int_map(nosubpix_map,nrecv_submap), stat=ierr )
       if ( ierr /= 0 ) stop 'No room to catenate messages'
       memory_all2all = memory_all2all + nosubpix_map*(nsend_submap + nrecv_submap)*4.
    end if

    allocate( subcc_send_map(nmap,nmap,nosubpix_map,nsend_submap), &
         subcc_recv_map(nmap,nmap,nosubpix_map,nrecv_submap), stat=ierr )
    if ( ierr /= 0 ) stop 'No room to catenate messages'
    memory_all2all = memory_all2all + nmap*nmap*nosubpix_map*(nsend_submap + nrecv_submap)*8.

    if ( nosubpix_cross /= nosubpix_map ) then
       allocate( submaps_send_cross(nmap,nosubpix_cross,nsend_submap), &
            submaps_recv_cross(nmap,nosubpix_cross,nrecv_submap), stat=ierr )
       if ( ierr /= 0 ) stop 'No room to catenate messages'
       memory_all2all = memory_all2all + nmap*nosubpix_cross*(nsend_submap + nrecv_submap)*8.
       if ( use_inmask .or. do_hits) then
          allocate( submaps_send_int_cross(nosubpix_cross,nsend_submap), &
               submaps_recv_int_cross(nosubpix_cross,nrecv_submap), stat=ierr )
          if ( ierr /= 0 ) stop 'No room to catenate messages'
          memory_all2all = memory_all2all + nosubpix_cross*(nsend_submap + nrecv_submap)*4.
       end if


       allocate( subcc_send_cross(nmap,nmap,nosubpix_cross,nsend_submap), &
            subcc_recv_cross(nmap,nmap,nosubpix_cross,nrecv_submap), stat=ierr )
       if ( ierr /= 0 ) stop 'No room to catenate messages'
       memory_all2all = memory_all2all + nmap*nmap*nosubpix_cross*(nsend_submap + nrecv_submap)*8.       
    else
       submaps_send_cross => submaps_send_map
       submaps_recv_cross => submaps_recv_map
       if ( use_inmask .or.do_hits) then
          submaps_send_int_cross => submaps_send_int_map
          submaps_recv_int_cross => submaps_recv_int_map
       end if
       subcc_send_cross => subcc_send_map
       subcc_recv_cross => subcc_recv_map
    end if
    
    allocate( sendcounts(0:ntasks-1), sendoffs(0:ntasks-1), recvcounts(0:ntasks-1), recvoffs(0:ntasks-1), stat=ierr )         
    if ( ierr /= 0 ) stop 'No room to catenate messages (2)'
    memory_all2all = memory_all2all + ntasks*4*4.0
       
    offset = 0
    ioffset = 0
    do itask = 0, ntasks-1
       nsend = count( ksubmap_table(:,id) .and. id_submap == itask ) ! * nmap * nosubpix_cross
       if ( nsend == 0 ) then
          sendcounts(itask) = 0
          sendoffs(itask) = offset
          cycle
       end if
       do i = 0,nosubmaps_tot-1
          if ( ksubmap_table(i,id) .and. id_submap(i) == itask ) then
             ioffset = ioffset + 1
             submaps_send_ind(ioffset) = count( ksubmap_table(0:i,id) ) - 1
          end if
       end do
       sendcounts(itask) = nsend
       sendoffs(itask) = offset
       offset = offset + nsend
    end do

    offset = 0
    ioffset = 0
    do itask = 0, ntasks-1
       nrecv = count( ksubmap_table(:,itask) .and. id_submap == id ) ! * nmap * nosubpix_cross
       if ( nrecv == 0 ) then
          recvcounts(itask) = 0
          recvoffs(itask) = offset             
          cycle
       end if
       do i = 0,nosubmaps_tot-1
          if ( ksubmap_table(i,itask) .and. id_submap(i) == id ) then
             ioffset = ioffset + 1
             submaps_recv_ind(ioffset) = count( id_submap(0:i) == id )
          end if
       end do
       recvcounts(itask) = nrecv
       recvoffs(itask) = offset
       offset = offset + nrecv
    end do

  END SUBROUTINE initialize_alltoallv


  !---------------------------------------------------------------------------


  SUBROUTINE assign_submaps( id_submap, nosubmaps, nopix_map, nopix_cross, nosubmaps_max )
    !
    ! Assign the submaps to minize communication
    !
    ! Updates the id_submap vector based on ksubmap_table
    ! Should be called before initialize_alltoallv
    !
    integer,intent(out) :: id_submap(0:nosubmaps_tot-1)
    integer,intent(out) :: nosubmaps, nopix_map, nopix_cross, nosubmaps_max

    logical, allocatable :: ksubmap(:,:)
    integer :: ierr, itask, ind, i, nosubmap_target, isubmap
    integer, allocatable :: nosubmaps_task(:)

    allocate( ksubmap(0:nosubmaps_tot-1,0:ntasks-1), nosubmaps_task(0:ntasks-1), stat=ierr )
    if ( ierr /= 0 ) STOP 'No room to assign submaps'
    ksubmap = ksubmap_table
    nosubmaps_task = 0

    nosubmap_target = ceiling( dble(nosubmaps_tot) / ntasks )

    id_submap = -1

    ! First assign submaps to processes with local data up to nosubmap_target submaps per process

    loop_target : do i = 1,nosubmap_target
       loop_task : do itask = 0, ntasks-1
          loop_submap : do isubmap = 0, nosubmaps_tot-1
             if ( ksubmap(isubmap,itask) ) then
                id_submap(isubmap) = itask
                ksubmap(isubmap,:) = .false.
                nosubmaps_task(itask) = nosubmaps_task(itask) + 1
                cycle loop_task
             end if
          end do loop_submap
          if ( nosubmaps_task(itask) == 0 ) then
             ! This process did not find any available submaps, pick the first available to have at least one
             loop_submap2 : do isubmap = 0, nosubmaps_tot-1
                if ( id_submap(isubmap) == -1 ) then
                   id_submap(isubmap) = itask
                   ksubmap(isubmap,:) = .false.
                   nosubmaps_task(itask) = nosubmaps_task(itask) + 1
                   cycle loop_task
                end if
             end do loop_submap2
          end if
       end do loop_task
    end do loop_target

    ! Then assign the rest of the maps. This time in a round robin fashion but 
    ! never more than nosubmap_target

    itask = 0
    do isubmap = 0, nosubmaps_tot-1
       if ( id_submap(isubmap) == -1 ) then
          ! Make sure the current task has free slots
          do while ( nosubmaps_task(itask) == nosubmap_target )
             itask = modulo( itask+1, ntasks )
          end do

          ! Assign the unassigned submap
          id_submap(isubmap) = itask
          nosubmaps_task(itask) = nosubmaps_task(itask) + 1

          ! Next submap is assigned to the next task
          itask = modulo( itask+1, ntasks )
       end if
    end do

    ! update the auxiliary information

    nosubmaps = count( id_submap == id )
    nopix_map = nosubmaps * nosubpix_map
    nopix_cross = nosubmaps * nosubpix_cross
    nosubmaps_max = nosubmaps
    call max_mpi( nosubmaps_max )

    deallocate( ksubmap )
    deallocate( nosubmaps_task )

  END SUBROUTINE assign_submaps


  !---------------------------------------------------------------------------


  SUBROUTINE free_locmaps

    if (allocated(locmap))  deallocate(locmap)
    if (allocated(loccc))   deallocate(loccc)
    if (allocated(locmask)) deallocate(locmask)
    if (allocated(lochits)) deallocate(lochits)
    if (allocated(ksubmap_table)) deallocate(ksubmap_table)

    if ( nosubpix_map /= nosubpix_cross ) then
       if (associated(submaps_send_cross)) deallocate(submaps_send_cross)
       if (associated(submaps_recv_cross)) deallocate(submaps_recv_cross)
       if ( use_inmask .or. do_hits ) then
          if (associated(submaps_send_int_cross)) deallocate(submaps_send_int_cross)
          if (associated(submaps_recv_int_cross)) deallocate(submaps_recv_int_cross)
       end if
       if (associated(subcc_send_cross)) deallocate(subcc_send_cross)
       if (associated(subcc_recv_cross)) deallocate(subcc_recv_cross)
    end if

    if (allocated(submaps_send_map)) deallocate(submaps_send_map)
    if (allocated(submaps_recv_map)) deallocate(submaps_recv_map)
    if ( use_inmask .or. do_hits ) then
       if (allocated(submaps_send_int_map)) deallocate(submaps_send_int_map)
       if (allocated(submaps_recv_int_map)) deallocate(submaps_recv_int_map)
    end if
    if (allocated(subcc_send_map)) deallocate(subcc_send_map)
    if (allocated(subcc_recv_map)) deallocate(subcc_recv_map)

    if (allocated(submaps_send_ind)) deallocate(submaps_send_ind)
    if (allocated(submaps_recv_ind)) deallocate(submaps_recv_ind)

    if (allocated(sendcounts)) deallocate(sendcounts)
    if (allocated(sendoffs)) deallocate(sendoffs)
    if (allocated(recvcounts)) deallocate(recvcounts)
    if (allocated(recvoffs)) deallocate(recvoffs)

    nsize_locmap = -1

  END SUBROUTINE free_locmaps


  !---------------------------------------------------------------------------


  SUBROUTINE collect_map(map,nosubpix)
    !
    ! locmap->map
    integer, intent(in) :: nosubpix
    real(dp),intent(inout) :: map(nmap,nosubpix,nosubmaps)
    integer :: i, j, k, m, mrecv, id_tod, id_map, ndegrade, nmap0
    real(dp) :: buffer(nmap,nosubpix)
    integer :: ierr, ind, id_thread, num_threads
    real(dp), pointer :: submaps_send(:,:,:), submaps_recv(:,:,:)

    ndegrade = nosubpix_max/nosubpix
    mrecv = 0
    m = 0

    if (concatenate_messages) then
       ! use alltoallv to reduce the local maps into global

       if ( nosubpix == nosubpix_map ) then
          submaps_send => submaps_send_map
          submaps_recv => submaps_recv_map
       else
          submaps_send => submaps_send_cross
          submaps_recv => submaps_recv_cross
       end if

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,k)
       if ( ndegrade == 1 ) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             submaps_send(:,:,i) = locmap(:,m:m+nosubpix-1)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix * ndegrade
             do k = 1,nosubpix
                submaps_send(:,k,i) = sum( locmap(:,m:m+ndegrade-1), 2 )
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL

       nmap0 = size(submaps_send,1) ! Workaround for unpolarized subsets

       call mpi_alltoallv( submaps_send, sendcounts*nmap0*nosubpix, sendoffs*nmap0*nosubpix, MPI_DOUBLE_PRECISION, &
            submaps_recv, recvcounts*nmap0*nosubpix, recvoffs*nmap0*nosubpix, MPI_DOUBLE_PRECISION, comm, ierr )

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(id_thread,num_threads,i,m)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       do i = 1, nrecv_submap
          m = submaps_recv_ind(i)
          if ( num_threads > 1 ) then ! don't do the costly modulo with one thread
             if ( modulo(m,num_threads) /= id_thread ) cycle
          end if
          map(1:nmap,:,m) = map(1:nmap,:,m) + submaps_recv(1:nmap,:,i)
       end do
       !$OMP END PARALLEL
       
    else
       do i = 0,nosubmaps_tot-1
          
          id_map = id_submap(i)
          if (ID==id_map) mrecv=mrecv+1

          do id_tod = 0,ntasks-1

             if (.not.ksubmap_table(i,id_tod)) cycle

             if (ID==id_tod) then  ! prepare send buffer

                buffer = 0.0
                do k = 1,nosubpix
                   do j = 1,ndegrade
                      buffer(:,k) = buffer(:,k)+locmap(:,m)
                      m = m+1
                   enddo
                enddo
             endif

             call send_mpi_vec_dp(buffer,nmap*nosubpix,id_tod,id_map)

             if (ID==id_map) map(:,:,mrecv)=map(:,:,mrecv)+buffer

          enddo
       enddo
    end if

  END SUBROUTINE collect_map


  !---------------------------------------------------------------------------


  SUBROUTINE collect_cc(cc, nosubpix)

    integer, intent(in)    :: nosubpix
    real(dp),intent(inout) :: cc(nmap,nmap,nosubpix,nosubmaps)
    integer :: i, j, k, m, mrecv, id_tod, id_map, ndegrade, nmap0
    real(dp) :: buffer(nmap,nmap,nosubpix)
    integer :: ierr, ind, id_thread, num_threads
    real(dp), pointer :: subcc_send(:,:,:,:), subcc_recv(:,:,:,:)

    ndegrade = nosubpix_max/nosubpix
    mrecv = 0
    m = 0

    if (concatenate_messages) then
       ! use alltoallv to reduce the local maps into global

       if ( nosubpix == nosubpix_map ) then
          subcc_send => subcc_send_map
          subcc_recv => subcc_recv_map
       else
          subcc_send => subcc_send_cross
          subcc_recv => subcc_recv_cross
       end if

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,k)
       if ( ndegrade == 1 ) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             subcc_send(:,:,:,i) = loccc(:,:,m:m+nosubpix-1)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix * ndegrade
             do k = 1,nosubpix
                subcc_send(:,:,k,i) = sum( loccc(:,:,m:m+ndegrade-1), 3 )
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL

       nmap0 = size(subcc_send,1) ! Workaround for unpolarized subsets

       call mpi_alltoallv( subcc_send, sendcounts*nmap0*nmap0*nosubpix, sendoffs*nmap0*nmap0*nosubpix, MPI_DOUBLE_PRECISION, &
            subcc_recv, recvcounts*nmap0*nmap0*nosubpix, recvoffs*nmap0*nmap0*nosubpix, MPI_DOUBLE_PRECISION, comm, ierr )

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(id_thread,num_threads,i,m)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       do i = 1, nrecv_submap
          m = submaps_recv_ind(i)
          if ( num_threads > 1 ) then ! don't do the costly modulo with one thread
             if ( modulo(m,num_threads) /= id_thread ) cycle
          end if
          cc(1:nmap,1:nmap,:,m) = cc(1:nmap,1:nmap,:,m) + subcc_recv(1:nmap,1:nmap,:,i)
       end do
       !$OMP END PARALLEL
       
    else
       do i = 0,nosubmaps_tot-1
          
          id_map = id_submap(i)
          if (ID==id_map) mrecv = mrecv + 1
          
          do id_tod = 0,ntasks-1
             
             if (.not.ksubmap_table(i,id_tod)) cycle
             
             if (ID==id_tod) then  ! prepare send buffer
                
                buffer = 0.0
                do k = 1,nosubpix
                   do j = 1,ndegrade
                      buffer(:,:,k) = buffer(:,:,k) + loccc(:,:,m)
                      m = m + 1
                   end do
                end do
             end if
             
             call send_mpi_vec_dp(buffer, nmap**2*nosubpix, id_tod, id_map)
             
             if (ID==id_map) cc(:,:,:,mrecv) = cc(:,:,:,mrecv)+buffer
             
          end do
       end do
    end if

  END SUBROUTINE collect_cc

  !---------------------------------------------------------------------------


  SUBROUTINE collect_hits(hits,nosubpix)

    integer, intent(in)    :: nosubpix
    integer, intent(inout) :: hits(nosubpix,nosubmaps)
    integer                :: i, j, k, m, mrecv, id_tod, id_map, ndegrade
    integer                :: buffer(nosubpix)
    integer, pointer :: submaps_send(:,:), submaps_recv(:,:)
    integer :: ierr, ind, id_thread, num_threads

    ndegrade = nosubpix_max/nosubpix
    mrecv = 0
    m = 0

    if (concatenate_messages) then
       ! use alltoallv to reduce the local maps into global

       if ( nosubpix == nosubpix_map ) then
          submaps_send => submaps_send_int_map
          submaps_recv => submaps_recv_int_map
       else
          submaps_send => submaps_send_int_cross
          submaps_recv => submaps_recv_int_cross
       end if

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,k)
       if ( ndegrade == 1 ) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             submaps_send(:,i) = lochits(m:m+nosubpix-1)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix * ndegrade
             do k = 1,nosubpix
                submaps_send(k,i) = sum( lochits(m:m+ndegrade-1) )
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL

       call mpi_alltoallv( submaps_send, sendcounts*nosubpix, sendoffs*nosubpix, MPI_INTEGER, &
            submaps_recv, recvcounts*nosubpix, recvoffs*nosubpix, MPI_INTEGER, comm, ierr )

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(id_thread,num_threads,i,m)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       do i = 1, nrecv_submap
          m = submaps_recv_ind(i)
          if ( num_threads > 1 ) then ! don't do the costly modulo with one thread
             if ( modulo(m,num_threads) /= id_thread ) cycle
          end if
          hits(:,m) = hits(:,m) + submaps_recv(:,i)
       end do
       !$OMP END PARALLEL
       
    else
       do i = 0,nosubmaps_tot-1

          id_map = id_submap(i)
          if (ID==id_map) mrecv=mrecv+1
          
          do id_tod = 0,ntasks-1
             
             if (.not.ksubmap_table(i,id_tod)) cycle
             
             if (ID==id_tod) then  ! prepare send buffer
                
                buffer = 0
                do k = 1,nosubpix
                   do j = 1,ndegrade
                      buffer(k) = buffer(k)+lochits(m)
                      m = m+1
                   enddo
                enddo
             endif
             
             call send_mpi(buffer,nosubpix,id_tod,id_map)

             if (ID==id_map) hits(:,mrecv)=hits(:,mrecv)+buffer
             
          enddo
       enddo
    end if

  END SUBROUTINE collect_hits


  !---------------------------------------------------------------------------


  SUBROUTINE scatter_map(map,nosubpix)

    integer, intent(in)  :: nosubpix
    real(dp),intent(in)  :: map(nmap,nosubpix,nosubmaps)
    integer              :: i, j, k, m, msend, id_tod, id_map, ndegrade
    real(dp) :: buffer(nmap,nosubpix)
    integer :: ierr, ind, id_thread, num_threads
    real(dp), pointer :: submaps_send(:,:,:), submaps_recv(:,:,:)

    ndegrade = nosubpix_max/nosubpix
    msend = 0
    m     = 0

    locmap = 0.0

    if (concatenate_messages .and. nosubpix==nosubpix_cross) then
       ! use alltoallv to scatter the global map into locals       
       ! This is the inverse operation of the alltoallv in collect_map so we swap the
       ! send and recv arrays

       if ( nosubpix == nosubpix_map ) then
          submaps_send => submaps_send_map
          submaps_recv => submaps_recv_map
       else
          submaps_send => submaps_send_cross
          submaps_recv => submaps_recv_cross
       end if

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(id_thread,num_threads,i,m,k,msend,id_tod)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       m = 0       
       do id_tod = 0, ntasks-1
          msend = 0
          do i = 0, nosubmaps_tot-1
             if ( id_submap(i) /= id ) cycle
             msend = msend + 1
             if ( .not. ksubmap_table(i,id_tod) ) cycle
             m = m + 1
             if ( num_threads > 1 ) then ! don't do the costly modulo with one thread
                if ( modulo(m,num_threads) /= id_thread ) cycle
             end if
             submaps_recv(:,:,m) = map(:,:,msend)
          end do
       end do
       !$OMP END PARALLEL

       call mpi_alltoallv( submaps_recv, recvcounts*nmap*nosubpix, recvoffs*nmap*nosubpix, MPI_DOUBLE_PRECISION, &
            submaps_send, sendcounts*nmap*nosubpix, sendoffs*nmap*nosubpix, MPI_DOUBLE_PRECISION, comm, ierr )

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,k)
       if ( ndegrade == 1 ) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             locmap(:,m:m+nosubpix-1) = submaps_send(:,:,i)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             m = m * ndegrade
             do k = 1, nosubpix
                locmap(:,m:m+ndegrade-1) = spread( submaps_send(:,k,i), 2, ndegrade )
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL
    else

       do i = 0,nosubmaps_tot-1
          id_map = id_submap(i)
          
          if (ID==id_map) then
             msend = msend+1
             buffer = map(:,:,msend)
          endif
          
          do id_tod = 0,ntasks-1
             
             if (.not.ksubmap_table(i,id_tod)) cycle

             call send_mpi_vec_dp(buffer,nmap*nosubpix,id_map,id_tod)
             
             if (ID==id_tod) then
                do k = 1,nosubpix
                   do j = 1,ndegrade
                      locmap(:,m) = buffer(:,k)
                      m = m+1
                   enddo
                enddo
             endif
             
          enddo
       enddo
       
    end if

  END SUBROUTINE scatter_map


  !---------------------------------------------------------------------------


  SUBROUTINE scatter_mask(mask,nosubpix)

    integer, intent(in)  :: nosubpix
    integer, intent(in)  :: mask(nosubpix,nosubmaps)
    integer              :: i, j, k, m, msend, id_tod, id_map, ndegrade
    integer              :: buffer(nosubpix)
    integer :: ierr, ind, id_thread, num_threads
    integer, pointer :: submaps_send(:,:), submaps_recv(:,:)


    ndegrade = nosubpix_max/nosubpix
    msend = 0
    m     = 0

    locmap = 0.0

    if (concatenate_messages .and. nosubpix==nosubpix_cross) then
       ! use alltoallv to scatter the global map into locals       
       ! This is the inverse operation of the alltoallv in collect_map so we swap the
       ! send and recv arrays

       if ( nosubpix == nosubpix_map ) then
          submaps_send => submaps_send_int_map
          submaps_recv => submaps_recv_int_map
       else
          submaps_send => submaps_send_int_cross
          submaps_recv => submaps_recv_int_cross
       end if

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(id_thread,num_threads,i,m,k,msend,id_tod)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       m = 0       
       do id_tod = 0, ntasks-1
          msend = 0
          do i = 0, nosubmaps_tot-1
             if ( id_submap(i) /= id ) cycle
             msend = msend + 1
             if ( .not. ksubmap_table(i,id_tod) ) cycle
             m = m + 1
             if ( num_threads > 1 ) then ! don't do the costly modulo with one thread
                if ( modulo(m,num_threads) /= id_thread ) cycle
             end if
             submaps_recv(:,m) = mask(:,msend)
          end do
       end do
       !$OMP END PARALLEL

       call mpi_alltoallv( submaps_recv, recvcounts*nosubpix, recvoffs*nosubpix, MPI_INTEGER, &
            submaps_send, sendcounts*nosubpix, sendoffs*nosubpix, MPI_INTEGER, comm, ierr )

       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,k)
       if ( ndegrade == 1 ) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             locmask(m:m+nosubpix-1) = submaps_send(:,i)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             m = m * ndegrade
             do k = 1, nosubpix
                locmask(m:m+ndegrade-1) = spread( submaps_send(k,i), 1, ndegrade )
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL
    else
       do i = 0,nosubmaps_tot-1
          id_map = id_submap(i)
          
          if (ID==id_map) then
             msend = msend+1
             buffer = mask(:,msend)
          endif
          
          do id_tod = 0,ntasks-1
             
             if (.not.ksubmap_table(i,id_tod)) cycle
             
             call send_mpi(buffer,nosubpix,id_map,id_tod)
             
             if (ID==id_tod) then
                do k = 1,nosubpix
                   do j = 1,ndegrade
                      locmask(m) = buffer(k)
                      m = m+1
                   enddo
                enddo
             endif
             
          enddo
       enddo
    end if

  END SUBROUTINE scatter_mask


  !---------------------------------------------------------------------------


END MODULE maptod_transfer
