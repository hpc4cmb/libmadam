MODULE maptod_transfer

  use commonparam
  use mpi_wrappers
  use memory_and_time, only : write_memory
  implicit none
  private

  real(dp), allocatable, public :: locmap(:, :)
  real(dp), allocatable, public :: loccc(:, :, :)

  real(dp), allocatable, target :: submaps_send_map(:, :, :)
  real(dp), allocatable, target :: submaps_recv_map(:, :, :)
  real(dp), pointer :: submaps_send_cross(:, :, :)
  real(dp), pointer :: submaps_recv_cross(:, :, :)
  integer, allocatable, target :: submaps_send_int_map(:, :)
  integer, allocatable, target :: submaps_recv_int_map(:, :)
  integer, pointer :: submaps_send_int_cross(:, :)
  integer, pointer :: submaps_recv_int_cross(:, :)
  integer, allocatable :: submaps_send_ind(:), submaps_recv_ind(:)
  integer, allocatable :: sendcounts(:), sendoffs(:)
  integer, allocatable :: recvcounts(:), recvoffs(:)
  integer, allocatable :: recvcounts_gather(:), displs_gather(:)
  integer :: nsend_submap, nrecv_submap, nsend_gather

  integer, allocatable, public :: locmask(:)
  integer, allocatable, public :: lochits(:)

  ! Flag ksubmap_table(i,k) tells if the TOD handled by process k
  !  contributes to submap i
  logical, allocatable, public :: ksubmap_table(:, :)
  integer, save, public :: nsize_locmap = -1

  character(len=40), parameter :: mstr='(x,a,t32,f9.1," MB")'
  character(len=40), parameter :: mstr3='(x,a,t32,3(f9.1," MB"))'

  public update_maptod_transfer, collect_map, collect_cc, collect_hits, &
       scatter_map, scatter_mask, free_locmaps, initialize_alltoallv, &
       assign_submaps

  real(dp), save, public :: memory_locmap = 0, memory_all2all = 0, &
       memory_ksubmap = 0

CONTAINS

  !---------------------------------------------------------------------------


  SUBROUTINE update_maptod_transfer(ksubmap)
    !
    ! Update parameters nolocmaps, nolocpix, ksubmap_table
    !   which control the MPI communication.
    ! Ksubmap tells which submaps are hit by the TOD handled by this process.
    !
    logical, allocatable, intent(in) :: ksubmap(:)
    integer(i4b) :: ierr
    real(sp) :: memsum, mem_min, mem_max
    integer :: nolocmaps_min, nolocmaps_max
    real(dp) :: nolocmaps_sum

    nolocmaps = count(ksubmap(0:nosubmaps_tot-1))
    nolocpix = nolocmaps * nosubpix_max

    nolocmaps_min = nolocmaps
    nolocmaps_max = nolocmaps
    nolocmaps_sum = nolocmaps
    call min_mpi(nolocmaps_min)
    call max_mpi(nolocmaps_max)
    call sum_mpi(nolocmaps_sum)
    if (ID == 0 .and. info > 0) then
       write(*, '(x,"Total submaps = ",i0,"  submap size = ",i0)') &
            nosubmaps_tot, nosubpix_max
       write(*, '(x,"Local submaps:  min = ",i0,"  max = ",i0,"  mean = ",f8.2)') &
            nolocmaps_min, nolocmaps_max, nolocmaps_sum / ntasks
    end if

    if (.not. allocated(ksubmap_table)) then
       ! Default Fortran logical is 32 bits (!)
       if (allreduce) then
          allocate(ksubmap_table(0:nosubmaps_tot - 1, 0:0), stat=ierr)
          memory_ksubmap = nosubmaps_tot*4
       else
          allocate(ksubmap_table(0:nosubmaps_tot - 1, 0:ntasks-1), stat=ierr)
          memory_ksubmap = nosubmaps_tot * ntasks * 4
       end if
       if (ierr /= 0) call abort_mpi('No room for ksubmap_table')
       ksubmap_table = .false.
       call write_memory('Submap table memory', memory_ksubmap)
    end if

    if (allreduce) then
       ! All processes share ksubmap
       ksubmap_table(:, 0) = ksubmap(0:nosubmaps_tot - 1)
    else
       call mpi_allgather(ksubmap, nosubmaps_tot, MPI_LOGICAL, &
            ksubmap_table, nosubmaps_tot, MPI_LOGICAL, comm, ierr)
       if (ierr /= 0) call abort_mpi('Gathering ksubmaps failed.')
    end if

    if (nsize_locmap < nolocpix) then
       if (allocated(locmap)) deallocate(locmap)
       if (allocated(loccc)) deallocate(loccc)
       if (allocated(locmask)) deallocate(locmask)
       if (allocated(lochits)) deallocate(lochits)

       nsize_locmap = nolocpix
       allocate(locmap(nmap, 0:nsize_locmap), &
            loccc(nmap, nmap, 0:nsize_locmap), locmask(0:nsize_locmap), &
            lochits(0:nsize_locmap), stat=ierr)

       if (ierr /= 0) call abort_mpi('Failed to allocate locmap')

       memory_locmap = (nsize_locmap + 1.) * (nmap * 8. + nmap ** 2 * 8. + 4)

    end if

    locmap = 0
    loccc = 0
    locmask = 0
    lochits = 0

    call write_memory('local maps memory', memory_locmap)

  END SUBROUTINE update_maptod_transfer


  !---------------------------------------------------------------------------


  SUBROUTINE initialize_alltoallv()
    !
    ! set up buffers and index tables to perform map scatter/gather using
    ! alltoallv rather than point-to-point communication
    !
    ! Should be called after last call to update_maptod_transfer, before any
    ! calls to collect_map/scatter_map
    ! User may want to update id_submap table calling assign submaps
    ! before this routine.
    !
    integer :: ierr
    integer :: nsend, nrecv, itask, offset, ioffset, i
    real(sp) :: memsum, mem_min, mem_max

    if (allreduce) then
       nsend_gather = 0
       do i = 0, nosubmaps_tot-1
          if ((id_submap(i) == id) .and. ksubmap_table(i, 0)) then
             nsend_gather = nsend_gather + 1
          end if
       end do

       if (allocated(recvcounts_gather)) deallocate(recvcounts_gather)
       if (allocated(displs_gather)) deallocate(displs_gather)
       allocate(recvcounts_gather(ntasks), displs_gather(ntasks), stat=ierr)
       if (ierr /= 0) stop 'No room for allgatherv counts'
       memory_all2all = ntasks*2*4.

       call mpi_allgather(nsend_gather, 1, MPI_INTEGER, recvcounts_gather, &
            1, MPI_INTEGER, comm, ierr)

       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to gather counts with allgather')

       displs_gather(1) = 0
       do itask = 2, ntasks
          displs_gather(itask) = displs_gather(itask-1) &
               + recvcounts_gather(itask-1)
       end do

       call write_memory("Allreduce memory", memory_all2all)
    end if

    ! allreduce implies concatenate_messages = .false.

    if (.not. concatenate_messages) return

    ! allocate memory for collective alltoallv operations during CG iteration

    nsend_submap = count(ksubmap_table(:, id))
    nrecv_submap = 0
    do i = 0, nosubmaps_tot-1
       if (id_submap(i) == id) then
          nrecv_submap = nrecv_submap + count(ksubmap_table(i, :))
       end if
    end do
    if (allocated(submaps_send_map)) then
       deallocate(submaps_send_map, submaps_recv_map)
       if (use_inmask .or. do_hits) then
          deallocate(submaps_send_int_map, submaps_recv_int_map)
          if (nosubpix_cross /= nosubpix_map) then
             deallocate(submaps_send_int_cross, submaps_recv_int_cross)
          end if
       end if
       if (nosubpix_cross /= nosubpix_map) then
          deallocate(submaps_send_cross, submaps_recv_cross)
       end if
       deallocate(submaps_send_ind, submaps_recv_ind)
       deallocate(sendcounts, sendoffs, recvcounts, recvoffs)
    end if

    allocate(submaps_send_ind(nsend_submap), submaps_recv_ind(nrecv_submap), &
         stat=ierr)
    if (ierr /= 0) stop 'No room to concatenate messages'
    memory_all2all = nmap*nosubpix_cross*(nsend_submap + nrecv_submap)*4.

    allocate(submaps_send_map(nmap, nosubpix_map, nsend_submap), &
         submaps_recv_map(nmap, nosubpix_map, nrecv_submap), stat=ierr)
    if (ierr /= 0) stop 'No room to concatenate messages'
    memory_all2all = memory_all2all &
         + nmap*nosubpix_map*(nsend_submap + nrecv_submap)*8.

    if (use_inmask .or. do_hits) then
       allocate(submaps_send_int_map(nosubpix_map, nsend_submap), &
            submaps_recv_int_map(nosubpix_map, nrecv_submap), stat=ierr)
       if (ierr /= 0) stop 'No room to concatenate messages'
       memory_all2all = memory_all2all &
            + nosubpix_map*(nsend_submap + nrecv_submap)*4.
    end if

    if (nosubpix_cross /= nosubpix_map) then
       allocate(submaps_send_cross(nmap, nosubpix_cross, nsend_submap), &
            submaps_recv_cross(nmap, nosubpix_cross, nrecv_submap), stat=ierr)
       if (ierr /= 0) stop 'No room to concatenate messages'
       memory_all2all = memory_all2all &
            + nmap*nosubpix_cross*(nsend_submap + nrecv_submap)*8.
       if (use_inmask .or. do_hits) then
          allocate(submaps_send_int_cross(nosubpix_cross, nsend_submap), &
               submaps_recv_int_cross(nosubpix_cross, nrecv_submap), stat=ierr)
          if (ierr /= 0) stop 'No room to concatenate messages'
          memory_all2all = memory_all2all &
               + nosubpix_cross*(nsend_submap + nrecv_submap)*4.
       end if

       if (ierr /= 0) stop 'No room to concatenate messages'
       memory_all2all = memory_all2all &
            + nmap*nmap*nosubpix_cross*(nsend_submap + nrecv_submap)*8.
    else
       submaps_send_cross => submaps_send_map
       submaps_recv_cross => submaps_recv_map
       if (use_inmask .or.do_hits) then
          submaps_send_int_cross => submaps_send_int_map
          submaps_recv_int_cross => submaps_recv_int_map
       end if
    end if

    allocate(sendcounts(0:ntasks-1), sendoffs(0:ntasks-1), &
         recvcounts(0:ntasks-1), recvoffs(0:ntasks-1), stat=ierr)
    if (ierr /= 0) stop 'No room to concatenate messages (2)'
    memory_all2all = memory_all2all + ntasks*4*4.

    offset = 0
    ioffset = 0
    do itask = 0, ntasks-1
       nsend = count(ksubmap_table(:, id) .and. id_submap == itask)
       if (nsend == 0) then
          sendcounts(itask) = 0
          sendoffs(itask) = offset
          cycle
       end if
       do i = 0, nosubmaps_tot-1
          if (ksubmap_table(i, id) .and. id_submap(i) == itask) then
             ioffset = ioffset + 1
             submaps_send_ind(ioffset) = count(ksubmap_table(0:i, id)) - 1
          end if
       end do
       sendcounts(itask) = nsend
       sendoffs(itask) = offset
       offset = offset + nsend
    end do

    offset = 0
    ioffset = 0
    do itask = 0, ntasks-1
       nrecv = count(ksubmap_table(:, itask) .and. id_submap == id)
       if (nrecv == 0) then
          recvcounts(itask) = 0
          recvoffs(itask) = offset
          cycle
       end if
       do i = 0, nosubmaps_tot-1
          if (ksubmap_table(i, itask) .and. id_submap(i) == id) then
             ioffset = ioffset + 1
             submaps_recv_ind(ioffset) = count(id_submap(0:i) == id)
          end if
       end do
       recvcounts(itask) = nrecv
       recvoffs(itask) = offset
       offset = offset + nrecv
    end do

    call write_memory("All2allv memory", memory_all2all)

  END SUBROUTINE initialize_alltoallv


  !---------------------------------------------------------------------------


  SUBROUTINE assign_submaps(id_submap, nosubmaps, nopix_map, nopix_cross, &
       nosubmaps_max)
    !
    ! Assign the submaps to minimize communication
    !
    ! Updates the id_submap vector based on ksubmap_table
    ! Should be called before initialize_alltoallv
    !
    integer, intent(out) :: id_submap(0:nosubmaps_tot-1)
    integer, intent(out) :: nosubmaps, nopix_map, nopix_cross, nosubmaps_max

    logical, allocatable :: ksubmap(:, :)
    integer :: ierr, itask, i, nosubmap_target, isubmap
    integer, allocatable :: nosubmaps_task(:)
    integer :: isubmap_start, isubmap_stop

    if (id == 0 .and. nosubmaps_tot > 100000) then
       write(*,'(a,i0,a,i0,a)') 'WARNING: You have a LOT of submaps (', &
            nosubmaps_tot, &
            '). Reassigning submaps will take time. Reduce nside_submap (', &
               nside_submap, ') to divide the map in larger chunks'
    end if

    nosubmap_target = ceiling(dble(nosubmaps_tot) / ntasks)
    id_submap = -1

    if (allreduce) then
       ! Assign the submaps in contiguous blocs.  This will allow
       ! fast MPI_allgather operations.
       isubmap_start = 0
       itask = 0
       do while (isubmap_start < nosubmaps_tot)
          isubmap_stop = min(isubmap_start+nosubmap_target-1, nosubmaps_tot-1)
          id_submap(isubmap_start:isubmap_stop) = itask
          isubmap_start = isubmap_stop + 1
          itask = itask + 1
       end do
    else
       ! Every process owns every submap in allreduce mode
       if (.not. allocated(ksubmap_table)) &
            call abort_mpi('assign_submaps: ksubmap_table not allocated')

       if (ID == 0) then
          allocate(ksubmap(0:nosubmaps_tot-1, 0:ntasks-1), &
               nosubmaps_task(0:ntasks-1), stat=ierr)
          if (ierr /= 0) call abort_mpi('No room to assign submaps')
          ksubmap = ksubmap_table
          nosubmaps_task = 0

          ! First assign submaps to processes with local data up to
          ! nosubmap_target submaps per process

          loop_target : do i = 1, nosubmap_target
             loop_task : do itask = 0, ntasks-1
                loop_submap : do isubmap = 0, nosubmaps_tot-1
                   if (ksubmap(isubmap, itask)) then
                      id_submap(isubmap) = itask
                      ksubmap(isubmap, :) = .false.
                      nosubmaps_task(itask) = nosubmaps_task(itask) + 1
                      cycle loop_task
                   end if
                end do loop_submap
                if (nosubmaps_task(itask) == 0) then
                   ! This process did not find any available submaps,
                   ! pick the first available to have at least one
                   loop_submap2 : do isubmap = 0, nosubmaps_tot-1
                      if (id_submap(isubmap) == -1) then
                         id_submap(isubmap) = itask
                         ksubmap(isubmap, :) = .false.
                         nosubmaps_task(itask) = nosubmaps_task(itask) + 1
                         cycle loop_task
                      end if
                   end do loop_submap2
                end if
             end do loop_task
          end do loop_target

          ! Then assign the rest of the maps. This time in a round robin fashion
          ! but never more than nosubmap_target

          itask = 0
          do isubmap = 0, nosubmaps_tot-1
             if (id_submap(isubmap) == -1) then
                ! Make sure the current task has free slots
                do while (nosubmaps_task(itask) == nosubmap_target)
                   itask = modulo(itask+1, ntasks)
                end do

                ! Assign the unassigned submap
                id_submap(isubmap) = itask
                nosubmaps_task(itask) = nosubmaps_task(itask) + 1

                ! Next submap is assigned to the next task
                itask = modulo(itask+1, ntasks)
             end if
          end do

          deallocate(ksubmap)
          deallocate(nosubmaps_task)
       end if
       call broadcast_mpi(id_submap, nosubmaps_tot, 0)
    end if

    ! update the auxiliary information

    nosubmaps = count(id_submap == id)
    nopix_map = nosubmaps * nosubpix_map
    nopix_cross = nosubmaps * nosubpix_cross
    nosubmaps_max = nosubmaps
    call max_mpi(nosubmaps_max)

  END SUBROUTINE assign_submaps


  !---------------------------------------------------------------------------


  SUBROUTINE free_locmaps

    if (allocated(locmap))  deallocate(locmap)
    if (allocated(loccc))   deallocate(loccc)
    if (allocated(locmask)) deallocate(locmask)
    if (allocated(lochits)) deallocate(lochits)
    if (allocated(ksubmap_table)) deallocate(ksubmap_table)

    if (nosubpix_map /= nosubpix_cross) then
       if (associated(submaps_send_cross)) deallocate(submaps_send_cross)
       if (associated(submaps_recv_cross)) deallocate(submaps_recv_cross)
       if (use_inmask .or. do_hits) then
          if (associated(submaps_send_int_cross)) &
               deallocate(submaps_send_int_cross)
          if (associated(submaps_recv_int_cross)) &
               deallocate(submaps_recv_int_cross)
       end if
    end if

    if (allocated(submaps_send_map)) deallocate(submaps_send_map)
    if (allocated(submaps_recv_map)) deallocate(submaps_recv_map)
    if (use_inmask .or. do_hits) then
       if (allocated(submaps_send_int_map)) deallocate(submaps_send_int_map)
       if (allocated(submaps_recv_int_map)) deallocate(submaps_recv_int_map)
    end if

    if (allocated(submaps_send_ind)) deallocate(submaps_send_ind)
    if (allocated(submaps_recv_ind)) deallocate(submaps_recv_ind)

    if (allocated(sendcounts)) deallocate(sendcounts)
    if (allocated(sendoffs)) deallocate(sendoffs)
    if (allocated(recvcounts)) deallocate(recvcounts)
    if (allocated(recvoffs)) deallocate(recvoffs)

    if (allocated(recvcounts_gather)) deallocate(recvcounts_gather)
    if (allocated(displs_gather)) deallocate(displs_gather)

    nsize_locmap = -1

  END SUBROUTINE free_locmaps


  !---------------------------------------------------------------------------


  SUBROUTINE collect_map(map, nosubpix, inplace)
    !
    ! locmap->map
    integer, intent(in) :: nosubpix
    real(dp), intent(inout) :: map(nmap, nosubpix, nosubmaps)
    logical, optional :: inplace
    integer :: i, j, k, m, n, mrecv, id_tod, id_map, ndegrade, nmap0
    real(dp) :: buffer(nmap, nosubpix)
    integer :: ierr, id_thread, num_threads, submap_type
    real(dp), pointer :: submaps_send(:, :, :), submaps_recv(:, :, :)

    ndegrade = nosubpix_max / nosubpix

    if (allreduce) then
       call sum_mpi(locmap)
       if (ndegrade == 1) then
          !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
          !$OMP     PRIVATE(i, j, k, m, id_thread) &
          !$OMP     SHARED(ksubmap_table, id_submap, id, nosubmaps_tot, &
          !$OMP         nosubpix, locmap, map, nthreads)
          id_thread = omp_get_thread_num()
          m = 0
          k = 0
          do i = 0, nosubmaps_tot-1
             if (ksubmap_table(i, 0)) k = k + 1
             if (id == id_submap(i)) m = m + 1
             if (nthreads > 1) then
                if (modulo(m, nthreads) /= id_thread) cycle
             end if
             if (ksubmap_table(i, 0) .and. id == id_submap(i)) then
                j = (k-1) * nosubpix
                map(:, :, m) = map(:, :, m) + locmap(:, j:j+nosubpix-1)
             end if
          end do
          !$OMP END PARALLEL
       else
          !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
          !$OMP     PRIVATE(i, j, k, m, n, id_thread) &
          !$OMP     SHARED(ksubmap_table, id_submap, id, nosubmaps_tot, &
          !$OMP         nosubpix, locmap, map, nthreads, ndegrade)
          id_thread = omp_get_thread_num()
          m = 0
          k = 0
          do i = 0, nosubmaps_tot-1
             if (ksubmap_table(i, 0)) k = k + 1
             if (id == id_submap(i)) m = m + 1
             if (nthreads > 1) then
                if (modulo(m, nthreads) /= id_thread) cycle
             end if
             if (ksubmap_table(i, 0) .and. id == id_submap(i)) then
                n = 0
                do j = (k-1)*nosubpix*ndegrade, k*nosubpix*ndegrade-1, ndegrade
                   n = n + 1
                   map(:, n, m) = map(:, n, m) + sum(locmap(:, j:j+ndegrade-1), 2)
                end do
             end if
          end do
          !$OMP END PARALLEL
       end if
       ! Replace locmap values with an average so subsequent calls won't fail
       if (.not. present(inplace)) then
          !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) &
          !$OMP     SHARED(locmap, ntasks, nsize_locmap)
          do i = 0, nsize_locmap
             locmap(:, i) = locmap(:, i) / ntasks
          end do
          !$OMP END PARALLEL DO
       end if
    else if (concatenate_messages) then
       ! use alltoallv to reduce the local maps into global

       if (nosubpix == nosubpix_map) then
          submaps_send => submaps_send_map
          submaps_recv => submaps_recv_map
       else
          submaps_send => submaps_send_cross
          submaps_recv => submaps_recv_cross
       end if

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, m, k) &
       !$OMP     SHARED(ndegrade, nsend_submap, submaps_send_ind, nosubpix, &
       !$OMP         submaps_send, locmap)
       if (ndegrade == 1) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             submaps_send(:, :, i) = locmap(:, m:m+nosubpix-1)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix * ndegrade
             do k = 1, nosubpix
                submaps_send(:, k, i) = sum(locmap(:, m:m+ndegrade-1), 2)
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL

       nmap0 = size(submaps_send, 1) ! Workaround for unpolarized subsets

       ! Manipulate the auxiliary arrays in place rather than have
       ! temporary arrays created

       call mpi_type_contiguous( &
            nmap0 * nosubpix, MPI_DOUBLE_PRECISION, submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to create mpi_type')

       call mpi_type_commit(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to commit mpi_type')

       call mpi_alltoallv( &
            submaps_send, sendcounts, sendoffs, submap_type, &
            submaps_recv, recvcounts, recvoffs, submap_type, &
            comm, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to collect map with alltoallv')

       call mpi_type_free(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to free mpi_type')

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(id_thread, num_threads, i, m) &
       !$OMP     SHARED(nrecv_submap, submaps_recv_ind, map, submaps_recv, nmap)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       do i = 1, nrecv_submap
          m = submaps_recv_ind(i)
          if (num_threads > 1) then
             ! don't do the costly modulo with one thread
             if (modulo(m, num_threads) /= id_thread) cycle
          end if
          map(1:nmap, :, m) = map(1:nmap, :, m) + submaps_recv(1:nmap, :, i)
       end do
       !$OMP END PARALLEL
    else
       mrecv = 0
       m = 0
       do i = 0, nosubmaps_tot-1
          id_map = id_submap(i)
          if (ID == id_map) mrecv = mrecv + 1
          do id_tod = 0, ntasks-1
             if (.not. ksubmap_table(i, id_tod)) cycle
             if (ID == id_tod) then ! prepare send buffer
                buffer = 0.0
                do k = 1, nosubpix
                   do j = 1, ndegrade
                      buffer(:, k) = buffer(:, k) + locmap(:, m)
                      m = m + 1
                   end do
                end do
             end if
             call send_mpi_vec_dp(buffer, nmap*nosubpix, id_tod, id_map)
             if (ID == id_map) map(:, :, mrecv) = map(:, :, mrecv) + buffer
          end do
       end do
    end if

  END SUBROUTINE collect_map


  !---------------------------------------------------------------------------


  SUBROUTINE collect_cc(cc, nosubpix)

    integer, intent(in) :: nosubpix
    real(dp), intent(inout) :: cc(nmap, nmap, nosubpix, nosubmaps)
    integer :: i, j, k, m, n, mrecv, id_tod, id_map, ndegrade, nmap0, col
    real(dp) :: buffer(nmap, nmap, nosubpix)
    integer :: ierr, id_thread, num_threads, submap_type
    real(dp), pointer :: submaps_send(:, :, :), submaps_recv(:, :, :)

    ndegrade = nosubpix_max / nosubpix

    if (allreduce) then
       call sum_mpi(loccc)
       if (ndegrade == 1) then
          !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
          !$OMP     PRIVATE(i, j, k, m, id_thread) &
          !$OMP     SHARED(ksubmap_table, id_submap, id, nosubmaps_tot, &
          !$OMP         nosubpix, loccc, cc, nthreads)
          id_thread = omp_get_thread_num()
          m = 0
          k = 0
          do i = 0, nosubmaps_tot-1
             if (ksubmap_table(i, 0)) k = k + 1
             if (id == id_submap(i)) m = m + 1
             if (nthreads > 1) then
                if (modulo(m, nthreads) /= id_thread) cycle
             end if
             if (ksubmap_table(i, 0) .and. id == id_submap(i)) then
                j = (k-1) * nosubpix
                cc(:, :, :, m) = cc(:, :, :, m) + loccc(:, :, j:j+nosubpix-1)
             end if
          end do
          !$OMP END PARALLEL
       else
          !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
          !$OMP     PRIVATE(i, j, k, m, n, col, id_thread) &
          !$OMP     SHARED(ksubmap_table, id_submap, id, nosubmaps_tot, &
          !$OMP         nosubpix, loccc, cc, nthreads, ndegrade, nmap)
          id_thread = omp_get_thread_num()
          m = 0
          k = 0
          do i = 0, nosubmaps_tot-1
             if (ksubmap_table(i, 0)) k = k + 1
             if (id == id_submap(i)) m = m + 1
             if (nthreads > 1) then
                if (modulo(m, nthreads) /= id_thread) cycle
             end if
             if (ksubmap_table(i, 0) .and. id == id_submap(i)) then
                do col = 1, nmap
                   n = 0
                   do j = (k-1)*nosubpix*ndegrade, k*nosubpix*ndegrade-1, &
                        ndegrade
                      n = n + 1
                      cc(:, col, n, m) = cc(:, col, n, m) &
                           + sum(loccc(:, col, j:j+ndegrade-1), 2)
                   end do
                end do
             end if
          end do
          !$OMP END PARALLEL
       end if
       ! scale back the value of loccc after allreduce
       !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i) &
       !$OMP     SHARED(loccc, ntasks, nsize_locmap)
       do i = 0, nsize_locmap
          loccc(:, :, i) = loccc(:, :, i) / ntasks
       end do
       !$OMP END PARALLEL DO
    else if (concatenate_messages) then
       ! use alltoallv to reduce the local maps into global

       if (nosubpix == nosubpix_map) then
          submaps_send => submaps_send_map
          submaps_recv => submaps_recv_map
       else
          submaps_send => submaps_send_cross
          submaps_recv => submaps_recv_cross
       end if

       nmap0 = size(submaps_send, 1) ! Workaround for unpolarized subsets

       do col = 1, nmap0

          !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, m, k) &
          !$OMP     SHARED(ndegrade, nsend_submap, submaps_send_ind, nosubpix, &
          !$OMP         submaps_send, loccc, col)
          if (ndegrade == 1) then
             !$OMP DO
             do i = 1, nsend_submap
                m = submaps_send_ind(i) * nosubpix
                submaps_send(:, :, i) = loccc(:, col, m:m+nosubpix-1)
             end do
             !$OMP END DO
          else
             !$OMP DO
             do i = 1, nsend_submap
                m = submaps_send_ind(i) * nosubpix * ndegrade
                do k = 1, nosubpix
                   submaps_send(:, k, i) = sum(loccc(:, col, m:m+ndegrade-1), 2)
                   m = m + ndegrade
                end do
             end do
             !$OMP END DO
          end if
          !$OMP END PARALLEL

          ! Manipulate the auxiliary arrays in place rather than have
          ! temporary arrays created

          call mpi_type_contiguous( &
               nmap0 * nosubpix, MPI_DOUBLE_PRECISION, submap_type, ierr)
          if (ierr /= MPI_SUCCESS) &
               call abort_mpi('Failed to create mpi_type')

          call mpi_type_commit(submap_type, ierr)
          if (ierr /= MPI_SUCCESS) &
               call abort_mpi('Failed to commit mpi_type')

          call mpi_alltoallv( &
               submaps_send, sendcounts, sendoffs, submap_type, &
               submaps_recv, recvcounts, recvoffs, submap_type, &
               comm, ierr)
          if (ierr /= MPI_SUCCESS) call abort_mpi('Failed to collect cc')

          call mpi_type_free(submap_type, ierr)
          if (ierr /= MPI_SUCCESS) &
               call abort_mpi('Failed to free mpi_type')

          !$OMP PARALLEL DEFAULT(NONE) PRIVATE(id_thread, num_threads, i, m) &
          !$OMP     SHARED(nrecv_submap, submaps_recv_ind, cc, submaps_recv, &
          !$OMP         nmap, col)
          id_thread = omp_get_thread_num()
          num_threads = omp_get_num_threads()

          do i = 1, nrecv_submap
             m = submaps_recv_ind(i)
             if (num_threads > 1) then
                ! don't do the costly modulo with one thread
                if (modulo(m, num_threads) /= id_thread) cycle
             end if
             cc(1:nmap, col, :, m) = cc(1:nmap, col, :, m) &
                  + submaps_recv(1:nmap, :, i)
          end do
          !$OMP END PARALLEL
       end do
    else
       mrecv = 0
       m = 0
       do i = 0, nosubmaps_tot-1
          id_map = id_submap(i)
          if (ID == id_map) mrecv = mrecv + 1
          do id_tod = 0, ntasks-1
             if (.not. ksubmap_table(i, id_tod)) cycle
             if (ID == id_tod) then  ! prepare send buffer
                buffer = 0.0
                do k = 1, nosubpix
                   do j = 1, ndegrade
                      buffer(:, :, k) = buffer(:, :, k) + loccc(:, :, m)
                      m = m + 1
                   end do
                end do
             end if
             call send_mpi_vec_dp(buffer, nmap**2*nosubpix, id_tod, id_map)
             if (ID == id_map) cc(:, :, :, mrecv) = cc(:, :, :, mrecv) + buffer
          end do
       end do
    end if

  END SUBROUTINE collect_cc

  !---------------------------------------------------------------------------


  SUBROUTINE collect_hits(hits, nosubpix)

    integer, intent(in) :: nosubpix
    integer, intent(inout) :: hits(nosubpix, nosubmaps)
    integer :: i, j, k, m, n, mrecv, id_tod, id_map, ndegrade
    integer :: buffer(nosubpix)
    integer, pointer :: submaps_send(:, :), submaps_recv(:, :)
    integer :: ierr, id_thread, num_threads, submap_type

    ndegrade = nosubpix_max / nosubpix

    if (allreduce) then
       call sum_mpi(lochits)
       m = 0
       k = 0
       if (ndegrade == 1) then
          do i = 0, nosubmaps_tot-1
             if (ksubmap_table(i, 0)) k = k + 1
             if (id == id_submap(i)) m = m + 1
             if (ksubmap_table(i, 0) .and. id == id_submap(i)) then
                j = (k-1) * nosubpix
                hits(:, m) = hits(:, m) + lochits(j:j+nosubpix-1)
             end if
          end do
       else
          do i = 0, nosubmaps_tot-1
             if (ksubmap_table(i, 0)) k = k + 1
             if (id == id_submap(i)) m = m + 1
             if (ksubmap_table(i, 0) .and. id == id_submap(i)) then
                n = 0
                do j = (k-1)*nosubpix*ndegrade, k*nosubpix*ndegrade-1, ndegrade
                   n = n + 1
                   hits(n, m) = hits(n, m) + sum(lochits(j:j+ndegrade-1))
                end do
             end if
          end do
       end if
       ! Collect hits will not be called again. No need to average
    else if (concatenate_messages) then
       ! use alltoallv to reduce the local maps into global

       if (nosubpix == nosubpix_map) then
          submaps_send => submaps_send_int_map
          submaps_recv => submaps_recv_int_map
       else
          submaps_send => submaps_send_int_cross
          submaps_recv => submaps_recv_int_cross
       end if

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, m, k) &
       !$OMP     SHARED(ndegrade, nsend_submap, submaps_send_ind, nosubpix, &
       !$OMP         submaps_send, lochits)
       if (ndegrade == 1) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             submaps_send(:, i) = lochits(m:m+nosubpix-1)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix * ndegrade
             do k = 1, nosubpix
                submaps_send(k, i) = sum(lochits(m:m+ndegrade-1))
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL

       ! Manipulate the auxiliary arrays in place rather than have
       ! temporary arrays created

       call mpi_type_contiguous( &
            nosubpix, MPI_INTEGER, submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to create mpi_type')

       call mpi_type_commit(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to commit mpi_type')

       call mpi_alltoallv( &
            submaps_send, sendcounts, sendoffs, submap_type, &
            submaps_recv, recvcounts, recvoffs, submap_type, comm, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to collect hits with alltoallv')

       call mpi_type_free(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to free mpi_type')

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(id_thread, num_threads, i, m) &
       !$OMP     SHARED(nrecv_submap, submaps_recv_ind, hits, submaps_recv)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       do i = 1, nrecv_submap
          m = submaps_recv_ind(i)
          if (num_threads > 1) then
             ! don't do the costly modulo with one thread
             if (modulo(m, num_threads) /= id_thread) cycle
          end if
          hits(:, m) = hits(:, m) + submaps_recv(:, i)
       end do
       !$OMP END PARALLEL

    else
       mrecv = 0
       m = 0
       do i = 0, nosubmaps_tot-1
          id_map = id_submap(i)
          if (ID == id_map) mrecv = mrecv + 1
          do id_tod = 0, ntasks-1
             if (.not. ksubmap_table(i, id_tod)) cycle
             if (ID == id_tod) then ! prepare send buffer
                buffer = 0
                do k = 1, nosubpix
                   do j = 1, ndegrade
                      buffer(k) = buffer(k) + lochits(m)
                      m = m + 1
                   end do
                end do
             end if
             call send_mpi(buffer, nosubpix, id_tod, id_map)
             if (ID == id_map) hits(:, mrecv) = hits(:, mrecv) + buffer
          end do
       end do
    end if

  END SUBROUTINE collect_hits


  !---------------------------------------------------------------------------


  SUBROUTINE scatter_map(map, nosubpix)

    integer, intent(in) :: nosubpix
    real(dp), intent(in) :: map(nmap, nosubpix, nosubmaps)
    integer :: i, j, k, m, msend, id_tod, id_map, ndegrade
    real(dp) :: buffer(nmap, nosubpix)
    integer :: ierr, id_thread, num_threads, submap_type
    real(dp), pointer :: submaps_send(:, :, :), submaps_recv(:, :, :)
    real(dp), allocatable :: recvbuf(:, :, :), sendbuf(:, :, :)

    ndegrade = nosubpix_max / nosubpix
    locmap = 0

    if (allreduce) then
       allocate(sendbuf(nmap, nosubpix, nsend_gather), &
            recvbuf(nmap, nosubpix, nolocmaps), stat=ierr)
       if (ierr /= 0) stop 'No room for allgatherv buffers'

       m = 0
       k = 0
       sendbuf = 0
       do i = 0, nosubmaps_tot-1
          if (id == id_submap(i)) then
             m = m + 1
             if (ksubmap_table(i, 0)) then
                k = k + 1
                sendbuf(:, :, k) = map(:, :, m)
             end if
          end if
       end do

       call mpi_allgatherv(sendbuf, nsend_gather*nmap*nosubpix, &
            MPI_DOUBLE_PRECISION, recvbuf, recvcounts_gather*nmap*nosubpix, &
            displs_gather*nmap*nosubpix, MPI_DOUBLE_PRECISION, comm, ierr)

       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to gather map with allgatherv')

       deallocate(sendbuf)

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, k, m) &
       !$OMP     SHARED(ndegrade, nolocmaps, nosubpix, recvbuf, locmap)
       if (ndegrade == 1) then
          !$OMP DO
          do i = 1, nolocmaps
             m = (i-1) * nosubpix
             locmap(:, m:m+nosubpix-1) = recvbuf(:, :, i)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nolocmaps
             m = (i-1) * nosubpix * ndegrade
             do k = 1, nosubpix
                locmap(:, m:m+ndegrade-1) = &
                     spread(recvbuf(:, k, i), 2, ndegrade)
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL

       deallocate(recvbuf)
    else if (concatenate_messages .and. nosubpix == nosubpix_cross) then
       ! use alltoallv to scatter the global map into locals
       ! This is the inverse operation of the alltoallv in collect_map so we
       ! swap the send and recv arrays

       if (nosubpix == nosubpix_map) then
          submaps_send => submaps_send_map
          submaps_recv => submaps_recv_map
       else
          submaps_send => submaps_send_cross
          submaps_recv => submaps_recv_cross
       end if

       !$OMP PARALLEL DEFAULT(NONE) &
       !$OMP    PRIVATE(id_thread, num_threads, i, m, k, msend, id_tod) &
       !$OMP    SHARED(ntasks, nosubmaps_tot, id_submap, ksubmap_table, &
       !$OMP        submaps_recv, map, id)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       m = 0
       do id_tod = 0, ntasks-1
          msend = 0
          do i = 0, nosubmaps_tot-1
             if (id_submap(i) /= id) cycle
             msend = msend + 1
             if (.not. ksubmap_table(i, id_tod)) cycle
             m = m + 1
             if (num_threads > 1) then
                ! don't do the costly modulo with one thread
                if (modulo(m, num_threads) /= id_thread) cycle
             end if
             submaps_recv(:, :, m) = map(:, :, msend)
          end do
       end do
       !$OMP END PARALLEL

       ! Manipulate the auxiliary arrays in place rather than have
       ! temporary arrays created

       call mpi_type_contiguous( &
            nmap * nosubpix, MPI_DOUBLE_PRECISION, submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to create mpi_type')

       call mpi_type_commit(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to commit mpi_type')

       call mpi_alltoallv( &
            submaps_recv, recvcounts, recvoffs, submap_type, &
            submaps_send, sendcounts, sendoffs, submap_type, comm, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to scatter map with alltoallv')

       call mpi_type_free(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to free mpi_type')

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, m, k) &
       !$OMP     SHARED(ndegrade, nsend_submap, submaps_send_ind, nosubpix, &
       !$OMP         locmap, submaps_send)
       if (ndegrade == 1) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             locmap(:, m:m+nosubpix-1) = submaps_send(:, :, i)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             m = m * ndegrade
             do k = 1, nosubpix
                locmap(:, m:m+ndegrade-1) = &
                     spread(submaps_send(:, k, i), 2, ndegrade)
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL
    else
       msend = 0
       m = 0
       do i = 0, nosubmaps_tot-1
          id_map = id_submap(i)
          if (ID == id_map) then
             msend = msend + 1
             buffer = map(:, :, msend)
          end if
          do id_tod = 0, ntasks-1
             if (.not. ksubmap_table(i, id_tod)) cycle
             call send_mpi_vec_dp(buffer, nmap*nosubpix, id_map, id_tod)
             if (ID == id_tod) then
                do k = 1, nosubpix
                   do j = 1, ndegrade
                      locmap(:, m) = buffer(:, k)
                      m = m + 1
                   end do
                end do
             end if
          end do
       end do
    end if

  END SUBROUTINE scatter_map


  !---------------------------------------------------------------------------


  SUBROUTINE scatter_mask(mask, nosubpix)

    integer, intent(in) :: nosubpix
    integer, intent(in) :: mask(nosubpix, nosubmaps)
    integer :: i, j, k, m, msend, id_tod, id_map, ndegrade
    integer :: buffer(nosubpix)
    integer :: ierr, id_thread, num_threads, submap_type
    integer, pointer :: submaps_send(:, :), submaps_recv(:, :)
    integer, allocatable :: recvbuf(:, :), sendbuf(:, :)

    ndegrade = nosubpix_max / nosubpix
    locmask = 0

    if (allreduce .and. nosubpix == nosubpix_cross) then
       allocate(sendbuf(nosubpix, nsend_gather), &
            recvbuf(nosubpix, nolocmaps), stat=ierr)
       if (ierr /= 0) stop 'No room for allgatherv buffers'

       m = 0
       k = 0
       sendbuf = 0
       do i = 0, nosubmaps_tot-1
          if (id == id_submap(i)) then
             m = m + 1
             if (ksubmap_table(i, 0)) then
                k = k + 1
                sendbuf(:, k) = mask(:, m)
             end if
          end if
       end do

       call mpi_allgatherv(sendbuf, nsend_gather*nosubpix, MPI_INTEGER, &
            recvbuf, recvcounts_gather*nosubpix, displs_gather*nosubpix, &
            MPI_INTEGER, comm, ierr)

       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to gather mask with allgatherv')

       deallocate(sendbuf)

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, k, m) &
       !$OMP     SHARED(ndegrade, nolocmaps, nosubpix, recvbuf, locmask)
       if (ndegrade == 1) then
          !$OMP DO
          do i = 1, nolocmaps
             m = (i-1) * nosubpix
             locmask(m:m+nosubpix-1) = recvbuf(:, i)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nolocmaps
             m = (i-1) * nosubpix * ndegrade
             do k = 1, nosubpix
                locmask(m:m+ndegrade-1) = spread(recvbuf(k, i), 1, ndegrade)
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL

       deallocate(recvbuf)
    else if (concatenate_messages .and. nosubpix == nosubpix_cross) then
       ! use alltoallv to scatter the global map into locals
       ! This is the inverse operation of the alltoallv in collect_map so we
       ! swap the send and recv arrays

       if (nosubpix == nosubpix_map) then
          submaps_send => submaps_send_int_map
          submaps_recv => submaps_recv_int_map
       else
          submaps_send => submaps_send_int_cross
          submaps_recv => submaps_recv_int_cross
       end if

       !$OMP PARALLEL DEFAULT(NONE) &
       !$OMP     PRIVATE(id_thread, num_threads, i, m, k, msend, id_tod) &
       !$OMP     SHARED(ntasks, nosubmaps_tot, id_submap, id, ksubmap_table, &
       !$OMP         submaps_recv, mask)
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       m = 0
       do id_tod = 0, ntasks-1
          msend = 0
          do i = 0, nosubmaps_tot-1
             if (id_submap(i) /= id) cycle
             msend = msend + 1
             if (.not. ksubmap_table(i, id_tod)) cycle
             m = m + 1
             if (num_threads > 1) then
                ! don't do the costly modulo with one thread
                if (modulo(m, num_threads) /= id_thread) cycle
             end if
             submaps_recv(:, m) = mask(:, msend)
          end do
       end do
       !$OMP END PARALLEL

       ! Manipulate the auxiliary arrays in place rather than have
       ! temporary arrays created

       call mpi_type_contiguous( &
            nosubpix, MPI_INTEGER, submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to create mpi_type')

       call mpi_type_commit(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to commit mpi_type')

       call mpi_alltoallv( &
            submaps_recv, recvcounts, recvoffs, submap_type, &
            submaps_send, sendcounts, sendoffs, submap_type, comm, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to scatter mask with alltoallv')

       call mpi_type_free(submap_type, ierr)
       if (ierr /= MPI_SUCCESS) &
            call abort_mpi('Failed to free mpi_type')

       !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, m, k) &
       !$OMP     SHARED(ndegrade, nsend_submap, submaps_send_ind, nosubpix, &
       !$OMP         locmask, submaps_send)
       if (ndegrade == 1) then
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             locmask(m:m+nosubpix-1) = submaps_send(:, i)
          end do
          !$OMP END DO
       else
          !$OMP DO
          do i = 1, nsend_submap
             m = submaps_send_ind(i) * nosubpix
             m = m * ndegrade
             do k = 1, nosubpix
                locmask(m:m+ndegrade-1) = &
                     spread(submaps_send(k, i), 1, ndegrade)
                m = m + ndegrade
             end do
          end do
          !$OMP END DO
       end if
       !$OMP END PARALLEL
    else
       msend = 0
       m = 0
       do i = 0, nosubmaps_tot-1
          id_map = id_submap(i)
          if (ID == id_map) then
             msend = msend + 1
             buffer = mask(:, msend)
          end if
          do id_tod = 0, ntasks-1
             if (.not. ksubmap_table(i, id_tod)) cycle
             call send_mpi(buffer, nosubpix, id_map, id_tod)
             if (ID == id_tod) then
                do k = 1, nosubpix
                   do j = 1, ndegrade
                      locmask(m) = buffer(k)
                      m = m + 1
                   end do
                end do
             end if
          end do
       end do
    end if

  END SUBROUTINE scatter_mask


  !---------------------------------------------------------------------------


END MODULE maptod_transfer
