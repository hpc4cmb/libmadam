MODULE output

  use commonparam
  use submap_transfer
  use mpi_wrappers
  use timing
  use maps_and_baselines
  use tod_storage
  use fitsmod2

  implicit none
  private

  real(sp), parameter :: Healpix_undef = -1.6375e30
  real(sp), save :: cputime_output = 0
  integer, save :: map_repcount = 1

  ! for binary output

  integer :: comm_bin, ntasks_bin, id_bin, writegroup, groupsize

  public init_output, write_matrix, &
       write_map, write_binmap, write_mask, write_leakmatrix, &
       write_hits, close_output, addi, addi2

CONTAINS

  !----------------------------------------------------------------------


  SUBROUTINE init_output


    integer :: ierr

    if (write_cut) then
       map_repcount = 1
    else
       map_repcount = nosubpix_map
    end if

    call addpath(path_output, file_map)
    call addpath(path_output, file_binmap)
    call addpath(path_output, file_hit)
    call addpath(path_output, file_base)
    call addpath(path_output, file_mask)
    call addpath(path_output, file_matrix)
    call addpath(path_output, file_leakmatrix)
    call addpath(path_output, file_wcov)
    call addpath(path_output, file_covmat)

    ! output flags
    do_binmap  = (len_trim(file_binmap) > 0)
    do_hits    = (len_trim(file_hit) > 0)
    do_mask    = (len_trim(file_mask) > 0)

    if (binary_output) then

       if (nwrite_binary > 1) then

          groupsize = ntasks / nwrite_binary
          writegroup = id / groupsize

          call mpi_comm_split(MPI_COMM_WORLD, writegroup, id, comm_bin, ierr)
          if (ierr /= 0) then
             call abort_mpi('Failed to split the communicator for output')
          end if

          call mpi_comm_size(comm_bin, ntasks_bin, ierr)
          call mpi_comm_rank(comm_bin , id_bin, ierr)

       else

          comm_bin = MPI_COMM_WORLD
          ntasks_bin = ntasks
          id_bin = id

       end if

    end if

  CONTAINS

    !--------------------------------------------------------------------

    SUBROUTINE addpath(path, filename)

      character(len=SLEN) :: path, filename
      integer :: n

      n = len_trim(path)
      if (n==0) return

      if (path(n:n).ne.'/') path=path(1:n)//'/'

      n = len_trim(filename)

      if (n > 0) then
         if (filename(1:1)=='!') then
            filename = '!' //trim(path) // filename(2:n)
         else
            filename = trim(path) // trim(filename)
         endif
      endif

    END SUBROUTINE addpath

  END SUBROUTINE init_output


  !----------------------------------------------------------------------


  SUBROUTINE write_map(map, mask)
    !
    ! Write the destriped output map into file.
    !
    real(dp), intent(in) :: map(nmap, nosubpix_map, nosubmaps)
    integer, intent(in) :: mask(nosubpix_map, nosubmaps)
    integer :: idwr, imap, isubmap, skycover, coloffset
    integer(idp) :: offset, nhit, writeoffset, i
    real(sp), allocatable :: sbuffer(:)
    integer, allocatable :: ibuffer(:)
    integer(i8b), allocatable :: i8buffer(:)
    type(fitshandle) :: out
    type(fitscolumn), pointer :: columns(:)

    real(sp), allocatable :: mapbuffer(:,:)
    character(len=SLEN) :: file_map_bin, file_map_txt
    INTEGER(i4b) :: status(MPI_STATUS_SIZE)
    INTEGER(i4b) :: fh, filemode, fileinfo, ierr, pix
    INTEGER(MPI_OFFSET_KIND) :: fileoffset

    INTEGER :: nsend, nrecv, mysubmap, rec_len, rec
    logical :: there
    INTEGER, allocatable :: recvcounts(:), displs(:)
    REAL(sp), allocatable :: sendbuffer(:), recvbuffer(:)

    if (len_trim(file_map).le.0) return

    call reset_time(10)

    skycover = count(mask > 0)
    call sum_mpi(skycover)

    idwr = 1
    if (ntasks==1) idwr=0

    if (ID == idwr) then
       if (info > 1) write(*,*) 'Writing destriped map...'

       inquire(file=trim(file_map), exist=there)

       if (there .and. binary_output .and. concatenate_binary) then
          print *,'Skipping existing header file and concatenating binary maps'
       else

          call check_outfile(file_map)

          call fits_create(out,file_map)

          call write_header_code(out)
          call write_header_simulation(out)
          call write_header_destriping(out)


          if (write_cut) then
             allocate(columns(nmap+1))
             columns(1)%repcount = map_repcount
             columns(1)%name = 'Pixel index'
             columns(1)%type = fits_int8
             coloffset = 1
          else
             allocate(columns(nmap))
             coloffset = 0
          end if

          do imap = 1,nmap
             columns(imap+coloffset)%repcount = map_repcount
             columns(imap+coloffset)%unit = trim(unit_tod)
             columns(imap+coloffset)%type = fits_real4
             if (nmap == 1 .or. nmap == 3) then
                ! healpix naming
                if (imap==1) columns(imap+coloffset)%name = 'TEMPERATURE'
                if (imap==2) columns(imap+coloffset)%name = 'Q-POLARISATION'
                if (imap==3) columns(imap+coloffset)%name = 'U-POLARISATION'
             else if (nmap == 2) then
                if (imap==1) columns(imap+coloffset)%name = 'Q-POLARISATION'
                if (imap==2) columns(imap+coloffset)%name = 'U-POLARISATION'
             else
                columns(imap+coloffset)%name = addi('COMP', imap)
             end if
          enddo

          call fits_insert_bintab(out, columns)
          deallocate(columns)

          call fits_add_key(out, 'objtype', 'madam.map', 'Object type')
          call write_header_healpix(out, nside_map)
          call write_header_sky(out, skycover)
          call fits_add_key(out, 'MONOPOLE', rm_monopole, &
               'Monopole removed (T map)')
       endif

    end if

    if (binary_output) then

       file_map_bin = trim(adjustl(file_map)) // '.bin'
       file_map_txt = trim(adjustl(file_map)) // '.txt'

       if (nwrite_binary > 1) then
          file_map_bin = trim(file_map_bin) // '_' // addi('', writegroup)
       end if

       !call mpi_info_create(fileinfo, ierr)
       !filemode = ior(MPI_MODE_WRONLY, MPI_MODE_CREATE)
       !call mpi_file_open(MPI_COMM_WORLD, file_map_bin, filemode, fileinfo, &
       !    fh, ierr)
       !if (ierr /= 0) call abort_mpi('Failed to open binary map file')

       !fileoffset = 0
       !call mpi_file_set_view(fh, fileoffset, MPI_REAL4, MPI_REAL4, 'native', &
       !    fileinfo, ierr)
       !if (ierr /= 0) call abort_mpi('Failed to set view to binary map file')

       ! Collect the maps to root process for writing

       allocate(recvcounts(ntasks_bin), displs(ntasks_bin), stat=ierr)
       if (ierr /= 0) call abort_mpi('No room for recvcounts')

       nsend = 0
       do isubmap = 0, nosubmaps_tot-1
          if (id_submap(isubmap) == id) nsend = nsend + nmap * nosubpix_map
       end do

       recvcounts = 0
       call mpi_gather(nsend, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, 0, &
            comm_bin, ierr)

       allocate(sendbuffer(nsend), stat=ierr)
       if (ierr /= 0) call abort_mpi('No room for sendbuffer')

       mysubmap = 0
       i = 0
       do isubmap = 0, nosubmaps_tot-1
          if (id_submap(isubmap) == id) then
             mysubmap = mysubmap + 1
             do pix = 1, nosubpix_map
                if (mask(pix,mysubmap) <= 0) then
                   do imap = 1, nmap
                      i = i + 1
                      sendbuffer(i) = Healpix_undef
                   end do
                else
                   do imap = 1, nmap
                      i = i + 1
                      sendbuffer(i) = map(imap,pix,mysubmap)
                   end do
                end if
             end do
          end if
       end do

       if (i /= nsend) call abort_mpi('Did not fill sendbuffer') ! redundant check

       if (id_bin == 0) then
          nrecv = sum(recvcounts)

          allocate(recvbuffer(nrecv), stat=ierr)
          if (ierr /= 0) call abort_mpi('No room for recvbuffer')

          displs = 0
          do i = 2, ntasks_bin
             displs(i) = recvcounts(i-1) + displs(i-1)
          end do
       end if

       call mpi_gatherv(sendbuffer, nsend, MPI_REAL4, recvbuffer, recvcounts, &
            displs, MPI_REAL4, 0, comm_bin, ierr)

       if (id_bin == 0) then

          !fileoffset = count(id_submap < id) * nosubpix_map * nmap
          !call mpi_file_write_at(fh, fileoffset, recvbuffer, nrecv, &
          !    MPI_REAL4, status, ierr)

          inquire(iolength=rec_len) recvbuffer
          open(unit=55, file=trim(file_map_bin), action='readwrite', &
               form='unformatted', access='direct', recl=rec_len)
          !inquire(unit=55, nextrec=rec)
          !print *,id,' : nextrec = ',rec ! DEBUG
          write(55, rec=record_number) recvbuffer
          close(unit=55)

       end if

       if (id == 0) then

          ! Finally, store the submap distribution so that the map
          ! may be assembled later

          inquire(file=trim(file_map_txt), exist=there)

          if (.not. there) then

             open(unit=55, file=trim(file_map_txt), status='replace', &
                  form='formatted')
             write(55, *) nosubmaps_tot, ntasks
             write(55, *) nosubpix_map, nmap
             do isubmap = 0, nosubmaps_tot-1
                write(55, *) isubmap, id_submap(isubmap)
             end do
             close(unit=55)

          end if

       end if

       deallocate(recvcounts, displs, sendbuffer)
       if (id_bin == 0) deallocate(recvbuffer)

       !call mpi_file_close(fh, ierr)
       !if (ierr /= 0) call abort_mpi('Failed to close binary file')

    else

       allocate(sbuffer(nosubpix_map))
       allocate(ibuffer(nosubpix_map))
       if (write_cut) allocate(i8buffer(nosubpix_map))

       offset = 0
       writeoffset = 0
       do isubmap = 0, nosubmaps_tot-1

          if (sync_output) call wait_mpi

          call get_submap(ibuffer, nosubpix_map, mask, isubmap, idwr)

          if (id == idwr .and. write_cut) then
             nhit = 0
             do i = 1, nosubpix_map
                if (ibuffer(i) > 0) then
                   nhit = nhit + 1
                   i8buffer(nhit) = offset + i - 1
                end if
             end do
             call fits_write_column(out, 1, i8buffer(:nhit), writeoffset)
          else
             nhit = nosubpix_map
          end if

          do imap = 1,nmap

             call get_submap(sbuffer, nosubpix_map, map, nmap, imap, &
                  isubmap, idwr)

             if (ID==idwr) then
                if (write_cut) then
                   ! Pack the observed pixels to the beginning
                   nhit = 0
                   do i = 1, nosubpix_map
                      if (ibuffer(i) > 0) then
                         nhit = nhit + 1
                         sbuffer(nhit) = sbuffer(i)
                      end if
                   end do
                else
                   ! Set unobserved pixels to undefined value
                   where (ibuffer.le.0) sbuffer=Healpix_undef
                end if

                if (nhit > 0) then
                   call fits_write_column(out, imap+coloffset, sbuffer(:nhit), &
                        writeoffset)
                end if
             end if
          end do

          offset = offset + nosubpix_map
          writeoffset = writeoffset + nhit
       end do

       deallocate(sbuffer, ibuffer)
       if (write_cut) deallocate(i8buffer)

    end if

    if (ID==idwr) then
       call fits_close(out)
       write(*,*) 'Map written in '//trim(file_map)
    endif

    cputime_output = cputime_output + get_time(10)

  END SUBROUTINE write_map


  !----------------------------------------------------------------------


  SUBROUTINE write_summap(map, mask)
    !
    ! Write the sum map into file.
    !
    real(dp),intent(in) :: map(nmap, nosubpix_map, nosubmaps)
    integer, intent(in) :: mask(nosubpix_map, nosubmaps)
    integer :: idwr, imap, isubmap, skycover, coloffset
    integer(idp) :: offset, nhit, writeoffset, i
    real(sp), allocatable :: sbuffer(:)
    integer, allocatable :: ibuffer(:)
    integer(i8b), allocatable :: i8buffer(:)
    type(fitshandle) :: out
    type(fitscolumn), pointer :: columns(:)

    if (len_trim(file_map).le.0) return

    call reset_time(10)

    skycover = count(mask > 0)
    call sum_mpi(skycover)

    idwr = 1
    if (ntasks == 1) idwr = 0

    if (ID == idwr) then
       if (info > 1) write(*,*) 'Writing destriped map...'

       call check_outfile(file_map)

       call fits_create(out,file_map)

       call write_header_code(out)
       call write_header_simulation(out)
       call write_header_destriping(out)


       if (write_cut) then
          allocate(columns(nmap+1))
          columns(1)%repcount = map_repcount
          columns(1)%name = 'Pixel index'
          columns(1)%type = fits_int8
          coloffset = 1
       else
          allocate(columns(nmap))
          coloffset = 0
       end if

       do imap = 1,nmap
          columns(imap+coloffset)%repcount = map_repcount
          columns(imap+coloffset)%unit = trim(unit_tod)
          columns(imap+coloffset)%type = fits_real4
          if (nmap == 1 .or. nmap == 3) then
             ! healpix naming
             if (imap==1) columns(imap+coloffset)%name = 'TEMPERATURE'
             if (imap==2) columns(imap+coloffset)%name = 'Q-POLARISATION'
             if (imap==3) columns(imap+coloffset)%name = 'U-POLARISATION'
          else if (nmap == 2) then
             if (imap==1) columns(imap+coloffset)%name = 'Q-POLARISATION'
             if (imap==2) columns(imap+coloffset)%name = 'U-POLARISATION'
          else
             columns(imap+coloffset)%name = addi('STOKES', imap)
          end if
       enddo

       call fits_insert_bintab(out, columns)
       deallocate(columns)

       call fits_add_key(out, 'objtype', 'madam.map', 'Object type')
       call write_header_healpix(out, nside_map)
       call write_header_sky(out, skycover)
       call fits_add_key(out, 'MONOPOLE', rm_monopole, &
            'Monopole removed (T map)')
    endif

    allocate(sbuffer(nosubpix_map))
    allocate(ibuffer(nosubpix_map))
    if (write_cut) allocate(i8buffer(nosubpix_map))

    offset = 0
    writeoffset = 0
    do isubmap = 0, nosubmaps_tot-1

       if (sync_output) call wait_mpi

       call get_submap(ibuffer, nosubpix_map, mask, isubmap, idwr)

       if (id == idwr .and. write_cut) then
          nhit = 0
          do i = 1, nosubpix_map
             if (ibuffer(i) > 0) then
                nhit = nhit + 1
                i8buffer(nhit) = offset + i - 1
             end if
          end do
          call fits_write_column(out, 1, i8buffer(:nhit), writeoffset)
       else
          nhit = nosubpix_map
       end if

       do imap = 1,nmap

          call get_submap(sbuffer, nosubpix_map, map, nmap, imap, isubmap, idwr)

          if (ID==idwr) then
             if (write_cut) then
                ! Pack the observed pixels to the beginning
                nhit = 0
                do i = 1, nosubpix_map
                   if (ibuffer(i) > 0) then
                      nhit = nhit + 1
                      sbuffer(nhit) = sbuffer(i)
                   end if
                end do
             else
                ! Set unobserved pixels to undefined value
                where (ibuffer.le.0) sbuffer=Healpix_undef
             end if

             if (nhit > 0) then
                call fits_write_column(out, imap+coloffset, sbuffer(:nhit), &
                     writeoffset)
             end if
          end if
       end do

       offset = offset + nosubpix_map
       writeoffset = writeoffset + nhit
    end do

    if (ID==idwr) then
       call fits_close(out)
       write(*,*) 'Map written in '//trim(file_map)
    endif

    deallocate(sbuffer,ibuffer)
    if (write_cut) deallocate(i8buffer)

    cputime_output = cputime_output + get_time(10)

  END SUBROUTINE write_summap


  !----------------------------------------------------------------------


  SUBROUTINE write_binmap(binmap, mask)
    !
    ! Write the binned map into file.
    !
    real(dp),intent(in) :: binmap(nmap, nosubpix_map, nosubmaps)
    integer, intent(in) :: mask(nosubpix_map, nosubmaps)
    integer :: idwr, imap, isubmap, skycover, coloffset
    integer(idp) :: offset, writeoffset, nhit, i
    real(sp), allocatable :: sbuffer(:)
    integer, allocatable :: ibuffer(:)
    integer(i8b), allocatable :: i8buffer(:)
    type(fitshandle) :: out
    type(fitscolumn), pointer :: columns(:)

    real(sp), allocatable :: mapbuffer(:, :)
    character(len=SLEN) :: file_map_bin, file_map_txt
    INTEGER(i4b) :: status(MPI_STATUS_SIZE)
    INTEGER(i4b) :: fh, filemode, fileinfo, ierr, pix
    INTEGER(MPI_OFFSET_KIND) :: fileoffset

    INTEGER :: nsend, nrecv, mysubmap, rec_len, rec
    logical :: there
    INTEGER, allocatable :: recvcounts(:), displs(:)
    REAL(sp), allocatable :: sendbuffer(:), recvbuffer(:)

    if (.not. do_binmap .or. len_trim(file_binmap) == 0) return

    call reset_time(10)

    skycover = count(mask > 0)
    call sum_mpi(skycover)

    idwr = 1
    if (ntasks == 1) idwr = 0

    if (ID == idwr) then
       if (info > 1) write(*,*) 'Writing binned map...'

       inquire(file=trim(file_binmap), exist=there)

       if (there .and. binary_output .and. concatenate_binary) then
          print *,'Skipping existing header file and concatenating binary maps'
       else

          call check_outfile(file_binmap)

          call fits_create(out,file_binmap)

          call write_header_code(out)
          call write_header_simulation(out)

          if (write_cut) then
             allocate(columns(nmap+1))
             columns(1)%repcount = map_repcount
             columns(1)%name = 'Pixel index'
             columns(1)%type = fits_int8
             coloffset = 1
          else
             allocate(columns(nmap))
             coloffset = 0
          end if

          do imap = 1,nmap
             columns(imap+coloffset)%repcount = map_repcount
             columns(imap+coloffset)%unit = trim(unit_tod)
             columns(imap+coloffset)%type = fits_real4
             if (nmap == 1 .or. nmap == 3) then
                ! healpix naming
                if (imap==1) columns(imap+coloffset)%name = 'TEMPERATURE'
                if (imap==2) columns(imap+coloffset)%name = 'Q-POLARISATION'
                if (imap==3) columns(imap+coloffset)%name = 'U-POLARISATION'
             else if (nmap == 2) then
                if (imap==1) columns(imap+coloffset)%name  ='Q-POLARISATION'
                if (imap==2) columns(imap+coloffset)%name  ='U-POLARISATION'
             else
                columns(imap+coloffset)%name = addi('COMP', imap)
             end if
          enddo

          call fits_insert_bintab(out, columns)
          deallocate(columns)

          call fits_add_key(out, "objtype", "madam.binmap", 'Object type')
          call write_header_healpix(out, nside_map)
          call write_header_sky(out, skycover)
          call fits_add_key(out, 'MONOPOLE', rm_monopole, &
               'Monopole removed (T map)')
       endif

    end if

    if (binary_output) then

       file_map_bin = trim(adjustl(file_binmap)) // '.bin'
       file_map_txt = trim(adjustl(file_binmap)) // '.txt'

       if (nwrite_binary > 1) then
          file_map_bin = trim(file_map_bin) // '_' // addi('', writegroup)
       end if

       !call mpi_info_create(fileinfo, ierr)
       !filemode = ior(MPI_MODE_WRONLY, MPI_MODE_CREATE)
       !call mpi_file_open(mpi_comm_world, file_map_bin, filemode, fileinfo, &
       !    fh, ierr)
       !if (ierr /= 0) call abort_mpi('Failed to open binary map file')

       !fileoffset = 0
       !call mpi_file_set_view(fh, fileoffset, MPI_REAL4, MPI_REAL4, 'native', &
       !    fileinfo, ierr)
       !if (ierr /= 0) call abort_mpi('Failed to set view to binary map file')

       ! Collect the maps to root process for writing

       allocate(recvcounts(ntasks_bin), displs(ntasks_bin), stat=ierr)
       if (ierr /= 0) call abort_mpi('No room for recvcounts')

       nsend = 0
       do isubmap = 0, nosubmaps_tot-1
          if (id_submap(isubmap) == id) nsend = nsend + nmap * nosubpix_map
       end do

       recvcounts = 0
       call mpi_gather(nsend, 1, MPI_INTEGER, recvcounts, 1, MPI_INTEGER, 0, &
            comm_bin, ierr)

       allocate(sendbuffer(nsend), stat=ierr)
       if (ierr /= 0) call abort_mpi('No room for sendbuffer')

       mysubmap = 0
       i = 0
       do isubmap = 0, nosubmaps_tot-1
          if (id_submap(isubmap) == id) then
             mysubmap = mysubmap + 1
             do pix = 1, nosubpix_map
                if (mask(pix,mysubmap) <= 0) then
                   do imap = 1, nmap
                      i = i + 1
                      sendbuffer(i) = Healpix_undef
                   end do
                else
                   do imap = 1, nmap
                      i = i + 1
                      sendbuffer(i) = binmap(imap,pix,mysubmap)
                   end do
                end if
             end do
          end if
       end do

       if (i /= nsend) call abort_mpi('Did not fill sendbuffer') ! redundant check

       if (id_bin == 0) then
          nrecv = sum(recvcounts)

          allocate(recvbuffer(nrecv), stat=ierr)
          if (ierr /= 0) call abort_mpi('No room for recvbuffer')

          displs = 0
          do i = 2, ntasks_bin
             displs(i) = recvcounts(i-1) + displs(i-1)
          end do
       end if

       call mpi_gatherv(sendbuffer, nsend, MPI_REAL4, recvbuffer, recvcounts, &
            displs, MPI_REAL4, 0, comm_bin, ierr)

       if (id_bin == 0) then

          !fileoffset = count(id_submap < id) * nosubpix_map * nmap
          !call mpi_file_write_at(fh, fileoffset, recvbuffer, nrecv, &
          !    MPI_REAL4, status, ierr)

          inquire(iolength=rec_len) recvbuffer
          open(unit=55, file=trim(file_map_bin), action='readwrite', &
               form='unformatted', access='direct', recl=rec_len)
          !inquire(unit=55, nextrec=rec)
          !print *,id,' : nextrec = ',rec ! DEBUG
          write(55, rec=record_number) recvbuffer
          close(unit=55)

       end if

       if (id == 0) then

          ! Finally, store the submap distribution so that the map may be
          ! assembled later

          inquire(file=trim(file_map_txt), exist=there)

          if (.not. there) then

             open(unit=55, file=trim(file_map_txt), status='replace', &
                  form='formatted')
             write(55, *) nosubmaps_tot, ntasks
             write(55, *) nosubpix_map, nmap
             do isubmap = 0, nosubmaps_tot-1
                write(55, *) isubmap, id_submap(isubmap)
             end do
             close(unit=55)

          end if

       end if

       deallocate(recvcounts, displs, sendbuffer)
       if (id_bin == 0) deallocate(recvbuffer)

       !call mpi_file_close(fh, ierr)
       !if (ierr /= 0) call abort_mpi('Failed to close binary file')

    else

       allocate(sbuffer(nosubpix_map))
       allocate(ibuffer(nosubpix_map))
       if (write_cut) allocate(i8buffer(nosubpix_map))

       offset = 0
       writeoffset = 0
       do isubmap = 0, nosubmaps_tot-1

          if (sync_output) call wait_mpi

          call get_submap(ibuffer, nosubpix_map, mask, isubmap, idwr)

          if (id == idwr .and. write_cut) then
             nhit = 0
             do i = 1, nosubpix_map
                if (ibuffer(i) > 0) then
                   nhit = nhit + 1
                   i8buffer(nhit) = offset + i - 1
                end if
             end do
             call fits_write_column(out, 1, i8buffer(:nhit), writeoffset)
          else
             nhit = nosubpix_map
          end if

          do imap = 1,nmap

             !print *,id,' : getting submap ',isubmap,imap ! debug
             call get_submap(sbuffer, nosubpix_map, binmap, nmap, imap, &
                  isubmap, idwr)

             if (ID==idwr) then
                if (write_cut) then
                   ! Pack the observed pixels to the beginning
                   nhit = 0
                   do i = 1, nosubpix_map
                      if (ibuffer(i) > 0) then
                         nhit = nhit + 1
                         sbuffer(nhit) = sbuffer(i)
                      end if
                   end do
                else
                   ! Set unobserved pixels to undefined value
                   where (ibuffer <= 0) sbuffer = Healpix_undef
                end if

                if (nhit > 0) then
                   call fits_write_column(out, imap+coloffset, sbuffer(:nhit), &
                        writeoffset)
                end if
             endif
          end do

          offset = offset + nosubpix_map
          writeoffset = writeoffset + nhit
       end do

       deallocate(sbuffer,ibuffer)
       if (write_cut) deallocate(i8buffer)

    end if

    if (ID == idwr) then
       call fits_close(out)
       write(*,*) 'Binned map written in ' // trim(file_binmap)
    endif

    cputime_output = cputime_output + get_time(10)

  END SUBROUTINE write_binmap


  !----------------------------------------------------------------------


  SUBROUTINE write_mask(mask,crit)
    !
    ! Write into file mask and criterion use to discard badly defined
    ! pixels.
    !
    integer, intent(in) :: mask(nosubpix_map, nosubmaps)
    real(sp),intent(in) :: crit(nosubpix_map, nosubmaps)
    integer :: idwr, isubmap, skycover, coloffset
    integer(idp) :: offset, writeoffset, nhit, i
    real(sp), allocatable :: sbuffer(:)
    integer, allocatable  :: ibuffer(:)
    integer(i8b),allocatable :: i8buffer(:)
    type(fitshandle) :: out
    type(fitscolumn), pointer :: columns(:)

    if (.not.do_mask .or. len_trim(file_mask) == 0) return

    call reset_time(10)

    skycover = count(mask > 0)
    call sum_mpi(skycover)

    idwr = 2
    if (ntasks < 3) idwr = 0

    if (ID == idwr) then
       if (info > 1) write(*,*) 'Writing mask...'

       call check_outfile(file_mask)

       call fits_create(out,file_mask)

       call write_header_code(out)
       call write_header_simulation(out)
       call write_header_destriping(out)

       allocate(columns(2))

       if (write_cut) then
          columns(1)%repcount = map_repcount
          columns(1)%type = fits_int8
          columns(1)%name = 'Pixel index'
       else
          ! The mask is redundant in cut sky format
          columns(1)%repcount = map_repcount
          columns(1)%type = fits_int4
          columns(1)%name = 'Mask'
       end if

       columns(2)%repcount = map_repcount
       columns(2)%type = fits_real4
       columns(2)%name = 'criterion'

       call fits_insert_bintab(out, columns)
       deallocate(columns)

       call fits_add_key(out, "objtype", "madam.mask", 'Object type')

       call write_header_healpix(out,nside_map)
       call write_header_sky(out,skycover)
    endif

    allocate(sbuffer(nosubpix_map))
    allocate(ibuffer(nosubpix_map))
    if (write_cut) allocate(i8buffer(nosubpix_map))

    offset = 0
    writeoffset = 0
    do isubmap = 0, nosubmaps_tot-1

       if (sync_output) call wait_mpi

       call get_submap(ibuffer, nosubpix_map, mask, isubmap, idwr)

       if (id == idwr .and. write_cut) then
          nhit = 0
          do i = 1, nosubpix_map
             if (ibuffer(i) > 0) then
                nhit = nhit + 1
                i8buffer(nhit) = offset + i - 1
             end if
          end do
          call fits_write_column(out, 1, i8buffer(:nhit), writeoffset)
       else
          nhit = nosubpix_map
       end if

       call get_submap(sbuffer, nosubpix_map, crit, isubmap, idwr)

       if (ID==idwr) then

          if (write_cut) then
             ! Pack the observed pixels to the beginning
             nhit = 0
             do i = 1, nosubpix_map
                if (ibuffer(i) > 0) then
                   nhit = nhit + 1
                   sbuffer(nhit) = sbuffer(i)
                end if
             end do
             call fits_write_column(out, 2, sbuffer(1:nhit), writeoffset)
          else
             call fits_write_column(out, 1, ibuffer, offset)
             call fits_write_column(out, 2, sbuffer, offset)
          end if
       end if

       offset = offset + nosubpix_map
       writeoffset = writeoffset + nhit
    end do

    if (ID==idwr) then
       call fits_close(out)
       write(*,*) 'Mask written in '//trim(file_mask)
    end if

    deallocate(sbuffer,ibuffer)

    cputime_output = cputime_output + get_time(10)

  END SUBROUTINE write_mask


  !----------------------------------------------------------------------


  SUBROUTINE write_hits(hits)
    !
    ! Write hit map into file.
    !
    integer, intent(in) :: hits(nosubpix_map, nosubmaps,*)
    integer :: idwr, isubmap, idet, coloffset
    integer(idp) :: offset, nhit, buflen, i
    integer, allocatable :: ibuffer(:), total(:)
    integer(i8b), allocatable :: i8buffer(:)
    type(fitshandle) :: out
    type(fitscolumn), pointer :: columns(:)
    integer :: npix, ierr, repeat
    integer(i8b) :: bufoffset, writeoffset

    if (.not. do_hits .or. len_trim(file_hit) == 0) return

    repeat = max(12*1024, nosubpix_map)
    npix = nosubmaps_tot * nosubpix_map
    repeat = min(repeat, npix)
    buflen = 2 * repeat

    if (write_cut) repeat = 1

    call reset_time(10)

    idwr = 1
    if (ntasks < 2) idwr = 0

    if (ID == idwr) then
       if (info > 1) write(*,*) 'Writing hits...'

       call check_outfile(file_hit)

       call fits_create(out,file_hit)

       call write_header_code(out)
       call write_header_simulation(out)

       if (write_cut) then
          allocate(columns(2))
          columns(1)%repcount = repeat
          columns(1)%name = 'Pixel index'
          columns(1)%type = fits_int8
          coloffset = 1
       else
          allocate(columns(1))
          coloffset = 0
       end if

       columns(1+coloffset)%repcount = repeat ! map_repcount
       columns(1+coloffset)%name='Total hits'
       columns(1+coloffset)%type=fits_int4

       call fits_insert_bintab(out, columns)
       deallocate(columns)

       call fits_add_key(out, 'objtype', 'madam.hits', 'Object type')
       call write_header_healpix(out, nside_map)

    endif

    allocate(ibuffer(nosubpix_map))
    if (write_cut) allocate(i8buffer(nosubpix_map))

    offset = 0
    writeoffset = 0
    do isubmap = 0, nosubmaps_tot-1

       if (sync_output) call wait_mpi ! clean the unexpected message buffers

       call get_submap_int(ibuffer, nosubpix_map, hits(1,1,1), 1, 1, isubmap, &
            idwr)

       if (ID == idwr) then
          if (write_cut) then
             nhit = 0
             do i = 1, nosubpix_map
                if (ibuffer(i) > 0) then
                   nhit = nhit + 1
                   i8buffer(nhit) = offset + i - 1
                   ibuffer(nhit) = ibuffer(i)
                end if
             end do
             if (nhit > 0) then
                call fits_write_column(out, 1, i8buffer(:nhit), writeoffset)
                call fits_write_column(out, 2, ibuffer(:nhit), writeoffset)
             end if
          else
             nhit = nosubpix_map
             call fits_write_column(out, 1, ibuffer(:nhit), writeoffset)
          end if

          writeoffset = writeoffset + nhit
       end if

       offset = offset + nosubpix_map
    end do

    if (ID == idwr) then
       call fits_close(out)

       write(*,*) 'Hit count written in ' // trim(file_hit)
    end if

    deallocate(ibuffer)

    cputime_output = cputime_output + get_time(10)

  END SUBROUTINE write_hits


  !----------------------------------------------------------------------


  SUBROUTINE write_matrix(cc, mask)
    !
    ! Write the pixel matrices into file
    !
    real(dp),intent(in) :: cc(nmap, nmap, nosubpix_map, nosubmaps)
    integer, intent(in), optional :: mask(nosubpix_map, nosubmaps)
    integer :: idwr, isubmap, imap, jmap, coloffset
    integer(idp) :: offset, writeoffset, i, j, nhit
    real(dp), allocatable :: dbuffer(:)
    integer, allocatable :: ibuffer(:)
    integer(i8b), allocatable :: i8buffer(:)
    type(fitshandle) :: out
    type(fitscolumn), pointer :: columns(:)
    character(len=SLEN) :: outfile

    outfile = file_matrix
    if (present(mask)) outfile = file_wcov

    if (len_trim(outfile).le.0) return

    call reset_time(10)

    idwr = 2
    if (ntasks < 3) idwr = 0

    if (ID==idwr) then
       if (info > 1) write(*,*) 'Writing pixel matrix...'

       call check_outfile(outfile)

       call fits_create(out, outfile)

       call write_header_code(out)
       call write_header_simulation(out)

       if (write_cut) then
          allocate(columns(ncc+1))
          columns(1)%repcount = map_repcount
          columns(1)%name = 'Pixel index'
          columns(1)%type = fits_int8
          coloffset = 1
       else
          allocate(columns(ncc))
          coloffset = 0
       end if

       do imap = 1,ncc
          columns(imap+coloffset)%repcount = map_repcount
          columns(imap+coloffset)%type = fits_real8
       enddo

       do i = 1,nmap
          do j = i,nmap
             coloffset = coloffset + 1
             write(columns(coloffset)%name, '("cc(",i0,",",i0,")")') i, j
          end do
       end do

       call fits_insert_bintab(out, columns)
       deallocate(columns)

       call fits_add_key(out, "objtype", "madam.matrix", 'Object type')
       call write_header_healpix(out,nside_map)
    endif

    allocate(dbuffer(nosubpix_map))
    if (write_cut) then
       allocate(ibuffer(nosubpix_map))
       allocate(i8buffer(nosubpix_map))
    end if

    offset = 0
    writeoffset = 0
    do isubmap = 0, nosubmaps_tot-1

       if (sync_output) call wait_mpi

       if (write_cut) then
          if (present(mask)) then
             call get_submap(ibuffer, nosubpix_map, mask, isubmap, idwr)
          else
             ibuffer = 0
             call get_submap(dbuffer, nosubpix_map, cc, nmap, 1, 1, isubmap, &
                  idwr)
             where(dbuffer /= 0) ibuffer = 1
          end if

          if (id == idwr) then
             nhit = 0
             do i = 1, nosubpix_map
                if (ibuffer(i) > 0) then
                   nhit = nhit + 1
                   i8buffer(nhit) = offset + i - 1
                end if
             end do
             if (nhit > 0) &
                  call fits_write_column(out, 1, i8buffer(:nhit), writeoffset)
          end if
       else
          nhit = nosubpix_map
       end if

       coloffset = coloffset - ncc
       do imap = 1,nmap
          do jmap = imap,nmap

             call get_submap(dbuffer, nosubpix_map, cc, nmap, imap, jmap, &
                  isubmap, idwr)

             if (ID == idwr) then
                if (write_cut) then
                   ! Pack the observed pixels to the beginning
                   nhit = 0
                   do i = 1, nosubpix_map
                      if (ibuffer(i) > 0) then
                         nhit = nhit + 1
                         dbuffer(nhit) = dbuffer(i)
                      end if
                   end do
                end if

                coloffset = coloffset + 1
                if (nhit > 0) then
                   call fits_write_column(out, coloffset, dbuffer(:nhit), &
                        writeoffset)
                end if
             end if
          end do
       end do

       offset = offset + nosubpix_map
       writeoffset = writeoffset + nhit
    end do

    if (ID==idwr) then
       call fits_close(out)
       write(*,*) 'Pixel matrix written in ' // trim(outfile)
    endif

    deallocate(dbuffer)
    if (write_cut) deallocate(ibuffer, i8buffer)

    cputime_output = cputime_output + get_time(10)

  END SUBROUTINE write_matrix


  !----------------------------------------------------------------------


  SUBROUTINE write_leakmatrix(cc, mask, detector_name)
    !
    ! Write the leakage matrices into file
    !
    real(dp),intent(in) :: cc(nmap,nmap, nosubpix_map, nosubmaps)
    integer, intent(in) :: mask(nosubpix_map, nosubmaps)
    character(len=*), intent(in) :: detector_name
    integer :: idwr, isubmap, imap, jmap, coloffset, ncol
    integer(idp):: offset, writeoffset, i, j, nhit
    real(dp), allocatable :: dbuffer(:)
    integer, allocatable  :: ibuffer(:)
    integer(i8b), allocatable :: i8buffer(:)
    type(fitshandle) :: out
    type(fitscolumn), pointer :: columns(:)
    character(len=SLEN) :: outfile

    if (len_trim(file_leakmatrix) <= 0) return

    outfile = trim(file_leakmatrix) // '_' // trim(detector_name) // '.fits'

    call reset_time(10)

    idwr = 2
    if (ntasks < 3) idwr=0

    ncol = nmap * nmap

    if (ID==idwr) then
       if (info > 1) then
          write(*,*) 'Writing leakage matrix for ', trim(detector_name), '...'
       end if

       call check_outfile(outfile)

       call fits_create(out, outfile)

       call write_header_code(out)
       call write_header_simulation(out)

       if (write_cut) then
          allocate(columns(ncol+1))
          columns(1)%repcount = map_repcount
          columns(1)%name = 'Pixel index'
          columns(1)%type = fits_int8
          coloffset = 1
       else
          allocate(columns(ncol))
          coloffset = 0
       end if

       do imap = 1,ncol
          columns(imap + coloffset)%repcount = map_repcount
          columns(imap + coloffset)%type = fits_real8
       enddo

       do i = 1,nmap
          do j = 1,nmap
             coloffset = coloffset + 1
             write(columns(coloffset)%name, '("gamma(",i0,",",i0,")")') i, j
          end do
       end do

       call fits_insert_bintab(out, columns)
       deallocate(columns)

       call fits_add_key(out, "objtype", "madam.leakmatrix", 'Object type')
       call write_header_healpix(out, nside_map)
    endif

    allocate(dbuffer(nosubpix_map))
    if (write_cut) then
       allocate(ibuffer(nosubpix_map))
       allocate(i8buffer(nosubpix_map))
    end if

    offset = 0
    writeoffset = 0
    do isubmap = 0, nosubmaps_tot-1

       if (sync_output) call wait_mpi

       if (write_cut) then
          call get_submap(ibuffer, nosubpix_map, mask, isubmap, idwr)

          if (id == idwr) then
             nhit = 0
             do i = 1, nosubpix_map
                if (ibuffer(i) > 0) then
                   nhit = nhit + 1
                   i8buffer(nhit) = offset + i - 1
                end if
             end do
             if (nhit > 0) then
                call fits_write_column(out, 1, i8buffer(:nhit), writeoffset)
             end if
          end if
       else
          nhit = nosubpix_map
       end if

       coloffset = coloffset - ncol
       do imap = 1,nmap
          do jmap = 1,nmap

             call get_submap(dbuffer, nosubpix_map, cc, nmap, imap, jmap, &
                  isubmap, idwr)

             if (ID == idwr) then
                if (write_cut) then
                   ! Pack the observed pixels to the beginning
                   nhit = 0
                   do i = 1, nosubpix_map
                      if (ibuffer(i) > 0) then
                         nhit = nhit + 1
                         dbuffer(nhit) = dbuffer(i)
                      end if
                   end do
                end if

                coloffset = coloffset + 1
                if (nhit > 0) then
                   call fits_write_column(out, coloffset, dbuffer(:nhit), &
                        writeoffset)
                end if
             end if
          end do
       end do

       offset = offset + nosubpix_map
       writeoffset = writeoffset + nhit
    end do

    if (ID==idwr) then
       call fits_close(out)
       write(*,*) 'Leakage matrix written in ' // trim(outfile)
    endif

    deallocate(dbuffer)
    if (write_cut) deallocate(ibuffer, i8buffer)

    cputime_output = cputime_output + get_time(10)

  END SUBROUTINE write_leakmatrix


  !------------------------------------------------------------------------------


  SUBROUTINE check_outfile(filename)
    !
    ! Create a new filename if file exists.
    !
    character(len=*),intent(inout) :: filename
    logical :: there
    integer :: k, n
    character(len=20) :: s, ending

    if (filename(1:1) == '!') return

    inquire(file=filename, exist=there)

    if (.not.there) return

    n = index(filename, '.', .true.)
    if (n <= 0) then
       ending = ''
       n = len_trim(filename) + 1
    else
       ending = filename(n:len_trim(filename))
    endif

    do k = 1,1000

       if (k==1000) then
          filename = ''
          return
       endif

       write(s,'("_",i3.3)') k
       filename = filename(1:n-1) // trim(s) // trim(ending)

       inquire(file=filename, exist=there)
       if (.not. there) exit
    enddo

  END SUBROUTINE check_outfile


  !------------------------------------------------------------------------------
  !
  !   Routines for writing groups of keywords
  !
  !------------------------------------------------------------------------------


  SUBROUTINE write_header_healpix(out, nside)

    type(fitshandle) :: out
    integer :: nside

    call fits_add_comment(out, '----------------------------------')
    call fits_add_comment(out, '   Healpix Map Specific Keywords  ')
    call fits_add_comment(out, '----------------------------------')

    call fits_add_key(out, 'PIXTYPE', 'HEALPIX','HEALPIX Pixelisation')
    call fits_add_key(out, 'ORDERING','NESTED', 'Pixel ordering scheme')
    call fits_add_key(out, 'NSIDE', nside, 'Resolution parameter for HEALPIX')
    call fits_add_key(out, 'FIRSTPIX', 0, 'First pixel # (0 based)')
    call fits_add_key(out, 'LASTPIX', 12*nside**2-1, 'Last pixel # (0 based)')
    call fits_add_key(out, 'POLAR', (nmap/=1),  &
         'Polarisation included (True/False)')
    call fits_add_key(out,'POLCONV','COSMO',   &
         'Coord. convention for polarisation (COSMO/IAU)')
    if (write_cut) then
       call fits_update_key(out, 'EXTNAME', "'CUT SKY MAP'")
       call fits_add_key(out, 'COMMENT', 'Cut sky data')
       call fits_add_key(out, 'OBJECT', 'PARTIAL')
       call fits_add_key(out, 'INDXSCHM', 'EXPLICIT', &
            ' Indexing : IMPLICIT or EXPLICIT')
       call fits_add_key(out, 'GRAIN', 1, ' Grain of pixel indexing')
    else
       call fits_update_key(out, 'EXTNAME', "'FULL SKY MAP'")
       call fits_add_key(out, 'COMMENT', 'Full sky data')
       call fits_add_key(out, 'OBJECT', 'FULLSKY')
       call fits_add_key(out, 'INDXSCHM', 'IMPLICIT', &
            ' Indexing : IMPLICIT or EXPLICIT')
       call fits_add_key(out, 'GRAIN', 0, ' Grain of pixel indexing')
    end if
    call fits_add_key(out, 'COMMENT', &
         'GRAIN=0 : no indexing of pixel data                         (IMPLICIT)')
    call fits_add_key(out, 'COMMENT', &
         'GRAIN=1 : 1 pixel index -> 1 pixel data                     (EXPLICIT)')
    call fits_add_key(out, 'COMMENT', &
         'GRAIN>1 : 1 pixel index -> data of GRAIN consecutive pixels (EXPLICIT)')

  END SUBROUTINE write_header_healpix


  !-----------------------------------------------------------------------------


  SUBROUTINE write_header_code(out)

    type(fitshandle) :: out

    call fits_add_comment(out,'----------------------------------')
    call fits_add_comment(out,'         CODE             ')
    call fits_add_comment(out,'----------------------------------')
    call fits_add_key(out, 'CREATOR', 'Madam', 'Software creating FITS file')
    call fits_add_key(out, 'VERSION', version, 'Software version')

    if (parallel_version) then
       call fits_add_comment(out, 'Parallelized version')
       call fits_add_key(out,'NOPROCS', ntasks, 'Number of processes')
    else
       call fits_add_comment(out,'Serial version')
    endif

    call fits_add_key(out, 'PARFILE', file_param, 'Parameter file')

  END SUBROUTINE write_header_code


  !------------------------------------------------------------------------------


  SUBROUTINE write_header_sky(out, skycover)

    type(fitshandle) :: out
    integer, intent(in) :: skycover

    call fits_add_comment(out, '----------------------------------')
    call fits_add_comment(out, '         Sky coverage             ')
    call fits_add_comment(out, '----------------------------------')
    call fits_add_key(out, 'PIXCOV', skycover, 'Recovered pixels')
    call fits_add_key(out, 'SKYCOV', real(skycover/(12.d0*nside_map**2)),  &
         'Sky coverage')
    call fits_add_key(out, 'PIXMMAP', pixmode_map,  &
         'Pixel removal criterion')
    call fits_add_key(out, 'PIXLMAP', real(pixlim_map), 'Pixel removal limit')

  END SUBROUTINE write_header_sky


  !---------------------------------------------------------------------------


  SUBROUTINE write_header_simulation(out)

    type(fitshandle) :: out
    real(dp) :: samples_to_d
    integer :: i, j, k

    samples_to_d = 1.d0/fsample/24.d0/3600.d0

    call fits_add_comment(out, '-------------------------------')
    call fits_add_comment(out, '          Simulation           ')
    call fits_add_comment(out, '-------------------------------')

    call fits_add_key(out, 'FSAMPLE', real(fsample),'Sampling frequency (Hz)')
    call fits_add_key(out, 'TMISSION', real(nosamples_tot*samples_to_d), &
         'Mission duration/days')
    call fits_add_key(out,'TSTART', real(istart_mission*samples_to_d), &
         'Mission start/days')
    call fits_add_key(out, 'IMISSION', nosamples_tot, &
         'Mission duration/samples')
    call fits_add_key(out,'ISTART', istart_mission, &
         'Mission start/samples')
    call fits_add_key(out, 'NSUBCHNK', nsubchunk,'Number of subchunks')
    call fits_add_key(out, 'ISUBCHNK', isubchunk,'Index of subchunk')
    call fits_add_comment(out, '---------------------------------')
    call fits_add_comment(out, '         Detector info           ')
    call fits_add_comment(out, '---------------------------------')

    call fits_add_key(out, 'NODET', nodetectors, 'Number of detectors')
    call fits_add_key(out, 'MDETW', mode_detweight, 'Detector weighting mode')

    do i = 1, nodetectors
       call fits_add_key(out, addi('DETNAM',i), detectors(i)%name, &
            'Detector name')
    enddo
    call fits_add_key(out, 'TEMPONLY', temperature_only, 'Temperature map only')

    call fits_add_comment(out, '----------------------------')
    call fits_add_comment(out, '          Pointing          ')
    call fits_add_comment(out, '----------------------------')

  END SUBROUTINE write_header_simulation


  !--------------------------------------------------------------------------


  SUBROUTINE write_header_destriping(out)

    type(fitshandle) :: out
    integer :: i

    call fits_add_comment(out, '---------------------------')
    call fits_add_comment(out, '     Data distribution     ')
    call fits_add_comment(out, '---------------------------')

    call fits_add_comment(out, 'Standard mode')

    call fits_add_key(out, 'NSIDESUB', nside_submap, 'Submap resolution')

    call fits_add_comment(out, '-------------------------------------------')
    call fits_add_comment(out, '        Destriping scheme                  ')
    call fits_add_comment(out, '-------------------------------------------')

    call fits_add_key(out, 'KFIRST', kfirst, 'First destriping phase on/off')

    if (kfirst) then
       call fits_add_key(out, 'NSHORT', nshort, 'Baseline length ')
       call fits_add_key(out, 'KFILTER', kfilter, 'Noise filter on/off')
       select case(basis_func)
       case (basis_poly)
          call fits_add_key(out, 'BASIS', 'POLYNOMIAL', &
               'Destriping function basis')
       case (basis_fourier)
          call fits_add_key(out, 'BASIS', 'FOURIER', &
               'Destriping function basis')
       case (basis_cheby)
          call fits_add_key(out, 'BASIS', 'CHEBYSHEV', &
               'Destriping function basis')
       case (basis_legendre)
          call fits_add_key(out, 'BASIS', 'LEGENDRE', &
               'Destriping function basis')
       case default
          call abort_mpi('Unknown function basis')
       end select
       call fits_add_key(out, 'BORDER', basis_order, &
            'Destriping function order')
    endif

    call fits_add_key(out, 'NSIDECR', nside_cross, 'Crossing point resolution')
    call fits_add_key(out, 'PIXMCR', pixmode_cross, 'Pixel removal criterion')
    call fits_add_key(out, 'PIXLCR', real(pixlim_cross), 'Pixel removal limit')
    call fits_add_key(out, 'CGLIMIT',real(cglimit), &
         'Conjugate gradient convergence limit')
    call fits_add_key(out,'ITERMAX',iter_max, &
         'Maximum number of iterations')

    if (use_inmask) call fits_add_key(out, 'INMASK', ftrim(file_inmask), &
         'input mask')

    if (kfilter) then

       call fits_add_comment(out, '---------------------------------------')
       call fits_add_comment(out, '           Noise filter                ')
       call fits_add_comment(out, '---------------------------------------')

       call fits_add_key(out, 'FILTTIM', real(filter_time), &
            'Noise filter length/seconds')
       call fits_add_key(out, 'TAILTIM', real(tail_time), 'Overlap/seconds')

       if (len_trim(file_spectrum) >= 0) then
          do i = 1, nodetectors
             call fits_add_comment(out, detectors(i)%name)
             call fits_add_key(out, addi('SIGMA', i), &
                  real(detectors(i)%sigmas(1)), 'White noise std')
             call fits_add_key(out, addi('SLOPE', i), &
                  real(detectors(i)%slope), 'Spectral slope')
             call fits_add_key(out, addi('FKNEE', i), &
                  real(detectors(i)%fknee), 'Knee frequency (Hz)')
             call fits_add_key(out, addi('FMIN', i), &
                  real(detectors(i)%fmin), 'Minimum frequency (Hz)')
          end do
       else
          call fits_add_key(out, 'FSPEC', trim(file_spectrum), &
               'Noise spectrum file')
       end if

    end if

  END SUBROUTINE write_header_destriping

  !------------------------------------------------------------------------------


  SUBROUTINE close_output

    integer :: ierr

    if (binary_output) then

       if (nwrite_binary > 1) then
          call mpi_comm_free(comm_bin, ierr)
          if (ierr /= 0) then
             call abort_mpi('Failed to close binary output communicator')
          end if
       end if

    end if

    return

  END SUBROUTINE close_output


  !------------------------------------------------------------------------------
  !
  ! Functions for handling header keywords
  !
  !------------------------------------------------------------------------------


  FUNCTION degs(angle)

    real(dp) :: degs
    real(dp) :: angle

    degs = angle / pi * 180.d0

  END FUNCTION degs


  !-----------------------------------------------------------------------


  FUNCTION ftrim(filename)
    ! remove path from a file name
    character(len=68) :: ftrim
    character(len=*)  :: filename
    integer :: n

    n = index(filename, '/', .true.)

    ftrim = filename(n+1:len_trim(filename))

  END FUNCTION ftrim


  !------------------------------------------------------------------------------


  FUNCTION addi(key_in, i) result(key)

    character(len=*) :: key_in
    character(len=8) :: key
    integer :: i, ilen

    write(key, '(i0)') i
    ilen = len_trim(key)

    if (len_trim(key_in) + ilen > 8) then
       write (key,'(a,i0)') key_in(:8-ilen), i
    else
       write (key,'(a,i0)') key_in, i
    end if

  END FUNCTION addi


  FUNCTION addi2(key_in, i) result(key)

    character(len=*) :: key_in
    character(len=24) :: key
    integer :: i, ilen

    write(key, '(i0)') i
    ilen = len_trim(key)

    if (len_trim(key_in) + ilen > 24) then
       write (key,'(a,i0)') key_in(:24-ilen), i
    else
       write (key,'(a,i0)') key_in, i
    end if

  END FUNCTION addi2


  !--------------------------------------------------

END MODULE output
