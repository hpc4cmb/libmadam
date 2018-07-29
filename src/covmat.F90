! Rework of covmat.F90
!
! Old module initiatied on 2007-05-31
!
! New module initiated on 2013-11-11
! Key changes:
!   - calculation is performed per pointing period. No correlations beyond the
!     pointing period boundaries
!   - allow for non-circulant middle matrix
!
! This module implements the formula (28) in astro-ph/0412517:
!
! C_m^-1 = P^T C_n^-1 P - P^T C_n^-1 F (C_a^-1 +  F^T C_n^-1 F)^-1 F^T C_n^-1 P
!        = P^T(C_n^-1 - C_n^-1 F (C_a^-1 +  F^T C_n^-1 F)^-1 F^T C_n^-1)P
!       := P^T F middle_matrix F^T P

MODULE covmat

  USE commonparam
  USE planck_config, ONLY : i4b, i8b, sp, dp
  USE fourier
  USE pointing
  USE covmat_util, ONLY : tic, toc
  USE mpi_wrappers
  USE noise_routines
  USE maps_and_baselines

  IMPLICIT NONE

  real(dp), save, public :: memory_ncm = 0, cputime_ncm=0, cputime_hitmap=0, &
       cputime_ptf=0, cputime_middlematrix=0, cputime_accumulate=0, &
       cputime_white=0, cputime_save_matrix=0, cputime_symmetrize=0, &
       cputime_invert_middlematrix=0, cputime_symmetrize_middlematrix=0

  real(dp), allocatable :: local_covmat(:,:,:,:), local_ptf(:,:,:), &
       middlematrix(:,:)
  integer(i4b), allocatable :: local_hitmap(:), local_basehitmap(:,:)
  integer :: wband



contains



  subroutine write_covmat(outroot)
    character(len=*), intent(in) :: outroot

    integer(i4b) :: ierr, ipix, ibase, jpix, jbase, i, j, kstart, noba, lag
    integer(i8b) :: idet, ipsd, ipsd_det, ichunk
    real(dp) :: detweight, mm, ptf1
    character(len=SLEN) :: outfile
    logical :: there

    call reset_time(5)

    ! initialize workspace

    allocate(local_covmat(nmap, nmap, 0:nolocpix, 0:nolocpix), &
         local_hitmap(0:nolocpix), stat=ierr)
    if (ierr /= 0) then
       print *, 'Failed to allocate ', &
            (nmap**2*(nolocpix+1)**2*8 + (nolocpix+1)*4)/2**20, &
            ' MB for local_covmat. nolocpix = ', nolocpix
       call abort_mpi('No room for local P^T F')
    end if
    if (ierr /= 0) call abort_mpi('No room for local covmat')

    memory_ncm = memory_ncm + nmap**2*(nolocpix+1)**2*8 + (nolocpix+1)*4

    ! accumulate the matrix

    loop_detector : do idet = 1, nodetectors

       outfile = trim(outroot) // '_' // trim(detectors(idet)%name) &
            // '_inv_bin.dat'

       INQUIRE(file=outfile, exist=there)
       IF (there) THEN
          IF (id == 0) WRITE (*,'(1x,a," EXISTS! skipping ...")') TRIM(outfile)
          cycle loop_detector
       ELSE
          IF (id == 0) WRITE (*,'(1x,"Computing ",a)') trim(outfile)
       END IF

       call tic

       if (id == 0) print *,'Processing detector ', trim(detectors(idet)%name)

       local_covmat = 0

       call tic(555)

       loop_chunk : do ichunk = 1, ninterval ! global indices
          noba = noba_short_pp(ichunk) ! baselines on this pointing period
          kstart = sum(noba_short_pp(1:ichunk-1)) ! first baseline, local index
          ipsd = psd_index(idet, baselines_short_time(kstart+1))
          ipsd_det = psd_index_det(idet, baselines_short_time(kstart+1))
          if (ipsd_det < 0) cycle loop_chunk
          ! ALWAYS use the optimal weights, even if the map has some
          ! other weighting scheme
          detweight = 1 / detectors(idet)%sigmas(ipsd_det)**2
          if (detweight < tiny(detweight)) cycle loop_chunk

          call get_ptf(idet, kstart, noba, detweight)

          call get_hitmap(idet, kstart, noba)

          call get_middlematrix(fcov, ipsd, detweight, noba)

          ! double loop over all non-zeros in P^T F

          call reset_time(10)

          loop_pix : do ipix = 0, nolocpix-1
             if (local_hitmap(ipix) == 0) cycle loop_pix
             loop_pix2 : do jpix = ipix, nolocpix-1
                if (local_hitmap(jpix) == 0) cycle loop_pix2
                loop_baseline : do i = 1, local_hitmap(ipix)
                   ibase = local_basehitmap(i, ipix)
                   loop_baseline2 : do j = 1, local_hitmap(jpix)
                      jbase = local_basehitmap(j, jpix)

                      lag = abs(jbase - ibase) + 1
                      mm = -middlematrix(lag, 1)

                      if (nmap == 1) then
                         local_covmat(1, 1, jpix, ipix) = &
                              local_covmat(1, 1, jpix, ipix) &
                              + mm * local_ptf(1, ibase, ipix) &
                              * local_ptf(1, jbase, jpix)
                      else
                         ptf1 = mm * local_ptf(1, ibase, ipix)
                         local_covmat(1:3, 1, jpix, ipix) = &
                              local_covmat(1:3, 1, jpix, ipix) &
                              + local_ptf(1:3, jbase, jpix) * ptf1

                         ptf1 = mm * local_ptf(2, ibase, ipix)
                         local_covmat(1:3, 2, jpix, ipix) = &
                              local_covmat(1:3, 2, jpix, ipix) &
                              + local_ptf(1:3, jbase, jpix) * ptf1

                         ptf1 = mm * local_ptf(3, ibase, ipix)
                         local_covmat(1:3, 3, jpix, ipix) = &
                              local_covmat(1:3, 3, jpix, ipix) &
                              + local_ptf(1:3, jbase, jpix) * ptf1
                      end if

                   end do loop_baseline2
                end do loop_baseline
             end do loop_pix2
          end do loop_pix

          cputime_accumulate = cputime_accumulate + get_time(10)

       end do loop_chunk

       ! Free workspace

       if (allocated(local_ptf)) deallocate(local_ptf)
       if (allocated(middlematrix)) deallocate(middlematrix)

       ! add the white noise contribution, use BLAS for the outer product

       call add_white(idet, detweight)

       call symmetrize()

       ! save the matrix

       call wait_mpi

       if (id == 0) call toc('Accumulate matrix', 555)

       if (id == 0) then
          print *,'Saving the matrix to ', trim(outroot) // '_' &
               // trim(detectors(idet)%name) // '_inv_bin.dat'
       end if

       call tic(555)

       call save_matrix(outfile)

       if (id == 0) call toc('Save matrix', 555)

       if (id == 0) call toc('Process detector')

    end do loop_detector

    ! free workspace

    if (allocated(local_hitmap)) deallocate(local_hitmap)
    if (allocated(local_basehitmap)) deallocate(local_basehitmap)
    if (allocated(local_covmat)) deallocate(local_covmat)

    cputime_ncm = cputime_ncm + get_time(5)

  end subroutine write_covmat



  subroutine add_white(idet, detweight)
    integer(i8b), intent(in) :: idet
    real(dp), intent(in) :: detweight

    integer(i4b) :: i
    integer(i8b) :: ip, k

    ! Add white noise to the diagonal

    call reset_time(10)

    do k = 1, noba_short
       loop_sample : do i = baselines_short_start(k), baselines_short_stop(k)
          ip = pixels(i, idet)
          if (ip == dummy_pixel) cycle loop_sample
          !call dsyr('U', nmap, detweight, weights(:, i, idet), 1, &
          !    local_covmat(:, :, ip, ip), nmap)
          call dger(nmap, nmap, detweight, weights(:, i, idet), 1, &
               weights(:, i, idet), 1, local_covmat(:, :, ip, ip), nmap)
          !do imap = 1,nmap
          !   do jmap = 1,nmap
          !      local_covmat(jmap, imap, ip, ip) = &
          !           local_covmat(jmap, imap, ip, ip) &
          !           + weights(imap, i, idet)
          !           * weights(jmap, i, idet) * detweight
          !   end do
          !end do
       end do loop_sample
    end do

    cputime_white = cputime_white + get_time(10)

  end subroutine add_white



  subroutine symmetrize()

    integer(i8b) :: ipix, jpix

    ! symmetrize the matrix since we only accumulated the lower diagonal

    call reset_time(10)

    do ipix = 0, nolocpix-1
       do jpix = ipix+1, nolocpix-1
          local_covmat(:, :, ipix, jpix) = &
               transpose(local_covmat(:, :, jpix, ipix))
       end do
    end do

    cputime_symmetrize = cputime_symmetrize + get_time(10)

  end subroutine symmetrize



  subroutine save_matrix(covmatfile)
    ! create a global pixel-pixel matrix, accumulate the local
    ! contributions and write the matrix out
    character(len=*), intent(in) :: covmatfile
    integer(i4b) :: ierr, mypix1, mypix2, isend, &
         ip1, ip2, evenodd, itarget, isource, nelem
    integer(i8b) :: firstpix, lastpix, npix, npix_proc, pix1, pix2, map1, &
         row, col
    real(dp), allocatable, target :: covmat1(:, :), covmat2(:, :)
    real(dp), pointer :: covmat_send(:, :), covmat_recv(:, :)
    INTEGER(i4b) :: filemode, fileinfo, outfilehandle, count, &
         status(MPI_STATUS_SIZE)
    INTEGER(MPI_OFFSET_KIND) :: fileoffset

    call reset_time(10)

    npix = 12 * nside_map**2
    npix_proc = ceiling(dble(npix) / ntasks)
    nelem = int(npix * nmap * npix_proc * nmap, i4b)

    allocate(covmat1(npix*nmap, npix_proc*nmap), &
         covmat2(npix*nmap, npix_proc*nmap), stat=ierr)
    if (ierr /= 0) then
       print *,'Failed to allocate ', npix*nmap*npix_proc*nmap*2*8/2**20, &
            ' MB for collecting matrix'
       call abort_mpi('No room to collect matrix')
    end if
    if (ierr /= 0) call abort_mpi('No room to collect covmat')

    memory_ncm = memory_ncm + npix*nmap*npix_proc*nmap*2*8

    covmat1 = 0
    covmat2 = 0

    firstpix = id * npix_proc
    ! It is ok to have lastpix > npix. No process will have entries to it
    lastpix = firstpix + npix_proc - 1

    cycle_matrices : DO isend = 0, ntasks - 1
       ! pass the pieces of the covariance matrix in a circle

       ! swap send and receive buffers
       IF (MODULO(isend, 2) == 0) THEN
          covmat_send => covmat1
          covmat_recv => covmat2
       ELSE
          covmat_send => covmat2
          covmat_recv => covmat1
       END IF

       loop_col : do mypix1 = 0, nolocpix-1
          ip1 = mypix1 / nosubpix_max
          pix1 = mypix1 + subtable2(ip1)*nosubpix_max
          if (pix1 < firstpix .or. pix1 > lastpix) cycle loop_col
          col = (pix1-firstpix)*nmap + 1

          loop_row : do mypix2 = 0, nolocpix-1
             ip2 = mypix2 / nosubpix_max
             pix2 = mypix2 + subtable2(ip2)*nosubpix_max
             row = pix2*nmap + 1

             if (local_covmat(1, 1, mypix2, mypix1) < tiny(1._dp)) then
                cycle loop_row
             end if

             covmat_send(row:row+nmap-1, col:col+nmap-1) = &
                  covmat_send(row:row+nmap-1, col:col+nmap-1) &
                  + local_covmat(:, :, mypix2, mypix1)
          end do loop_row
       end do loop_col

       ! cycle the local covariance matrices
       IF (ntasks == 1) THEN
          covmat_recv = covmat_send
       ELSE
          DO evenodd = 0, 1
             IF (MODULO(ID, 2) == evenodd) THEN
                ! send
                itarget = MODULO(ID+1, ntasks)
                CALL MPI_send(covmat_send, nelem, MPI_DOUBLE_PRECISION, &
                     itarget, ID, comm, ierr)
             ELSE
                !receive
                isource = MODULO(ID-1, ntasks)
                CALL MPI_recv(covmat_recv, nelem, MPI_DOUBLE_PRECISION, &
                     isource, isource, comm, status, ierr)
             END IF
          END DO

          firstpix = firstpix - npix_proc
          if (firstpix < 0) firstpix = (ntasks-1) * npix_proc
          lastpix = firstpix + npix_proc - 1
       END IF

    END DO cycle_matrices

    if (lastpix > npix - 1) lastpix = npix - 1

    ! add the sentinel value to unobserved pixels

    do pix1 = firstpix, lastpix
       if (covmat_recv(pix1*nmap+1, (pix1-firstpix)*nmap+1) < tiny(1._dp)) then
          do map1 = 1, nmap
             covmat_recv(pix1*nmap+map1, (pix1-firstpix)*nmap+map1) = -1
          end do
       end if
    end do


    ! Write to disk

    call wait_mpi

    IF (id == 0) WRITE (*,*) 'Writing ' // TRIM(covmatfile)
    CALL mpi_info_create(fileinfo, ierr)

    filemode = IOR(MPI_MODE_WRONLY, MPI_MODE_CREATE)
    CALL mpi_file_open(comm, TRIM(covmatfile), filemode, fileinfo, &
         outfilehandle, ierr)
    IF (ierr /= 0) CALL abort_mpi('Unable to open covmat file')

    if (firstpix < npix) then
       fileoffset = firstpix * nmap * npix * nmap * 8
    else
       fileoffset = 0
    end if

    CALL mpi_file_set_view(outfilehandle, fileoffset, MPI_REAL8, MPI_REAL8, &
         'native', fileinfo, ierr)
    IF (ierr /= 0) CALL abort_mpi('Unable to establish view to covmat file')

    count = int(npix * nmap * (lastpix - firstpix + 1) * nmap, i4b)
    if (count > 0) then
       CALL mpi_file_write(outfilehandle, covmat_recv, count, MPI_REAL8, &
            status, ierr)
       IF (ierr /= 0) CALL abort_mpi('Failed to write covmat file')
    end if

    call wait_mpi

    CALL mpi_file_close(outfilehandle, ierr)
    IF (ierr /= 0) CALL abort_mpi('Failed to close covmat file')


    deallocate(covmat1, covmat2)

    cputime_save_matrix = cputime_save_matrix + get_time(10)

  end subroutine save_matrix



  subroutine get_hitmap(idet, kstart, noba)
    ! construct a hitmap for the current pointing period
    integer(i8b), intent(in) :: idet
    integer(i4b), intent(in) :: kstart, noba
    integer(i4b) :: ierr, k, i, hitmax, nhit
    integer(i8b) :: ip
    integer(i4b), save :: hitmax_max = 0

    call reset_time(10)

    ! measure how much space we need for the basehitmap

    hitmax = 0
    do ip = 0, nolocpix
       nhit = count(abs(local_ptf(1, :, ip)) >= tiny(1._dp))
       hitmax = max(nhit, hitmax)
    end do

    if (allocated(local_basehitmap)) deallocate(local_basehitmap)
    allocate(local_basehitmap(hitmax,0:nolocpix), stat=ierr)
    if (ierr /= 0) then
       print *, 'Failed to allocate ', hitmax*(nolocpix+1)*4/2**20, &
            ' MB for hitmap'
       call abort_mpi('No room for local hitmap')
    end if


    if (hitmax_max < hitmax) then
       memory_ncm = memory_ncm - hitmax_max*(nolocpix+1)*8
       memory_ncm = memory_ncm + hitmax*(nolocpix+1)*8
       hitmax_max = hitmax
    end if

    ! local_basehitmap will have a list of all baselines that hit
    ! a particular pixel

    local_hitmap = 0
    local_basehitmap = 0

    do k = 1, noba
       do i = baselines_short_start(k+kstart), baselines_short_stop(k+kstart)
          ip = pixels(i, idet)
          if (ip == dummy_pixel) cycle

          nhit = local_hitmap(ip)
          if (nhit == 0) then
             local_hitmap(ip) = 1
             local_basehitmap(1, ip) = k
          else
             if (local_basehitmap(nhit, ip) /= k) then
                nhit = nhit + 1
                if (nhit > hitmax) then
                   print *, id, ' : extra hits! ', k, i, nhit, hitmax, &
                        ip, nolocpix
                end if
                local_hitmap(ip) = nhit
                local_basehitmap(nhit, ip) = k
             end if
          end if

       end do
    end do

    cputime_hitmap = cputime_hitmap + get_time(10)

  end subroutine get_hitmap



  subroutine get_ptf(idet, kstart, noba, detweight)
    ! Accumulate P^T F
    integer(i8b), intent(in) :: idet
    integer(i4b), intent(in) :: kstart, noba
    real(dp) :: detweight
    integer(i4b) :: ierr, k, i
    integer(i8b) :: ip
    integer(i4b), save :: nobamax = 0

    call reset_time(10)

    if (allocated(local_ptf)) deallocate(local_ptf)

    allocate(local_ptf(nmap, noba, 0:nolocpix), stat=ierr)
    if (ierr /= 0) then
       print *,'Failed to allocate ', nmap*noba*(nolocpix+1)*8/2**20, &
            ' MB for PTF. noba = ', noba,', nolocpix = ', nolocpix
       call abort_mpi('No room for local P^T F')
    end if

    if (nobamax < noba) then
       memory_ncm = memory_ncm - nmap*(nolocpix+1)*nobamax*8
       memory_ncm = memory_ncm + nmap*(nolocpix+1)*noba*8
       nobamax = noba
    end if

    local_ptf = 0

    do k = 1, noba
       do i = baselines_short_start(k+kstart), baselines_short_stop(k+kstart)
          ip = pixels(i, idet)
          if (ip == dummy_pixel) cycle
          local_ptf(:, k, ip) = local_ptf(:, k, ip) + weights(:, i, idet)
       end do
    end do

    local_ptf = local_ptf * detweight

    cputime_ptf = cputime_ptf + get_time(10)

  end subroutine get_ptf



  subroutine get_middlematrix(fcov, ipsd, detweight, noba)
    ! compute the noba x noba middle matrix block:
    ! M = (C_a^-1 + F^T C_w^-1 F)^-1

    complex(dp), intent(inout) :: fcov(nof/2+1,*)
    integer(i8b), intent(in) :: ipsd
    real(dp), intent(in) :: detweight
    integer(i4b), intent(in) :: noba
    integer(i4b) :: ierr

    integer(i8b), save :: last_ipsd = -1, last_noba = -1, noba_max = 0

    if (ipsd == last_ipsd .and. noba == last_noba) return

    call reset_time(10)

    if (nof < 2 * noba) then
       call abort_mpi( &
            'ERROR: filter is too short to reliably estimate the middle matrix')
    end if

    if (allocated(middlematrix)) deallocate(middlematrix)

    allocate(middlematrix(noba, 1), stat=ierr)
    if (ierr /= 0) then
       print *,'Failed to allocate ', noba*8/2**20, &
            ' MB for middlematrix. noba = ', noba
       call abort_mpi('No room for local middlematrix')
    end if
    if (noba > noba_max) then
       memory_ncm = memory_ncm - noba_max*8
       memory_ncm = memory_ncm + noba*8
       noba_max = noba
    end if

    middlematrix = 0

    ! Adding the white noise term in Fourier domain means adding
    ! a constant offset
    call dfftinv(xx, 1/(fcov(:, ipsd) + nshort * detweight))
    middlematrix(:,1) = xx(1:noba)

    last_ipsd = ipsd
    last_noba = noba

    cputime_middlematrix = cputime_middlematrix + get_time(10)

  end subroutine get_middlematrix


end MODULE covmat
