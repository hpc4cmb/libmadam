! Small utility module for covariance matrices
! 2007-09-07 RK

MODULE covmat_util

  USE planck_config, ONLY : i4b, i8b, sp, dp
  !USE healpix_types
#ifdef HFI
  USE scalapack_tools
#endif

#ifdef HFI
  use piolibkind
  use piolib
#endif

  use mpi_wrappers, only : mpi_wtime

  IMPLICIT NONE

  INTEGER, PARAMETER, private :: covmat_util_unit = 613 ! output unit
  REAL(dp), private :: time_start, time_end, tic_times(1000)
  INTEGER, private :: ID = 0

  INTERFACE eigen_invert_matrix
     MODULE PROCEDURE eigen_invert_symmetric_matrix, &
          eigen_invert_hermitian_matrix
  END INTERFACE



CONTAINS



#ifdef HFI
  subroutine hfi_dbncm2dat(object_in, file_out, dir, nside_out, nstokes_out)
    character(len=*) :: object_in
    character(len=filenamelen), optional :: file_out
    character(len=*), optional :: dir
    integer(i4b), optional :: nside_out, nstokes_out

    integer(i4b) :: nside, nstokes
    character(len=filenamelen) :: file
    character(len=DMCPIOSTRINGMAXLEN) :: command
    integer(pioint) :: pioerr
    integer(piolong) :: index_begin, index_end, read_start, read_end, nread
    integer(i4b) :: bufflen, npix, start_name, rec_len=0, ibuffer, ierr
    real(dp), pointer :: buffer(:)

    integer(MPI_OFFSET_KIND) :: fileoffset
    integer(i4b) :: status(MPI_STATUS_SIZE), filemode, fileinfo, filehandle

    index_begin = piogetbeginobjectidx(object_in)!, group)
    index_end = piogetendobjectidx(object_in)!, group)

    npix = SQRT(DBLE(index_end-index_begin+1))
    if (present(nside_out) .or. present(nstokes_out)) then
       nstokes = 1
       IF (modulo(npix,9) == 0) nstokes = 3
       nside = nint(sqrt(npix/nstokes/12.0))
       if (present(nside_out)) nside_out = nside
       if (present(nstokes_out)) nstokes_out = nstokes
    end if
    bufflen = 10000000

    start_name = index(object_in, ':', .true.) + 1
    file = object_in(start_name:)
    if (present(dir)) file = trim(dir) // trim(file)

    !call mpi_info_create(fileinfo, ierr)
    fileinfo = MPI_INFO_NULL
    filemode = MPI_MODE_WRONLY + MPI_MODE_CREATE
    fileoffset = 0
    call mpi_file_open(MPI_COMM_SELF, TRIM(file), filemode, &
         fileinfo, filehandle, ierr)
    call mpi_file_set_view(filehandle, fileoffset, MPI_REAL8, MPI_REAL8, &
         'native', fileinfo, ierr)

    read_start = index_begin
    ibuffer = 1
    do
       read_end = read_start + bufflen - 1
       if (read_end > index_end) read_end = index_end
       write(command, '("begin=",i0,";end=",i0)') read_start, read_end

       nread = pioreadvectobject(buffer, object_in, command)!, group)

       !call mpi_file_write_at(filehandle, fileoffset, buffer, &
       !     read_end-read_start+1, MPI_REAL8, status, ierr)
       call mpi_file_write(filehandle, buffer, &
            nread, MPI_REAL8, MPI_STATUS_IGNORE, ierr)
       pioerr =  piodeletevecttable(buffer)

       read_start = read_start + bufflen
       if (read_start > index_end) exit

       !fileoffset = fileoffset + bufflen
       ibuffer = ibuffer + 1
    end do

    call mpi_file_close(filehandle, ierr)

    if (present(file_out)) file_out = file

  end subroutine hfi_dbncm2dat



  subroutine hfi_datncm2db(file_in, object_out, npix, npix2_in)
    character(len=*) :: file_in
    character(len=*) :: object_out
    integer(i4b), optional :: npix
    integer(i4b), optional :: npix2_in

    integer (i4b) :: npix2, nread
    character(len=DMCPIOSTRINGMAXLEN) :: groupname, command
    integer(pioint) :: pioerr
    integer(piolong) :: group, write_start, write_end, nwrote
    integer(i4b) :: bufflen, rec_len, ibuffer, ierr
    real(dp), pointer :: buffer(:)

    integer(MPI_OFFSET_KIND) :: fileoffset, filesize
    integer(i4b) :: status(MPI_STATUS_SIZE), filemode, fileinfo, filehandle
    integer(i8b) :: nelem

    npix2 = npix
    if (present(npix2_in)) npix2 = npix2_in

    pioerr = piogetgrpname(groupname, object_out)
    group = pioopenvectgrp(groupname, 'w')
    pioerr = piocreatevectobject(object_out, 'PIODOUBLE', group)

    fileinfo = MPI_INFO_NULL
    filemode = MPI_MODE_RDONLY
    fileoffset = 0
    call mpi_file_open(MPI_COMM_SELF, TRIM(file_in), filemode, &
         fileinfo, filehandle, ierr)
    call mpi_file_set_view(filehandle, fileoffset, MPI_REAL8, MPI_REAL8, &
         'native', fileinfo, ierr)

    bufflen = 10000000
    allocate(buffer(bufflen), stat=ierr)
    if (ierr /= 0) stop 'No room for writebuffer'

    write_start = 0
    do
       call mpi_file_read(filehandle, buffer, bufflen, MPI_REAL8, &
            status, ierr)       
       call mpi_get_count(status, MPI_REAL8, nread, ierr)
       !write (*,'("Read ",i0," elements from ",a)') nread,trim(file_in)
       if (nread < 1) exit

       write_end = write_start + nread - 1
       write(command, '("begin=",i0,";end=",i0)') write_start, write_end

       nwrote = piowritevectobject(buffer, object_out, command, group)
       !write (*,'("Wrote ",i0," elements to ",a)') nwrote,trim(object_out)

       if (nwrote /= nread) stop 'ERROR: failed to write all elements'

       write_start = write_start + nread
    end do

    call mpi_file_close(filehandle, ierr)

    pioerr = pioclosevectgrp(group)

  end subroutine hfi_datncm2db




  subroutine hfi_dbvec2ascii(object_in, file_out, dir)
    character(len=*) :: object_in
    character(len=filenamelen), optional :: file_out
    character(len=*), optional :: dir

    character(len=DMCPIOSTRINGMAXLEN), save :: groupname
    integer(piolong), save :: group

    character(len=filenamelen) :: file
    character(len=DMCPIOSTRINGMAXLEN) :: command
    integer(pioint) :: pioerr
    integer(piolong) :: index_begin, index_end, read_start, read_end, nread
    integer(i4b) :: bufflen, npix, start_name, rec_len, ibuffer, ierr
    real(dp), pointer :: buffer(:)

    if (index(object_in, 'OPEN') > 0 .or. index(object_in, 'open') > 0 ) then
       !pioerr = piogetgrpname(groupname, object_in)
       !group = pioopenvectgrp(groupname, 'r')
       return
    end if

    if (index(object_in, 'CLOSE') > 0 .or. index(object_in, 'close') > 0 ) then
       !pioerr = pioclosevectgrp(group)
       return
    end if

    index_begin = piogetbeginobjectidx(object_in, group)
    index_end = piogetendobjectidx(object_in, group)
    npix = index_end-index_begin+1
    
    allocate(buffer(npix), stat=ierr)
    if (ierr /= 0) stop 'No room for readbuffer'

    start_name = index(object_in, ':', .true.) + 1
    file = object_in(start_name:)
    if (present(dir)) file = trim(dir) // trim(file)
    inquire(iolength=rec_len) buffer

    write(command, '("begin=",i0,";end=",i0)') index_begin, index_end
    nread = pioreadvectobject(buffer, object_in, command, group)

    call save_vector(buffer, trim(file))

    if (present(file_out)) file_out = file

  end subroutine hfi_dbvec2ascii



  subroutine hfi_asciivec2db(file_in, object_out, npix)
    character(len=*) :: file_in
    character(len=*) :: object_out
    integer(i4b) :: npix

    character(len=DMCPIOSTRINGMAXLEN) :: groupname
    integer(piolong) :: group

    character(len=filenamelen) :: file
    character(len=DMCPIOSTRINGMAXLEN) :: command
    integer(pioint) :: pioerr
    integer(piolong) :: nwrite
    integer(i4b) :: ierr
    real(dp), pointer :: buffer(:)

    pioerr = piogetgrpname(groupname, object_out)
    group = pioopenvectgrp(groupname, 'w')
    pioerr = piocreatevectobject(object_out, 'PIODOUBLE', group)

    allocate(buffer(npix), stat=ierr)
    if (ierr /= 0) stop 'No room for writebuffer'

    open(unit=55, file=file_in, form='formatted', status='old')

    read (55, *, iostat=ierr) buffer
    if (ierr /= 0) then
       write (*,*) 'ERROR reading ' // trim(file_in) // ' -- RETURNING'
       return
    end if

    write(command, '("begin=",i0,";end=",i0)') 0, npix-1
    nwrite = piowritevectobject(buffer, object_out, command, group)

    close(55)

    pioerr = pioclosevectgrp(group)

  end subroutine hfi_asciivec2db
#endif


  FUNCTION wtime()
    REAL(dp) :: wtime
    INTEGER :: count, rate
    !call cpu_time(wtime)
    !call system_clock(count, rate)
    !wtime = dble(count)/rate

    wtime =  mpi_wtime()
  END FUNCTION wtime



  SUBROUTINE reduce_matrix(matrix, total_hits)
    ! This routine removes zero lines and columns from the matrix for 
    ! inversion
    REAL(dp), POINTER :: matrix(:,:)
    INTEGER, POINTER  :: total_hits(:) ! hit mask map
    REAL(dp), POINTER :: matrix_temp(:,:)
    INTEGER           :: npix, sigs, pixel, ierr, hits, &
         a, b, row, col, reduced_row, reduced_col

    npix = SIZE(total_hits, 1)
    sigs = SIZE(matrix, 1) / npix
    hits  = COUNT(total_hits >= 0)

    ALLOCATE(matrix_temp(0:sigs*npix-1, 0:sigs*npix-1), stat=ierr)
    IF (ierr/=0) THEN
       WRITE (*,*) 'ERROR: unable to allocate in reduce_matrix'
       RETURN
    END IF
    matrix_temp = matrix
    matrix = 0.0

    reduced_col = 0
    DO col = 0, npix-1
       IF (total_hits(col) == 0) CYCLE
       reduced_row = 0
       DO row = 0, npix-1
         IF (total_hits(row) == 0) CYCLE
         DO a = 0,sigs-1
            DO b = 0,sigs-1
               matrix(reduced_row+a*hits, reduced_col+b*hits) = &
                    matrix_temp(row+a*npix, col+b*npix)
            END DO
         END DO
         reduced_row = reduced_row + 1
       END DO
       reduced_col = reduced_col + 1
    END DO

    DEALLOCATE(matrix_temp)

  END SUBROUTINE reduce_matrix



  SUBROUTINE expand_matrix(matrix, total_hits)
    ! This routine restores original matrix size after call to
    ! reduce_matrix() 
    REAL(dp), POINTER :: matrix(:,:)
    INTEGER, POINTER  :: total_hits(:)
    REAL(dp), POINTER :: matrix_temp(:,:)
    INTEGER           :: npix, sigs, pixel, ierr, hits, &
         a, b, row, col, reduced_row, reduced_col

    npix = SIZE(total_hits, 1)
    sigs = SIZE(matrix, 1) / npix
    hits  = COUNT(total_hits >= 0)

    ALLOCATE(matrix_temp(0:sigs*npix-1, 0:sigs*npix-1), stat=ierr)
    IF (ierr/=0) THEN
       WRITE (*,*) 'ERROR: unable to allocate in reduce_matrix'
       RETURN
    END IF
    matrix_temp = matrix
    matrix = 0.0

    col = -1
    DO reduced_col = 0, hits-1
       DO ! find the real column corresponding to the reduced one
          col = col + 1
          IF (total_hits(col) > 0) EXIT
       END DO
       
       row = -1
       DO reduced_row = 0, hits-1
          DO ! find the real row corresponding to the reduced one
             row = row + 1
             IF (total_hits(row) > 0) EXIT
          END DO
         DO a = 0,sigs-1
            DO b = 0,sigs-1
               matrix(row+a*npix, col+b*npix) = &
                    matrix_temp(reduced_row+a*hits, reduced_col+b*hits)
            END DO
         END DO
       END DO
    END DO

    DEALLOCATE(matrix_temp)

  END SUBROUTINE expand_matrix



  SUBROUTINE convert_to_interlaced(matrix)
    ! This routine takes a block formed polarization map
    ! and returns it in interlaced form where each
    ! pixel is represented by a 3x3 matrix
    REAL(dp), POINTER :: matrix(:,:)
    REAL(dp), POINTER :: matrix_temp(:,:)
    INTEGER           :: npix, ierr, a, b

    IF (LBOUND(matrix,1) /= 0 .OR. MODULO(UBOUND(matrix,1)+1,3) /= 0) THEN
       WRITE (*,*) 'ERROR in convert_to_interlaced, matrix does not '//&
            'include polarization'
       RETURN
    END IF

    npix = (UBOUND(matrix,1)-LBOUND(matrix,1)+1)/3

    ALLOCATE(matrix_temp(0:npix*3-1, 0:npix*3-1), stat=ierr)

    IF (ierr/=0) THEN
       WRITE (*,*) 'ERROR : no room to convert into interlaced'
       RETURN
    END IF

    matrix_temp(:,:) = matrix(:,:)

    DO a = 0,2
       DO b = 0,2
          matrix(a:3*npix-1:3, b:3*npix-1:3) = &
               matrix_temp(a*npix:(a+1)*npix-1,b*npix:(b+1)*npix-1)
       END DO
    END DO

    DEALLOCATE(matrix_temp)
    
  END SUBROUTINE convert_to_interlaced



  SUBROUTINE convert_to_block_form(matrix)
    ! This routine takes an interlaced polarization map
    ! and returns it in block form where correlations are in
    ! II IQ IU
    ! QU QQ QU
    ! UI UQ UU form
    REAL(dp), POINTER :: matrix(:,:)
    REAL(dp), POINTER :: matrix_temp(:,:)
    INTEGER           :: npix, ierr, a, b

    IF (LBOUND(matrix,1) /= 0 .OR. MODULO(UBOUND(matrix,1)+1,3) /= 0) THEN
       WRITE (*,*) 'ERROR in convert_to_block_form, matrix does not contain'//&
            ' polarization part'
       RETURN
    END IF

    npix = (UBOUND(matrix,1)-LBOUND(matrix,1)+1)/3

    ALLOCATE(matrix_temp(0:npix*3-1, 0:npix*3-1), stat=ierr)

    IF (ierr/=0) THEN
       WRITE (*,*) 'ERROR : no room to convert into block form'
       RETURN
    END IF

    matrix_temp(:,:) = matrix(:,:)

    DO a = 0,2
       DO b = 0,2
          matrix(a*npix:(a+1)*npix-1,b*npix:(b+1)*npix-1) = &
          matrix_temp(a:3*npix-1:3, b:3*npix-1:3)
       END DO
    END DO

    DEALLOCATE(matrix_temp)
    
  END SUBROUTINE convert_to_block_form



  SUBROUTINE save_matrix(matrix, rows, cols, filename, unit)
    REAL(dp), POINTER     :: matrix(:,:)
    INTEGER(dp), OPTIONAL :: rows, cols
    CHARACTER(len=*)      :: filename
    INTEGER, OPTIONAL     :: unit
    INTEGER               :: out_unit=covmat_util_unit, i, j, elems

    IF (PRESENT(unit)) out_unit = unit

    elems = SIZE(matrix,1)

    OPEN(unit=out_unit, file=TRIM(filename), status='replace', recl=elems*15)
    DO i = LBOUND(matrix, 1), UBOUND(matrix, 1)
       DO j = LBOUND(matrix, 2), UBOUND(matrix, 2)
          WRITE (out_unit, '(ES15.5)', advance='no') &
               matrix(i, j)
       END DO
       WRITE (out_unit,*)
    END DO
    CLOSE(out_unit)

  END SUBROUTINE save_matrix



  SUBROUTINE save_matrix_bin(matrix, rows, cols, filename, unit)
    REAL(dp), POINTER     :: matrix(:,:)
    INTEGER(dp), OPTIONAL :: rows, cols
    CHARACTER(len=*)      :: filename
    INTEGER, OPTIONAL     :: unit
    INTEGER               :: out_unit=covmat_util_unit, i, j

    IF (PRESENT(unit)) out_unit = unit

    OPEN(unit=out_unit,file=TRIM(filename),status='replace',form='unformatted')
    WRITE(out_unit) matrix
    CLOSE(out_unit)

  END SUBROUTINE save_matrix_bin



  SUBROUTINE save_vector(vector, filename, unit)
    REAL(dp), POINTER     :: vector(:)
    CHARACTER(len=*)      :: filename
    INTEGER, OPTIONAL     :: unit
    INTEGER               :: out_unit=covmat_util_unit, i

    IF (PRESENT(unit)) out_unit = unit

    OPEN(unit=out_unit, file=TRIM(filename), status='replace')
    WRITE (out_unit, '(ES20.10)') vector
    CLOSE(out_unit)

  END SUBROUTINE save_vector



  SUBROUTINE save_ivector(vector, filename, unit)
    INTEGER, POINTER      :: vector(:)
    CHARACTER(len=*)      :: filename
    INTEGER, OPTIONAL     :: unit
    INTEGER               :: out_unit=covmat_util_unit, i

    IF (PRESENT(unit)) out_unit = unit

    OPEN(unit=out_unit, file=TRIM(filename), status='replace')
    WRITE (out_unit, '(I12)') vector
    CLOSE(out_unit)

  END SUBROUTINE save_ivector



  SUBROUTINE eigen_invert_symmetric_matrix(matrix, limit)
    ! computes the eigenvalues of given matrix and inverts
    ! it using the eigenvalues and vectors
    REAL(dp), POINTER :: matrix(:,:)
    REAL, OPTIONAL    :: limit
    REAL              :: rcond_limit = 1E-30
    INTEGER(dp)       :: row, col, workspace_length
    INTEGER           :: ierr, N, good, i
    REAL(dp), POINTER :: workspace(:), eigenvectors(:,:), eigenvectorsT(:,:),&
         eigenvalues(:)
    INTEGER, POINTER  :: pivotings(:)

    IF (PRESENT(limit)) rcond_limit = limit

    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    WRITE (*,'(a,i8)') 'eigen_invert : N == ', N
    ALLOCATE(eigenvalues(N), eigenvectors(N,N), &
         eigenvectorsT(N,N), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') &
            ID, ' : Sorry! no room for eigenvalue inverting matrix'
       RETURN
    END IF

    workspace_length = MAX(5*N, 10000)
    ALLOCATE(workspace(workspace_length), stat=ierr)
    IF (ierr /=0 ) THEN
       WRITE  (*,'(i3,a)') &
            ID, ' : Sorry! no room for workspace in eigenvalue inversion'
       RETURN
    END IF

    ! use dsyev to compute the eigenvalues
    !SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    eigenvectors = matrix
    CALL tic(634)
    CALL DSYEV('V', 'U', N, eigenvectors, N, eigenvalues, &
         workspace, workspace_length, ierr)
    CALL toc('compute eigenvalues', 634)

    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a,i10)') ID, ' : PROBLEM with eigen_invert, info ==',ierr
       RETURN
    END IF

    ! save eigenvalues and eigenvectors
    !CALL save_vector(eigenvalues, 'eigenvalues.dat')
    OPEN(unit=98, file='eigenvalues.dat', status='replace', form='formatted')
    WRITE(98, '(1000000(ES30.20,/))') eigenvalues
    CLOSE(unit=98)
    OPEN(unit=99, file='evecs.dat', status='replace', &
         form='unformatted', access='direct', recl=N**2*8)
    WRITE(99, rec=1) eigenvectors
    CLOSE(unit=99)

    ! construct the inverted matrix    
    matrix = 0
    eigenvectors(:,1) = 0
    eigenvectorsT = TRANSPOSE(eigenvectors)
    DO i = 0, N ! ignore worst eigenvalue
       IF ( ABS(eigenvalues(i)/eigenvalues(N)) > rcond_limit ) &
            eigenvectors(:,i) = eigenvalues(i)**(-1)*eigenvectors(:,i)
    END DO
    matrix = MATMUL(eigenvectors, eigenvectorsT)

  END SUBROUTINE eigen_invert_symmetric_matrix



  SUBROUTINE eigen_invert_hermitian_matrix(matrix, limit)
    ! computes the eigenvalues of given matrix and inverts
    ! it using the eigenvalues and vectors
    COMPLEX(dp), POINTER :: matrix(:,:)
    REAL, OPTIONAL        :: limit
    REAL                  :: rcond_limit = 1E-6
    INTEGER(dp)           :: row, col, workspace_length
    INTEGER               :: ierr, N, good, i
    COMPLEX(dp), POINTER :: workspace(:)
    COMPLEX(dp), POINTER :: eigenvectors(:,:), eigenvectorsT(:,:)
    REAL(dp), POINTER     :: eigenvalues(:), workspace2(:)
    INTEGER, POINTER      :: pivotings(:)

    IF (PRESENT(limit)) rcond_limit = limit

    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    WRITE (*,'(a,i8)') 'eigen_invert : N == ', N
    ALLOCATE(eigenvalues(N), eigenvectors(N,N), &
         eigenvectorsT(N,N), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') &
            ID, ' : Sorry! no room for eigenvalue inverting matrix'
       RETURN
    END IF

    workspace_length = MAX(5*N, 10000)
    ALLOCATE(workspace(workspace_length), workspace2(3*N-2), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE  (*,'(i3,a)') &
            ID, ' : Sorry! no room for workspace in eigenvalue inversion'
       RETURN
    END IF

    ! use zheev to compute the eigenvalues
    !SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO )
    eigenvectors = matrix
    CALL tic(634)
    CALL ZHEEV('V', 'U', N, eigenvectors, N, eigenvalues, &
         workspace, workspace_length, workspace2, ierr)
    CALL toc('compute eigenvalues', 634)

    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a,i10)') ID, ' : PROBLEM with eigen_invert, info ==',ierr
       RETURN
    END IF

    ! save eigenvalues
    CALL save_vector(eigenvalues, 'eigenvalues.dat')

    ! save eigenvectors
    open(unit=145, file='eigenvectors.dat', status='replace', &
         form='unformatted', access='direct', recl=N*8*2)
    do i = 1, N
       write (unit=145, rec=i) eigenvectors(:,i)
    end do

    ! construct the inverted matrix    
    matrix = 0
    eigenvectorsT = CONJG(TRANSPOSE(eigenvectors))
    DO i = 1, N
       IF (ABS(eigenvalues(i)/eigenvalues(N)) > rcond_limit) THEN
          eigenvectors(:,i) = eigenvalues(i)**(-1)*eigenvectors(:,i)
       ELSE
          eigenvectors( :,i) = 0
          eigenvectorsT(i,:) = 0
          WRITE (*,'(a,i6)') ' Skipping eigenmode #', i
       END IF
    END DO
    matrix = MATMUL(eigenvectors, eigenvectorsT)

  END SUBROUTINE eigen_invert_hermitian_matrix




  SUBROUTINE svd_invert_matrix(matrix, limit)
    ! Computes the pseudoinverse of a matrix
    REAL(dp), POINTER :: matrix(:,:)
    REAL, OPTIONAL    :: limit
    REAL              :: rcond_limit = 1E-3
    INTEGER(dp)       :: row, col, workspace_length
    INTEGER           :: ierr, N, good
    REAL(dp), POINTER :: workspace(:), U(:,:), SIGMA(:), VT(:,:)
    INTEGER, POINTER  :: pivotings(:)

    IF (PRESENT(limit)) rcond_limit = limit

    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    ALLOCATE(SIGMA(N), U(N,N), VT(N,N), stat=ierr)
    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') &
            ID, ' : Sorry! no room for SVD inverting matrix'
       RETURN
    END IF

    workspace_length = MAX(10*N, 10000)
    ALLOCATE(workspace(workspace_length), stat=ierr)
    IF (ierr /=0 ) THEN
       WRITE  (*,'(i3,a)') &
            ID, ' : Sorry! no room for workspace in svd inversion'
       RETURN
    END IF

    ! decompose:
    ! This routine returns the singular value decomposition
    ! of matrix A : A = U * SIGMA * transpose(V), where
    ! U and V are orthogonal matrices and SIGMA is a diagonal matrix
    ! with the singular values in decreasing order on the diagonal
    !
    !SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
    !                   WORK, LWORK, INFO )
    CALL DGESVD('A', 'A', N, N, matrix, N, SIGMA, U, N, VT, N, &
         workspace, workspace_length, ierr)

    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') ID, ' : PROBLEM with SVD'
       RETURN
    END IF

    CALL save_vector(SIGMA, 'singular_values.dat')

    ! Check how many good singular values we have. Omit bad ones from
    ! the pseudoinverse
    good = 2
    DO 
       IF (good > N) EXIT
       IF (SIGMA(good)/SIGMA(1) < rcond_limit) EXIT
       good = good + 1
    END DO
    good = good - 1

    !good = 189 ! remove after testing

    WRITE (*,'(i3,a,i6,a,i6)') ID, ' : Number of good singular values : ', &
         good, '/', N

    DO col = 1,good
       U(:, col) = U(:, col)/SIGMA(col)
    END DO

    ! The pseudoinverse is now:
    ! A^-1 = V * SIGMA^inv * transpose(U), where
    ! SIGMA^inv is the original diagonal SIGMA-matrix, with the nonzero
    ! singular values inverted

    matrix = TRANSPOSE(MATMUL(U(:,1:good),VT(1:good,:)))

    DEALLOCATE(SIGMA, U, VT, workspace)

  END SUBROUTINE svd_invert_matrix
    


  SUBROUTINE invert_matrix(matrix)
    REAL(dp), POINTER :: matrix(:,:)
    INTEGER(dp)       :: row, col, workspace_length, N
    INTEGER           :: ierr
    REAL(dp), POINTER :: workspace(:)
    INTEGER, POINTER  :: pivotings(:)
    
    ! factorize, dsytrf assumes covmat_inv to be in upper triangular form
    ! SUBROUTINE DSYTRF( UPLO, N, A, LDA, IPIV, WORK, LWORK, INFO )
    ! LDA=leading dimension, tells the routine how many elements there
    ! are between two consecutive columns on the same row
    ierr = 0
    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    workspace_length = N**2 * 5
    ALLOCATE(workspace(workspace_length), pivotings(N), stat=ierr)
    IF (ierr/=0) THEN
       WRITE (*,'(i3,a)') &
            ID, ' : Sorry, no room for inverting matrix'
    END IF

    ! Factorize
    CALL DSYTRF('upper', N, matrix, N, pivotings, workspace, &
         workspace_length, ierr)
    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') &
            ID,' : PROBLEM with factorizing the covariance matrix'
    END IF

    ! Invert
    CALL DSYTRI('upper', N, matrix, N, pivotings, workspace, &
         workspace_length, ierr)
    IF (ierr /= 0) THEN
       WRITE (*,'(i3,a)') &
            ID, ' : PROBLEM with inverting the covariance matrix'
    END IF

    ! Finally just copy the upper triangle to lower half as well
    DO col = 1, N-1
       DO row = 0, col-1
          matrix(col, row) = matrix(row, col)
       END DO
    END DO

    DEALLOCATE(workspace, pivotings)

  END SUBROUTINE invert_matrix



  SUBROUTINE tic(index)
    INTEGER, OPTIONAL :: index
    
    IF (PRESENT(index)) THEN
       tic_times(index) = wtime()
    ELSE
       time_start = wtime()
    END IF

  END SUBROUTINE tic



  SUBROUTINE toc(label, index, threshold)
    CHARACTER(len=*), OPTIONAL :: label
    INTEGER, OPTIONAL :: index
    REAL(dp), OPTIONAL :: threshold
    CHARACTER(len=512) :: middle_string
    REAL(dp) :: elapsed, time_now
    REAL(dp), SAVE :: reporting_threshold=1e-3

    IF (PRESENT(threshold)) reporting_threshold = threshold

    time_now = wtime()

    IF (PRESENT(index)) THEN
       elapsed = time_now - tic_times(index)
    ELSE
       elapsed = time_now - time_start
    END IF

    IF (elapsed < reporting_threshold) RETURN

    IF (PRESENT(label)) THEN
       middle_string = ' : '//TRIM(label)//' completed in '
    ELSE
       middle_string = ' : elapsed time : '
    END IF

    WRITE (*,'(i3,a,f8.3,a)') &
         ID, TRIM(middle_string), elapsed, ' s'

  END SUBROUTINE toc


END MODULE covmat_util
