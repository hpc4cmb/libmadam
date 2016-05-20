program test_covmat
  
  use commonparam
  use inputparam
  use simulation
  use parameter_control
  use pointing
  use maptod_transfer
  use submap_transfer
  use compression
  use fourier
  use noise_routines
  use madam_routines
  use map_routines
  use read_data
  use output
  use mpi_wrappers
  use maps_and_baselines
  use tod_storage
  use memory_and_time
  use timing
  use covmat
  use covmat_util, only : tic, toc

  implicit none

  integer :: fftlen, n=256, blen=4, i, j, lim, offset, ierr, npix, rec_len, ibase
  real(dp), allocatable, dimension(:) :: invcov, mmvec
  real(dp), allocatable, dimension(:,:) :: p, f, cn, ca, ptemp, ftemp, ctemp, mm
  real(dp), allocatable, dimension(:,:) :: cf, ptf, ptp, ppcov, covmat_in, diff, ptfa
  real(dp) :: mindiff, maxdiff
  character(len=1024) :: file_in

  call init_mpi(ntasks, id)

  print *,' ID == ', id

  allocate(detectors(1))
  detectors(1)%name = 'testdet1'
  detectors(1)%sigma = 1.0

  nosamples_proc = n
  nosamples_tot = ntasks * n
  nshort = blen
  noba_short_tot = nosamples_tot / blen
  noba_short = n / blen
  kshort_start = id * noba_short
  nostokes = 1
  nside_map = 8
  npix = 12 * nside_map**2
  nthreads = 1
  id_thread = 0
  fsample = 1.0
  allocate(pixels(n,1), qw(n,1), uw(n,1))
  pixels(:,1) = (/ (modulo(n*id+i-1, npix), i=1,n) /)
  if (id == 1) pixels = -1
  if (nostokes == 3) then
     qw = 1.0
     uw = -1.0
  end if

  if (id == 0) then
     print *,' Number of samples ', nosamples_tot, '(', nosamples_proc, ')'
     print *,' Number of baselines ', noba_short_tot, '(', noba_short, ')'
  end if

  ! construct the baseline covariance

  fftlen = noba_short_tot
  allocate(fcov(fftlen/2+1,1), invcov(fftlen))
  fcov = 0
  fcov(1,1) = 1
  fcov(2,1) = .5
  fcov(fftlen/2+1,1) = .5
  call init_fourier(fftlen)
  call dfftinv(invcov, fcov)

  call write_covmat(.false., .false., 'testcovmat')

  ! produce the same matrix using dense matrix operations
  allocate(p(nosamples_tot, 0:npix*nostokes-1), f(nosamples_tot, noba_short_tot), &
       cn(nosamples_tot, nosamples_tot), ca(noba_short_tot, noba_short_tot), &
       ptemp(nosamples_tot, 0:npix*nostokes-1), ftemp(nosamples_tot, noba_short_tot), &
       ctemp(nosamples_tot, nosamples_tot))
  p = 0
  f = 0
  cn = 0
  offset = kshort_start*blen
  do i = 1,nosamples_proc
     if (pixels(i,1) >= 0) then
        ! pointing in the covmat module is interlaced
        p(i+offset, pixels(i,1)*nostokes + 0) = 1
        if (nostokes == 3) then
           p(i+offset, pixels(i,1)*nostokes + 1) = qw(i,1)
           p(i+offset, pixels(i,1)*nostokes + 2) = uw(i,1)
        end if
     end if

     ibase = (i-1) / nshort
     f(i+offset, 1+ibase+kshort_start) = 1

     cn(i+offset,i+offset) = 1 / detectors(1)%sigma**2
  end do

  call mpi_allreduce(p, ptemp, nosamples_tot*npix*nostokes, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  p = ptemp

  call mpi_allreduce(f, ftemp, nosamples_tot*noba_short_tot, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  f = ftemp

  call mpi_allreduce(cn, ctemp, nosamples_tot**2, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_comm_world, ierr)
  cn = ctemp

  ca = 0
  do i = 1,noba_short_tot
     ca(i, :) = cshift(invcov, 1-i)
  end do

  if (id == 0) then
     !print *,'P'
     !do i = 1,nosamples_tot
     !   write (*,'(100(f5.2))') p(i,:)
     !end do

     !print *,'F'
     !do i = 1,nosamples_tot
     !   write (*,'(100(f5.2))') f(i,:)
     !end do

     !print *,'Ca'
     !do i = 1,noba_short_tot
     !   write (*,'(100(f5.2))') ca(i,:)
     !end do

     !print *,'Cn'
     !do i = 1,nosamples_tot
     !   write (*,'(100(f5.2))') cn(i,:)
     !end do

     allocate(mm(noba_short_tot, noba_short_tot), ptf(npix*nostokes, noba_short_tot), &
          ptp(npix*nostokes, npix*nostokes), ppcov(npix*nostokes, npix*nostokes), &
          covmat_in(npix*nostokes, npix*nostokes), diff(npix*nostokes, npix*nostokes), &
          ptfa(npix*nostokes, noba_short_tot), stat=ierr)
     if (ierr /= 0) stop 'No room for arrays'

     ! perform the multiplications for the matrix

     mm = ca + matmul(matmul(transpose(f), cn), f)
     call eigen_invert_symmetric_matrix(mm)
     print *,'First row of middle matrix'
     write (*,'(100(f15.5))') mm(1,:)
     do lim = noba_short_tot/2, 1, -1
        if (mm(lim,1) / mm(1,1) > zero_limit) exit
     end do
     if (lim > noba_short) lim = noba_short
     print *,'Truncating middle matrix at ', lim
     do i = 1,noba_short_tot
        do j = 1,noba_short_tot
           if (abs(i-j) >= lim) mm(i,j) = 0
        end do
     end do

     !print *,'middle matrix'
     !do i = 1,noba_short_tot
     !   write (*,'(100(f10.4))') mm(i,:)
     !end do

     ptf = matmul(transpose(p), matmul(cn, f))
     ptp = matmul(transpose(p), matmul(cn, p))
     ptfa = matmul(ptf, mm)
     ppcov = ptp - matmul(ptfa, transpose(ptf))

     ! For true comparison, empty elements along the diagonal need to be substituted by -1
     do i = 1,npix
        if (ppcov(i,i) == 0) ppcov(i,i) = -1
     end do

     !print *,'PTFA: '
     !do i = 1,npix*nostokes
     !   write (*,'(100(f10.4))') ptfa(i,:)
     !end do     

     !print *,'Right term: '
     !do i = 1,npix*nostokes
     !   write (*,'(100(f15.5))') ptp(i,:) - ppcov(i,:)
     !end do     

     ! read in the one written from write_covmat

     file_in = 'testcovmat_inv_bin.dat'
     inquire(iolength=rec_len) covmat_in
     OPEN(unit=55, file=file_in, status='old', form='unformatted', access='direct', recl=rec_len)          
     READ (55, rec=1) covmat_in
     CLOSE(unit=55)

     ! compare
     print *,'        dense,          sparse'
     print *,'Minima : ', minval(ppcov), minval(covmat_in)
     print *,'Maxima : ', maxval(ppcov), maxval(covmat_in)
     print *,'Diff : ', minval(covmat_in-ppcov), maxval(covmat_in-ppcov)

     if (npix*nostokes < 37) then
        print *,'Dense results'
        do i = 1,npix*nostokes
           write (*,'(100(f10.4))') ppcov(:,i)
        end do

        print *,'Sparse results'
        do i = 1,npix*nostokes
           write (*,'(100(f10.4))') covmat_in(:,i)
        end do

        print *,'Sparse - dense'
        diff = covmat_in - ppcov
        where(abs(diff) < 1e-10) diff = 0
        do i = 1,npix*nostokes
           write (*,'(100(f10.4))') diff(:,i)
        end do
     end if

     mindiff = minval(covmat_in - ppcov)
     maxdiff = maxval(covmat_in - ppcov)

     if (mindiff > -1e-10 .and. maxdiff < 1e-10) then
        print *,'Test PASSED'
     else
        print *,'Test FAILED'
     end if

     deallocate(mm, ptf, ptp, ppcov, covmat_in, diff, ptfa)
  end if

  deallocate(detectors, pixels, qw, uw, fcov, invcov, p, f, cn, ca, ptemp, ftemp, ctemp)
  
  call close_mpi()



CONTAINS



  SUBROUTINE eigen_invert_symmetric_matrix(matrix, limit)
    ! computes the eigenvalues of given matrix and inverts
    ! it using the eigenvalues and vectors
    REAL(dp), ALLOCATABLE :: matrix(:,:)
    REAL, OPTIONAL    :: limit
    REAL              :: rcond_limit = 1E-30
    INTEGER(dp)       :: row, col, workspace_length
    INTEGER           :: ierr, N, good, i
    REAL(dp), POINTER :: workspace(:), eigenvectors(:,:), eigenvectorsT(:,:),&
         eigenvalues(:)
    INTEGER, POINTER  :: pivotings(:)
    CHARACTER(len=*), parameter :: covmat_file='middlematrix'

    IF (PRESENT(limit)) rcond_limit = limit

    N = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    WRITE (*,'(a,i8)') 'eigen_invert : N == ', N
    ALLOCATE(eigenvalues(N), eigenvectors(N,N), &
         eigenvectorsT(N,N), stat=ierr)
    IF (ierr /= 0) THEN
       stop ' Sorry! no room for eigenvalue inverting matrix'
    END IF

    workspace_length = MAX(5*N, 10000)
    ALLOCATE(workspace(workspace_length), stat=ierr)
    IF (ierr /=0 ) THEN
       stop ' Sorry! no room for workspace in eigenvalue inversion'
    END IF

    ! use dsyev to compute the eigenvalues
    !SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
    eigenvectors = matrix
    CALL tic(634)
    CALL DSYEV('V', 'U', N, eigenvectors, N, eigenvalues, &
         workspace, workspace_length, ierr)
    CALL toc('compute eigenvalues', 634)

    IF (ierr /= 0) THEN
       WRITE (*,'(a,i10)') ' PROBLEM with eigen_invert, info ==',ierr
       stop
    END IF

    ! save eigenvalues and eigenvectors
    OPEN(unit=98, file=trim(covmat_file)//'.eigen', status='replace', &
         form='formatted')
    WRITE(98, '(1000000(ES30.20,/))') eigenvalues
    CLOSE(unit=98)

    inquire(iolength=rec_len) eigenvectors
    OPEN(unit=99, file=trim(covmat_file)//'-evecs', status='replace', &
         form='unformatted', access='direct', recl=rec_len)
    WRITE(99, rec=1) eigenvectors
    CLOSE(unit=99)

    ! construct the inverted matrix    
    matrix = 0
    eigenvectorsT = TRANSPOSE(eigenvectors)
    DO i = 1, N
       IF ( eigenvalues(i)/eigenvalues(N) > rcond_limit ) &
            eigenvectors(:,i) = eigenvalues(i)**(-1)*eigenvectors(:,i)
    END DO
    matrix = MATMUL(eigenvectors, eigenvectorsT)

    deallocate(eigenvalues, eigenvectors, eigenvectorsT, workspace)

  END SUBROUTINE eigen_invert_symmetric_matrix



end program test_covmat
