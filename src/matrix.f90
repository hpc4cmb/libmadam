MODULE matrix
  !
  ! Routines for handling of symmetric matrices
  !

  use commonparam, only : incomplete_matrices, i4b, i8b, sp, dp

  implicit none
  private

  public invert_eig, invert_eig_s, invert_lu

CONTAINS

  SUBROUTINE invert_eig(cc, cdet, rcond, sing, n)
    !
    ! Invert a symmetric positive-definite matrix through eigenvalue analysis.
    ! Also returns CDET=determinant and RCOND=reciprocal condition number.
    ! Dlag SING is true if matrix is singular.
    !
    integer(i8b) :: n
    double precision, intent(inout) :: cc(n,n)
    double precision, intent(out) :: cdet
    double precision, intent(out) :: rcond
    logical, intent(out) :: sing
    double precision, allocatable :: eigenvalues(:), eigenvectors(:, :), &
         eigenvectorsT(:, :), workspace(:)
    logical :: found(n)
    integer :: info, ierr
    integer(i8b) :: i, j, nfound
    integer(i8b), allocatable :: good(:)

    if (incomplete_matrices) then
       found = .false.
       do i = 1, n
          if (abs(cc(i, i)) > 1e-30) then
             found(i) = .true.
          end if
       end do
       nfound = count(found)
    else
       found = .true.
       nfound = n
    end if

    allocate(eigenvalues(nfound), eigenvectors(nfound,nfound), &
         eigenvectorsT(nfound,nfound), &
         workspace(10*nfound), good(nfound), stat=ierr)
    if (ierr /= 0) stop 'No room to invert matrices'

    j = 1
    do i = 1, n
       if (found(i)) then
          good(j) = i
          j = j + 1
       end if
    end do

    eigenvectors = cc(good,good)

    call dsyev(&
         'V', 'U', nfound, eigenvectors, nfound, eigenvalues, workspace, &
         10*nfound, info)

    cdet = product(eigenvalues)
    rcond = minval(eigenvalues) / maxval(eigenvalues)

    if (info /= 0 .or. maxval(eigenvalues) < 1e-30 .or. rcond < 1e-12) then
       sing = .true.
       cc = 0
       cdet = 0
       rcond = 0
    else
       sing = .false.
       eigenvectorsT = transpose(eigenvectors)
       do i = 1,nfound
          eigenvectors(:,i) = eigenvectors(:,i) / eigenvalues(i)
       end do
       cc(good,good) = matmul(eigenvectors, eigenvectorsT)
    end if

    deallocate(eigenvalues, eigenvectors, eigenvectorsT, workspace, good)

    return

  END SUBROUTINE invert_eig


  !----------------------------------------------------------------------------


  SUBROUTINE invert_eig_s(cc, cdet, rcond, noeig, plim, n)
    !
    ! Invert a symmetric positive-definite 3x3 matrix through
    ! eigenvalue analysis.
    ! Special treatment for nearly singular matrices:
    ! Only the part corresponding to eigenvalues below PLIM included
    ! in solution.
    !
    integer(i8b) :: n
    double precision,intent(inout) :: cc(n, n)
    double precision,intent(out) :: cdet
    double precision,intent(out) :: rcond
    integer, intent(out) :: noeig ! number of accepted eigenvalues
    double precision,intent(in) :: plim  ! lowest accepted eigenvalue
    double precision :: eigenvalues(n), eigenvectors(n, n), &
         eigenvectorsT(n, n), workspace(10*n)
    integer :: info
    integer(i8b) :: i

    eigenvectors = cc
    call dsyev('V', 'U', n, eigenvectors, n, eigenvalues, workspace, 10*n, info)

    cdet = product(eigenvalues)
    rcond = minval(eigenvalues) / maxval(eigenvalues)

    !if (info /= 0 .or. maxval(eigenvalues) < 1e-30 .or. rcond < 1e-12) then
    if (info /= 0 .or. maxval(eigenvalues) < 1e-30) then
       cc = 0
       cdet = 0
       rcond = 0
       noeig = 0
       return
    end if

    noeig = count(eigenvalues >= abs(plim*maxval(eigenvalues)))

    where(eigenvalues < abs(plim*maxval(eigenvalues)))
       eigenvalues = 0
    elsewhere
       eigenvalues = 1 / eigenvalues
    end where

    eigenvectorsT = transpose(eigenvectors)

    do i = 1, n
       eigenvectors(:, i) = eigenvectors(:, i) * eigenvalues(i)
    end do

    cc = matmul(eigenvectors, eigenvectorsT)

  END SUBROUTINE invert_eig_s


  !-----------------------------------------------------------------------


  SUBROUTINE invert_LU(cc, cdet, sing, n)
    !
    ! Invert a symmetric matrix by LU decomposition
    !
    integer(i8b) :: n
    double precision, intent(inout) :: cc(n, n)
    double precision, intent(out) :: cdet
    logical, intent(out) :: sing
    integer :: info
    integer(i8b) :: i, j

    call dpotrf('U', n, cc, n, info)

    if (info /= 0) then
       sing = .true.
       return
    end if

    cdet = 1
    do i = 1, n
       cdet = cdet * cc(i, i)
    end do

    call dpotri('U', n, cc, n, info)

    if (info /= 0) then
       cc = 0
       sing = .true.
       return
    end if

    ! copy to lower triangle

    do i = 1, n
       do j = i+1, n
          cc(j, i) = cc(i, j)
       end do
    end do

    sing = .false.

  END SUBROUTINE invert_LU

  !-----------------------------------------------------------------------

END MODULE matrix
