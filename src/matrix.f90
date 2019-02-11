MODULE matrix
  !
  ! Routines for handling of symmetric matrices
  !

  use commonparam, only : incomplete_matrices, allow_decoupling

  implicit none
  private

  public invert_eig, invert_eig_s, invert_lu, test_decoupling

CONTAINS

  subroutine test_decoupling(cc, n, istart)
    !
    ! This routine checks if the matrix has a perfectly decoupled
    ! intensity block.  If one is found, it is inverted and
    ! istart = 2 is returned.  Otherwise cc is unchanged and istart = 1
    !
    integer, intent(in) :: n
    double precision, intent(inout) :: cc(n, n)
    integer, intent(out) :: istart
    integer :: i

    if (.not. allow_decoupling) then
       istart = 1
    else
       istart = 2
       do i = 1, n
          if (cc(i, i) > 1e-30) then
             if (any(abs(cc(1, 2:)) / cc(i, i) > 1e-10)) then
                istart = 1
             end if
          end if
       end do
       if (any(abs(cc(2:, 1)) > 1e-10)) then
          istart = 1
       else
          ! Decoupled case
          istart = 2
       end if

       if (istart == 2 .and. cc(1, 1) > 1e-30) then
          cc(1, 1) = 1 / cc(1, 1)
          cc(2:, 1) = 0
          cc(1, 2:) = 0
       end if
    end if

  end subroutine test_decoupling

  SUBROUTINE invert_eig(cc, cdet, rcond, sing, n)
    !
    ! Invert a symmetric positive-definite matrix through eigenvalue analysis.
    ! Also returns CDET=determinant and RCOND=reciprocal condition number.
    ! Dlag SING is true if matrix is singular.
    !
    integer, intent(in) :: n
    double precision, intent(inout) :: cc(n, n)
    double precision, intent(out) :: cdet
    double precision, intent(out) :: rcond
    logical, intent(out) :: sing
    double precision, allocatable :: eigenvalues(:), eigenvectors(:, :), &
         eigenvectorsT(:, :), workspace(:)
    logical :: found(n)
    integer :: info, i, j, nfound, ierr
    integer, allocatable :: good(:)

    if (incomplete_matrices) then
       found = .false.
       do i = 1, n
          if (abs(cc(i, i)) > 1e-30) then
             found(i) = .true.
          end if
       end do
       nfound = count(found)
       if (nfound == 0) then ! DEBUG
          print *,'WARNING: nfound == 0: cc = ', cc ! DEBUG
          nfound = 1 ! DEBUG
       end if ! DEBUG
    else
       found = .true.
       nfound = n
    end if

    allocate(eigenvalues(nfound), eigenvectors(nfound, nfound), &
         eigenvectorsT(nfound, nfound), &
         workspace(10 * nfound), good(nfound), stat=ierr)
    if (ierr /= 0) stop 'No room to invert matrices'

    j = 1
    do i = 1, n
       if (found(i)) then
          good(j) = i
          j = j + 1
       end if
    end do

    eigenvectors = cc(good, good)

    call dsyev('V', 'U', nfound, eigenvectors, nfound, eigenvalues, workspace, &
         10 * nfound, info)

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
       do i = 1, nfound
          eigenvectors(:, i) = eigenvectors(:, i) / eigenvalues(i)
       end do
       cc(good, good) = matmul(eigenvectors, eigenvectorsT)
    end if

    deallocate(eigenvalues, eigenvectors, eigenvectorsT, workspace, good)

    return

  END SUBROUTINE invert_eig


  !----------------------------------------------------------------------------


  SUBROUTINE invert_eig_s(cc, cdet, rcond, noeig, plim, n)
    !
    ! Invert a symmetric positive-definite 3x3 matrix through eigenvalue analysis.
    ! Special treatment for nearly singular matrices:
    ! Only the part corresponding to eigenvalues below PLIM included in solution.
    !
    integer :: n
    double precision, intent(inout) :: cc(n, n)
    double precision, intent(out) :: cdet
    double precision, intent(out) :: rcond
    integer, intent(out) :: noeig ! number of accepted eigenvalues
    double precision, intent(in) :: plim  ! lowest accepted eigenvalue
    double precision :: eigenvalues(n), eigenvectors(n, n)
    double precision :: eigenvectorsT(n, n), workspace(10 * n)
    integer :: info, i

    eigenvectors = cc
    call dsyev('V', 'U', n, eigenvectors, n, eigenvalues, workspace, 10 * n, info)

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

    where(eigenvalues < abs(plim * maxval(eigenvalues)))
       eigenvalues = 0
    elsewhere
       eigenvalues = 1 / eigenvalues
    end where

    noeig = count(eigenvalues /= 0)

    eigenvectorsT = transpose(eigenvectors)

    do i = 1, n
       eigenvectors(:, i) = eigenvectors(:, i) * eigenvalues(i)
    end do

    cc = matmul(eigenvectors, eigenvectorsT)

  END SUBROUTINE invert_eig_s


  !-----------------------------------------------------------------------


  SUBROUTINE Jacobi(a, n, d, v, nrot)
    !
    ! Eigenvalues d and eigenvectors v of a matrix a.
    ! Jacobi algorithm from Numerical Recipes
    !
    integer, intent(in) :: n
    integer, intent(out) :: nrot
    double precision, intent(inout) :: a(n, n)
    double precision, intent(out) :: d(n), v(n, n)
    integer :: i, j, ip, iq
    double precision :: c, g, h, s, sm, t, tau, theta, tresh
    double precision :: b(n), z(n)

    v = 0
    do ip = 1, n
       v(ip, ip) = 1
       d(ip) = a(ip, ip)
    end do
    b = d
    z = 0

    nrot = 0
    !do i = 1,50
    do i = 1, 1000 ! -RK

       sm = 0
       do ip = 1, n - 1
          do iq = ip + 1, n
             sm = sm + abs(a(ip, iq))
          end do
       end do
       if (sm == 0) return
       if (i < 4) then
          tresh = 0.2 * sm / (n * n)
       else
          tresh = 0
       endif

       do ip = 1, n - 1
          do iq = ip + 1, n
             g = 100 * abs(a(ip, iq))
             if ((i > 4) .and. (abs(d(ip)) + g == abs(d(ip))) .and. &
                  (abs(d(iq)) + g == abs(d(iq)))) then
                a(ip, iq) = 0
             else if (abs(a(ip, iq)) > tresh) then
                h = d(iq) - d(ip)
                if (abs(h) + g == abs(h)) then
                   t = a(ip, iq) / h
                else
                   theta = 0.5d0 * h / a(ip, iq)
                   t = 1 / (abs(theta) + sqrt(1 + theta * theta))
                   if (theta < 0) t = -t
                end if
                c = 1 / sqrt(1 + t * t)
                s = t * c
                tau = s / (1 + c)
                h = t * a(ip, iq)
                z(ip) = z(ip) - h
                z(iq) = z(iq) + h
                d(ip) = d(ip) - h
                d(iq) = d(iq) + h
                a(ip, iq) = 0
                do j = 1, ip - 1
                   g = a(j, ip)
                   h = a(j, iq)
                   a(j, ip) = g - s*(h + g * tau)
                   a(j, iq) = h + s*(g - h * tau)
                end do
                do j = ip + 1, iq - 1
                   g = a(ip, j)
                   h = a(j, iq)
                   a(ip, j) = g - s * (h + g * tau)
                   a(j, iq) = h + s * (g - h * tau)
                end do
                do j = iq + 1, n
                   g = a(ip, j)
                   h = a(iq, j)
                   a(ip, j) = g - s * (h + g * tau)
                   a(iq, j) = h + s * (g - h * tau)
                end do
                do j = 1, n
                   g = v(j, ip)
                   h = v(j, iq)
                   v(j, ip) = g - s * (h + g * tau)
                   v(j, iq) = h + s * (g - h * tau)
                end do
                nrot = nrot + 1
             end if
          end do
       end do
       b = b + z
       d = b
       z = 0
    end do
    write(*,*) 'Too many iterations'
    stop

  END SUBROUTINE Jacobi

  !------------------------------------------------------------------------


  SUBROUTINE invert_LU(cc, cdet, sing, n)
    !
    ! Invert a symmetric matrix by LU decomposition
    !
    integer :: n
    double precision, intent(inout) :: cc(n, n)
    double precision, intent(out) :: cdet
    logical, intent(out) :: sing
    integer :: info, i, j

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
       do j = i + 1, n
          cc(j, i) = cc(i, j)
       end do
    end do

    sing = .false.

  END SUBROUTINE invert_LU

  !-----------------------------------------------------------------------

END MODULE matrix
