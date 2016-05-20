MODULE matrix
  !
  ! Routines for handling of symmetric matrices
  !
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
    integer :: n
    double precision,intent(inout) :: cc(n,n)
    double precision,intent(out) :: cdet
    double precision,intent(out) :: rcond
    logical,         intent(out) :: sing
    double precision :: eigenvalues(n), eigenvectors(n,n), eigenvectorsT(n,n), workspace(10*n)
    integer :: info, i

    eigenvectors = cc
    call dsyev( 'V', 'U', n, eigenvectors, n, eigenvalues, workspace, 10*n, info )

    cdet = product( eigenvalues )
    rcond = minval( eigenvalues ) / maxval( eigenvalues )

    if (info /= 0 .or. maxval(eigenvalues) < 1e-30 .or. rcond < 1e-12) then
       sing = .true.
       cc = 0
       cdet = 0
       rcond = 0
       return
    else
       sing = .false.
    end if

    eigenvectorsT = transpose(eigenvectors)
    do i = 1,n
       eigenvectors(:,i) = eigenvectors(:,i) / eigenvalues(i)
    end do

    cc = matmul( eigenvectors, eigenvectorsT )

  END SUBROUTINE invert_eig


  !----------------------------------------------------------------------------


  SUBROUTINE invert_eig_s(cc, cdet, rcond, noeig, plim, n)
    !
    ! Invert a symmetric positive-definite 3x3 matrix through eigenvalue analysis.
    ! Special treatment for nearly singular matrices:
    ! Only the part corresponding to eigenvalues below PLIM included in solution.
    !
    integer :: n
    double precision,intent(inout) :: cc(n,n)
    double precision,intent(out)   :: cdet
    double precision,intent(out)   :: rcond
    integer,         intent(out)   :: noeig ! number of accepted eigenvalues
    double precision,intent(in)    :: plim  ! lowest accepted eigenvalue
    double precision :: eigenvalues(n), eigenvectors(n,n), eigenvectorsT(n,n), workspace(10*n)
    integer :: info, i

    eigenvectors = cc
    call dsyev( 'V', 'U', n, eigenvectors, n, eigenvalues, workspace, 10*n, info )

    cdet = product( eigenvalues )
    rcond = minval( eigenvalues ) / maxval( eigenvalues )

    !if (info /= 0 .or. maxval(eigenvalues) < 1e-30 .or. rcond < 1e-12) then
    if (info /= 0 .or. maxval(eigenvalues) < 1e-30) then
       cc = 0
       cdet = 0
       rcond = 0
       noeig = 0
       return
    end if

    where( eigenvalues < abs(plim*maxval(eigenvalues)) )
       eigenvalues = 0
    elsewhere
       eigenvalues = 1 / eigenvalues 
    end where

    noeig = count( eigenvalues /= 0 )

    eigenvectorsT = transpose(eigenvectors)

    do i = 1,n       
       eigenvectors(:,i) = eigenvectors(:,i) * eigenvalues(i)
    end do

    cc = matmul( eigenvectors, eigenvectorsT )

  END SUBROUTINE invert_eig_s


  !-----------------------------------------------------------------------


  SUBROUTINE Jacobi(a,n,d,v,nrot)
    !
    ! Eigenvalues d and eigenvectors v of a matrix a.
    ! Jacobi algorithm from Numerical Recipes
    !
    integer,         intent(in)    :: n
    integer,         intent(out)   :: nrot
    double precision,intent(inout) :: a(n,n)
    double precision,intent(out)   :: d(n), v(n,n)
    integer                        :: i, j, ip, iq
    double precision               :: c, g, h, s, sm, t, tau, theta, tresh
    double precision               :: b(n), z(n)

    v = 0.0
    do ip = 1,n
       v(ip,ip) = 1.d0
       d(ip) = a(ip,ip)
    enddo
    b = d
    z = 0.0

    nrot = 0
    !do i = 1,50
    do i = 1,1000 ! -RK

       sm = 0.0
       do ip = 1,n-1
          do iq = ip+1,n
             sm = sm+abs(a(ip,iq))
          enddo
       enddo
       if (sm.eq.0) return
       if (i.lt.4) then
          tresh = 0.2*sm/(n*n)
       else
          tresh = 0.0
       endif

       do ip = 1,n-1
          do iq = ip+1,n
             g = 100.*abs(a(ip,iq))
             if ((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.  &
                  (abs(d(iq))+g.eq.abs(d(iq)))) then
                a(ip,iq) = 0.
             elseif (abs(a(ip,iq)).gt.tresh) then
                h = d(iq)-d(ip)
                if (abs(h)+g.eq.abs(h)) then
                   t = a(ip,iq)/h
                else
                   theta = 0.5d0*h/a(ip,iq)
                   t = 1/(abs(theta)+sqrt(1+theta*theta))
                   if (theta.lt.0) t=-t
                endif
                c = 1/sqrt(1+t*t)
                s = t*c
                tau = s/(1+c)
                h = t*a(ip,iq)
                z(ip) = z(ip)-h
                z(iq) = z(iq)+h
                d(ip) = d(ip)-h
                d(iq) = d(iq)+h
                a(ip,iq) = 0.
                do j = 1,ip-1
                   g = a(j,ip)
                   h = a(j,iq)
                   a(j,ip) = g-s*(h+g*tau)
                   a(j,iq) = h+s*(g-h*tau)
                enddo
                do j = ip+1,iq-1
                   g = a(ip,j)
                   h = a(j,iq)
                   a(ip,j) = g-s*(h+g*tau)
                   a(j,iq) = h+s*(g-h*tau)
                enddo
                do j = iq+1,n
                   g = a(ip,j)
                   h = a(iq,j)
                   a(ip,j) = g-s*(h+g*tau)
                   a(iq,j) = h+s*(g-h*tau)
                enddo
                do j = 1,n
                   g = v(j,ip)
                   h = v(j,iq)
                   v(j,ip) = g-s*(h+g*tau)
                   v(j,iq) = h+s*(g-h*tau)
                enddo
                nrot = nrot+1
             endif
          enddo
       enddo
       b = b+z
       d = b
       z = 0.
    enddo
    write(*,*) 'Too many iterations'
    stop

  END SUBROUTINE Jacobi

  !------------------------------------------------------------------------


  SUBROUTINE invert_LU(cc, cdet, sing, n)
    !
    ! Invert a symmetric matrix by LU decomposition
    !
    integer :: n 
    double precision,intent(inout) :: cc(n,n)
    double precision,intent(out) :: cdet
    logical,         intent(out) :: sing
    integer :: info, i, j

    call dpotrf( 'U', n, cc, n, info )

    if (info /= 0) then
       sing = .true.
       return
    end if

    cdet = 1
    do i = 1,n
       cdet = cdet * cc(i,i)
    end do

    call dpotri( 'U', n, cc, n, info )

    if (info /= 0) then
       cc = 0
       sing = .true.
       return
    end if

    ! copy to lower triangle

    do i = 1, n
       do j = i+1, n
          cc(j,i) = cc(i,j)
       end do
    end do

    sing = .false.

  END SUBROUTINE invert_LU

  !-----------------------------------------------------------------------

END MODULE matrix
