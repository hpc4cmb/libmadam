MODULE fourier

   implicit none
   private

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
   real(dp),parameter :: pi = 3.14159265358979323846d0

   integer,save :: nof
   complex(dp),allocatable,save :: datab(:)

   public init_fourier, close_fourier, dfft, dfftinv

CONTAINS


 SUBROUTINE init_fourier(nof_in)

   integer :: nof_in

   nof = nof_in
   if (nof.gt.0) allocate(datab(nof))

 END SUBROUTINE


 !-----------------------------------------------------------------------


 SUBROUTINE close_fourier

   if (allocated(datab)) deallocate(datab)

 END SUBROUTINE


 !-----------------------------------------------------------------------


 SUBROUTINE dfft(fyy,yy)
! Fourier transform of a real vector.

    complex(dp),intent(out) :: fyy(nof/2+1)
    real(dp),   intent(in)  :: yy(nof)

    datab = yy

    call fourc(datab,nof,1)

    fyy = datab(1:nof/2+1)

 END SUBROUTINE


!-----------------------------------------------------------------------


 SUBROUTINE dfftinv(yy,fyy)
! Inverse Fourier transform complex -> real.

    real(dp),   intent(out) :: yy(nof)
    complex(dp),intent(in)  :: fyy(nof/2+1)
    integer                 :: i

    datab(1:nof/2+1) = fyy
    do i = 2,nof/2
       datab(nof-i+2) = conjg(fyy(i))
    enddo

    call fourc(datab,nof,-1)

    yy = real(datab)/nof

 END SUBROUTINE


!-----------------------------------------------------------------------


 SUBROUTINE fourc(datab,n,isign)

   integer     :: isign,n
   complex(dp) :: datab(n)
   integer     :: i,istep,j,m,mmax
   complex(dp) :: temp,w,wp
   real(dp)    :: theta,wpr,wpi

   j = 1
   do i = 1,n
      if (j.gt.i) then
         temp = datab(j)
         datab(j) = datab(i)
         datab(i) = temp
      endif
      m = n/2
      do while ((m.ge.2).and.(j.gt.m))
         j = j-m
         m = m/2
      enddo
      j = j+m
   enddo

   mmax = 1
   do while (n.gt.mmax)
      istep = 2*mmax
      theta = pi/(isign*mmax)
      wpr = -2*sin(0.5d0*theta)**2
      wpi = sin(theta)
      wp = cmplx(wpr,wpi,dp)
      w = 1.d0
      do m = 1,mmax
         do i = m,n,istep
            j = i+mmax
            temp = w*datab(j)
            datab(j) = datab(i)-temp
            datab(i) = datab(i)+temp
         enddo
         w = w*wp+w
      enddo
      mmax = istep
   enddo

 END SUBROUTINE

END MODULE