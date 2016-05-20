
MODULE healpix_routines

! Routines related to Healpix pixelization

   implicit none
   private

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: isp = selected_int_kind(9)
   integer, parameter :: idp = selected_int_kind(18)

   real(dp),parameter :: pi = 3.14159265358979323846264d0
   real(dp),parameter :: twopi = 2*pi

   real(dp),parameter :: sqrt6  = 2.4494897427831780982d0
   real(dp),parameter :: twoopi = 2.d0/pi

   integer :: nside_nested = -1
   integer :: nside_ring   = -1

   integer,allocatable :: npixtab_nested(:)
   integer,allocatable :: npixtab_ring(:)


   public angle_to_nested,      &
          cartesian_to_nested,  &
          cartesian_to_ring

CONTAINS

!------------------------------------------------------------------------------


   SUBROUTINE angle_to_nested(npix,theta,phii,n,nside,idum_in)
!
! Transforms spherical coordinates (theta,phii) into nested healpix
! number NPIX at resolution NSIDE.
! Theta must be in range [0,pi] and phi in range [-8*pi,8*pi]
! If not, a dummy pixel value is returned

      integer, intent(out)         :: npix(n)          ! Pixel number
      real(dp),intent(in)          :: phii(n),theta(n) ! Spherical coordinates
      integer, intent(in)          :: n, nside
      integer, intent(in),optional :: idum_in          ! dummy pixel value

      integer  :: face, zone, ix, iy, k, idum
      real(dp) :: sx, sy, phi, z, x, y, a, p

      idum = 12*nside*nside
      if (present(idum_in)) idum=idum_in

      call set_pixtab_nested(nside)

      do k = 1,n

         if (theta(k).lt.0.or.theta(k).gt.pi.or.  &
             phii(k).lt.-8*pi.or.phii(k).gt.8*pi) then
             npix(k) = idum
             cycle
         endif

         if (theta(k).lt..5d0*pi) then
            a = sqrt6*sin(.5d0*theta(k))
            z = 0.75d0 -0.25d0*a*a
         else
            a = sqrt6*cos(.5d0*theta(k))
            z = 0.25d0*a*a -0.75d0
         endif

         phi = phii(k)
         do
            if (phi.lt.twopi) exit
            phi = phi-twopi
         enddo
         do
            if (phi.ge.0) exit
            phi = phi+twopi
         enddo

         p = phi*twoopi
         zone = int(p)
         p = p-zone-.5d0

         if (z.gt..5d0) then  ! north polar zone
            face = zone
            sy = 1.d0-.5d0*a -a*p
            sx = 1.d0-.5d0*a +a*p
         elseif (z.lt.-.5d0) then! south polar zone
            face = zone+8
            sx = .5d0*a +a*p
            sy = .5d0*a -a*p
         else
            y = z-p
            x = z+p
            if (x.ge.0.and.y.ge.0) then! northern equatorial zone
               face = zone
               sy = y
               sx = x
            elseif (x.ge.0) then ! equatorial zone
               if (zone==3) then
                  face = 4
               else
                  face = zone+5
               endif
               sy = y +1.d0
               sx = x
            elseif (y.ge.0) then ! equatorial zone
               face = zone+4
               sy = y
               sx = x +1.d0
            else ! southern equatorial zone
               face = zone+8
               sy = 1.d0 +y
               sx = 1.d0 +x
            endif
         endif

         ix = sx*nside
         iy = sy*nside

         npix(k) = face*nside**2 +npixtab_nested(ix)+2*npixtab_nested(iy)
      enddo

   END SUBROUTINE angle_to_nested


!------------------------------------------------------------------------------


   SUBROUTINE cartesian_to_nested(npix,xyz,n,nside,idum_in)
!
! Transforms cartesian coordibates (x,y,z) into nested healpix number
!   at resolution NSIDE.
! Call init_healpix first to speed up computation.

      integer, intent(out)         :: npix(n)    ! Pixel number
      real(dp),intent(in)          :: xyz(3,n)   ! Location in cartesian coordinates
      integer, intent(in)          :: n, nside
      integer, intent(in),optional :: idum_in ! dummy pixel value

      integer              :: face, zone, ix, iy, k, idum
      real(dp)             :: sx, sy, phi, z, x, y, a, p
      real(dp),allocatable :: phitab(:)

      call set_pixtab_nested(nside)

      allocate(phitab(n))

      call atan_approx(phitab,xyz,n)

      idum = 12*nside*nside
      if (present(idum_in)) idum=idum_in

      do k = 1,n

         if (all(xyz(:,k)==0.0).or.   &
             any(xyz(:,k).gt.1.d0).or. &
             any(xyz(:,k).lt.-1.d0)) then
            npix(k) = idum
            cycle
         endif

         z = 0.75d0*xyz(3,k)

!         phi = atan2(xyz(2,k),xyz(1,k))
         phi = phitab(k)

         if (phi.lt.0) phi=phi+twopi
         p = phi*twoopi

         zone = int(p)
         p = p-zone-.5d0

         if (z.gt..5d0) then  ! north polar zone
            face = zone
            a = sqrt(3.d0-4.d0*z)
            sy = 1.d0-.5d0*a -a*p
            sx = 1.d0-.5d0*a +a*p
         elseif (z.lt.-.5d0) then ! south polar zone
            face = zone+8
            a = sqrt(3.d0+4.d0*z)
            sx = .5d0*a +a*p
            sy = .5d0*a -a*p
         else
            y = z-p
            x = z+p
            if (x.ge.0.and.y.ge.0) then! northern equatorial zone
               face = zone
               sy = y
               sx = x
            elseif (x.ge.0) then ! equatorial zone
               if (zone==3) then
                  face = 4
               else
                  face = zone+5
               endif
               sy = y +1.d0
               sx = x
            elseif (y.ge.0) then ! equatorial zone
               face = zone+4
               sy = y
               sx = x +1.d0
            else  ! southern equatorial zone
               face = zone+8
               sy = 1.d0 +y
               sx = 1.d0 +x
            endif
         endif

         ix = sx*nside
         iy = sy*nside

         npix(k) = face*nside**2 +npixtab_nested(ix)+2*npixtab_nested(iy)

      enddo

      deallocate(phitab)

   END SUBROUTINE cartesian_to_nested


!------------------------------------------------------------------------------


   SUBROUTINE angle_to_ring(npix,phii,theta,n,nside,idum_in)
!
! Transforms spherical coordinates (theta,phi) into ring-scheme healpix number
!   at resolution NSIDE.

      integer, intent(out)         :: npix(n)          ! Pixel number
      real(dp),intent(in)          :: phii(n),theta(n) ! Spherical coordinates
      integer, intent(in)          :: n, nside
      integer, intent(in),optional :: idum_in  ! dummy pixel value

      integer  :: zone, ix, iy, k, ktheta, kphi, idum
      real(dp) :: sx, sy, phi, z, a, p

      idum = 12*nside*nside
      if (present(idum_in)) idum=idum_in

      call set_pixtab_ring(nside)

      do k = 1,n

         if (theta(k).lt.0.or.theta(k).gt.pi.or.  &
             phii(k).lt.-8*pi.or.phii(k).gt.8*pi) then
             npix(k) = idum
             cycle
         endif

         if (theta(k).lt..5d0*pi) then
            a = sqrt6*sin(.5d0*theta(k))
            z = 0.75d0 -0.25d0*a*a
         else
            a = sqrt6*cos(.5d0*theta(k))
            z = 0.25d0*a*a -0.75d0
         endif

         phi = phii(k)
         do
            if (phi.lt.twopi) exit
            phi = phi-twopi
         enddo
         do
            if (phi.ge.0) exit
            phi = phi+twopi
         enddo

         p = phi*twoopi
         zone = int(p)

         if (z.gt..5d0) then  ! north polar zone
            p = p-zone-.5d0

            sy = 1.d0-.5d0*a -a*p
            sx = 1.d0-.5d0*a +a*p

            ix = sx*nside
            iy = sy*nside

            ktheta = 2*nside-1-ix-iy
            kphi = zone*ktheta +nside-iy-1

         elseif (z.lt.-.5d0) then  ! south polar zone
            p = p-zone-.5d0

            sx = .5d0*a +a*p
            sy = .5d0*a -a*p

            ix = sx*nside
            iy = sy*nside

            ktheta = ix+iy+1
            kphi = zone*ktheta +ix

            ktheta = 4*nside-ktheta
         else   ! equatorial zone
            iy = (z-p)*nside
            ix = (z+p)*nside

            ktheta = 2*nside-ix-iy
            kphi = (ix-iy)/2

         endif

         npix(k) = npixtab_ring(ktheta) +kphi

      enddo

    END SUBROUTINE angle_to_ring


!------------------------------------------------------------------------------


   SUBROUTINE cartesian_to_ring(npix,xyz,n,nside,idum_in)
!
! Transforms cartesian coordibates (x,y,z) into ring-scheme healpix number
!   at resolution NSIDE.

      integer, intent(out)         :: npix(n)  ! Pixel number
      real(dp),intent(in)          :: xyz(3,n) ! Location in cartesian coord.
      integer, intent(in)          :: n, nside
      integer, intent(in),optional :: idum_in

      integer              :: zone, ix, iy, k, ktheta, kphi, idum
      real(dp)             :: sx, sy, phi, z, a, p
      real(dp),allocatable :: phitab(:)

      call set_pixtab_ring(nside)

      allocate(phitab(n))

      call atan_approx(phitab,xyz,n)

      idum = 12*nside*nside
      if (present(idum_in)) idum=idum_in

      do k = 1,n

         if (all(xyz(:,k)==0.0).or.   &
             any(xyz(:,k).gt.1.d0).or. &
             any(xyz(:,k).lt.-1.d0)) then
            npix(k) = idum
            cycle
         endif

         z = 0.75d0*xyz(3,k)

!         phi = atan2(xyz(2,k),xyz(1,k))
         phi = phitab(k)
         if (phi.lt.0) phi=phi+twopi

         p = phi*twoopi
         zone = int(p)

         if (z.gt..5d0) then  ! north polar zone
            a = sqrt(3.d0-4.d0*z)
            p = p-zone-.5d0

            sy = 1.d0-.5d0*a -a*p
            sx = 1.d0-.5d0*a +a*p

            ix = sx*nside
            iy = sy*nside

            ktheta = 2*nside-1-ix-iy
            kphi = zone*ktheta +nside-iy-1

         elseif (z.lt.-.5d0) then  ! south polar zone
            a = sqrt(3.d0+4.d0*z)
            p = p-zone-.5d0

            sx = .5d0*a +a*p
            sy = .5d0*a -a*p

            ix = sx*nside
            iy = sy*nside

            ktheta = ix+iy+1
            kphi = zone*ktheta +ix

            ktheta = 4*nside-ktheta
         else   ! equatorial zone
            iy = (z-p)*nside
            ix = (z+p)*nside

            ktheta = 2*nside-ix-iy
            kphi = (ix-iy)/2

         endif

         npix(k) = npixtab_ring(ktheta) +kphi

      enddo

      deallocate(phitab)

   END SUBROUTINE cartesian_to_ring


!------------------------------------------------------------------------------


   SUBROUTINE set_pixtab_nested(nside)

! In NESTED Healpix pixeling scheme, the index of a pixel identified by
! cartesian (integer) coordinates ix,iy (in range 0...nside-1) inside 
! a base pixel, can be written as
! ipix = face*nside**2 +npixtab_nested(ix)+2*npixtab_nested(iy)
! where FACE identifies the base pixel, and NPIXTAB_NESTED is a function of 
! one coordinate only.
! This routine computes and stores table NPIXTAB_NESTED for given resolution.
! The precomputed table is used to speed up routines
!  angle_to_nested and cartesian_to_nested.

      integer,intent(in) :: nside
      integer            :: ipix, dx, ix, i

      if (nside==nside_nested) return

      if (allocated(npixtab_nested)) deallocate(npixtab_nested)

      if (nside.lt.0) then
          nside_nested = -1
          return
      endif

      nside_nested = nside
      allocate(npixtab_nested(0:nside_nested))

      do i = 0,nside

         ix = i
         ipix = 0
         dx = nside

         do
            dx = dx/2
            if (dx==0) exit

            ipix = 4*ipix
            if (ix.ge.dx) then
               ipix = ipix+1
               ix = ix-dx
            endif

         enddo

         npixtab_nested(i) = ipix
      enddo

   END SUBROUTINE


!-------------------------------------------------------------------------


   SUBROUTINE set_pixtab_ring(nside)

! Compute the pixel number of the first pixel in each row at given resolution.
! The pixel numbers are stored in table NPIXTAB_RING. The last element gives
! the total number of pixels.
! The precomputed table will be used to speed up routines angle_to_ring
! and cartesian_to_ring

      integer,intent(in) :: nside
      integer            :: k, nrow, ncum

      if (nside==nside_ring) return

      if (allocated(npixtab_ring)) deallocate(npixtab_ring)

      if (nside.lt.0) then
          nside_ring = -1
          return
      endif

      nside_ring = nside
      allocate(npixtab_ring(4*nside))

      ncum = 0
      nrow = 0

      do k = 1,nside
         ncum = ncum+nrow
         nrow = nrow+4
         npixtab_ring(k) = ncum
      enddo

      do k = nside+1,3*nside
         ncum = ncum+nrow
         npixtab_ring(k) = ncum
      enddo

      do k = 3*nside+1,4*nside
         ncum = ncum+nrow
         nrow = nrow-4
         npixtab_ring(k) = ncum
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE downgrade_nested(ipix,n,nside1,nside2)

      integer,intent(in)    :: n, nside1, nside2
      integer,intent(inout) :: ipix(n)
      integer               :: i, nside

      if (powertwo(nside1).lt.0.or.powertwo(nside2).lt.0) then
         write(*,*) 'ERROR in downgrade_nested: '
         write(*,*) 'nside1, nside2 must be powers of two.'
         write(*,*) 'nside1, nside2 =',nside1,nside2
         stop
      elseif (nside2.gt.nside1) then
         write(*,*) 'ERROR in downgrade_nested: '
         write(*,*) 'nside2 must be equal to or lower than nside1.'
         write(*,*) 'nside1, nside2 =',nside1,nside2
         stop
      endif

      nside = nside1
      do
         if (nside==nside2) exit

         do i = 1,n
            ipix(i) = ipix(i)/4
         enddo
         nside = nside/2
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE downgrade_ring(ipix,n,nside1,nside2)

      integer,intent(in)    :: n, nside1, nside2
      integer,intent(inout) :: ipix(n)
      integer               :: i, k, nside, ir, nopix, nzone, izone, ie
      integer,allocatable   :: irow(:)

      if (powertwo(nside1).lt.0.or.powertwo(nside2).lt.0) then
         write(*,*) 'ERROR in downgrade_ring: '
         write(*,*) 'nside1, nside2 must be powers of two.'
         write(*,*) 'nside1, nside2 =',nside1,nside2
         stop
      elseif (nside2.gt.nside1) then
         write(*,*) 'ERROR in downgrade_ring: '
         write(*,*) 'nside2 must be equal to or lower than nside1.'
         write(*,*) 'nside1, nside2 =',nside1,nside2
         stop
      endif

      allocate(irow(n))

      nside = nside1
      do
         if (nside==nside2) exit

         nopix = 12*nside*nside

         call set_pixtab_ring(nside)

         ir = 1
         do i = 1,n
            if (ipix(i).lt.0.or.ipix(i).ge.nopix) then
               write(*,*) 'ERROR in downgrade_ring:'
               write(*,*) 'Illegal pixel number',ipix(i)
               stop
            endif

! Find the row where the pixel belongs
            do
               if (ipix(i).ge.npixtab_ring(ir+1)) then
                  ir = ir+1
               else
                  exit
               endif
            enddo
            do
               if (ipix(i).lt.npixtab_ring(ir)) then
                   ir = ir-1
               else
                  exit
               endif
            enddo
            ie = ipix(i)-npixtab_ring(ir)

            if (ir.le.nside) then
               if (mod(ir,2)==0) then
                  irow(i) = ir/2
                  ie = ie/2
               else
                  nzone = ir
                  izone = ie/nzone
                  ie = mod(ie,nzone)
                  if (mod(ie,2)==0) then
                     irow(i) = (ir+1)/2
                     ie = ie/2 +izone*(nzone+1)/2
                  else
                     irow(i) = (ir-1)/2
                     ie = (ie-1)/2 +izone*(nzone-1)/2
                  endif
               endif
            elseif (ir.le.3*nside) then
               k = mod(ir-nside,4)
               if (k==0) then
                  irow(i) = ir/2
                  ie = ie/2
               elseif (k==1) then
                  if (mod(ie,2)==0) then
                     irow(i) = (ir+1)/2
                  else
                     irow(i) = (ir-1)/2
                  endif
                  ie = ie/2
               elseif (k==2) then
                  irow(i) = ir/2
                  ie = (ie+1)/2
                  if (ie==2*nside) ie=0
               elseif (k==3) then
                  if (mod(ie,2)==0) then
                     irow(i) = (ir-1)/2
                  else
                     irow(i) = (ir+1)/2
                  endif
                  ie = ie/2
               endif
            else
               if (mod(ir,2)==0) then
                  irow(i) = ir/2
                  ie = ie/2
               else
                  nzone = 4*nside-ir
                  izone = ie/nzone
                  ie = mod(ie,nzone)
                  if (mod(ie,2)==0) then
                     irow(i) = (ir-1)/2
                     ie = ie/2 +izone*(nzone+1)/2
                  else
                     irow(i) = (ir+1)/2
                     ie = (ie-1)/2 +izone*(nzone-1)/2
                  endif
               endif
            endif

            ipix(i) = ie
         enddo

         nside = nside/2
         call set_pixtab_ring(nside)

         do i = 1,n
            ipix(i) = ipix(i)+npixtab_ring(irow(i))
         enddo

      enddo

      deallocate(irow)

   END SUBROUTINE


!------------------------------------------------------------------------------


   integer FUNCTION powertwo(n)

      integer :: n, k, m

      powertwo = -1

      m = 1
      do k = 0,31
         if (n==m) then
             powertwo = k
             exit
         else
            m = 2*m
         endif
      enddo

   END FUNCTION


!------------------------------------------------------------------------------


   SUBROUTINE get_costheta(xtheta,nside)
! Find cos(theta) corresponding to the 4*nside-1 healpix rows
      integer, intent(in)  :: nside
      real(dp),intent(out) :: xtheta(4*nside-1)
      integer              :: i, npix, iarea
      real(dp)             :: area

      npix = 0  ! pixels/row
      iarea = 0
      do i = 1,4*nside-1
         iarea = iarea+npix/2

         if (i.le.nside) then
             npix = npix+4
         elseif (i.gt.3*nside) then
             npix = npix-4
         endif
         iarea = iarea+npix/2
         area = iarea/(12.d0*nside*nside)

         xtheta(i) = 1.d0-2.d0*area
      enddo

   END SUBROUTINE


!-------------------------------------------------------------------------------


   SUBROUTINE ringmap_to_nested(map,nopix,nomaps)
! Transforms a map from RING to NESTED scheme.

      real,   intent(inout) :: map(0:nopix-1,nomaps)
      integer,intent(in)    :: nopix,nomaps
      integer               :: ns,nside,i,k,imax,ipix,i4,dx,ix,iy,face,n
      integer,allocatable   :: nring(:),iv(:),ih(:),row(:),col(:),nrow(:),ind(:)
      real,   allocatable   :: buffer(:)

      nside = sqrt(nopix/12+1.)
      ns = 0
      k = nside
      do while (k>1)
         k = k/2
         ns = ns+1
      enddo

      n = nside*nside
      allocate(iv(n),ih(n),row(n),col(n),nring(n),nrow(n))
      allocate(ind(0:nopix-1))

      do i = 0,n-1
        ipix = i
        dx = nside
        imax = nside*nside
        ix = 0
        iy = 0
        do k = 1,ns
            dx = dx/2
            imax = imax/4
            i4 = ipix/imax
            ipix = MOD(ipix,imax)
            if (i4==1.OR.i4==3)  ix=ix+dx
            if (i4==2.OR.i4==3)  iy=iy+dx
         enddo
         iv(i+1) = ix+iy
         ih(i+1) = ix-iy
      enddo

      face = 0
      row = 2*nside-1-iv
      where (row.lt.nside)
         nring = 2*row*(row-1)+(row-1+ih)/2
      elsewhere
         col = nside+ih
         nring = 2*nside*(nside-1)+(row-nside)*4*nside +col/2
      end where
      ind(0:n-1) = nring

      do face=1,3
         where (row.lt.nside)
            nring = nring+row
         elsewhere
            nring = nring+nside
         end where
         ind(face*n:(face+1)*n-1) = nring
      enddo

      row = row+nside

      face = 4
      col = ih
      where (col<0) col=col+8*nside
      nring = 2*nside*(nside-1)+(row-nside)*4*nside+col/2
      ind(4*n:5*n-1) = nring

      face = 5
      col = 2*nside+ih
      nring = 2*nside*(nside-1)+(row-nside)*4*nside +col/2
      ind(5*n:6*n-1) = nring

      face = 6
      nring = nring+nside
      ind(6*n:7*n-1) = nring

      face = 7
      nring = nring+nside
      ind(7*n:8*n-1) = nring

      row = row+nside

      face = 8
      where (row.le.3*nside)
         col = nside+ih
         nring = 2*nside*(nside-1)+(row-nside)*4*nside+col/2
      elsewhere
         nrow = 4*nside-row
         nring = 2*nside*(nside-1)+8*nside*nside
         nring = nring+2*nside*(nside+1)-2*nrow*(nrow+1)
         nring = nring+(nrow-1+ih)/2
      end where
      ind(8*n:9*n-1) = nring

      do face = 9,11
         where (row.le.3*nside)
            nring = nring+nside
         elsewhere
            nring = nring+nrow
         end where
         ind(face*n:(face+1)*n-1)=nring
      enddo

      deallocate(iv,ih,row,col,nrow,nring)

      allocate(buffer(0:nopix-1))

      do k = 1,nomaps
         do i = 0,nopix-1
           buffer(i) = map(ind(i),k)
           if (ind(i).lt.0.or.ind(i).gt.nopix-1) write(*,*) i,ind(i)
         enddo
         do i = 0,nopix-1
            map(i,k) = buffer(i)
         enddo
      enddo

      deallocate(ind,buffer)

   END SUBROUTINE ringmap_to_nested


!-------------------------------------------------------------------------

   SUBROUTINE atan_approx(a,xyz,n)
!  Approximate atan2(y,x) by a rational function
!  A routine obtained from Ted Kisner.

      integer, intent(in)  :: n
      real(dp),intent(in)  :: xyz(3,n)
      real(dp),intent(out) :: a(n)
      integer              :: i, sector
      real(dp)             :: r, r2, x, y, theta
      logical              :: mirror, add_pio6

      real(dp), parameter :: cheb1 = 48.70107004404898384
      real(dp), parameter :: cheb2 = 49.5326263772254345
      real(dp), parameter :: cheb3 =  9.40604244231624
      real(dp), parameter :: cheb4 = 48.70107004404996166
      real(dp), parameter :: cheb5 = 65.7663163908956299
      real(dp), parameter :: cheb6 = 21.587934067020262

      real(dp),parameter :: pio2 = pi/2
      real(dp),parameter :: pio6 = pi/6

      real(dp),parameter :: tan_pio12 = 0.267949192431123
      real(dp),parameter :: tan_pio6  = 0.577350269189626

      do i = 1,n

         x = xyz(1,i)
         y = xyz(2,i)

         if (x.lt.0) then
            x = -x
            if (y.lt.0) then
               sector = 3
               y = -y
            else
               sector = 2
            endif
         else
            if (y.lt.0) then
               sector = 4
               y = -y
            else
               sector = 1
            endif
         endif

         if (x==0.0.and.y==0.0) then
            r = 0.0
            mirror = .false.
         elseif (y.gt.x) then
            r = x/y
            mirror = .true.
         else
            r = y/x
            mirror = .false.
         endif

         if (r.gt.tan_pio12) then
            r = (r-tan_pio6)/(1+tan_pio6*r)
            add_pio6 = .true.
         else
            add_pio6 = .false.
         endif

         r2 = r*r
         theta = r*(cheb1+(cheb2+cheb3*r2)*r2)/(cheb4 +(cheb5+r2*(cheb6+r2))*r2)

         if (add_pio6) theta = theta+pio6
         if (mirror)   theta = pio2-theta

         if (sector==2) then
            theta = pi-theta
         elseif (sector==3) then
            theta = theta+pi
         elseif (sector==4) then
            theta = twopi-theta
         endif

         a(i) = theta
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------

END MODULE healpix_routines
