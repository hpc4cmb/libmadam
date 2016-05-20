MODULE coordinate_conversion

! Routines for transform between coordinate systems.
!
! A detector's position may be given either
! (1) as three orthogonal vectors COORD(3,3) or
! (2) as angles (theta,phi,psi), where theta and phi define a point
!     in the celestial sphere and phi defines detector's orientation.
! This module provides coordinate trasform routines for both.
!
! Definitions from Trafo (LevelS)
! -----------------------------------------------------------------
!        trtype: specifies type of transform
!                0     No conversion
!                1     RA-Dec (2000) -> Galactic
!                2     Galactic      -> RA-Dec
!                3     RA-Dec        -> Ecliptic
!                4     Ecliptic      -> RA-Dec
!                5     Ecliptic      -> Galactic
!                6     Galactic      -> Ecliptic
! Celestial coordinates (RA, Dec) should be given in equinox
! J2000 unless ifk is present and equal to 1.
! -----------------------------------------------------------------

! Euler angles copied from Trafo module from LevelS

   implicit none
   private

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)

   real(dp), dimension(6,2), parameter :: gammat = reshape( &
     (/ 0.57595865315_dp, 4.92619181360_dp, 0.00000000000_dp, &
        0.00000000000_dp, 0.11129056012_dp, 4.70053728340_dp, &
        0.57477043300_dp, 4.93682924650_dp, 0.00000000000_dp, &
        0.00000000000_dp, 0.11142137093_dp, 4.71279419371_dp /), (/6,2/) )

  real(dp), dimension(6,2), parameter :: betat = reshape( &
     (/ 1.09257611175_dp,-1.09257611175_dp, 0.40920620808_dp, &
       -0.40920620808_dp, 1.05047957320_dp,-1.05047957320_dp, &
        1.09731904398_dp,-1.09731904398_dp, 0.40909280422_dp, &
       -0.40909280422_dp, 1.05048856981_dp,-1.05048856981_dp /), (/6,2/) )

   real(dp), dimension(6,2), parameter :: alphat = reshape( &
     (/ 4.92619181360_dp, 0.57595865315_dp, 0.00000000000_dp, &
        0.00000000000_dp, 4.70053728340_dp, 0.11129056012_dp, &
        4.93682924650_dp, 0.57477043300_dp, 0.00000000000_dp, &
        0.00000000000_dp, 4.71279419371_dp, 0.11142137093_dp /) , (/6,2/) )

   real(dp),save :: sin_alpha=0.0, cos_alpha=1.0,  &
                    sin_beta =0.0, cos_beta =1.0,  &
                    sin_gamma=0.0, cos_gamma=1.0

   real(dp),save :: alpha=0.0,beta=0.0,gamma=0.0

   real(dp),save :: trot(3,3) = 0.0
   integer, save :: ifk_stored = 2
   integer, save :: trtype_stored = -1

   public convert_coord, convert_angle

CONTAINS

!------------------------------------------------------------------------------


   SUBROUTINE convert_coord(coord,n,trtype,ifk_in)
!
! Given a detector's positition as 3 vectors (COORD) in one coordinate system,
!  transform it into another coordinate system.
!
      integer, intent(in)          :: n, trtype
      real(dp),intent(inout)       :: coord(3,3,n)
      integer, intent(in),optional :: ifk_in
      integer                      :: ifk

      ifk = 2
      if (present(ifk_in)) ifk=ifk_in

      call init_coordinate_transform(trtype,ifk)

      call rotate_coord(coord,n)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE convert_angle(theta,phi,psi,n,trtype,ifk_in)
!
! Given a detector's position as angles (theta,phi,psi) in one coordinate
!  system, transform it into another.
!
      integer, intent(in)          :: n, trtype
      integer, intent(in),optional :: ifk_in
      real(dp),intent(inout)       :: theta(n), phi(n), psi(n)
      integer                      :: ifk, k
      real(dp)                     :: x, y, z
      real(dp)                     :: sin_theta, cos_theta, sin_phi, cos_phi

      ifk = 2
      if (present(ifk_in)) ifk=ifk_in

      call init_coordinate_transform(trtype,ifk)

      do k = 1,n
         sin_theta = sin(theta(k))
         cos_theta = cos(theta(k))
         sin_phi   = sin(phi(k)-alpha)
         cos_phi   = cos(phi(k)-alpha)

         z = -sin_beta*sin_theta*sin_phi +cos_beta*cos_theta
         theta(k) = acos(z)

         y = cos_beta*sin_theta*sin_phi +sin_beta*cos_theta
         x = sin_theta*cos_phi
         phi(k) = atan2(y,x) +gamma

         y = -cos_phi*sin_beta
         x = sin_theta*cos_beta +cos_theta*sin_phi*sin_beta
         psi(k) = psi(k) +atan2(y,x)

      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE convert_angle2coord(coord,theta,phi,psi,n,trtype,ifk_in)

      integer, intent(in)          :: n, trtype
      real(dp),intent(out)         :: coord(3,3,n)
      integer, intent(in),optional :: ifk_in
      real(dp),intent(in)          :: theta(n), phi(n), psi(n)
      integer                      :: ifk, k
      real(dp)                     :: ephi(3), etheta(3)
      real(dp)                     :: sin_theta, cos_theta, &
                                      sin_phi, cos_phi, sin_psi, cos_psi

      ifk = 2
      if (present(ifk_in)) ifk=ifk_in

      call init_coordinate_transform(trtype,ifk)

      do k = 1,n
         sin_theta = sin(theta(k))
         cos_theta = cos(theta(k))

         sin_phi = sin(phi(k))
         cos_phi = cos(phi(k))

         sin_psi = sin(psi(k))
         cos_psi = cos(psi(k))

         coord(1,3,k) = sin_theta*cos_phi
         coord(2,3,k) = sin_theta*sin_phi
         coord(3,3,k) = cos_theta

         ephi(1) = -sin_phi
         ephi(2) =  cos_phi
         ephi(3) =  0.0

         etheta(1) =  cos_theta*cos_phi
         etheta(2) =  cos_theta*sin_phi
         etheta(3) = -sin_theta

         coord(:,1,k) = -cos_psi*etheta+sin_psi*ephi
         coord(:,2,k) =  cos_psi*ephi  +sin_psi*etheta

      enddo

      call rotate_coord(coord,n)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE convert_coord2angle(theta,phi,psi,coord,n,trtype,ifk_in)

      integer, intent(in)          :: n, trtype
      real(dp),intent(out)         :: theta(n), phi(n), psi(n)
      real(dp),intent(in)          :: coord(3,3,n)
      integer, intent(in),optional :: ifk_in
      integer                      :: ifk, k
      real(dp)                     :: x, y, z, c(3,3), ephi(3), cpsi, spsi

      ifk = 2
      if (present(ifk_in)) ifk=ifk_in

      call init_coordinate_transform(trtype,ifk)

      do k = 1,n

         c = coord(:,:,k)
         call rotate_coord(c,1)

         x = c(1,3)
         y = c(2,3)
         z = c(3,3)

         theta(k) = acos(z)
         phi(k) = atan2(y,x)

         ephi(1) = -y
         ephi(2) =  x
         ephi(3) =  0.0

         spsi = c(1,1)*ephi(1) +c(2,1)*ephi(2)
         cpsi = c(1,2)*ephi(2) +c(2,2)*ephi(2)

         psi(k) = atan2(spsi,cpsi)
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE coord2quweight(qw,uw,coord,n)
!
! Given detector's position COORD, compute Stokes weights
!   qw=cos(2*psi), uw=sin(2*psi).
!
      integer, intent(in)  :: n
      real(dp),intent(out) :: qw(n), uw(n)
      real(dp),intent(in)  :: coord(3,3,n)
      real(dp)             :: ephi(3), spsi, cpsi, cf
      integer              :: k

      do k = 1,n
         ephi(1) = -coord(2,3,k)
         ephi(2) =  coord(1,3,k)
         ephi(3) =  0.0

         spsi = coord(1,1,k)*ephi(1) +coord(2,1,k)*ephi(2)
         cpsi = coord(1,2,k)*ephi(2) +coord(2,2,k)*ephi(2)

         cf =  1.d0/(cpsi*cpsi +spsi*spsi)
         qw(k) = (cpsi*cpsi -spsi*spsi)*cf  ! cos(2*psi)
         uw(k) = 2.d0*cpsi*spsi*cf          ! sin(2*psi)
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE init_coordinate_transform(trtype,ifk)
!
! Store the rotation matrix for given type of coordinate transform.
!
      integer,intent(in)          :: trtype
      integer,intent(in),optional :: ifk

      if (ifk==ifk_stored.and.trtype==trtype_stored) return

      if (trtype.gt.6) then
         write(*,*) 'ERROR: Illegal value of trtype =',trtype
         stop
      elseif (ifk.gt.2.or.ifk.le.0) then
         write(*,*) 'ERROR: Illegal value of ifk =',ifk
         stop
      elseif (trtype==0) then
         alpha = 0.0
         beta = 0.0
         gamma = 0.0
      else
         alpha = alphat(trtype,ifk)
         beta  = betat(trtype,ifk)
         gamma = gammat(trtype,ifk)
      endif

      trtype_stored = trtype
      ifk_stored = ifk

      cos_alpha = cos(alpha)
      sin_alpha = sin(alpha)
      cos_beta  = cos(beta)
      sin_beta  = sin(beta)
      cos_gamma = cos(gamma)
      sin_gamma = sin(gamma)

      trot(1,1) =  cos_alpha*cos_gamma +sin_alpha*cos_beta*sin_gamma
      trot(2,1) =  sin_alpha*cos_gamma -cos_alpha*cos_beta*sin_gamma
      trot(3,1) = -sin_beta*sin_gamma

      trot(1,2) =  cos_alpha*sin_gamma -sin_alpha*cos_beta*cos_gamma
      trot(2,2) =  sin_alpha*sin_gamma +cos_alpha*cos_beta*cos_gamma
      trot(3,2) =  sin_beta*cos_gamma

      trot(1,3) =  sin_alpha*sin_beta
      trot(2,3) = -cos_alpha*sin_beta
      trot(3,3) =  cos_beta

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE rotate_coord(coord,n)

      integer, intent(in)    :: n
      real(dp),intent(inout) :: coord(3,3,n)
      integer                :: j, k
      real(dp)               :: x, y, z

      do k = 1, n
         do j = 1,3
            x = coord(1,j,k)
            y = coord(2,j,k)
            z = coord(3,j,k)
            coord(1,j,k) = x*trot(1,1) +y*trot(2,1) +z*trot(3,1)
            coord(2,j,k) = x*trot(1,2) +y*trot(2,2) +z*trot(3,2)
            coord(3,j,k) = x*trot(1,3) +y*trot(2,3) +z*trot(3,3)
         enddo
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------

END MODULE