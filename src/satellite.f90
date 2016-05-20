MODULE satellite
!
! satellite pointing -> detector pointing
!
! How to use this module:
!
! Call Read_satellite_position to read and store satellite pointing data
!  for an interval STARTTIME-ENDTIME.
! Select a detector by calling SELECT_DETECTOR_SATPOINT
! You may then repeatedly call DETECTOR_POINTING_IQU or DETECTOR_POINTING_PIX
!   to compute detector pointing for an interval inside the interval previously
!   defined in the call of READ_SATELLITE_POSITION
!
   use healpix_routines
   use planck_config
   use fitsmod2
   use coordinate_conversion

   implicit none
   private

   real(dp),parameter :: deg2rad = pi/180.d0

 ! Time interval in satellite pointing
 ! and time of the first sample
   real(dp),parameter :: sat_step   = 1.d0
   real(dp),parameter :: sat_start = 0.5d0

   character(len=*),parameter :: idf  = '(i3,": ",a,x,i0,2x,i0)'
   character(len=*),parameter :: idfe = '(i3,": ",a,x,2es12.4)'

   logical, save,public :: verbose_satellite = .false.
   integer, save,public :: ID_sat = 0

   real(sp),save,public :: memory_satellite = 0.0

   real(dp),save :: detrot(3,3) =  0.0
   real(dp),save :: fsample_det = -1.0
   logical, save :: polar_det   = .false.

   integer, save,public :: nosteps_satellite_file = -1
   integer, save        :: nosteps_allocated = -1
   integer, save        :: nosteps_read = -1
   integer, save        :: firststep_read, laststep_read
   real(dp),save        :: starttime_read, endtime_read
   real(dp),allocatable :: satpos(:,:,:)

   integer, save,public        :: noperiods = -1
   real(dp),allocatable,public :: starttime_period(:)
   real(dp),allocatable,public :: endtime_period(:)
   integer, allocatable,public :: firststep_period(:)
   integer, allocatable,public :: laststep_period(:)

   public read_satellite_position, read_pointing_periods, close_satellite, &
          select_detector_satpoint, detector_pointing_pix, detector_pointing_iqu

 CONTAINS

!------------------------------------------------------------------------------


   SUBROUTINE select_detector_satpoint(theta,phi,psi,theta_b,fsample,  &
                                       kpolar,unit)
!
! Set detector's position on the focal plane and sampling frequency
! This routine must be called before calling DETECTOR_POINTING_IQU/PIX.

     ! Detector's position on the focal plane
      real(dp),        intent(in)          :: theta, phi, psi !detector position
      real(dp),        intent(in)          :: theta_b ! opening angle
      real(dp),        intent(in)          :: fsample ! sampling frequency (Hz)
      logical,         intent(in),optional :: kpolar ! default: T
      character(len=*),intent(in),optional :: unit ! deg/rad (default:rad)

      if (verbose_satellite.and.ID_sat==0) then
         write(*,*) 'Selecting detector for satellite pointing conversion:'
         write(*,'(x,a,f10.5)') 'fsample =',fsample
         write(*,'(x,a,f10.5)') 'theta   =',theta
         write(*,'(x,a,f10.5)') 'phi     =',phi
         write(*,'(x,a,f10.5)') 'theta_b =',theta_b
         if (present(kpolar)) write(*,*) 'kpolar  =',kpolar
         if (present(unit))   write(*,*) 'unit    =',unit
      endif

      if (fsample.le.0) then
         write(*,idf)  ID_sat,'ERROR in select_detector_satpoint:'
         write(*,idfe) ID_sat,'Illegal value of fsample =',fsample
         call exit_with_status(1)
      endif

      fsample_det = fsample

      polar_det = .true.
      if (present(kpolar)) polar_det=kpolar

      if (present(unit)) then

         if (unit=='RAD'.or.unit=='rad') then
            detrot = detector_rotation(theta,phi,psi,theta_b)

         elseif (unit=='DEG'.or.unit=='deg') then
            detrot = detector_rotation(theta*deg2rad,phi*deg2rad,  &
                                        psi*deg2rad,theta_b*deg2rad)
         else
            write(*,*) 'ERROR in select_detector_satpoint: unknown unit'
            call exit_with_status(1)
         endif

      else  ! default: radians
         detrot = detector_rotation(theta,phi,psi,theta_b)

      endif

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE read_satellite_position(starttime,endtime,file_satellite,trtype)

! Read satellite pointing to cover the time interval STARTTIME-ENDTIME.
! Store in table SATPOS.
! If ENDTIME<0 read all data.

      real(dp),        intent(in)          :: starttime, endtime
      character(len=*),intent(in)          :: file_satellite
      integer,         intent(in),optional :: trtype !coordinate conversion

      type(fitshandle)     :: inhandle
      integer(i8b)         :: offset
      integer              :: nosteps_left, i, m, n, iperiod
      real(dp)             :: total_time

      integer, parameter   :: nbuffer = 1000
      real(dp)             :: x_theta(nbuffer),x_phi(nbuffer),   &
                              y_theta(nbuffer),y_phi(nbuffer),   &
                              z_theta(nbuffer),z_phi(nbuffer)

      if (verbose_satellite) then
         write(*,idf)  ID_sat,'  Reading satellite position'
         write(*,idfe) ID_sat,'  Start time =',starttime
         write(*,idfe) ID_sat,'  End time =',endtime
      endif

      if (noperiods.le.0) call read_pointing_periods(file_satellite)

      total_time = nosteps_satellite_file*sat_step

      starttime_read = starttime
      endtime_read = endtime

      if (endtime.lt.0) endtime_read=total_time

      if (endtime.gt.total_time) then
         if (verbose_satellite) &
             write(*,idf) ID_sat,'  Warning: Cutting end time to',total_time
         endtime_read = total_time
      endif

      if (starttime_read.lt.0.or.endtime_read.lt.starttime_read) then
         write(*,idf)  ID_sat,'ERROR in read_satellite_position:'
         write(*,idf)  ID_sat,'Time out of range'
         write(*,idfe) ID_sat,'Start/end time =',starttime_read,endtime_read
         write(*,idfe) ID_sat,'Allowed range  =',0.d0,total_time
         call exit_with_status(1)
      endif

      iperiod = 1
      do
         if (starttime_period(iperiod).gt.starttime_read) then ! between periods
            firststep_read = (starttime_read-sat_start)/sat_step+1
            starttime_read = (firststep_read-1)*sat_step+sat_start
            exit
         elseif (endtime_period(iperiod).gt.starttime_read) then
! Starttime inside this pointing period: Read the whole period
            firststep_read = firststep_period(iperiod)
            starttime_read = starttime_period(iperiod)
            exit
         else
            iperiod = iperiod+1

            if (iperiod==noperiods+1) then
               firststep_read = (starttime_read-sat_start)/sat_step+1
               starttime_read = (firststep_read-1)*sat_step+sat_start
               exit
            endif
         endif
      enddo

      if (firststep_read.le.0) then
         firststep_read = 1
         starttime_read = 0.0
      endif

      iperiod = 1
      do
         if (starttime_period(iperiod).gt.endtime_read) then
            laststep_read =(endtime_read-sat_start)/sat_step+2
            endtime_read = sat_start+(laststep_read-1)*sat_step
            exit
         elseif (endtime_period(iperiod).gt.endtime_read) then
            laststep_read = laststep_period(iperiod)
            endtime_read = endtime_period(iperiod)
            exit
         else
            iperiod = iperiod+1

            if (iperiod==noperiods+1) then
               laststep_read =(endtime_read-sat_start)/sat_step+2
               endtime_read = sat_start+(laststep_read-1)*sat_step
               exit
            endif
         endif
      enddo

      if (verbose_satellite) then
          write(*,idfe) ID_sat,'starttime_read',starttime_read
          write(*,idfe) ID_sat,'endtime_read  ',endtime_read
          write(*,idf)  ID_sat,'firststep_read',firststep_read
          write(*,idf)  ID_sat,'laststep_read ',laststep_read
      endif

      if (laststep_read.gt.nosteps_satellite_file) then
          laststep_read = nosteps_satellite_file
          endtime_read = total_time
      endif

      nosteps_read = laststep_read-firststep_read+1
      offset = firststep_read-1

      if (nosteps_read.gt.nosteps_allocated) then

         if (nosteps_allocated.gt.0) deallocate(satpos)

         nosteps_allocated = nosteps_read+10
! Allocate a few extra bytes to avoid unnecessary reallocation
         allocate(satpos(3,3,nosteps_allocated))

         memory_satellite = 72.*nosteps_allocated
      endif

      nosteps_left = nosteps_read

      call fits_open(inhandle,trim(file_satellite),fits_readonly,3)
 
      m = 0
      do
         n = min(nosteps_left,nbuffer)

         call fits_read_column(inhandle,1,x_theta(1:n),offset)
         call fits_read_column(inhandle,2,x_phi(1:n),offset)
         call fits_read_column(inhandle,3,y_theta(1:n),offset)
         call fits_read_column(inhandle,4,y_phi(1:n),offset)
         call fits_read_column(inhandle,5,z_theta(1:n),offset)
         call fits_read_column(inhandle,6,z_phi(1:n),offset)

         do i = 1,n
            m = m+1
            satpos(:,1,m) = get_cartesian(x_theta(i),x_phi(i))
            satpos(:,2,m) = get_cartesian(y_theta(i),y_phi(i))
            satpos(:,3,m) = get_cartesian(z_theta(i),z_phi(i))
         enddo

         nosteps_left = nosteps_left-n
         offset = offset+n

         if (nosteps_left==0) exit
      enddo

      if (present(trtype)) call convert_coord(satpos,nosteps_read,trtype)

      call fits_close(inhandle)

      if (verbose_satellite) write(*,idf) ID_sat,'  Done'

   END SUBROUTINE



!------------------------------------------------------------------------------


   SUBROUTINE read_pointing_periods(file_satellite)

      character(len=*)     :: file_satellite
      integer              :: iperiod
      real(dp)             :: starttime_int, starttime_frac
      type(fitshandle)     :: inhandle
      real(dp),allocatable :: dbuffer(:)

! Read and store starting/ending times of all pointing periods

      if (verbose_satellite) write(*,idf) ID_sat,'  Reading pointing periods'

      call fits_open(inhandle,trim(file_satellite),fits_readonly)

      noperiods = inhandle%nrows*inhandle%columns(1)%repcount

      if (allocated(starttime_period)) deallocate(starttime_period)
      if (allocated(endtime_period))   deallocate(endtime_period)

      allocate(starttime_period(noperiods))
      allocate(endtime_period(noperiods))
      allocate(firststep_period(noperiods))
      allocate(laststep_period(noperiods))
      allocate(dbuffer(noperiods))

      call fits_read_column(inhandle,1,dbuffer)
      starttime_int = dbuffer(1)
      starttime_period = dbuffer-starttime_int

      call fits_read_column(inhandle,2,dbuffer)
      starttime_frac = dbuffer(1)
      starttime_period = starttime_period +dbuffer-starttime_frac

      call fits_read_column(inhandle,13,dbuffer)
      endtime_period = dbuffer-starttime_int

      call fits_read_column(inhandle,14,dbuffer)
      endtime_period = endtime_period +dbuffer-starttime_frac

      deallocate(dbuffer)

      do iperiod = 1,noperiods
         firststep_period(iperiod) = &
             (starttime_period(iperiod)-sat_start)/sat_step+2
         laststep_period(iperiod) =  &
              (endtime_period(iperiod)-sat_start)/sat_step+1
      enddo

! Read the number of satellite pointing steps
      call fits_goto_hdu(inhandle,3)

      nosteps_satellite_file = inhandle%nrows*inhandle%columns(1)%repcount

      call fits_close(inhandle)

      if (verbose_satellite)   &
           write(*,idf) ID_sat,'  Done'

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE close_satellite()

      if (allocated(satpos))           deallocate(satpos)
      if (allocated(starttime_period)) deallocate(starttime_period)
      if (allocated(endtime_period))   deallocate(endtime_period)
      if (allocated(firststep_period)) deallocate(firststep_period)
      if (allocated(laststep_period))  deallocate(laststep_period)

      nosteps_allocated = -1
      nosteps_read = -1
      nosteps_satellite_file = -1
      noperiods = -1

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE detector_pointing_iqu(pixels,qw,uw,nosamples,start_time, &
                                    nside,scheme,shift_in)
! Compute detector pointings, transform into pixel number
!  and factors cos(2*psi),sin(2*psi), and store into tables PIXELS, QW, UW.
! NOSAMPLES samples, starting at START_TIME.

      integer,         intent(in)          :: nosamples, nside
      integer,         intent(out)         :: pixels(nosamples)
      real(sp),        intent(out)         :: qw(nosamples), uw(nosamples)
      real(dp),        intent(in)          :: start_time
      character(len=*),intent(in),optional :: scheme
      real(dp),        intent(in),optional :: shift_in

      integer  :: i, n, nosamples_left, ir, nosamples_done
      logical  :: do_nested
      real(dp) :: time, cf, cpsi, spsi
      real(dp) :: ephi(3), shift

      integer, parameter   :: nbuff=10000
      real(dp),allocatable :: detpos(:,:,:), detdir(:,:)
      logical, allocatable :: flag(:)

      if (present(scheme)) then
         if (scheme=='nested'.or.scheme=='NESTED') then
            do_nested = .true.
         elseif (scheme=='ring'.or.scheme=='RING') then
            do_nested = .false.
         else
            write(*,idf) ID_sat,'ERROR in detector_pointing_iqu:'
            write(*,idf) ID_sat,'Unknown scheme ',scheme
            call exit_with_status(1)
         endif
      else
         do_nested = .true.  !default
      endif

      shift = 0.0
      if (present(shift_in)) shift=shift_in

      allocate(detpos(3,3,nbuff))
      allocate(detdir(3,nbuff))
      allocate(flag(nbuff))

      nosamples_left = nosamples
      nosamples_done = 0
      time = start_time
      ir = 0
      do
         n = min(nosamples_left,nbuff)
         if (n==0) exit

         call detector_position(detpos,flag,n,time,shift)

         do i = 1,n
            detdir(:,i) = detpos(:,3,i)
         enddo

         if (do_nested) then
            call cartesian_to_nested(pixels(ir+1),detdir,n,nside)
         else
            call cartesian_to_ring(pixels(ir+1),detdir,n,nside)
         endif

         do i = 1,n
            ir = ir+1

            if (.not.flag(i)) then
               pixels(ir) = 12*nside*nside
               qw(ir) = 0.0
               uw(ir) = 0.0
            else

               ephi(1) = -detpos(2,3,i)
               ephi(2) =  detpos(1,3,i)
               ephi(3) =  0.0

               spsi = detpos(1,1,i)*ephi(1) +detpos(2,1,i)*ephi(2)
               cpsi = detpos(1,2,i)*ephi(1) +detpos(2,2,i)*ephi(2)

               cf =  1.d0/(cpsi*cpsi +spsi*spsi)
               qw(ir) = (cpsi*cpsi -spsi*spsi)*cf  ! cos(2*psi)
               uw(ir) = 2.d0*cpsi*spsi*cf      ! sin(2*psi)
            endif
         enddo

         nosamples_left = nosamples_left-n
         nosamples_done = nosamples_done+n

         time = start_time +nosamples_done/fsample_det
      enddo

      deallocate(detpos,detdir,flag)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE detector_pointing_pix(pixels,nosamples,start_time,nside,  &
                                    scheme,shift_in)

! Compute detector pointings, transform into pixel number,
! and store them into table PIXELS
! NOSAMPLES samples, start at START_TIME

      integer,         intent(in)          :: nosamples, nside
      integer,         intent(out)         :: pixels(nosamples)
      real(dp),        intent(in)          :: start_time
      character(len=*),intent(in),optional :: scheme
      real(dp),        intent(in),optional :: shift_in

      integer  :: i, n, nosamples_left, nosamples_done, ir
      logical  :: do_nested
      real(dp) :: time, shift

      integer, parameter   :: nbuff=10000
      real(dp),allocatable :: detpos(:,:,:), detdir(:,:)
      logical, allocatable :: flag(:)

      if (present(scheme)) then
         if (scheme=='nested'.or.scheme=='NESTED') then
            do_nested = .true.
         elseif (scheme=='ring'.or.scheme=='RING') then
            do_nested = .false.
         else
            write(*,idf) ID_sat,': ERROR in detector_pointing_pix:'
            write(*,idf) ID_sat,': Unknown scheme '//scheme
            call exit_with_status(1)
         endif
      else
         do_nested = .true.  !default
      endif

      shift = 0.0
      if (present(shift_in)) shift=shift_in

      allocate(detpos(3,3,nbuff))
      allocate(detdir(3,nbuff))
      allocate(flag(nbuff))

      nosamples_left = nosamples
      nosamples_done = 0
      time = start_time
      ir = 0
      do
         n = min(nosamples_left,nbuff)
         if (n==0) exit

         call detector_position(detpos,flag,n,time,shift)

         do i = 1,n
            detdir(:,i) = detpos(:,3,i)
         enddo

         if (do_nested) then
            call cartesian_to_nested(pixels(ir+1),detdir,n,nside)
         else
            call cartesian_to_ring(pixels(ir+1),detdir,n,nside)
         endif

         do i = 1,n
            ir = ir+1
            if (.not.flag(i)) pixels(ir)=12*nside*nside
         enddo

         nosamples_left = nosamples_left-n
         nosamples_done = nosamples_done+n

         time = start_time +nosamples_done/fsample_det

      enddo

      deallocate(detpos,detdir,flag)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE detector_position(detpos,flag,nosamples,start_time,shift_in)

! This routine computes the actual conversion.
! Compute detector pointings:
! NOSAMPLES samples starting at time START_TIME

      integer, intent(in)          :: nosamples
      real(dp),intent(out)         :: detpos(3,3,nosamples)
      logical, intent(out)         :: flag(nosamples)
      real(dp),intent(in)          :: start_time
      real(dp),intent(in),optional :: shift_in

      integer  :: k, iperiod, ip, ind, nosteps_period
      real(dp) :: time, subtime, period_time, invsat_step
      real(dp) :: period_duration, end_time
      real(dp) :: da, cosstep, sinstep, cosrot, sinrot, swap
      real(dp) :: axis(3), angle, dtime, per_start, shift
      logical  :: newstep, newperiod
      real(dp) :: tsatpos1(3,3), satpos2(3,3)
      real(dp) :: tsatrotc(3,3), tsatrots(3,3), satrot(3,3)
      real(dp) :: detposc(3,3), detposs(3,3), detpos0(3,3), unit3(3,3)

      if (fsample_det.le.0) then
         write(*,idf) ID_sat,'ERROR in satellite-detector conversion:'
         write(*,idf) ID_sat,'No detector selected.'
         call exit_with_status(1)
      endif

      if (nosteps_read.lt.0) then
         write(*,idf) ID_sat,'ERROR in satellite-detector conversion:'
         write(*,idf) ID_sat,'No satellite data read.'
         call exit_with_status(1)
      endif

      shift = 0.0
      if (present(shift_in)) shift=shift_in

      dtime = 1.d0/fsample_det
      invsat_step = 1.d0/sat_step

      end_time = start_time +(nosamples-1)*dtime

      if (start_time.lt.starttime_read) then
         write(*,idf)  ID_sat,'ERROR in satellite-detector conversion:'
         write(*,idf)  ID_sat,'Start time out of range.'
         write(*,idfe) ID_sat,'Start_time, starttime_read =',  &
                                 start_time, starttime_read
         call exit_with_status(1)
      elseif (end_time.gt.endtime_read) then
         write(*,idf)  ID_sat,'ERROR in satellite-detector conversion:'
         write(*,idf)  ID_sat,'End time out of range.'
         write(*,idfe) ID_sat,'End_time, endtime_read =',end_time,endtime_read
         call exit_with_status(1)
      endif

! Find first pointing period before the first sample
      iperiod = noperiods
      do
         if (endtime_period(iperiod).le.start_time) exit

         iperiod = iperiod-1
         if (iperiod==0) exit
      enddo

      unit3 = 0.0
      unit3(1,1) = 1.d0
      unit3(2,2) = 1.d0
      unit3(3,3) = 1.d0

      newperiod = .true.

      time = start_time

      do k = 1,nosamples

         if (newperiod) then

! Between periods: no data
            if (iperiod==noperiods) then
               flag(k)       = .false.
               detpos(:,:,k) = unit3
               time          = time+dtime
               cycle

            elseif (starttime_period(iperiod+1).gt.time) then
               flag(k)       = .false.
               detpos(:,:,k) = unit3
               time          = time+dtime
               cycle
            endif

! Entering a new pointing period
            iperiod = iperiod+1

            period_duration = endtime_period(iperiod)-starttime_period(iperiod)
            nosteps_period = laststep_period(iperiod)  &
                             -firststep_period(iperiod)+1

! Time from beginning of the pointing period
! to the first satellite step inside it
            per_start = (firststep_period(iperiod)-1)*sat_step +sat_start  &
                          +shift/fsample_det  -starttime_period(iperiod)

!Location inside the pointing period
            period_time = time-starttime_period(iperiod)
            ip = (period_time-per_start)/sat_step+1

            if (ip.ge.nosteps_period) then
               ip = nosteps_period-1
            elseif (ip.le.0) then
               ip = 1
            endif

            ind = ip +firststep_period(iperiod)-firststep_read

            subtime = period_time -per_start-(ip-1)*sat_step

            if (ind.lt.1.or.ind.ge.nosteps_read) then
               write(*,idf) ID_sat,'ERROR in satellite-detector conversion'
               write(*,idf) ID_sat,'Index exceeds allowed range'
               write(*,idf) ID_sat,'IND, nosteps_read =',ind,nosteps_read
               call exit_with_status(1)
            endif

            satpos2 = satpos(:,:,ind)

            newperiod = .false.
            newstep = .true.

         endif

         if (newstep) then

            tsatpos1 = transpose(satpos2)
            satpos2 = satpos(:,:,ind+1)

            satrot = matmul(satpos2,tsatpos1)
            call get_axis(axis,angle,satrot)

            call satrotation(tsatrots,tsatrotc,axis)

            detpos0 = matmul3(tsatpos1,detrot)
            detposs = matmul3(tsatrots,detpos0)
            detposc = matmul3(tsatrotc,detpos0)

            da = angle*subtime*invsat_step
            cosrot = cos(da)
            sinrot = sin(da)

            da = angle*dtime*invsat_step
            cosstep = cos(da)
            sinstep = sin(da)

            newstep = .false.

         else

            swap   = cosrot*cosstep-sinrot*sinstep
            sinrot = sinrot*cosstep+cosrot*sinstep
            cosrot = swap
         endif

         detpos(:,:,k) = detpos0 +sinrot*detposs +(1.d0-cosrot)*detposc
         flag(k) = .true.

!avoid invrementing TIME
         period_time = period_time+dtime
         subtime     = subtime+dtime

         if (period_time.ge.period_duration) then
            time = starttime_period(iperiod)+period_time
            newperiod = .true.

         elseif (subtime.ge.sat_step) then
            do
               if (ip==nosteps_period-1) exit

               subtime = subtime-sat_step
               ip = ip+1
               ind = ind+1
               newstep = .true.

               if (subtime.lt.sat_step) exit

            enddo
         endif

      enddo

   CONTAINS

      FUNCTION matmul3(tmat1,mat2) result(mat3)
!   mat3 = transpose(tmat1)*mat2
!
         real(dp) :: tmat1(3,3),mat2(3,3),mat3(3,3)
         integer  :: j

         do j = 1,3
            mat3(1,j) = sum(tmat1(:,1)*mat2(:,j))
            mat3(2,j) = sum(tmat1(:,2)*mat2(:,j))
            mat3(3,j) = sum(tmat1(:,3)*mat2(:,j))
         enddo

      END FUNCTION

   END SUBROUTINE


!------------------------------------------------------------------------------
!
! Internal routines
!
!------------------------------------------------------------------------------


   FUNCTION get_cartesian(theta,phi) result(coord)
!
      real(dp) :: coord(3)
      real(dp) :: phi,theta ! rad

      coord(1) = sin(theta)*cos(phi)
      coord(2) = sin(theta)*sin(phi)
      coord(3) = cos(theta)

   END FUNCTION


!------------------------------------------------------------------------------


   FUNCTION detector_rotation(theta,phi,psi,theta_b) result(mat)
!
! Find rotation matrix that rotates to detector position on the focal plane
!
      real(dp) :: mat(3,3)
      real(dp) :: phi,theta,psi,theta_b   ! rad
      real(dp) :: axis(3), m(3,3)

! Rotation around y-axis by 90-theta_b
      axis(1) = 0.0
      axis(2) = 1.d0
      axis(3) = 0.0
      mat = rotmatrix(axis,.5d0*pi-theta_b)

! Rotation around an axis in the xy-plane by theta_uv
      axis(1) = -sin(phi)
      axis(2) =  cos(phi)
      axis(3) = 0.0
      m = rotmatrix(axis,theta)
      mat = matmul(mat,m)

! Rotation aroud z-axis by psi_uv
      axis(1) = 0.0
      axis(2) = 0.0
      axis(3) = 1.d0
      m = rotmatrix(axis,psi)
      mat = matmul(mat,m)

   END FUNCTION


!------------------------------------------------------------------------------


   FUNCTION rotmatrix(axis,angle) result(mat)
!
! Construct rotation matrix
! Rotation by theta around axis
!
      real(dp)  :: mat(3,3)
      real(dp)  :: axis(3),angle
      real(dp)  :: ca, sa

      ca = cos(angle)
      sa = sin(angle)

      mat(1,1) = (1.d0-ca)*axis(1)*axis(1) +ca
      mat(2,1) = (1.d0-ca)*axis(1)*axis(2) +sa*axis(3)
      mat(3,1) = (1.d0-ca)*axis(1)*axis(3) -sa*axis(2)

      mat(1,2) = (1.d0-ca)*axis(2)*axis(1) -sa*axis(3)
      mat(2,2) = (1.d0-ca)*axis(2)*axis(2) +ca
      mat(3,2) = (1.d0-ca)*axis(2)*axis(3) +sa*axis(1)

      mat(1,3) = (1.d0-ca)*axis(3)*axis(1) +sa*axis(2)
      mat(2,3) = (1.d0-ca)*axis(3)*axis(2) -sa*axis(1)
      mat(3,3) = (1.d0-ca)*axis(3)*axis(3) +ca

   END FUNCTION


!------------------------------------------------------------------------------


   SUBROUTINE satrotation(tsatrots,tsatrotc,axis)
!
! Rotation around axis by angle phi:
! transpose(Rotation matrix) = 1+sin(phi)*tsatrots+(1-cos(phi))*tsatrotc
!

      real(dp),intent(out) :: tsatrots(3,3)
      real(dp),intent(out) :: tsatrotc(3,3)
      real(dp),intent(in)  :: axis(3)

      tsatrots(1,1) =  0.0
      tsatrots(2,1) = -axis(3)
      tsatrots(3,1) =  axis(2)

      tsatrots(1,2) =  axis(3)
      tsatrots(2,2) =  0.0
      tsatrots(3,2) = -axis(1)

      tsatrots(1,3) = -axis(2)
      tsatrots(2,3) =  axis(1)
      tsatrots(3,3) =  0.0

      tsatrotc(1,1) = axis(1)*axis(1) -1.d0
      tsatrotc(2,1) = axis(2)*axis(1)
      tsatrotc(3,1) = axis(3)*axis(1)

      tsatrotc(1,2) = axis(1)*axis(2)
      tsatrotc(2,2) = axis(2)*axis(2) -1.d0
      tsatrotc(3,2) = axis(3)*axis(2)

      tsatrotc(1,3) = axis(1)*axis(3)
      tsatrotc(2,3) = axis(2)*axis(3)
      tsatrotc(3,3) = axis(3)*axis(3) -1.d0

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE get_axis(axis,angle,mat)
!
! Find axis and rotation angle corresponding to rotation matrix mat
! Angle is limited in range [0,pi[
      real(dp),intent(out) :: axis(3),angle
      real(dp),intent(in)  :: mat(3,3)
      real(dp)             :: ca, sa, tr, cf

      tr = mat(1,1)+mat(2,2)+mat(3,3)
      ca = 0.5d0*(tr-1.d0)
      angle = acos(ca)

      if (angle.gt.1.e-20) then
         sa = sin(angle)
         cf = 0.5d0/sa
         axis(1) = (mat(3,2)-mat(2,3))*cf
         axis(2) = (mat(1,3)-mat(3,1))*cf
         axis(3) = (mat(2,1)-mat(1,2))*cf
      else ! angle=0
         axis(1) = 0.0
         axis(2) = 0.0
         axis(3) = 1.d0
         angle = 0.0
      endif
! I expect case angle = pi will never happen
   END SUBROUTINE


!------------------------------------------------------------------------------

END MODULE satellite
