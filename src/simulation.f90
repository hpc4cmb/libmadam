MODULE simulation
  !
  ! Routines for reading a simulation file and storing the information in it.
  ! This module may be used by a number of tools, therefore must not use 
  !   other Madam modules.

  implicit none
  public

  integer, parameter, private  :: SLEN  = 150
  integer, parameter, private  :: sp = kind(1.0)
  integer, parameter, private  :: dp = kind(1.0d0)
  integer, parameter, private  :: idp = selected_int_kind(18)

  real(dp),parameter, private  :: pi = 3.14159265358979323846264d0
  real(dp),parameter, private  :: deg2rad = pi/180.d0
  integer, parameter, private  :: MAXNODET = 99

  integer, private   :: ID = 0

  TYPE detector_data
     integer           :: idet   =  0
     integer           :: ipoint =  0
     real(dp)          :: slope  = .0
     real(dp)          :: fknee  = .0
     real(dp)          :: fmin   = .0
     !real(dp)          :: sigma  = .0 replaced by sigmas
     real(dp)          :: psipol = .0
     !real(dp)          :: weight = .0 replaced by weights
     character(len=20) :: name   = ''
     logical           :: kpolar = .false.
     !type(toast_psd), allocatable :: psds(:) ! -RK
     integer :: npsd
     real(dp), allocatable :: psdstarts(:)
     real(dp), allocatable :: psdfreqs(:)
     real(dp), allocatable :: psds(:,:)
     real(dp), allocatable :: sigmas(:), weights(:) ! -RK
     !real(dp), allocatable :: psd_freq(:), psd_val(:) ! -RK
  END TYPE detector_data

  TYPE pointing_data
     integer           :: ipoint =  0
     real(dp)          :: theta  = .0
     real(dp)          :: phi    = .0
     real(dp)          :: psi    = .0
     real(dp)          :: psipol = .0
  END TYPE pointing_data

  TYPE tod_component
     real(dp)          :: weight = 1.d0
     character(len=20) :: name   = ''
     character(len=20) :: key    = ''
  END TYPE tod_component

  !Parameters read directly from the simulation file

  type(detector_data),allocatable :: detectors_simu(:)
  type(pointing_data),allocatable :: pointings_simu(:)
  type(tod_component),allocatable :: tods_simu(:)

  character(len=SLEN),allocatable :: file_tod_simu(:,:,:)
  character(len=SLEN),allocatable :: file_point_simu(:,:)

  character(len=SLEN),allocatable, private :: path_tod(:)
  character(len=SLEN),             private :: path_point  = ''

  integer :: nodetectors_simu = 0   ! Number of detectors
  integer :: nopointings_simu = 0   ! Number of independent pointings (=horns)
  integer :: notods_simu      = 0   ! Number of TOD components (cmb,noise...)
  integer :: nofiles_simu     = 0   ! Number of files/objects per timeline

  ! If file naming follows the LevelS convention, it is enough to provide
  !  a name template with the hour label set to zero, 
  ! The hour label of Kth file will be label_start+(k-1)*label_step.
  ! As default all file names are listed in the simulation file
  integer :: label_start = -1   !default: list all files
  integer :: label_step  =  1

  real(dp) :: fsample      = 1.0   ! Sampling frequency
  real(dp) :: fobservation = 0.0   ! Observation frequency, needed for Tcmb-Tantenna conversion

  character(len=68)   :: simcomment = ''
  character(len=30)   :: unit_tod = 'unknown'
  character(len=SLEN) :: file_satellite = ''   ! Satellite pointing file

  ! The following parameters are nor read directly from the simulation file.
  ! They are defined here, but filled by the read_data module

  ! File/object sizes as a number of samples
  integer(idp),allocatable :: nosamples_file(:)          ! Samples/file

END MODULE simulation
