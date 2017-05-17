MODULE commonparam
  ! Common parameters and definitions

  ! TOAST additions -RK
  !use toast_mpi
  !use toast
  use iso_c_binding
  ! TOAST additions end -RK

  use planck_config
  use simulation, only : detector_data, pointing_data, tod_component

  implicit none
  public

  integer, parameter :: SLEN  = 150

  integer, parameter :: idp = i8b

  real(dp), parameter :: deg2rad = pi/180.d0

  integer, parameter :: UNIT_POINTING_PERIOD=-5

  ! OpenMP
  integer :: nthreads_max=1, nthreads=1, id_thread=0
  integer, external :: omp_get_num_procs, omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

  ! TOAST additions -RK
  integer(dp), allocatable :: chunk_offsets(:), chunk_sizes(:)
  integer(i4b) :: n_chunkset, n_chunkset_tot, n_channel, n_chunk, n_chunk_tot
  integer(i4b) :: n_chunk_max, n_sky_nonzero
  integer(i4b) :: first_chunk, last_chunk ! starts from 1
  logical, allocatable :: use_toast_channel(:)
  character(len=256) :: telescope_name, sky_name
  logical :: flag_by_horn=.false., force_pol=.false.
  logical :: concatenate_messages = .true.
  logical :: reassign_submaps = .true.
  logical :: noise_weights_from_psd = .false. ! integrate noise weights internally
  logical :: radiometers = .true. ! Assume well-behaved noise spectrum without low pass filtering
  integer(i8b) :: read_buffer_len = 1e6, psd_downsample=10
  integer (i8b) :: psdlen=1e6
  ! Enable sub ring map making
  integer(i2b) :: nsubchunk=0
  integer(i2b) :: isubchunk=0
  real(dp) :: fnoise_max=1000 ! When measuring noise variance, use this limit
  character(len=SLEN) :: file_profile = ''
  character(len=SLEN) :: file_intermediate_profile = ''
  logical :: checknan=.false. ! Can cost time
  logical :: sync_output=.true., skip_existing=.false.
  logical :: write_cut=.false.
  logical :: filter_mean=.false.
  logical :: tod_is_clean = .false.
  logical :: binary_output=.false., concatenate_binary=.false.
  integer :: record_number=1 ! Used for concatenate_binary when storing multiple MC maps
  integer(i4b) :: nwrite_binary = 10 ! Number of independent groups of processes writing binary maps
  integer(i4b), parameter :: basis_poly=1, basis_fourier=2, basis_cheby=3, basis_legendre=4
  integer(i4b) :: basis_func=basis_legendre, basis_order=0
  type :: basis_function_type
     integer(i8b) :: nsamp
     logical :: copy
     real(dp), pointer :: arr(:,:)
  end type basis_function_type
  type(basis_function_type), allocatable :: basis_functions(:)
  real(dp), pointer :: basis_function(:,:)

  integer, parameter :: NDETMAX=1000
  integer, parameter :: NDETSETMAX=1000
  type detset_type
     character(len=SLEN) :: name
     character(len=SLEN) :: detectors(NDETMAX)
     integer(i4b) :: ndet
     logical :: nopol
  end type detset_type
  type(detset_type) :: detsets(0:NDETSETMAX)
  integer(i4b) :: ndetset
  logical :: detflags(NDETMAX)

  integer, parameter :: NSPANMAX=1000
  integer, parameter :: NSURVEYMAX=1000
  type survey_type
     character(len=SLEN) :: name
     real(dp) :: starts(NSPANMAX), stops(NSPANMAX)
     integer(i4b) :: nspan
  end type survey_type  
  type(survey_type) :: surveys(0:NSURVEYMAX)
  integer(i4b) :: nsurvey
  logical, allocatable :: surveyflags(:)

  logical :: bin_subsets = .true.
  logical :: mcmode = .false., cached = .false.

  ! TOAST additions end -RK

  ! other additions -RK
  integer :: nread_concurrent=-1
  real(dp) :: good_baseline_fraction=0 ! default acceps all baselines
  ! monte Carlo mode
  integer(idp) :: mc_increment=1e7, mc_loops=1, mc_id=0, rng_base=0
  logical :: run_submap_test=.true.
  logical :: incomplete_matrices = .false.
  ! other additions end -RK

  integer :: ID = 0
  integer :: ntasks = 1

  character(len=*), parameter :: version = '3.7'
  character(len=*), parameter :: idf = '(i3,": ",a,1x,i0,2x,i0)'

  ! Input parameters
  integer :: info=2

  integer :: nside_map=512, nside_cross=-1, nside_submap=16

  real(dp) :: time_unit=-5.d0
  real(dp) :: mission_time=0.0
  integer :: nofiles=-1

  integer :: pixmode_map=2, pixmode_cross=2
  real(dp) :: pixlim_map=1e-6, pixlim_cross=1e-3

  real(dp) :: dnshort=-1
  integer :: nlong=-1, nshort=-1
  logical :: kfirst=.true., kfilter=.false.
  real(dp) :: filter_time=3600, tail_time=600

  real(dp) :: cglimit=1.d-12
  integer :: iter_min=3, iter_max=1000
  integer :: precond_width=10

  integer :: mode_detweight=0

  logical :: rm_monopole=.false., temperature_only=.false.       

  ! Input files
  character(len=SLEN) :: file_param='', file_simulation='', &
       file_inmask='', file_spectrum='', file_gap='',  &
       file_pntperiod='', file_objectsize='', &
       file_fpdb_supplement=''

  ! Output files
  character(len=SLEN) :: file_root='madam'
  character(len=SLEN) :: file_map='', file_hit='', file_base=''
  character(len=SLEN) :: file_matrix='', file_mask='', file_binmap=''
  character(len=SLEN) :: file_wcov='', file_leakmatrix=''
  character(len=SLEN) :: file_gap_out='', file_mc='',  path_output=''

  logical :: write_tod=.false.

  ! LFI specific keywords
  character(len=80) :: instrument  = ''

  ! NCVM specific parameters
  logical :: kwrite_covmat=.false., ksplit_covmat=.false., bfinvert=.false.
  logical :: kcompress_pixels=.false.
  character(len=SLEN) :: file_covmat = ''

  type(detector_data), allocatable :: detectors(:)

  ! Derived directly from input parameters
  !integer            :: nostokes, ncc, nside_max
  integer            :: nmap=0, ncc=0, nside_max
  integer            :: nodetectors=-1

  ! Pixels
  integer            :: nopix_map, nopix_cross
  integer            :: nosubpix_map, nosubpix_cross, nosubpix_max
  integer            :: nosubmaps_tot, nosubmaps, nosubmaps_max
  integer            :: nolocmaps, nolocpix

  ! Number of samples
  integer(idp)       :: nosamples_tot
  integer(idp)       :: nosamples_proc_max, nosamples_proc

  integer(idp)       :: istart_mission
  integer(idp)       :: istart_proc

  ! Baselines
  integer(idp)        :: noba_short_tot
  integer(idp)        :: noba_short_max

  integer(idp)        :: noba_short

  integer(idp)        :: kshort_start

  integer, allocatable :: baselines_short(:) ! short baselines per process
  integer, allocatable :: baselines_short_start(:)
  integer, allocatable :: baselines_short_stop(:)
  real(dp), allocatable :: baselines_short_time(:)
  integer, allocatable :: base_pntid_short(:)

  ! Number of pointing periods and their duration as a number of samples
  integer(i4b)              :: nopntperiods = -1
  integer(i8b), allocatable :: pntperiods(:)
  integer(i4b), allocatable :: pntperiod_id(:)

  ! Baselines per pointing period (only used in case use_pntperiods=T)
  integer,allocatable :: noba_short_pp(:)

  integer,allocatable :: id_submap(:)
  integer             :: id_next, id_prev

  logical             :: do_map=.true., do_binmap=.false., do_hits=.false.
  logical             :: do_mask=.false., do_matrix=.false., do_wcov=.false.
  logical             :: do_base=.false., do_leakmatrix=.false.
  logical             :: do_wnmap=.false., do_dethits=.false.
  logical             :: use_inmask
  logical             :: kinputcheck = .false.

  integer             :: noiter = 0

  !--------------------------------------------------------------------------

END MODULE commonparam
