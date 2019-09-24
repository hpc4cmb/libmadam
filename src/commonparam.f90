MODULE commonparam
  ! Common parameters and definitions

  use iso_c_binding

  use planck_config

  implicit none
  public

  integer, parameter :: SLEN  = 150

  integer, parameter :: idp = i8b
  integer, parameter :: byte = 1 ! 1-byte logical kind

  real(dp), parameter :: deg2rad = pi/180.d0

  integer, parameter :: UNIT_POINTING_PERIOD=-5

  ! From simulation.f90

  TYPE detector_data
     integer :: idet = 0
     real(dp) :: fknee = 0
     character(len=20) :: name = ''
     integer :: npsd
     real(dp), allocatable :: psdstarts(:)
     real(dp), allocatable :: psdfreqs(:)
     real(dp), allocatable :: psds(:,:)
     real(dp), allocatable :: sigmas(:) ! "white" noise sigma
     real(dp), allocatable :: weights(:) ! detector noise weight
     ! plateau is a constant offset subtracted from the PSD to get
     ! the 1/f (baseline) part
     real(dp), allocatable :: plateaus(:)
  END TYPE detector_data

  real(dp) :: fsample ! Sampling frequency
  character(len=30) :: unit_tod

  ! OpenMP
  integer :: nthreads_max, nthreads, id_thread
  integer, external :: omp_get_num_procs, omp_get_max_threads, &
       omp_get_thread_num, omp_get_num_threads

  logical :: flag_by_horn, force_pol
  logical :: concatenate_messages, allreduce, reassign_submaps
  logical :: noise_weights_from_psd ! integrate noise weights internally
  ! Assume well-behaved noise spectrum without low pass filtering
  logical :: radiometers
  integer(i8b) :: psd_downsample, psdlen
  ! Enable sub ring map making
  integer(i2b) :: nsubchunk, isubchunk
  real(dp) :: fnoise_max
  character(len=SLEN) :: file_profile = ''
  character(len=SLEN) :: file_intermediate_profile = ''
  logical :: checknan
  real(dp) :: diagfilter
  logical :: sync_output, skip_existing
  logical :: write_cut
  logical :: tod_is_clean
  logical :: binary_output, concatenate_binary
  ! Used for concatenate_binary when storing multiple MC maps
  integer :: record_number
  ! Number of independent groups of processes writing binary maps
  integer(i4b) :: nwrite_binary
  integer(i4b), parameter :: basis_poly=1, basis_fourier=2, basis_cheby=3, &
       basis_legendre=4
  integer(i4b) :: basis_func, basis_order
  type :: basis_function_type
     integer(i8b) :: nsamp
     logical :: copy
     real(dp), pointer :: arr(:, :)
  end type basis_function_type
  type(basis_function_type), allocatable :: basis_functions(:)
  real(dp), pointer :: basis_function(:, :)

  integer, parameter :: NDETMAX=10000
  integer, parameter :: NDETSETMAX=1000
  type detset_type
     character(len=SLEN) :: name
     character(len=SLEN), allocatable :: detectors(:)
     integer(i4b) :: ndet
     logical :: nopol
  end type detset_type
  type(detset_type) :: detsets(0:NDETSETMAX)
  integer(i4b) :: ndetset
  logical, allocatable :: detflags(:)

  integer, parameter :: NSPANMAX=1000
  integer, parameter :: NSURVEYMAX=1000
  type survey_type
     character(len=SLEN) :: name
     real(dp) :: starts(NSPANMAX), stops(NSPANMAX)
     integer(i4b) :: nspan
  end type survey_type
  type(survey_type) :: surveys(0:NSURVEYMAX)
  integer(i4b) :: nsurvey
  logical(byte), allocatable :: surveyflags(:)

  logical :: bin_subsets
  logical :: mcmode, cached

  real(dp) :: good_baseline_fraction
  ! monte Carlo mode
  integer(idp) :: mc_increment, mc_loops, mc_id, rng_base
  logical :: incomplete_matrices

  integer :: ID
  integer :: ntasks

  character(len=*), parameter :: version = '3.7'
  character(len=*), parameter :: idf = '(i4,": ",a,1x,i0,2x,i0)'

  ! Input parameters
  integer :: info

  integer :: nside_map, nside_cross, nside_submap

  real(dp) :: time_unit
  real(dp) :: mission_time
  integer :: nofiles

  integer :: pixmode_map, pixmode_cross
  real(dp) :: pixlim_map, pixlim_cross
  logical :: allow_decoupling

  real(dp) :: dnshort
  integer :: nlong, nshort
  logical :: kfirst, kfilter
  logical :: unaligned_fft

  real(dp) :: cglimit
  integer :: iter_min, iter_max
  integer :: precond_width_min, precond_width_max
  logical :: use_fprecond, use_cgprecond

  integer :: mode_detweight

  logical :: rm_monopole, temperature_only

  ! Input files
  character(len=SLEN) :: file_param, file_inmask, file_spectrum, file_gap

  ! Output files
  character(len=SLEN) :: file_root, file_map, file_hit, file_base, &
       file_matrix, file_mask, file_binmap, file_wcov, file_leakmatrix, &
       file_gap_out, file_mc, path_output

  logical :: write_tod

  ! LFI specific keywords
  character(len=80) :: instrument

  ! NCVM specific parameters
  logical :: kwrite_covmat
  character(len=SLEN) :: file_covmat

  type(detector_data), allocatable, target :: detectors(:)

  ! Derived directly from input parameters
  integer :: nmap, ncc, nside_max, nodetectors

  ! Pixels
  integer :: nopix_map, nopix_cross
  integer :: nosubpix_map, nosubpix_cross, nosubpix_max
  integer :: nosubmaps_tot, nosubmaps, nosubmaps_max
  integer :: nolocmaps, nolocpix

  ! Number of samples
  integer(i8b) :: nosamples_tot, nosamples_proc_max, nosamples_proc
  integer(i8b) :: istart_mission, istart_proc

  ! Baselines
  integer(i8b) :: noba_short_tot, noba_short_max, noba_short, noba_short_pp_max
  integer(i8b) :: kshort_start

  integer(i4b), allocatable :: baselines_short(:) ! short baselines per process
  integer(i4b), allocatable :: baselines_short_start(:)
  integer(i4b), allocatable :: baselines_short_stop(:)
  real(i8b), allocatable :: baselines_short_time(:)

  ! Number of pointing periods and their duration as a number of samples
  integer(i4b) :: ninterval, ninterval_tot
  integer(i8b), allocatable :: intervals(:)
  integer(i4b), allocatable :: interval_id(:)

  ! Baselines per pointing period (only used in case use_pntperiods=T)
  integer(i4b), allocatable :: noba_short_pp(:)

  integer(i4b), allocatable :: id_submap(:)
  integer(i4b) :: id_next, id_prev

  logical :: do_map, do_binmap, do_hits, do_mask, do_matrix, do_wcov, do_base, &
       do_leakmatrix, do_wnmap, use_inmask
  integer :: noiter

  !--------------------------------------------------------------------------

contains

  subroutine set_parameter_defaults()
    ! Set or restore commonparam values to defaults
    fsample = 1
    unit_tod = "unknown"
    nthreads_max = 1
    nthreads = 1
    id_thread = 1
    flag_by_horn = .false.
    force_pol = .false.
    concatenate_messages = .true.
    allreduce = .false.
    reassign_submaps = .true.
    noise_weights_from_psd = .false.
    radiometers = .true.
    psd_downsample = 10
    psdlen = 1e6
    nsubchunk = 0
    isubchunk = 0
    fnoise_max = 1000  ! When measuring noise variance, use this limit
    checknan = .false.  ! Can cost time
    diagfilter = 0
    sync_output = .true.
    skip_existing = .false.
    write_cut = .false.
    tod_is_clean = .false.
    binary_output = .false.
    concatenate_binary = .false.
    record_number = 1
    nwrite_binary = 10
    basis_func = basis_legendre
    basis_order = 0
    bin_subsets = .false.
    mcmode = .false.
    cached = .false.
    good_baseline_fraction = 0  ! default accepts all baselines
    mc_increment = 1e7
    mc_loops = 1
    mc_id = 0
    rng_base = 0
    incomplete_matrices = .false.
    ID = 0
    ntasks = 1
    info = 2
    nside_map = 512
    nside_cross = -1
    nside_submap = 16
    time_unit = -5
    mission_time = 0
    nofiles = -1
    pixmode_map = 2
    pixmode_cross = 2
    pixlim_map = 1e-6
    pixlim_cross = 1e-3
    allow_decoupling = .false.
    dnshort = -1
    nlong = -1
    nshort = -1
    kfirst = .true.
    kfilter = .false.
    unaligned_fft = .false.
    cglimit = 1e-12
    iter_min = 3
    iter_max = 1000
    precond_width_min = 10
    precond_width_max = 100
    use_fprecond = .false.
    use_cgprecond = .false.
    mode_detweight = 0
    rm_monopole = .false.
    temperature_only = .false.
    file_param = ""
    file_inmask = ""
    file_spectrum = ""
    file_gap = ""
    file_root = "madam"
    file_map = ""
    file_hit = ""
    file_base = ""
    file_matrix = ""
    file_mask = ""
    file_binmap = ""
    file_wcov = ""
    file_leakmatrix = ""
    file_gap_out = ""
    file_mc  = ""
    path_output = ""
    write_tod = .false.
    instrument  = ""
    kwrite_covmat = .false.
    file_covmat = ""
    nmap = 0
    ncc = 0
    nodetectors = -1
    ninterval = -1
    ninterval_tot = -1
    do_map = .true.
    do_binmap = .false.
    do_hits = .false.
    do_mask = .false.
    do_matrix = .false.
    do_wcov = .false.
    do_base = .false.
    do_leakmatrix = .false.
    do_wnmap = .false.
    noiter = 0
  end subroutine set_parameter_defaults

END MODULE commonparam
