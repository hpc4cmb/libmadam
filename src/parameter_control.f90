MODULE parameter_control

  use iso_c_binding

  use commonparam
  use mpi_wrappers
  use noise_routines, only : interpolate_psd, measure_noise_weights
  use pointing, only : subchunk
  use maps_and_baselines, only : memory_baselines, memory_maps

  implicit none
  private

  public init_parameters, init_parallelization, start_timeloop, &
       write_parameters, baseline_times

  real(sp), save, public :: memory_basis_functions = 0
  character(len=40), parameter :: mstr='(x,a,t32,f9.1," MB")'
  character(len=40), parameter :: mstr3='(x,a,t32,3(f9.1," MB"))'

CONTAINS

  !-------------------------------------------------------------------------


  SUBROUTINE init_parameters()

    ! extensively rewritten for TOAST: use_pntperiods is now
    ! enabled by default and looping over pointing, time and buffer
    ! are disabled

    integer(i4b) :: idet, ndet, i, m, n, nleft, ns, ierr, idet1, idet2, &
         horn1, horn2, ipsd
    real(i8b) :: weight
    logical :: use_all_data
    character(len=SLEN) :: subsetname
    real(dp), allocatable :: weights(:)
    real(dp) :: psdmin, psdmax

    real(dp), pointer :: psd(:)

    ! Make sure path_output contains the separator

    n = len_trim(path_output)
    if (n /= 0) then
       if (path_output(n:n) /= '/') path_output = path_output(1:n) // '/'
    end if

    ! expand the file names

    subsetname = trim(file_root)
    if (len_trim(detsets(0)%name) /= 0) &
         subsetname = trim(subsetname) // '_' // trim(detsets(0)%name)
    if (len_trim(surveys(0)%name) /= 0) &
         subsetname = trim(subsetname) // '_' // trim(surveys(0)%name)

    file_map = ''
    file_binmap = ''
    file_hit = ''
    file_matrix = ''
    file_leakmatrix = ''
    file_wcov = ''
    file_base = ''
    file_mask = ''

    if (do_map) file_map = trim(subsetname) // '_map.fits'
    if (do_binmap) file_binmap = trim(subsetname) // '_bmap.fits'
    if (do_hits) file_hit = trim(subsetname) // '_hmap.fits'
    if (do_matrix) file_matrix = trim(subsetname) // '_wcov_inv.fits'
    if (do_leakmatrix) file_leakmatrix = trim(subsetname) // '_leakmatrix'
    if (do_wcov) file_wcov = trim(subsetname) // '_wcov.fits'
    if (do_base) file_base = trim(subsetname) // '_base.fits'
    if (do_mask) file_mask = trim(subsetname) // '_mask.fits'

    ! in a Monte Carlo, the binary maps are appended into a single file
    ! rather than storing each realization separately

    if (binary_output .and. (mc_loops > 1 .or. mc_id > 0)) then
       concatenate_binary = .true.
    else
       concatenate_binary = .false.
    end if

    ! Set common parameters

    if (id == 0 .and. info > 0) then
       write(*,*)
       write(*,*) 'Initializing parameters'
    end if

    if (.not. radiometers) then

       ! Regularize the provided PSDs by adding power to the
       ! low frequency modes

       do idet = 1, nodetectors
          do ipsd = 1, detectors(idet)%npsd
             psd => detectors(idet)%psds(:, ipsd)
             i = maxloc(psd, 1)
             psdmax = psd(i)
             psd(1:i-1) = psdmax
          end do
       end do

    end if

    if (noise_weights_from_psd) then
       ! Improve noise weighting in case we have full PSDs from TOAST
       if (id == 0 .and. info > 0) then
          write (*,*) 'Adjusting noise weights using noise spectra '
       end if
       call measure_noise_weights(radiometers)
    end if

    ! Checkings

    if (bin_subsets .or. do_leakmatrix) then
       ! These options assume that all processes have the same detectors
       call detector_check()
    end if

    if (nsubchunk < 2) then
       nsubchunk = 0
       isubchunk = 0
    end if

    if (isubchunk > nsubchunk) call abort_mpi('ERROR: isubchunk > nsubchunk')

    if (len_trim(file_covmat) > 0) kwrite_covmat = .true.
    if (kwrite_covmat) then
       kfilter = .true.
       kfirst = .true.
    end if

    if (nthreads > nthreads_max) nthreads = nthreads_max

    if (good_baseline_fraction < 0 .or. good_baseline_fraction > 1) &
         call abort_mpi('ERROR: good_baseline_fraction must be in 0..1')

    if (nside_cross < 0) nside_cross = nside_map

    call nside_check(nside_map, 'nside_map')
    call nside_check(nside_submap, 'nside_submap')
    call nside_check(nside_cross, 'nside_cross')

    if (nside_map < nside_cross) &
         call abort_mpi('nside_cross > nside_map not supported by Madam/TOAST')

    if (pixlim_cross < 0) pixlim_cross = pixlim_map

    if (pixmode_map < 0 .or. pixmode_map > 2) then
       if (id == 0) write(*,*) 'ERROR in input parameters:',  &
            ' Pixmode_map must be in range [0-2]'
       call exit_with_status(1)
    elseif (pixmode_cross < 0 .or. pixmode_cross > 4) then
       if (id == 0) write(*,*) 'ERROR in input parameters:',  &
            ' Pixmode_cross must be in range [0-4]'
       call exit_with_status(1)
    end if

    use_inmask = (len_trim(file_inmask).gt.0)

    if (nmap == 1) then
       if (id == 0 .and. info > 0) then
          if (nmap == 1) &
               write(*,*) 'Only unpolarized detector weights.'
          write(*,*) '  Will produce only temperature map'
          if (.not. temperature_only) &
               write(*,*) '  Setting temperature_only = T'
       end if

       temperature_only = .true.
       nmap = 1
    elseif (temperature_only) then
       if (id == 0 .and. info > 0) then
          write(*,*) 'Temperature_only = T:' // &
               ' Only temperature map will be produced'
       end if
       nmap = 1
    else
       if (id == 0 .and. info > 0) then
          write(*,*) 'Polarized detectors present: ' // &
               ' Will produce polarization maps'
       end if
    end if

    ncc = nmap * (nmap + 1) / 2

    ! Detector weighting

    if (info > 0) then
       if (id == 0) print *,'noise_weights_from_psd = ', noise_weights_from_psd
       if (id == 0) print *,'radiometers = ', radiometers
       if (id == 0) print *,'mode_detweight = ', mode_detweight
    end if

    if (mode_detweight == 0) then ! Detector weights from RMS

       do idet = 1,nodetectors
          where (detectors(idet)%sigmas == 0)
             detectors(idet)%weights = 0
          elsewhere
             detectors(idet)%weights = 1.d0 / detectors(idet)%plateaus / fsample
          end where
       end do

    elseif (mode_detweight == 1) then ! uniform weights

       do idet = 1, nodetectors
          detectors(idet)%weights = 1
       end do

    elseif (mode_detweight == 2) then ! one weight per horn
       ! each detector has its own pointing so figuring out which
       ! detectors share a horn is a little more involved
       do idet1 = 1,nodetectors-1
          horn1 = get_horn(detectors(idet1)%name)
          do idet2 = idet1+1,nodetectors
             horn2 = get_horn(detectors(idet2)%name)
             if (horn1 == horn2) then
                if (detectors(idet1)%npsd /= detectors(idet2)%npsd) then
                   call abort_mpi('Detectors in the same horn do not ' &
                        // 'have the same number of PSDs')
                end if

                allocate(weights(detectors(idet1)%npsd))

                where (detectors(idet1)%plateaus == 0 &
                     .or. detectors(idet2)%plateaus == 0)
                   weights = 0
                elsewhere
                   weights = 2 / (detectors(idet1)%plateaus &
                        + detectors(idet2)%plateaus) / fsample
                end where
                if (id == 0) write (*,'(a,g13.5,a,g13.5,a,g13.5)') &
                     'Assigning horn weight, ', weights(1), ', to ' // &
                     trim(detectors(idet1)%name) // ' and ' // &
                     trim(detectors(idet2)%name)  // ' for the FIRST period' &
                     // '. Original weights were ', detectors(idet1)%weights(1), &
                     ' and ', detectors(idet2)%weights(1)
                detectors(idet1)%weights = weights
                detectors(idet2)%weights = weights

                deallocate(weights)
             end if
          end do
       end do
    else
       print *, ' ERROR: bad mode_detweight : ', mode_detweight
       call abort_mpi('bad mode_detweight')
    end if

    nside_max = max(nside_map, nside_cross)
    nside_submap = min(nside_submap, nside_map, nside_cross)

    nosubmaps_tot = 12*nside_submap**2

    use_all_data = .true.

    istart_mission = 0
    use_all_data = .true.

    if (id == 0 .and. info > 2) then
       write(*,*)
       write(*,'(x,a,t24,"= ",i12)') 'istart_mission', istart_mission
       write(*,'(x,a,t24,"= ",i12)') 'nosamples_tot', nosamples_tot
    end if

    ! Non-integer baseline length allowed in standard mode + use_pntperiods=T

    nshort = int(dnshort + .5)

    ! Baseline length

    if (.not. kfirst) then
       dnshort = 1e6
       nshort = int(dnshort  +.5)
       precond_width = 0
       kfilter = .false.
    end if

    ! Data divided by pointing periods

    if (ninterval_tot <= 0) then
       if (id == 0) write(*, *) 'ERROR: intervals not known.'
       call exit_with_status(1)
    end if

    ! Store the number of baselines per pointing period
    ! Most loops now iterate over baselines, even when no destriping is done
    allocate(noba_short_pp(ninterval), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for noba_short_pp')
    memory_baselines = memory_baselines + ninterval*4
    noba_short_pp = 0

    do i = 1, ninterval
       if (id == 0 .and. info > 3) then
          if (intervals(i) < dnshort) &
               write (*,'(a,i0,a,i0,a,i0,a)') &
               'WARNING: interval ', i,' is ', &
               intervals(i), ' samples long but the baseline length is ', &
               int(dnshort, i8b), ' samples. Baseline truncated.'
       endif
       noba_short_pp(i) = (intervals(i)-1) / dnshort + 1
    end do
    noba_short = sum(noba_short_pp)

    ! Resolution
    nosubpix_max = nside_max**2 / nside_submap**2
    nosubpix_map = nside_map**2 / nside_submap**2
    nosubpix_cross = nside_cross**2 / nside_submap**2

    ! Message passing strategy
    if (allreduce) then
       concatenate_messages = .false.
       reassign_submaps = .true.
    end if

    if (checknan) then
       do idet = 1,nodetectors
          if (any(isnan(detectors(idet)%weights)) &
               .or. any(isnan(detectors(idet)%sigmas))) then
             print *, id, &
                  ' : Bad detector params (NaN in sigmas or weights) : ', &
                  detectors(idet)%name
          end if
       end do
    end if



  CONTAINS


    !------------------------------------------------------------------------


    SUBROUTINE nside_check(nside, varname)

      integer,intent(in) :: nside
      character(len=*)   :: varname
      integer            :: n, k

      k = 1
      n = nside
      do
         if (n==1) exit
         n = n/2
         k = k*2
      end do

      if (n*k.ne.nside) then
         write(*,*) 'ERROR: Illegal parameter value.'
         write(*,*) varname,' must be a power of two.'
         call exit_with_status(1)
      end if

    END SUBROUTINE nside_check


    !------------------------------------------------------------------------


    SUBROUTINE detector_check()

      ! Check to see if all the processes have the same detectors

      integer :: nodetectors_root, idet
      character(len=SLEN) :: detname_root

      nodetectors_root = nodetectors
      call broadcast_mpi(nodetectors_root, 0)
      if (nodetectors /= nodetectors_root) &
           call abort_mpi('Processes have different numbers of detectors')

      do idet = 1,nodetectors
         detname_root = trim(detectors(idet)%name)
         call broadcast_mpi(detname_root, 0)
         if (trim(detectors(idet)%name) /= trim(detname_root)) &
              call abort_mpi('Processes have different detectors')
      end do

    END SUBROUTINE detector_check


    !------------------------------------------------------------------------


    function get_horn(name)
      integer :: get_horn
      character(len=*) :: name

      integer :: ierr

      if (name(1:3) == 'LFI') then
         read(name(4:5), '(i2)', iostat=ierr) get_horn
      else
         read(name(5:5), '(i1)', iostat=ierr) get_horn
      end if

      if (ierr /= 0) &
           call abort_mpi('ERROR : Failed to extract horn from '//trim(name))

    end function get_horn


  END SUBROUTINE init_parameters


  !-------------------------------------------------------------------------


  SUBROUTINE init_parallelization

    ! Substantially rewritten: TOAST does now the data distribution.
    ! Looping over time, pointings and buffer are disabled.

    character(len=30) :: fi

    integer(i8b) :: nleft, kstart, n, nmax
    integer(i4b) :: i, k, iloop, nba, ierr

    fi = '(x,a,t24,"= ",i12,   2x,a)'

    if (id == 0 .and. info > 2) write(*,*) 'Initializing parallelization'
    if (info > 4) write(*,idf) id, 'Initializing parallelization'

    ! Compute the total number of short baselines (first destriping)
    ! and the maximum number per process.

    ! most loops now iterate over baselines regardless of kfirst
    noba_short_max = 0
    noba_short_tot = 0

    noba_short_max = noba_short
    call max_mpi(noba_short_max)

    noba_short_pp_max = maxval(noba_short_pp)
    call max_mpi(noba_short_pp_max)

    noba_short_tot = noba_short
    call sum_mpi(noba_short_tot)

    ! Distribute submaps

    if (allocated(id_submap)) then
       deallocate(id_submap)
       memory_maps = memory_maps - nosubmaps_tot*4
    end if

    allocate(id_submap(0:nosubmaps_tot-1), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for id_submap')
    memory_maps = memory_maps + nosubmaps_tot*4

    nosubmaps_max = (nosubmaps_tot-1) / ntasks + 1

    nosubmaps = 0
    do i = 0, nosubmaps_tot-1
       k = mod(i, ntasks)
       id_submap(i) = k
       if (id == k) nosubmaps = nosubmaps + 1
    end do

    nopix_map = nosubmaps * nosubpix_map
    nopix_cross = nosubmaps * nosubpix_cross

    if (id == 0 .and. info > 2) then
       write(*,fi) 'ntasks', ntasks, 'Number or processes'
       write(*,fi) 'nosamples_tot',     nosamples_tot,     'Total samples'
       write(*,fi) 'nosamples_proc_max',nosamples_proc_max,'Samples/process'

       if (kfirst) then
          write(*,fi) 'noba_short_tot', noba_short_tot, 'Total baselines'
          write(*,fi) 'noba_short_pp_max', noba_short_pp_max, &
               'Longest interval in baselines'
          write(*,fi) 'noba_short_max', noba_short_max, 'Baselines/process'
       end if

       write(*,*)
    end if

  END SUBROUTINE init_parallelization


  !--------------------------------------------------------------------------


  SUBROUTINE write_parameters()

    integer :: idet, itod, i, j
    character(len=30) :: fi, fe, ff, fk, fs
    real(dp) :: sigma

    fi = '(x,a,t24,"= ",i12,   t40,a)'
    fe = '(x,a,t24,"= ",es12.5,t40,a)'
    ff = '(x,a,t24,"= ",f12.4, t40,a)'
    fk = '(x,a,t24,"= ",l12,   t40,a)'
    fs = '(x,a,t24,"= ",a,     t40,a)'

    if (ntasks > 1 .and. id /= 0) return
    if (info == 0) return

    write (*,fk) 'MCMode', mcmode
    write (*,fk) 'write_cut', write_cut
    select case (basis_func)
    case (basis_poly)
       write (*,fs) 'basis_func', 'Polynomial', 'Destriping function basis'
    case (basis_fourier)
       write (*,fs) 'basis_func', 'Fourier', 'Destriping function basis'
    case (basis_cheby)
       write (*,fs) 'basis_func', 'Chebyshev', 'Destriping function basis'
    case (basis_legendre)
       write (*,fs) 'basis_func', 'Legendre', 'Destriping function basis'
    case default
       call abort_mpi('Unknown function basis')
    end select
    write (*,fi) 'basis_order', basis_order, 'Destriping function order'

    write (*, *)
    write (*, fk) 'bin_subsets', bin_subsets
    if (bin_subsets) then
       write (*, '(x,a,t24,"= ")') 'Surveys'
       do i = 1, nsurvey
          write (*, '(i4,4x,a15,"  :  ")', advance='no') i, trim(surveys(i)%name)
          do j = 1, surveys(i)%nspan
             write (*, '("( ",f15.0," -- ",f15.0," )")', advance='no') &
                  surveys(i)%starts(j), surveys(i)%stops(j)
          end do
          write (*, *)
       end do
       write (*, '(x,a,t24,"= ")') 'Detector sets'
       do i = 1, ndetset
          write (*, '(i4,4x,a15,"  :  ")', advance='no') i, trim(detsets(i)%name)
          do j = 1, detsets(i)%ndet
             write (*, '(a,x)', advance='no') trim(detsets(i)%detectors(j))
          end do
          write (*, *)
       end do
    end if

    write (*,*)
    write (*,fi) 'ntasks',ntasks,'Number of processes'
    write (*,fi) 'nthreads',nthreads,'Number of threads per process'
    if (nthreads < nthreads_max) &
         write (*,fi) 'nthreads_max',nthreads_max, &
         'Maximum number of threads per process'
    write (*,fi) 'info',info,'Screen output level'
    write (*,ff) 'fsample',fsample,'Sampling frequency (Hz)'

    if (nmap /= 1) then
       write (*,fi) 'nmap',nmap,'Polarization included'
    else
       write (*,fi) 'nmap',nmap,  &
            'Temperature-only treatment (no polarization)'
    end if
    write (*,fi) 'ncc',ncc,'Independent wcov elements'

    if (nsubchunk > 1) then
       write (*,*)
       write (*,*) 'Mapping SUBCHUNKS'
       write (*,fi) 'nsubchunk', nsubchunk, &
            'split chunks into subchunks and map all separately'
       if (isubchunk > 0) write (*,fi) 'isubchunk', isubchunk, &
            'Begin from a preset subchunk instead of all data'
       write (*,*)
    end if

    if (mc_loops > 1) then
       write (*,*)
       write (*,*) 'Monte Carlo mode ON:'
       write (*,fi) 'mc_loops', mc_loops,  &
            'number of MC maps to produce'
       write (*,fi) 'mc_increment', mc_increment,  &
            'random stream increment between MC loops'
       write (*,*)
    end if

    if (skip_existing) write (*,*) 'Will SKIP existing files'

    write (*,*) 'Input files:'
    if (len_trim(file_gap) > 0) write (*,fs) 'file_gap',trim(file_gap)

    if (len_trim(file_inmask) > 0) then
       write (*,fs) 'file_inmask',trim(file_inmask)
       use_inmask = .true.
    end if

    write (*,*)
    write (*,fi) 'nside_map', nside_map, 'Healpix resolution (output map)'
    write (*,fi) 'nside_cross', nside_cross, 'Healpix resolution (destriping)'
    if (info > 1) &
         write (*,fi) 'nside_submap', nside_submap, 'Submap resolution'
    write (*,fk) 'concatenate_messages', concatenate_messages, &
         'use mpi_alltoallv to communicate'
    write (*,fk) 'allreduce', allreduce, &
         'use allreduce to communicate'
    write (*,fk) 'reassign_submaps', reassign_submaps, &
         'minimize communication by reassigning submaps'

    if (info > 1) then
       write (*,fi) 'pixmode_map',  pixmode_map,    &
            'Pixel rejection criterion (output map)'
       write (*,fi) 'pixmode_cross',pixmode_cross,  &
            'Pixel rejection criterion (destriping)'
       write (*,fe)  'pixlim_map',   pixlim_map,    &
            'Pixel rejection limit (output map)'
       write (*,fe)  'pixlim_cross', pixlim_cross,    &
            'Pixel rejection limit (destriping)'
    end if

    write (*,*)
    write (*,*) 'Standard mode'

    if (basis_order > 0 .and. kfilter) &
         call abort_mpi('Filter only implemented for basis_order==0')

    if (noise_weights_from_psd .or. kfilter) then
       write (*,fi) 'psdlen',psdlen,'Length of requested noise PSD'
       write (*,fi) 'psd_downsample',psd_downsample,'PSD downsampling factor'
    end if

    if (kfirst) then
       write (*,*)
       write (*,fk) 'kfirst',kfirst, 'First destriping ON'
       write (*,ff) 'dnshort',dnshort,'Baseline length (samples)'
       write (*,ff) ' ',dnshort/fsample,'seconds'

       if (kfilter) then
          write (*,fk) 'kfilter',kfilter,'Noise filter ON'
       else
          write (*,fk) 'kfilter', kfilter, 'Noise filter OFF'
          write (*,fe) 'diagfilter', diagfilter, 'diagonal baseline filter'
          write (*,ff) 'good_baseline_fraction', good_baseline_fraction, &
               'fraction of samples needed to use baseline'
       end if
    else
       write (*,fk) 'kfirst',kfirst, 'First destriping OFF'
    end if

    write (*,*)
    if (info > 1) then
       write (*,fe) 'cglimit',cglimit, 'Iteration convergence limit'
       write (*,fi) 'iter_min',iter_min, 'Minimum number of iterations'
       write (*,fi) 'iter_max',iter_max, 'Maximum number of iterations'
    end if

    write(*,fi) 'precond_width', precond_width, &
         'Width of the preconditioner band matrix'
    if (precond_width==0) then
       write(*,*) 'No preconditioning'
    else
       write (*,fk) 'use_fprecond', use_fprecond, 'use C_a preconditioner'
    end if

    if (flag_by_horn) then
       write (*,fk) 'flag_by_horn',flag_by_horn,'Combining flags within horns'
    else
       write (*,fk) 'flag_by_horn',flag_by_horn,'Flags are independent'
    end if

    if (mode_detweight == 0) then
       write (*,fi) 'mode_detweight',mode_detweight,  &
            'Detector weighting mode: sigma from simulation file'
    elseif (mode_detweight == 1) then
       write (*,fi) 'mode_detweight',mode_detweight,  &
            'Detector weighting mode: uniform'
    elseif (mode_detweight == 2) then
       write (*,fi) 'mode_detweight',mode_detweight,  &
            'Detector weighting mode: one weight per horn'
    end if

    if (kwrite_covmat) then
       write (*,fk) 'kwrite_covmat', kwrite_covmat, 'Covariance matrix written'
       write (*,fs) 'file_covmat', trim(file_covmat)
    end if

    write (*,*)
    if (mc_loops > 1) then
       write (*,fi) 'mc_increment', mc_increment, &
            'Random number increment for iterations'
       write (*,fi) 'mc_loops', mc_loops, &
            'Number of Monte Carlo iterations'
       write (*,fi) 'mc_id', mc_id, &
            'Starting iteration identifier'
    end if

    write (*,*)
    if (time_unit >= 0) then
       write (*,fe) 'time_unit',time_unit,'Time unit/samples'
    else
       write (*,fs) 'time_unit','          pp','Time unit = pointing period'
    end if
    write (*,fi) 'mission_time', ninterval_tot, 'Mission length in time units'
    write (*,fi) 'nosamples_tot', nosamples_tot, 'Total samples'
    write (*,ff) '', nosamples_tot/fsample/3600., 'hours'

    write (*,*)
    write (*,'(x,a)') 'Detectors available on the FIRST process and noise ' &
         // 'according to the FIRST period'
    write (*,'(x,a12,3a15)') 'detector    ', 'sigma', 'weight', '1/sqrt(weight)'
    do idet = 1, nodetectors
       sigma = 0
       if (detectors(idet)%weights(1) > 0) then
          sigma = 1 / sqrt(detectors(idet)%weights(1))
       end if
       write (*,'(x,a12,3g15.5)') detectors(idet)%name, &
            detectors(idet)%sigmas(1), detectors(idet)%weights(1), sigma
    end do

    if (info > 1) then
       write (*,*)
       write (*,*) 'Output files'
       write (*,*) 'file_root    = ',trim(file_root)
       if (len_trim(file_map) > 0)      &
            write (*,*) 'file_map     = ',trim(file_map)
       if (len_trim(file_binmap) > 0)   &
            write (*,*) 'file_binmap  = ',trim(file_binmap)
       if (len_trim(file_mask) > 0)     &
            write (*,*) 'file_mask    = ',trim(file_mask)
       if (len_trim(file_hit) > 0)      &
            write (*,*) 'file_hit     = ',trim(file_hit)
       if (len_trim(file_matrix) > 0)   &
            write (*,*) 'file_matrix  = ',trim(file_matrix)
       if (len_trim(file_wcov) > 0)   &
            write (*,*) 'file_wcov    = ',trim(file_wcov)
       if (len_trim(file_leakmatrix) > 0)   &
            write (*,*) 'file_leakmatrix  = ',trim(file_leakmatrix)
       if (len_trim(file_base) > 0)     &
            write (*,*) 'file_base    = ',trim(file_base)
       if (len_trim(file_gap_out) > 0)  &
            write (*,*) 'file_gap_out = ',trim(file_gap_out)
       if (len_trim(file_mc) > 0)       &
            write (*,*) 'file_mc      = ',trim(file_mc)
       write (*,*) 'binary_output      = ',binary_output
       write (*,*) 'concatenate_binary = ',concatenate_binary
    end if
    write (*,*)

  END SUBROUTINE write_parameters


  !--------------------------------------------------------------------------


  SUBROUTINE start_timeloop()
    ! Initialize the next split-mode step.
    ! Find the TOD chunk handled in this step. In standard mode = all TOD.
    ! Also initialize parameters defining the division of TOD to processes.
    integer :: i, j, k, m, n0, n, ierr
    integer(i8b) :: nsamp, order, my_offset
    integer(i8b) :: sublen, suboffset, sub_start, sub_end, isub
    real(dp) :: dn, dn0, r, dr, rstart, ninv
    real(dp), pointer :: basis_function(:,:)
    real(sp) :: memsum, mem_min, mem_max

    call wait_mpi

    ! Samples handled by the current process

    istart_proc = 1

    ! Short baselines handled by the current process.

    ! most loops iterate over baselines when there is no destriping
    kshort_start = 0

    ! Store the baseline lengths for first destriping

    allocate(baselines_short(noba_short), &
         basis_functions(noba_short), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for baselines_short')

    basis_functions%copy = .true.
    basis_functions%nsamp = 0

    memory_baselines = memory_baselines + noba_short*(4+4+17)

    !Split the interval into short baselines of length int(dnshort)
    ! or int(dnshort+1).
    m = 0
    baselines_short = 0
    do i = 1, ninterval
       dn = 0.0
       n = 0
       do k = 1, noba_short_pp(i) - 1 ! Loop over all but the last
          dn0 = dn
          n0 = n
          dn = dn0 + dnshort
          n = int(dn + .5)
          m = m + 1
          baselines_short(m) = n - n0
       end do
       m = m + 1
       baselines_short(m) = intervals(i) - n ! Last baseline takes the rest

       if (kfirst) then
          if (baselines_short(m) < 0 &
               .or. baselines_short(m) > int(dnshort, i8b) + 1) then
             write (*,*) id, ' : ERROR in start_timeloop: ' &
                  // 'Lengths do not match., dnshort = ', dnshort
             write (*,*) id, ' : i, intervals(i), noba_short =', i, &
                  intervals(i), noba_short
             write (*,*) id, ' : last baseline =', baselines_short(m)
             print *,id, ' : noba_short_pp(i) = ', noba_short_pp(i)
             do k=1,noba_short_pp(i)
                print *, k, baselines_short(m-k+1)
             end do
             call exit_with_status(1)
          end if
       end if
    end do

    ! Split the interval into subchunks

    if (nsubchunk > 0) then
       my_offset = 0
       do i = 1, ninterval
          n = intervals(i)
          sublen = n / nsubchunk
          do isub = 1, nsubchunk
             suboffset = sublen * (isub - 1)
             sub_start = my_offset + 1 + suboffset
             sub_end = sub_start + sublen - 1
             if (isub == nsubchunk) sub_end = my_offset + n
             subchunk(sub_start:sub_end) = isub
          end do
          my_offset = my_offset + n
       end do
    end if

    ! Auxiliary arrays for OpenMP threading and output

    allocate(baselines_short_start(noba_short), &
         baselines_short_stop(noba_short), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for baselines_short_start')

    memory_baselines = memory_baselines + noba_short*(4+4)

    baselines_short_start = 1
    baselines_short_stop = 1

    do k = 1, noba_short
       if (k > 1) baselines_short_start(k) = baselines_short_start(k-1) &
            + baselines_short(k-1)
       baselines_short_stop(k) = baselines_short_start(k) &
            + baselines_short(k) - 1
    end do

    if (kfirst) then
       ! Set up the basis function arrays.

       do k = 1,noba_short
          basis_functions(k)%copy = .false.
          basis_functions(k)%nsamp = baselines_short(k)
          basis_functions(k)%arr => NULL()
          loop_previous: do j = 1, k-1
             if (basis_functions(j)%nsamp == basis_functions(k)%nsamp) then
                ! we already have a basis function of this length stored.
                basis_functions(k)%arr => basis_functions(j)%arr
                basis_functions(k)%copy = .true.
                exit loop_previous
             end if
          end do loop_previous
          if (.not. basis_functions(k)%copy) then
             ! Allocate and initialize the basis function array for
             ! this baseline length
             nsamp = basis_functions(k)%nsamp
             allocate(basis_functions(k)%arr(0:basis_order, 0:nsamp-1), stat=ierr)
             if (ierr /= 0) stop 'No room for basis function'
             memory_basis_functions = memory_basis_functions &
                  + (basis_order+1)*nsamp*8
             basis_function => basis_functions(k)%arr
             if (nsamp < 1) cycle
             dr = 2. / nsamp
             rstart = 0.5*dr - 1
             select case (basis_func)
             case (basis_poly)
                ! Simple polynomial basis
                do i = 0,nsamp-1
                   r = rstart + i*dr
                   do order = 0,basis_order
                      basis_function(order, i) = r**order
                   end do
                end do
             case (basis_fourier)
                ! Real Fourier basis
                ninv = pi / nsamp
                do i = 0,nsamp-1
                   r = i * ninv
                   do order = 0,basis_order
                      if (modulo(order, 2) == 0) then
                         basis_function(order, i) = cos(order * r)
                      else
                         basis_function(order, i) = sin((order + 1) * r)
                      end if
                   end do
                end do
             case (basis_cheby)
                ! use recursive formula for Chebyshev polynomials
                ! of the second kind
                do i = 0,nsamp-1
                   r = rstart + i*dr
                   do order = 0,basis_order
                      if (order == 0) basis_function(order, i) = 1
                      if (order == 1) basis_function(order, i) = 2 * r
                      if (order > 1) basis_function(order, i) = &
                           2*r*basis_function(order-1, i) &
                           - basis_function(order-2, i)
                   end do
                   ! normalize, Chebyshev polynomials are not orthogonal
                   ! wrt the L2 norm
                   basis_function(:, i) = basis_function(:, i) * (1 - r**2)**.25
                end do
             case (basis_legendre)
                ! use recursive formula for Legendre polynomials
                do i = 0,nsamp-1
                   r = rstart + i*dr
                   do order = 0,basis_order
                      if (order == 0) basis_function(order, i) = 1
                      if (order == 1) basis_function(order, i) = r
                      if (order > 1) basis_function(order, i) = &
                           ((2*order-1)*r*basis_function(order-1, i) &
                           - (order-1)*basis_function(order-2, i)) / order
                   end do
                end do
             case default
                call abort_mpi('Unknown function basis')
             end select
          end if
       end do
    end if ! if (kfirst)

    memsum = memory_basis_functions / 2**20
    mem_min = memsum
    mem_max = memsum
    call min_mpi(mem_min)
    call max_mpi(mem_max)
    call sum_mpi(memsum)
    if (ID == 0 .and. info > 0) then
       write(*,mstr3) 'Allocated memory for basis_functions:', &
            memsum, mem_min, mem_max
    end if

    !end if

    if (info > 3) then
       write (*,'(i3,a,i11)') ID,': nosamples_proc =', nosamples_proc
       write (*,'(i3,a,i11)') ID,': istart_proc    =', istart_proc
       if (kfirst) write (*,'(i3,a,i11)') ID,': noba_short     =', noba_short
    end if

  END SUBROUTINE start_timeloop


  !------------------------------------------------------------------------


  subroutine baseline_times(short_times, sample_times)

    ! new utility routine for madam/TOAST

    real(dp), allocatable, intent(out) :: short_times(:)
    real(c_double), pointer, intent(in) :: sample_times(:)

    integer(i8b) :: i, j
    integer :: ierr

    ! insert the baseline start times
    ! local array of short baseline start times
    allocate(short_times(noba_short), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for short_times')

    memory_baselines = memory_baselines + noba_short*8

    short_times = sample_times(baselines_short_start)

  end subroutine baseline_times

END MODULE parameter_control
