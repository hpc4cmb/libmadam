module smadam

  ! Program MADAM, version 3.7
  ! Elina Keihanen, September 2010
  !

  use iso_c_binding

  use commonparam
  use inputparam
  use simulation
  use parameter_control
  use pointing
  use maptod_transfer
  use submap_transfer
  !use compression
  use fourier
  use noise_routines
  use madam_routines
  use map_routines
  use read_data
  use output
  use mpi_wrappers
  use maps_and_baselines
  use tod_storage
  use memory_and_time
  use timing
  use covmat
  use covmat_util, only : tic, toc

  implicit none

  integer :: idet

  real(sp) :: cputime_init=.0, cputime_final=.0, cputime_total=.0, &
       cputime_wait=.0, cputime_read=.0

  real(dp), pointer :: temparr(:,:)=>null()
  character(len=256) :: outfile
  character(len=80) :: header(1000)
  integer :: ierr, best_nside_submap, test_nside_submap
  integer(idp) :: isample

  ! covmat
  integer(i4b), pointer :: ppix(:)
  real(dp), pointer :: piw(:), pqw(:), puw(:)
  integer(i4b), pointer :: pflag(:)
  integer(i4b) :: flagmask, detbit
  complex(dp), pointer :: pfcov(:)

  ! openmp
  integer :: nprocs

  integer(i4b) :: subchunk_start, i
  character(len=SLEN) :: subchunk_file_map, subchunk_file_base, subchunk_file_binmap
  character(len=SLEN) :: subchunk_file_hit, subchunk_file_mask, subchunk_file_matrix, subchunk_file_wcov
  ! subset mapping
  real(sp) :: cputime_flag_subset=.0, cputime_subset=.0, cputime_write_subset=.0

contains

  subroutine destripe( comm, parstring, ndet, detstring, detweights, &
       nsamp, nnz, timestamps, pix, pixweights, signal, &
       nperiod, periods, &
       npsd, npsdtot, psdstarts, npsdbin, psdfreqs, npsdval, psdvals ) bind( c, name='destripe' )

    integer(c_int), intent(in), value :: comm ! MPI communicator

    character(kind=c_char), intent(in) :: parstring(*)

    integer(c_long), intent(in), value :: ndet
    character(kind=c_char), intent(in) :: detstring(*)
    real(c_double), intent(in) :: detweights(ndet)

    integer(c_long), intent(in), value :: nsamp
    integer(c_long), intent(in), value :: nnz

    ! TOD

    type(c_ptr), intent(in), value :: timestamps
    type(c_ptr), intent(in), value :: pix
    type(c_ptr), intent(in), value :: pixweights
    type(c_ptr), intent(in), value :: signal

    ! Pointing periods

    integer(c_long), intent(in), value :: nperiod
    integer(c_long), intent(in) :: periods(nperiod)

    ! Noise PSDs

    integer(c_long), intent(in) :: npsd(ndet)
    integer(c_long), intent(in), value :: npsdtot
    real(c_double), intent(in) :: psdstarts(npsdtot)
    integer(c_long), intent(in), value :: npsdbin
    real(c_double), intent(in) :: psdfreqs(npsdbin)
    integer(c_long), intent(in), value :: npsdval
    real(c_double), intent(in) :: psdvals(npsdval)

    integer :: ierr, idet, i, pixmin, pixmax

    ! set up MPI

    call init_mpi( comm, ntasks, id )

    ! set up OpenMP

    nprocs = omp_get_num_procs()
    nthreads_max = omp_get_max_threads()
    nthreads = nthreads_max
    if (id == 0) then
       write (*,'("OMP: ",i0," tasks with ",i0," procs per node, ",i0, &
            " threads per task.")') ntasks, nprocs, nthreads
    end if

    call reset_time()

    if (id == 0) then
       write (*,*)
       write (*,*) 'Program MADAM'
       write (*,*) 'Destriping of CMB data with a noise filter'
       write (*,*) 'Version ',version
       write (*,*)
    endif

    ! Then parse the parameter string

    call read_parameters( parstring )

    if (mcmode .and. cached) call abort_mpi('Destripe called while caches are not empty.')
    if (mcmode .and. nsubchunk > 1) call abort_mpi('MCMode is not compatible with nsubchunk > 1.')
    
    nmap = nnz

    call read_detectors( detstring, ndet, detweights, npsd, npsdtot, psdstarts, npsdbin, psdfreqs, npsdval, psdvals )

    call c_f_pointer(timestamps, sampletime, (/nsamp/) )
    call c_f_pointer(pix, pixels, (/nsamp, ndet/) )
    call c_f_pointer(pixweights, weights, (/nnz, nsamp, ndet/) )
    call c_f_pointer(signal, tod_stored, (/nsamp, ndet/) )

    ! Do a basic check for pointing

    pixmin = minval( pixels )
    pixmax = maxval( pixels )

    nside_max = max( nside_map, nside_cross )

    if (pixmin < -1) then
       print *,'ERROR: smallest provided pixel number is ',pixmin
       call abort_mpi('Too small pixel number')
    end if
    if (pixmax > 12*nside_max**2-1 ) then
       print *,'ERROR: largest provided pixel number is ',pixmax, '>', 12*nside_max**2, ', nside = ',nside_max
       call abort_mpi('Too large pixel number')
    end if

    call check_files( nperiod, periods, nsamp )

    ! set output level for timing based on info level
    call tic
    select case(info)
    case (:1)
       call toc(threshold=1000.0_dp)
    case (2)
       call toc(threshold=1.0_dp)
    case (3)
       call toc(threshold=0.1_dp)
    case (4:)
       call toc(threshold=0.0_dp)
    end select

    call init_parameters
    call init_parallelization

    subchunk_start = isubchunk

    loop_subchunk : do isubchunk = subchunk_start, nsubchunk

       subchunk_file_map = file_map; subchunk_file_binmap = file_binmap
       subchunk_file_base = file_base; subchunk_file_hit = file_hit; subchunk_file_mask = file_mask
       subchunk_file_matrix = file_matrix; subchunk_file_wcov = file_wcov

       tod_is_clean = .false.

       detflags = .true.

       if (nsubchunk > 1) then
          baseline_open = .false. ! from output.f90
          if (id == 0) then
             write (*,'(/," ********** SUBCHUNK == ",i0,/)') isubchunk
          endif
          call add_subchunk_id(file_map, isubchunk, nsubchunk)
          call add_subchunk_id(file_hit, isubchunk, nsubchunk)
          call add_subchunk_id(file_mask, isubchunk, nsubchunk)
          call add_subchunk_id(file_matrix, isubchunk, nsubchunk)
          call add_subchunk_id(file_wcov, isubchunk, nsubchunk)
          call add_subchunk_id(file_binmap, isubchunk, nsubchunk)
          call add_subchunk_id(file_base, isubchunk, nsubchunk)
       end if

       if (skip_existing) then
          if ( file_exists(file_map) ) then
             file_map = subchunk_file_map; file_binmap = subchunk_file_binmap
             file_base = subchunk_file_base; file_hit = subchunk_file_hit; file_mask = subchunk_file_mask
             file_matrix = subchunk_file_matrix; file_wcov = subchunk_file_wcov
             call reset_timers
             cycle loop_subchunk
          end if
       end if

       call init_output

       call write_parameters
       call wait_mpi

       if ( isubchunk == subchunk_start ) then
          call allocate_tod

          call allocate_baselines

          call tic
          call init_pointing
          where ( pixels < 0 ) pixels = dummy_pixel
          if (id == 0) call toc('init_pointing')

          call wait_mpi

          cputime_init = cputime_init + get_time(1)

          call tic
          call build_filter
          if (id == 0) call toc('build_filter')

          call reset_time(1)

          call start_timeloop
          call time_stamp

          call baseline_times(baselines_short_time, sampletime)
       else
          if ( kfirst ) then
             aa = 0
             yba = 0
             nna = 0
             nna_inv = 0
          end if
       end if

       subchunkpp => subchunk(1:nosamples_proc)

       call tic
       call reduce_pixels_buff
       if (id == 0) call toc('reduce_pixels_buff')

       call tic
       call update_maptod_transfer(ksubmap)
       if (id == 0) call toc('update_maptod_transfer')

       if ( reassign_submaps ) then
          call tic
          call assign_submaps( id_submap, nosubmaps, nopix_map, nopix_cross, nosubmaps_max )
          if (id == 0) call toc('assign_submaps')
       end if

       call allocate_maps

       call tic
       call initialize_alltoallv()
       if (id == 0) call toc('initialize_alltoallv')

       dummy_pixel = 12*nside_max**2

       if (kfirst) then
          call tic
          call flag_bad_baselines
          if (id == 0) call toc('flag_bad_baselines')
       end if

       call tic
       call pixel_matrix(cca, cc)
       if (id == 0) call toc('pixel_matrix')

       call tic
       call count_hits(nohits)
       if (id == 0) call toc('count_hits')

       call tic
       call bin_tod(map, binmap, wamap, tod_stored)
       if (id == 0) call toc('bin_tod')

       call time_stamp

       call tic
       call restore_pixels_buff
       if (id == 0) call toc('restore_pixels_buff')

       !  First destriping

       if (kfirst) then

          if (id == 0) then
             write(*,*)
             write(*,*) 'First destriping phase'
          endif
          call time_stamp

          call tic
          call reduce_pixels_a
          if (id == 0) call toc('reduce_pixels_a')

          !call tic
          !call update_maptod_transfer(ksubmap)
          !if (id == 0) call toc('update_maptod_transfer')

          if (use_inmask) then
             if ( isubchunk == subchunk_start ) then
                call tic
                call read_inmask(inmask)
                if (id == 0) call toc('read_inmask')
             end if
             call tic
             call scatter_mask(inmask, nosubpix_cross)
             if (id == 0) call toc('scatter_mask')
          end if

          call tic
          call invert_pixelmatrix_cross(cca,inmask)
          if (id == 0) call toc('invert_pixelmatrix_cross')

          call tic
          call initialize_a(yba, nna, wamap, cca, tod_stored)
          if (id == 0) call toc('initialize_a')

          if (basis_order == 0) then
             call tic
             call construct_preconditioner(nna(0,0,:,:))
             if (id == 0) call toc('construct_preconditioner')
          end if

          if (kwrite_covmat) then

             call write_covmat(file_covmat) ! pixels must be reduced

          else

             call tic
             call iterate_a(aa, yba, nna, wamap, cca, tod_stored)
             if (id == 0) call toc('iterate_a')

             call tic
             call subtract_baselines_a(map, aa)
             if (id == 0) call toc('subtract_baselines_a')

             call tic
             call write_baselines_a(aa)
             if (id == 0) call toc('write_baselines_a')
             
          end if

       endif

       call wait_mpi
       if (id == 0) then
          write(*,*)
          write(*,*) 'Finalization begins'
          if (info.ge.2) write(*,*)
       endif
       call time_stamp
       call reset_time(1)

       call tic
       call write_matrix(cc)
       if (id == 0) call toc('Write matrix')
       !
       !  Binning of the final map
       !
       call tic
       call invert_pixelmatrix_map(cc, outmask, crit)
       if (id == 0) call toc('Invert pixelmatrix_map')

       if (do_binmap) call makemaps(binmap, cc, outmask)
       call makemaps(map, cc, outmask)
       call map_analysis(map, outmask)
       !
       ! output files
       !
       call tic
       call write_matrix(cc, outmask)
       if (id == 0) call toc('Write matrix')
       call tic
       call write_map(map, outmask)
       if (id == 0) call toc('Write map')
       call tic
       call write_binmap(binmap, outmask)
       if (id == 0) call toc('Write binmap')
       call tic
       call write_mask(outmask, crit)
       if (id == 0) call toc('Write mask')
       call tic
       call write_hits(nohits)
       if (id == 0) call toc('Write hits')
       call tic

       call run_subsets()

       if ( isubchunk == nsubchunk .and. kfirst ) then
          ! subtract the baselines to return the destriped TOD.
          ! If the baselines were already subtracted, the call has no effect.
          call tic
          call clean_tod(tod_stored, aa)
          if (id == 0) call toc('clean_tod')
       end if

       if ( kfirst ) then
          call restore_pixels_a
       end if

       if ( .not. mcmode ) then
          if (isubchunk == nsubchunk) then
             call free_baselines
          end if
          call free_maps
          call free_locmaps
       else
          cached = .true.
       end if

       call wait_mpi
       call time_stamp

       cputime_final = cputime_final + get_time(1)
       cputime_total = get_time(0)

       ! Memory

       if (info > 0) then
          if (id == 0) write(*,*)
          if (id == 0) write(*,*) 'MEMORY (MB):'

          call write_memory('Detector pointing',  memory_pointing)
          call write_memory('TOD buffer',         memory_tod)
          call write_memory('Maps',               memory_maps)
          call write_memory('Baselines',          memory_baselines)
          call write_memory('Basis functions',    memory_basis_functions)
          call write_memory('Noise filter',       memory_filter)
          call write_memory('Preconditioner',     memory_precond)
          call write_memory('Temporary maps',     memory_locmap)
          call write_memory('All2All buffers',    memory_all2all)
          call write_memory('CG work space',      memory_cg)
          call write_memory('NCM', memory_ncm)
          call write_memory('Total')
       endif

       ! Timing

       if (info > 0) then
          if (id == 0) write(*,*) 'WALL-CLOCK TIME (s):'

          call write_time('Initialization',     cputime_init)

          cputime_read = cputime_read_pointing_periods + cputime_read_detpointing + cputime_read_tod + cputime_read_timestamps

          call write_time('I/O', cputime_read)
          call write_time('- Pointing Period boundaries', cputime_read_pointing_periods)
          call write_time('- Pointing I/O',  cputime_read_detpointing)
          call write_time('- Time stamp I/O', cputime_read_timestamps)
          call write_time('- TOD I/O', cputime_read_tod)
          call write_time('- Other')

          call write_time('Waiting',                cputime_wait)
          call write_time('Building pixel matrices',cputime_build_matrix)
          call write_time('Sending pixel matrices',cputime_send_matrix)
          call write_time('Inverting pixel matrices',cputime_inv)
          call write_time('Binning TOD',            cputime_bin_maps)
          call write_time('Sending binned TOD',     cputime_send_maps)
          call write_time('Counting hits',          cputime_count_hits)
          call write_time('Building preconditioner',cputime_prec_construct)
          call write_time('Subtract/add baselines', cputime_clean_tod+cputime_unclean_tod)

          call write_time('Initialization (1. phase)', cputime_cga_init)
          call write_time('CG iteration',          cputime_cga)
          call write_time('- TOD - map',           cputime_cga_1)
          call write_time('- CG MPI',  cputime_cga_mpi)
          call write_time('- CG ccmultiply',  cputime_cga_cc)
          call write_time('- Map - TOD',           cputime_cga_2)
          call write_time('- Filtering',           cputime_filter)
          call write_time('- Preconditioning',     cputime_precond)
          call write_time('- Other')

          call write_time('NCM', cputime_ncm)
          call write_time('- P^T F', cputime_ptf)
          call write_time('- middlematrix', cputime_middlematrix)
          call write_time('- invert middlematrix', cputime_invert_middlematrix)
          call write_time('- symmetrize middlematrix', cputime_symmetrize_middlematrix)
          call write_time('- accumulate', cputime_accumulate)
          call write_time('- add white', cputime_white)
          call write_time('- symmetrize', cputime_symmetrize)
          call write_time('- save', cputime_save_matrix)
          call write_time('- Other')

          call write_time('Subset maps', cputime_subset )
          call write_time('- Flag subset', cputime_flag_subset )
          call write_time('- Write subset', cputime_write_subset )

          call write_time('Finalization and output',cputime_final)
          call write_time('Other',cputime_total-time_cum)
       endif

       call write_time('Total',cputime_total)
       if (id == 0) write(*,*)

       file_map = subchunk_file_map; file_binmap = subchunk_file_binmap
       file_base = subchunk_file_base; file_hit = subchunk_file_hit; file_mask = subchunk_file_mask
       file_matrix = subchunk_file_matrix; file_wcov = subchunk_file_wcov

       if (isubchunk == nsubchunk) exit loop_subchunk

       call reset_timers

    end do loop_subchunk

    isubchunk = subchunk_start

    call reset_timers

    if (.not. mcmode) then

       if (allocated(pntperiods)) deallocate( pntperiods, pntperiod_id, noba_short_pp )
       if (allocated(baselines_short)) deallocate( baselines_short, base_pntid_short, basis_functions )
       if (allocated(baselines_short_start)) deallocate( baselines_short_start, baselines_short_stop )
       if (allocated(prec_diag)) deallocate( prec_diag )
       if (allocated(bandprec)) deallocate( bandprec )

       call close_filter()
       call close_output()
       call close_pointing()
       call free_mask()
       call free_tod()
    end if

    call close_mpi()
    
    return

  end subroutine destripe


  subroutine destripe_with_cache( comm, ndet, nsamp, nnz, &
       timestamps, pix, pixweights, signal ) bind( c, name='destripe_with_cache' )

    integer(c_int), intent(in), value :: comm ! MPI communicator

    integer(c_long), intent(in), value :: ndet
    integer(c_long), intent(in), value :: nsamp
    integer(c_long), intent(in), value :: nnz

    ! TOD

    type(c_ptr), intent(in), value :: timestamps
    type(c_ptr), intent(in), value :: pix
    type(c_ptr), intent(in), value :: pixweights
    type(c_ptr), intent(in), value :: signal

    integer :: ierr, idet, i, pixmin, pixmax

    ! set up MPI

    call init_mpi( comm, ntasks, id )

    ! set up OpenMP

    nprocs = omp_get_num_procs()
    nthreads_max = omp_get_max_threads()
    nthreads = nthreads_max
    if (id == 0) then
       write (*,'("OMP: ",i0," tasks with ",i0," procs per node, ",i0, &
            " threads per task.")') ntasks, nprocs, nthreads
    end if

    call reset_time()

    if (id == 0) then
       write (*,*)
       write (*,*) 'Program MADAM (with cache)'
       write (*,*) 'Destriping of CMB data with a noise filter'
       write (*,*) 'Version ',version
       write (*,*)
    endif

    if (.not. cached) call abort_mpi('destripe_with_cache called with empty caches.')

    nmap = nnz

    call c_f_pointer(timestamps, sampletime, (/nsamp/) )
    call c_f_pointer(pix, pixels, (/nsamp, ndet/) )
    call c_f_pointer(pixweights, weights, (/nnz, nsamp, ndet/) )
    call c_f_pointer(signal, tod_stored, (/nsamp, ndet/) )

    subchunk_file_map = file_map; subchunk_file_binmap = file_binmap
    subchunk_file_base = file_base; subchunk_file_hit = file_hit; subchunk_file_mask = file_mask
    subchunk_file_matrix = file_matrix; subchunk_file_wcov = file_wcov

    tod_is_clean = .false.

    detflags = .true.

    if (nsubchunk > 1) then
       baseline_open = .false. ! from output.f90
       if (id == 0) then
          write (*,'(/," ********** SUBCHUNK == ",i0,/)') isubchunk
       endif
       call add_subchunk_id(file_map, isubchunk, nsubchunk)
       call add_subchunk_id(file_hit, isubchunk, nsubchunk)
       call add_subchunk_id(file_mask, isubchunk, nsubchunk)
       call add_subchunk_id(file_matrix, isubchunk, nsubchunk)
       call add_subchunk_id(file_wcov, isubchunk, nsubchunk)
       call add_subchunk_id(file_binmap, isubchunk, nsubchunk)
       call add_subchunk_id(file_base, isubchunk, nsubchunk)
    end if

    if (skip_existing) then
       if ( file_exists(file_map) ) then
          return
       end if
    end if

    call init_output

    call write_parameters
    call wait_mpi

    if (kfirst) then
       aa = 0.0
       yba = 0.0
       nna = 0
       nna_inv = 0
    end if

    call tic
    where ( pixels < 0 ) pixels = dummy_pixel
    if (id == 0) call toc('init_pointing')

    call wait_mpi

    cputime_init = cputime_init + get_time(1)

    call reset_time(1)

    call time_stamp

    subchunkpp => subchunk(1:nosamples_proc)

    call tic
    call reduce_pixels_buff
    if (id == 0) call toc('reduce_pixels_buff')

    call tic
    call update_maptod_transfer(ksubmap)
    if (id == 0) call toc('update_maptod_transfer')

    if ( reassign_submaps ) then
       call tic
       call assign_submaps( id_submap, nosubmaps, nopix_map, nopix_cross, nosubmaps_max )
       if (id == 0) call toc('assign_submaps')
    end if

    map = 0
    binmap = 0

    dummy_pixel = 12*nside_max**2

    if (kfirst) then
       call tic
       call flag_bad_baselines
       if (id == 0) call toc('flag_bad_baselines')
    end if

    call tic
    call bin_tod(map, binmap, wamap, tod_stored)
    if (id == 0) call toc('bin_tod')

    call time_stamp

    call tic
    call restore_pixels_buff
    if (id == 0) call toc('restore_pixels_buff')

    !  First destriping

    if (kfirst) then

       if (id == 0) then
          write(*,*)
          write(*,*) 'First destriping phase'
       endif
       call time_stamp

       call tic
       call reduce_pixels_a
       if (id == 0) call toc('reduce_pixels_a')

       call tic
       call initialize_a(yba, nna, wamap, cca, tod_stored)
       if (id == 0) call toc('initialize_a')

       call tic
       call iterate_a(aa, yba, nna, wamap, cca, tod_stored)
       if (id == 0) call toc('iterate_a')

       call tic
       call subtract_baselines_a(map, aa)
       if (id == 0) call toc('subtract_baselines_a')

       call tic
       call write_baselines_a(aa)
       if (id == 0) call toc('write_baselines_a')

    endif

    call wait_mpi
    if (id == 0) then
       write(*,*)
       write(*,*) 'Finalization begins'
       if (info.ge.2) write(*,*)
    endif
    call time_stamp
    call reset_time(1)

    !
    !  Binning of the final map
    !
    if (do_binmap) call makemaps(binmap, cc, outmask)
    call makemaps(map, cc, outmask)
    call map_analysis(map, outmask)
    !
    ! output files
    !
    call tic
    call write_map(map, outmask)
    if (id == 0) call toc('Write map')
    call tic
    call write_binmap(binmap, outmask)
    if (id == 0) call toc('Write binmap')
    call tic

    call run_subsets()

    if ( isubchunk == nsubchunk .and. kfirst ) then
       ! subtract the baselines to return the destriped TOD.
       ! If the baselines were already subtracted, the call has no effect.
       call tic
       call clean_tod(tod_stored, aa)
       if (id == 0) call toc('clean_tod')
    end if

    call restore_pixels_a

    call wait_mpi
    call time_stamp

    cputime_final = cputime_final + get_time(1)
    cputime_total = get_time(0)

    ! Memory

    if (info > 0) then
       if (id == 0) write(*,*)
       if (id == 0) write(*,*) 'MEMORY (MB):'

       call write_memory('Detector pointing',  memory_pointing)
       call write_memory('TOD buffer',         memory_tod)
       call write_memory('Maps',               memory_maps)
       call write_memory('Baselines',          memory_baselines)
       call write_memory('Basis functions',    memory_basis_functions)
       call write_memory('Noise filter',       memory_filter)
       call write_memory('Preconditioner',     memory_precond)
       call write_memory('Temporary maps',     memory_locmap)
       call write_memory('All2All buffers',    memory_all2all)
       call write_memory('CG work space',      memory_cg)
       call write_memory('NCM', memory_ncm)
       call write_memory('Total')
    endif

    ! Timing

    if (info > 0) then
       if (id == 0) write(*,*) 'WALL-CLOCK TIME (s):'

       call write_time('Initialization',     cputime_init)

       cputime_read = cputime_read_pointing_periods + cputime_read_detpointing + cputime_read_tod + cputime_read_timestamps

       call write_time('I/O', cputime_read)
       call write_time('- Pointing Period boundaries', cputime_read_pointing_periods)
       call write_time('- Pointing I/O',  cputime_read_detpointing)
       call write_time('- Time stamp I/O', cputime_read_timestamps)
       call write_time('- TOD I/O', cputime_read_tod)
       call write_time('- Other')

       call write_time('Waiting',                cputime_wait)
       call write_time('Building pixel matrices',cputime_build_matrix)
       call write_time('Sending pixel matrices',cputime_send_matrix)
       call write_time('Inverting pixel matrices',cputime_inv)
       call write_time('Binning TOD',            cputime_bin_maps)
       call write_time('Sending binned TOD',     cputime_send_maps)
       call write_time('Counting hits',          cputime_count_hits)
       call write_time('Building preconditioner',cputime_prec_construct)
       call write_time('Subtract/add baselines', cputime_clean_tod+cputime_unclean_tod)

       call write_time('Initialization (1. phase)', cputime_cga_init)
       call write_time('CG iteration',          cputime_cga)
       call write_time('- TOD - map',           cputime_cga_1)
       call write_time('- CG MPI',  cputime_cga_mpi)
       call write_time('- CG ccmultiply',  cputime_cga_cc)
       call write_time('- Map - TOD',           cputime_cga_2)
       call write_time('- Filtering',           cputime_filter)
       call write_time('- Preconditioning',     cputime_precond)
       call write_time('- Other')

       call write_time('Subset maps', cputime_subset )
       call write_time('- Flag subset', cputime_flag_subset )
       call write_time('- Write subset', cputime_write_subset )

       call write_time('Finalization and output',cputime_final)
       call write_time('Other',cputime_total-time_cum)
    endif

    call write_time('Total',cputime_total)
    if (id == 0) write(*,*)

    file_map = subchunk_file_map; file_binmap = subchunk_file_binmap
    file_base = subchunk_file_base; file_hit = subchunk_file_hit; file_mask = subchunk_file_mask
    file_matrix = subchunk_file_matrix; file_wcov = subchunk_file_wcov

    call reset_timers

    call close_mpi()

    return

  end subroutine destripe_with_cache


  subroutine clear_caches() bind( c, name='clear_caches' )

    if (.not. cached .and. id == 0) write (*,*) 'WARNING: Madam caches are already empty.'
    
    call free_baselines
    call free_maps
    call free_locmaps
    
    if (allocated(pntperiods)) deallocate( pntperiods, pntperiod_id, noba_short_pp )
    if (allocated(baselines_short)) deallocate( baselines_short, base_pntid_short, basis_functions )
    if (allocated(baselines_short_start)) deallocate( baselines_short_start, baselines_short_stop )
    if (allocated(prec_diag)) deallocate( prec_diag )
    if (allocated(bandprec)) deallocate( bandprec )

    call close_filter()
    call close_output()
    call close_pointing()
    call free_mask()
    call free_tod()

    cached = .false.

  end subroutine clear_caches


  subroutine run_subsets()

    type(detset_type) :: detset
    type(survey_type) :: survey
    integer(i8b) :: nhit_det, nhit_survey, idetset, idet, jdet, isurvey, i, j, nmap_save, ncc_save
    logical :: do_binmap_save, kfirst_save, temperature_only_save, concatenate_messages_save
    character(len=SLEN) :: file_binmap_save, file_hit_save, file_matrix_save, file_wcov_save, file_mask_save, subsetname
    character(len=SLEN) :: detsetname, surveyname
    integer(i4b) :: pass, npass

    if ( .not. bin_subsets ) return

    cputime_final = cputime_final + get_time_and_reset(1)

    ! bin and write all requested subset maps with full data destriping.

    kfirst_save = kfirst
    do_binmap_save = do_binmap
    file_binmap_save = file_binmap
    file_hit_save = file_hit
    file_matrix_save = file_matrix
    file_wcov_save = file_wcov
    file_mask_save = file_mask
    temperature_only_save = temperature_only
    nmap_save = nmap
    ncc_save = ncc
    concatenate_messages_save = concatenate_messages

    npass = 1
    if ( do_binmap .and. do_map ) npass = 2

    kfirst = .false.
    if ( do_map ) do_binmap = .true.

    loop_pass : do pass = 1,npass

       ! First pass produces the binned map, hit map, mask and the noise matrices
       ! then the baselines are subtracted and the second pass produces the destriped maps
       ! only if both binned and destriped maps are required are two passes performed

       if ( pass == npass .and. do_map ) then
          ! subtract the baselines
          call tic
          call clean_tod(tod_stored, aa)
          if (id == 0) call toc('clean_tod')
       end if

       loop_survey : do isurvey = 0,nsurvey

          !cputime_final = cputime_final + get_time_and_reset(1)
          call reset_time( 10 )

          if ( isurvey /= 0 ) then
             survey = surveys(isurvey)
             surveyname = survey%name

             nhit_survey = 0
             surveyflags = .false.
             do i = 1, nosamples_proc
                do j = 1,survey%nspan
                   if ( survey%starts(j) <= sampletime(i) .and. sampletime(i) <= survey%stops(j) ) then
                      surveyflags(i) = .true.
                   end if
                end do
             end do
          else
             surveyflags = .true.
             surveyname = surveys(0)%name
          end if

          nhit_survey = count(surveyflags)
          call sum_mpi(nhit_survey)

          cputime_flag_subset = cputime_flag_subset + get_time_and_reset(10)

          if ( nhit_survey == 0 ) cycle

          loop_detset : do idetset = 0,ndetset

             if ( idetset == 0 .and. isurvey == 0 ) cycle ! This case is already processed

             call reset_time( 10 )

             if ( idetset /= 0 ) then
                detset = detsets(idetset)
                detsetname = detset%name

                detflags = .false.
                do idet = 1,nodetectors
                   do jdet = 1,detset%ndet
                      if ( trim(detectors(idet)%name) == trim(detset%detectors(jdet)) ) then
                         detflags(idet) = .true.
                      end if
                   end do
                end do
                temperature_only = ( detset%nopol .and. .not. force_pol)
             else
                detflags = .true.
                detsetname = detsets(0)%name
                temperature_only = temperature_only_save
             end if

             nhit_det = count(detflags)
             call sum_mpi(nhit_det)

             ! FIXME: we could easily remove all detsets and surveys that do not overlap with present data.

             cputime_flag_subset = cputime_flag_subset + get_time_and_reset(10)

             if ( nhit_det == 0 ) cycle loop_detset

             if ( temperature_only ) then
                nmap = 1
                ncc = 1
                !if ( nmap_save /= nmap ) concatenate_messages = .false.
             else
                nmap = nmap_save
                ncc = ncc_save
                concatenate_messages = concatenate_messages_save
             end if

             if ( id == 0 ) then
                print *,''
                print *,' *** Processing detset=',trim(detsetname),', survey=',trim(surveyname)
                print *,''
             end if

             ! construct file names

             subsetname = trim(path_output) // trim(file_root)
             if ( len_trim(detsetname) /= 0 ) subsetname = trim(subsetname) // '_' // trim( detsetname )
             if ( len_trim(surveyname) /= 0 ) subsetname = trim(subsetname) // '_' // trim( surveyname )

             file_binmap = ''
             file_hit = ''
             file_matrix = ''
             file_wcov = ''
             file_mask = ''

             if ( pass == npass .and. do_map ) then
                file_binmap = trim(subsetname) // '_map.fits'
             else
                if ( do_binmap ) file_binmap = trim(subsetname) // '_bmap.fits'
             end if
             if ( do_hits )   file_hit    = trim(subsetname) // '_hmap.fits'
             if ( do_matrix ) file_matrix = trim(subsetname) // '_wcov_inv.fits'
             if ( do_wcov )   file_wcov   = trim(subsetname) // '_wcov.fits'
             if ( do_mask )   file_mask   = trim(subsetname) // '_mask.fits'

             if ( nsubchunk > 1 ) then
                call add_subchunk_id(file_hit, isubchunk, nsubchunk)
                call add_subchunk_id(file_mask, isubchunk, nsubchunk)
                call add_subchunk_id(file_matrix, isubchunk, nsubchunk)
                call add_subchunk_id(file_wcov, isubchunk, nsubchunk)
                call add_subchunk_id(file_binmap, isubchunk, nsubchunk)
             end if

             !if ( skip_existing ) then
             if ( file_exists(file_binmap) .and. .not. concatenate_binary ) file_binmap = ''
             if ( file_exists(file_hit) )    file_hit = ''
             if ( file_exists(file_matrix) ) file_matrix = ''
             if ( file_exists(file_wcov) )   file_wcov = ''
             if ( file_exists(file_mask) )   file_mask = ''
             !end if

             if ( len_trim(file_binmap) == 0 .and. len_trim(file_hit) == 0 &
                  .and. len_trim(file_matrix) == 0 .and. len_trim(file_wcov) == 0 &
                  .and. len_trim(file_mask) == 0 ) cycle

             cca = 0
             cc = 0
             binmap = 0
             nohits = 0
             outmask = 0
             crit = 0

             call tic

             ! This is where we actually build the observation matrices and bin the subset data

             ! These calls required some tricks to accomodate unpolarized subsets

             call reset_time(10)
             call pixel_matrix(cca(1:nmap,1:nmap,:), cc(1:nmap,1:nmap,:))
             cputime_write_subset = cputime_write_subset + get_time(10)

             if ( pass == 1 ) then
                call count_hits(nohits)

                call reset_time(10)
                call write_hits(nohits)
                cputime_write_subset = cputime_write_subset + get_time(10)
             end if

             call bin_tod(map(1:nmap,:), binmap(1:nmap,:), wamap(1:nmap,:), tod_stored)

             if ( pass == 1 ) then
                call reset_time(10)
                call write_matrix(cc(1:nmap,1:nmap,:))
                cputime_write_subset = cputime_write_subset + get_time(10)
             end if

             call invert_pixelmatrix_map(cc(1:nmap,1:nmap,:), outmask, crit)

             if ( pass == 1 ) then
                call reset_time(10)
                call write_matrix(cc(1:nmap,1:nmap,:), outmask)
                cputime_write_subset = cputime_write_subset + get_time(10)
             end if

             call makemaps(binmap(1:nmap,:), cc(1:nmap,1:nmap,:), outmask)
             call map_analysis(binmap(1:nmap,:), outmask)

             call reset_time(10)
             call write_binmap(binmap(1:nmap,:), outmask)
             if ( pass == 1 ) call write_mask(outmask, crit)                
             cputime_write_subset = cputime_write_subset + get_time(10)

             if (id == 0) call toc('subset map')

          end do loop_detset
       end do loop_survey

       if ( pass == npass .and. do_map .and. nsubchunk > 1 ) then
          ! add the baselines back
          call tic
          call unclean_tod(tod_stored, aa)
          if (id == 0) call toc('unclean_tod')
       end if

    end do loop_pass

    file_binmap = file_binmap_save
    file_hit = file_hit_save
    file_matrix = file_matrix_save
    file_wcov = file_wcov_save
    file_mask = file_mask_save

    kfirst = kfirst_save
    do_binmap = do_binmap_save
    temperature_only = temperature_only_save
    nmap = nmap_save
    ncc = ncc_save
    concatenate_messages = concatenate_messages_save

    cputime_subset = cputime_subset + get_time_and_reset(1)

  end subroutine run_subsets



  function file_exists( filename )
    character(len=*) :: filename
    logical :: file_exists

    if (len_trim(filename) == 0) then
       file_exists = .false.
    else
       inquire( file=trim(path_output) // trim(filename), exist=file_exists)
       if ( .not. file_exists ) inquire( file=filename, exist=file_exists)
    end if

    if (file_exists .and. id == 0) write (*,*) trim(filename) // ' exists!'
    if (.not. file_exists .and. id == 0) write (*,*) trim(path_output) // trim(filename) // ' does not exist!'
    if (.not. file_exists .and. id == 0) write (*,*) trim(filename) // ' does not exist!'

  end function file_exists



  subroutine add_subchunk_id(filename, subchunk_id, nsubchunk)
    character(len=*) :: filename
    integer(i2b) :: subchunk_id, nsubchunk

    integer :: i
    character(len=SLEN) :: stemp

    if (len_trim(filename) == 0 .or. subchunk_id == 0) return

    stemp = trim(filename)
    i = index(stemp, '.', .true.)
    if (i > 0) then
       if (nsubchunk < 10) then
          write (filename,'(a,"_sub",i1.1,"of",i1.1,a)') stemp(:i-1), subchunk_id, nsubchunk, trim(stemp(i:))
       else
          write (filename,'(a,"_sub",i2.2,"of",i2.2,a)') stemp(:i-1), subchunk_id, nsubchunk, trim(stemp(i:))
       end if
    else
       if (nsubchunk < 10) then
          write (filename,'(a,"_sub",i1.1,"of",i1.1)') trim(stemp), subchunk_id, nsubchunk
       else
          write (filename,'(a,"_sub",i2.2,"of",i2.2)') trim(stemp), subchunk_id, nsubchunk
       end if
    end if

  end subroutine add_subchunk_id



  subroutine test_mpi_speed(nside_submap_test, best_nside_submap)
    integer :: nside_submap_test, best_nside_submap

    real(dp), allocatable :: tempmap(:,:,:)
    real(dp), allocatable :: buffer(:,:)
    logical, allocatable :: ksubmap_table_test(:,:), ksubmap_test(:)
    logical, allocatable :: kbuffer(:)
    integer, allocatable :: id_submap_test(:)

    integer :: ierr, nosubpix_test, nosubmaps_test, nosubmaps_tot_test
    integer :: nosubpix_max_test, nosubmaps_max_test
    integer :: i, j, k, m, mrecv, id_tod, id_map, ndegrade
    integer :: id_send, msend, idet
    real(dp) :: time1, time2, time
    real(dp), save :: best_time=1e9

    if (id == 0) write (*,*) ' Running nside_submap test ...'

    ! initialize submap parameters for nside_submap_test
    nosubmaps_tot_test = 12*nside_submap_test**2
    if (id == 0) write (*,'(a,"==",i0)') 'nosubmaps_tot_test', nosubmaps_tot_test
    nosubmaps_max_test = (nosubmaps_tot_test-1)/ntasks+1
    nosubpix_test = nside_cross**2/nside_submap_test**2
    nosubpix_max_test = nside_max**2/nside_submap_test**2
    ndegrade = nosubpix_max_test/nosubpix_test

    ! simulate update_maptod_transfer to get ksubmap_table for the test
    allocate(ksubmap_test(0:nosubmaps_tot_test))
    ksubmap_test = .false.
    do idet = 1, nodetectors
       do i = 1,nosamples_proc
          ksubmap_test(pixels(i, idet)/nosubpix_max_test) = .true.
       end do
    end do

    allocate(ksubmap_table_test(0:nosubmaps_tot_test-1,0:ntasks-1), &
         kbuffer(0:nosubmaps_tot_test-1))
    do id_send = 0,ntasks-1
       if (ID == id_send) kbuffer = ksubmap_test(0:nosubmaps_tot_test-1)
       call broadcast_mpi(kbuffer,nosubmaps_tot_test,id_send)
       ksubmap_table_test(0:nosubmaps_tot_test-1,id_send) = kbuffer
    end do

    allocate(id_submap_test(0:nosubmaps_tot_test-1))
    nosubmaps_test = 0
    do i = 0,nosubmaps_tot_test-1
       k = mod(i,ntasks)
       id_submap_test(i) = k
       if (ID == k) nosubmaps_test = nosubmaps_test+1
    enddo

    ! Simulate collect_map
    call wait_mpi
    time1 = mpi_wtime()
    allocate(tempmap(nmap,nosubpix_test,nosubmaps_test), &
         buffer(nmap,nosubpix_test), stat=ierr)
    mrecv = 0
    m = 0
    loop_submap : do i = 0,nosubmaps_tot_test-1
       id_map = id_submap_test(i)
       if (ID == id_map) mrecv = mrecv+1

       loop_recv : do id_tod = 0,ntasks-1
          if (.not. ksubmap_table_test(i,id_tod)) cycle
          if (ID == id_tod) then  ! prepare send buffer
             buffer = 0.0
             loop_subpix : do k = 1,nosubpix_test
                do j = 1,ndegrade
                   buffer(:,k) = buffer(:,k) + 1.0 ! locmap(:,m)
                   m = m + 1
                end do
             end do loop_subpix
          end if

          call send_mpi_vec_dp(buffer,nmap*nosubpix_test,id_tod,id_map)
          if (ID == id_map) tempmap(:,:,mrecv) = tempmap(:,:,mrecv) + buffer
       end do loop_recv
    end do loop_submap
    !call wait_mpi
    time2 = mpi_wtime()
    if (id == 0) write (*,'(a,i0,a,f6.3)') &
         'collect_map@nside_submap == ',nside_submap_test,' : ',time2-time1

    ! simulate scatter_map
    msend = 0
    m     = 0
    do i = 0,nosubmaps_tot_test-1
       id_map = id_submap_test(i)
       if (id == id_map) then
          msend = msend+1
          buffer = 1.0 !map(:,:,msend)
       endif
       do id_tod = 0,ntasks-1
          if (.not.ksubmap_table_test(i,id_tod)) cycle
          call send_mpi_vec_dp(buffer,nmap*nosubpix_test,id_map,id_tod)
          if (id == id_tod) then
             do k = 1,nosubpix_test
                do j = 1,ndegrade
                   !locmap(:,m) = buffer(:,k)
                   buffer(:,k) = buffer(:,k) + 1.0
                   m = m + 1
                end do
             end do
          end if
       enddo
    enddo
    call wait_mpi
    time = mpi_wtime()
    if (id == 0) write (*,'(a,i0,a,f6.3)') &
         'scatter_map@nside_submap == ',nside_submap_test,' : ',time-time2

    time = time-time1
    if (id == 0 .and. time < best_time) then
       write (*,'(a,i0)') &
            ' -- New OPTIMAL nside_submap is ',nside_submap_test
       best_nside_submap = nside_submap_test
       best_time = time
    end if

    call broadcast_mpi(best_nside_submap, 0) ! synchronize nside_submap

    deallocate( ksubmap_table_test, id_submap_test, tempmap, buffer )

    if (id == 0) write (*,*) ' ... nside_submap test completed'

  end subroutine test_mpi_speed


  !---------------------------------------------------------------------------


  SUBROUTINE time_stamp()

    call wait_mpi

    if (id == 0) write(*,'(" Clock =",f10.3," s")') get_time(0)

  END SUBROUTINE time_stamp

  !---------------------------------------------------------------------------


  SUBROUTINE barrier

    call reset_time(5)

    call wait_mpi

    cputime_wait = cputime_wait +get_time(5)

  END SUBROUTINE barrier


  !---------------------------------------------------------------------------


  subroutine reset_timers()
    ! Reset all timers

    integer(i4b) :: i

    cputime_init = 0
    cputime_read_pointing_periods = 0
    cputime_read_detpointing = 0
    cputime_read_timestamps = 0
    cputime_read_tod = 0
    cputime_postprocess_tod = 0
    cputime_read = 0
    cputime_wait = 0
    cputime_build_matrix = 0
    cputime_send_matrix = 0
    cputime_inv = 0
    cputime_bin_maps = 0
    cputime_send_maps = 0
    cputime_count_hits = 0
    cputime_prec_construct = 0

    cputime_cga_init = 0
    cputime_cga = 0
    cputime_cga_1 = 0
    cputime_cga_mpi = 0
    cputime_cga_cc = 0
    cputime_cga_2 = 0

    cputime_clean_tod = 0
    cputime_unclean_tod = 0

    cputime_filter = 0
    cputime_precond = 0

    cputime_final = 0

    cputime_subset = 0
    cputime_flag_subset = 0
    cputime_write_subset = 0

    cputime_total = 0
    time_cum = 0
    time2_sum = 0

    memory_total = 0

    do i = 0,99
       call reset_time(i)
    end do

  end subroutine reset_timers

end module smadam
