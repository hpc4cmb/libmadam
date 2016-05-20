MODULE inputparam
  ! Read input parameters from file

  use commonparam
  use simulation
  use mpi_wrappers

  implicit none
  private

  integer, parameter :: MAXNODET = 9999

  character(len=40),allocatable :: detc(:)

  integer :: nodetc = 0

  public read_parameters, read_detectors

CONTAINS


  !--------------------------------------------------------------------------


  SUBROUTINE read_parameters( parstring )

    ! Completely rewritten to parse the parameter string representation

    character(kind=c_char), intent(in) :: parstring(*)

    integer :: n, istart, istop, i
    character(len=SLEN) :: line

    detsets(0)%name = ''
    surveys(0)%name = ''

    n = 1
    do while ( parstring(n) /= C_NULL_CHAR )
       n = n + 1
    end do
    
    istart = 1
    do while ( istart < n )
       
       ! Isolate the next definition. Using scan intrinsic does not work.
       
       istop = istart
       line = ''
       do while( parstring(istop) /= ';' )
          i = istop - istart + 1
          line(i:i) = parstring(istop)
          istop = istop + 1
          if ( istop == n ) exit
       end do

       call read_line( trim( adjustl(line) ) )
       
       istart = istop + 1       
       
    end do

    dnshort = nint(dnshort * fsample) ! baseline length is in seconds

    call check_missing_parameters

  END SUBROUTINE read_parameters


  !---------------------------------------------------------------------------


  SUBROUTINE read_detectors( detstring, ndet, detweights, &
       npsd, npsdtot, psdstarts, npsdbin, psdfreqs, npsdval, psdvals )

    ! Completely rewritten to parse the parameter string representation

    character(kind=c_char), intent(in) :: detstring(*)
    integer(c_long), intent(in) :: ndet
    real(c_double), intent(in) :: detweights(ndet)

    integer(c_long), intent(in) :: npsd(ndet)
    integer(c_long), intent(in), value :: npsdtot
    real(c_double), intent(in) :: psdstarts(npsdtot)
    integer(c_long), intent(in), value :: npsdbin
    real(c_double), intent(in) :: psdfreqs(npsdbin)
    integer(c_long), intent(in), value :: npsdval
    real(c_double), intent(in) :: psdvals(npsdval)

    integer :: n, istart, istop, i, idet, nint, psdoffset, ierr, ipsd
    real(dp) :: detweight
    character(len=SLEN) :: line

    nodetectors = ndet
    allocate( detectors(nodetectors) )

    n = 1
    do while ( detstring(n) /= C_NULL_CHAR )
       n = n + 1
    end do
    
    istart = 1
    psdoffset = 1
    ipsd = 1
    idet = 0
    do while ( istart < n )
       
       ! Isolate the next definition. Using scan intrinsic does not work.
       
       istop = istart
       line = ''
       do while( detstring(istop) /= ';' )
          i = istop - istart + 1
          line(i:i) = detstring(istop)
          istop = istop + 1
          if ( istop == n ) exit
       end do

       idet = idet + 1
       if ( idet > ndet ) stop 'Too many detectors'
       
       detectors(idet)%name = trim(adjustl(line))
       detectors(idet)%idet = idet

       ! first noise weights

       nint = npsd( idet )
       if ( nint == 0 ) nint = 1
       detectors(idet)%npsd = nint
       allocate( detectors(idet)%weights(nint), detectors(idet)%sigmas(nint), detectors(idet)%psdstarts(nint) )
       detweight = detweights( idet )
       detectors(idet)%weights = detweight
       if ( detweight > 0._dp ) then
          detectors(idet)%sigmas = 1._dp / SQRT( detweight )
       else
          detectors(idet)%sigmas = 0._dp
       end if
       detectors(idet)%psdstarts = 0._dp
       detectors(idet)%kpolar = ( nmap /= 1 )

       ! Then the noise PSDs

       nint = npsd( idet )
       if ( nint > 0 ) then
          allocate( detectors(idet)%psdfreqs(npsdbin), detectors(idet)%psds(npsdbin,nint), stat=ierr )
          if ( ierr /= 0 ) stop 'No room for detector PSDs'
          detectors(idet)%psdfreqs = psdfreqs
          do i = 1,nint
             detectors(idet)%psdstarts(i) = psdstarts(ipsd)
             detectors(idet)%psds(:,i) = psdvals(psdoffset:psdoffset+npsdbin-1)
             psdoffset = psdoffset + npsdbin
             ipsd = ipsd + 1
          end do
       end if
       
       istart = istop + 1
       
    end do

  END SUBROUTINE read_detectors


  !---------------------------------------------------------------------------


  SUBROUTINE read_line(line)
    !
    ! Extract one line from the input file
    character(len=*) :: line
    character(len=40) :: key
    character(len=200) :: value
    integer :: i, ierr

    if (len_trim(line) == 0) return

    i = scan(line,'=')
    if (line(1:1) == '#') return
    if (i == 0) call abort_mpi('Syntax error in parameter definition: ' // trim(line))

    key = trim( adjustl( line(:i-1) ) )
    value = trim( adjustl( line(i+1:) ) )

    i = scan(value, '#')
    if (i > 0) value=value(:i-1)

    ! Special cases

    if (key(1:8) == 'detector') key = 'detector'
    if (len_trim(value) == 0 .and. key(1:4) /= 'file') return

    ! Keyword parsing

    ierr = 0

    select case (key)
    case ('info')
       read(value,*,iostat=ierr) info
    case ('runfile')
       file_simulation = trim(value)
    case ('read_buffer_len')
       read(value,*,iostat=ierr) read_buffer_len
    case ('nthreads')
       read(value,*,iostat=ierr) nthreads
    case ('nsubchunk')
       read(value,*,iostat=ierr) nsubchunk
    case ('isubchunk')
       read(value,*,iostat=ierr) isubchunk
    case ('time_unit')
       if (value(1:1) == 'p' .or. value(1:1) == 'P') then
          time_unit = -5.d0 ! pointing periods
       else
          read(value,*,iostat=ierr) time_unit
       end if
    case ('base_first','nshort')
       read(value,*,iostat=ierr) dnshort
    case ('nside_map')
       read(value,*,iostat=ierr) nside_map
    case ('nside_cross')
       read(value,*,iostat=ierr) nside_cross
    case ('nside_submap')
       read(value,*,iostat=ierr) nside_submap
    case ('nread_concurrent')
       read(value,*,iostat=ierr) nread_concurrent
    case ('good_baseline_fraction')
       read(value,*,iostat=ierr) good_baseline_fraction
    case ('run_submap_test')
       read(value,*,iostat=ierr) run_submap_test
    case ('concatenate_messages')
       read(value,*,iostat=ierr) concatenate_messages
    case ('reassign_submaps')
       read(value,*,iostat=ierr) reassign_submaps
    case ('pixmode_map')
       read(value,*,iostat=ierr) pixmode_map
    case ('pixmode_cross')
       read(value,*,iostat=ierr) pixmode_cross
    case ('pixlim_map')
       read(value,*,iostat=ierr) pixlim_map
    case ('pixlim_cross')
       read(value,*,iostat=ierr) pixlim_cross
    case ('kfirst')
       read(value,*,iostat=ierr) kfirst
    case ('basis_func')
       select case (value(1:4))
       case ('poly')
          basis_func = basis_poly
       case ('four')
          basis_func = basis_fourier
       case ('cheb')
          basis_func = basis_cheby
       case ('lege')
          basis_func = basis_legendre
       case default
          call abort_mpi('Unknown function basis : ' // trim(value))
       end select
    case ('basis_order')
       read(value,*,iostat=ierr) basis_order
       if (basis_order < 0) call abort_mpi('basis order must be nonnegative')
    case ('filter_mean')
       read(value,*,iostat=ierr) filter_mean
    case ('iter_min')
       read(value,*,iostat=ierr) iter_min
    case ('iter_max')
       read(value,*,iostat=ierr) iter_max
    case ('cglimit')
       read(value,*,iostat=ierr) cglimit
    case ('fsample')
       read(value,*,iostat=ierr) fsample
    case ('initialize_by_regression')
       read(value,*,iostat=ierr) initialize_by_regression
    case ('mode_detweight')
       read(value,*,iostat=ierr) mode_detweight    
    case ('flag_by_horn')
       read(value,*,iostat=ierr) flag_by_horn
    case ('dist_by_obs')
       read(value,*,iostat=ierr) dist_by_obs
    case ('write_cut')
       read(value,*,iostat=ierr) write_cut
    case ('checknan')
       read(value,*,iostat=ierr) checknan
    case ('sync_output')
       read(value,*,iostat=ierr) sync_output
    case ('skip_existing')
       read(value,*,iostat=ierr) skip_existing
    case ('temperature_only')
       read(value,*,iostat=ierr) temperature_only
    case ('force_pol')
       read(value,*,iostat=ierr) force_pol
    case ('noise_weights_from_psd')
       read(value,*,iostat=ierr) noise_weights_from_psd
    case ('radiometers')
       read(value,*,iostat=ierr) radiometers
    case ('psdlen')
       read(value,*,iostat=ierr) psdlen
    case ('psd_down','psd_downsample')
       read(value,*,iostat=ierr) psd_downsample
    case ('kfilter')
       read(value,*,iostat=ierr) kfilter
    case ('precond_width')
       read(value,*,iostat=ierr) precond_width
    case ('filter_time')
       read(value,*,iostat=ierr) filter_time
    case ('tail_time')
       read(value,*,iostat=ierr) tail_time
    case ('rm_monopole')
       read(value,*,iostat=ierr) rm_monopole

    case ('path_output')
       path_output = trim(value)
    case ('file_root')
       file_root = trim(value)
    case ('write_map')
       read(value,*,iostat=ierr) do_map
    case ('write_binmap')
       read(value,*,iostat=ierr) do_binmap
    case ('write_hits')
       read(value,*,iostat=ierr) do_hits
    case ('write_matrix')
       read(value,*,iostat=ierr) do_matrix
    case ('write_wcov')
       read(value,*,iostat=ierr) do_wcov
    case ('write_base')
       read(value,*,iostat=ierr) do_base
    case ('write_mask')
       read(value,*,iostat=ierr) do_mask
    case ('unit_tod')
       unit_tod = trim(value)

    case ('file_gap_out')
       file_gap_out = trim(value)
    case ('file_mc')
       file_mc = trim(value)
    case ('write_tod')
       read(value,*,iostat=ierr) write_tod
    case ('file_inmask')
       file_inmask = trim(value)
    case ('file_spectrum')
       file_spectrum = trim(value)
    case ('file_gap')
       file_gap = trim(value)
    case ('file_pntperiod')
       file_pntperiod = trim(value)
    case ('file_objectsize')
       file_objectsize = trim(value)
    case ('file_fpdb_supplement')
       file_fpdb_supplement = trim(value)
    case('binary_output')
       read(value,*,iostat=ierr) binary_output
    case('nwrite_binary')
       read(value,*,iostat=ierr) nwrite_binary
    !NCVM
    case ('kcompress_pixels')
       read(value,*,iostat=ierr) kcompress_pixels
    case ('ksplit_covmat')
       read(value,*,iostat=ierr) ksplit_covmat
    case ('bfinvert')
       read(value,*,iostat=ierr) bfinvert
    case ('file_covmat')
       file_covmat = trim(value)

    case('detset')
       call parse_detset(value)
    case('detset_nopol')
       call parse_detset(value,nopol=.true.)
    case('survey')
       call parse_survey(value)
    case('bin_subsets')
       read(value,*,iostat=ierr) bin_subsets

    case default
       call abort_mpi('Unknown parameter: ' // trim(key))
    end select

    if (ierr /= 0) call abort_mpi('Failed to parse ' // trim(key) // ' from ' // trim(value))


  END SUBROUTINE read_line



  subroutine parse_detset(line, nopol)
    character(len=*) :: line
    logical, optional :: nopol
    integer(i4b) :: i, ndet
    character(len=SLEN) :: detname

    i = index( line, ':' )
    if ( i == 0 ) call abort_mpi('Failed to parse: ' // trim(line) // ' for valid detector set name' )

    ndetset = ndetset + 1
    if ( ndetset > NDETSETMAX ) call abort_mpi('Number of detector sets exceeds NDETSETMAX')
    detsets( ndetset )%name = trim( adjustl( line(:i-1) ) )
    if ( len(trim(detsets( ndetset )%name)) == 0 ) call abort_mpi('Empty detset name')
    if ( id == 0 ) print *,'Adding detset = ',trim(detsets(ndetset)%name) ! debug
    line = trim(adjustl(line(i+1:)))

    ndet = 0
    do
       i = index( line, ',' )
       ndet = ndet + 1
       if ( i /= 0 ) then
          detname = trim( adjustl( line(:i-1) ) )
          line = trim(adjustl(line(i+1:)))
       else
          detname = trim( adjustl( line ) )
       end if
       if ( len_trim(detname) == 0 .or. detname == 'all' .or. detname == 'All' .or. detname == 'ALL' ) then
          ! allow user to name the full detector set
          detsets(0)%name = detsets( ndetset )%name
          ndetset = ndetset - 1
          if ( id == 0 ) print *,'  Adding detector ALL' ! debug
          return
       end if
       if ( id == 0 ) print *,'  Adding detector ' // trim(detname) ! debug
       detsets(ndetset)%detectors(ndet) = detname
       if ( i == 0 ) exit
    end do

    detsets(ndetset)%ndet = ndet    
    if ( present( nopol ) ) then
       detsets(ndetset)%nopol = nopol
    else
       detsets(ndetset)%nopol = .false.
    end if

  end subroutine parse_detset



  subroutine parse_survey(line)
    character(len=*) :: line
    integer(i4b) :: i, nspan, ierr
    real(dp) :: sstart, sstop

    i = index( line, ':' )
    if ( i == 0 ) call abort_mpi('Failed to parse: ' // trim(line) // ' for valid survey name' )

    nsurvey = nsurvey + 1
    if ( nsurvey > NSURVEYMAX ) call abort_mpi('Number of detector sets exceeds NSURVEYMAX')
    surveys( nsurvey )%name = trim( adjustl( line(:i-1) ) )
    if ( len(trim(surveys( nsurvey )%name)) == 0 ) call abort_mpi('Empty survey name')
    if ( id == 0 ) print *,'Adding survey = ',trim(surveys(nsurvey)%name) ! debug
    line = trim(adjustl(line(i+1:)))

    nspan = 0
    do
       nspan = nspan + 1
       i = index( line, '-' )
       if ( i == 0 ) then
          if ( len_trim(line) == 0 .or. line == 'all' .or. line == 'All' .or. line == 'ALL') then
             ! allow user to name the full data span
             surveys(0)%name = surveys( nsurvey )%name
             nsurvey = nsurvey - 1
             if ( id == 0 ) print *,'  Adding span ALL' ! debug
             return
          end if
          call abort_mpi('Failed to parse: ' // trim(line) // ' for valid survey span marker (-)' )
       end if
       read( line(:i-1), *, iostat=ierr ) sstart
       if ( ierr /= 0 ) call abort_mpi('Failed to parse: ' // trim(line) // ' for valid survey start' )
       line = trim(adjustl(line(i+1:)))
       read( line, *, iostat=ierr ) sstop
       if ( ierr /= 0 ) call abort_mpi('Failed to parse: ' // trim(line) // ' for valid survey stop' )
       
       surveys(nsurvey)%starts(nspan) = sstart
       surveys(nsurvey)%stops(nspan) = sstop
       if ( id == 0 ) print *,'  Adding span ',sstart,' - ',sstop ! debug
       i = index( line, ',' )
       if ( i /= 0 ) then
          line = trim(adjustl(line(i+1:)))
       else
          exit
       end if
    end do

    surveys(nsurvey)%nspan = nspan

  end subroutine parse_survey


  !-------------------------------------------------------------------------


  SUBROUTINE check_missing_parameters

    ! disabled obsolete checks -RK

    ! By default use a very long baseline that is truncated to fit each interval
    if ( kfirst .and. dnshort < 0 ) dnshort = 1e10

  END SUBROUTINE check_missing_parameters


  !--------------------------------------------------------------------------


  SUBROUTINE check_includes
    ! This routine is completely rewritten -RK

    ! dropped ability to use same pointing
    ! for multiple detectors

    integer :: idet, idet2
    logical :: found

  END SUBROUTINE check_includes

  !--------------------------------------------------------------------------



END MODULE inputparam
