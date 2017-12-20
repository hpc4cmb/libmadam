!
MODULE read_data

  ! Routines for reading input data:
  ! TOD, pointing, input mask

  use simulation
  use commonparam
  use submap_transfer
  use fitsmod2
  use pointing
  use maps_and_baselines
  !use coordinate_conversion
  use healpix_routines
  use satellite
  !use tod_storage
  use timing
  use mpi_wrappers

  implicit none
  private

  real(sp), public :: cputime_read_tod=0,  cputime_read_detpointing=0, &
       cputime_read_timestamps=0, cputime_postprocess_tod, memory_sat=0.0, &
       cputime_read_pointing_periods=0

  public init_read_data, check_files, read_inmask, close_read_data, &
       toast_read_pointing_and_tod, toast_read_mc_tod, flag_bad_baselines

CONTAINS


  subroutine toast_read_pointing_and_tod(pixels, weights, tod, time, tod_only, flush)
    !
    ! Step through all available data and concatenate it into the supplied
    ! arrays.
    !

    integer, allocatable :: pixels(:,:)
    real(dp), allocatable :: tod(:,:)
    real(dp), allocatable :: weights(:,:,:)
    real(dp), allocatable :: time(:)
    logical :: tod_only, flush

    stop 'I/O should now be handled externally'

  end subroutine toast_read_pointing_and_tod



  subroutine toast_read_mc_tod( tod )
    !
    ! Step through all available data and concatenate it into the supplied
    ! array.
    !
    real(dp), allocatable :: tod(:,:)

    stop 'I/O should now be handled externally'

  end subroutine toast_read_mc_tod



  !---------------------------------------------------------------------------


  subroutine flag_bad_baselines

    use maptod_transfer, only : locmask

    integer(i8b) :: n_flagged, first, last, idet, ibase, len, n_bad, n_good, ipix, i

    ! include only baselines that have more than good_baseline_fraction
    ! of data unflagged and unmasked

    if (good_baseline_fraction == 0 .or. kfilter) return

    n_flagged = 0
    n_bad = 0

    do idet = 1,nodetectors
       do ibase = 1, noba_short
          n_good = 0
          first = baselines_short_start(ibase)
          last = baselines_short_stop(ibase)
          do i = first, last
             if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
             ipix = pixels(i,idet)
             if (ipix == dummy_pixel) cycle
             if (use_inmask .and. locmask(ipix) == 0) cycle

             n_good = n_good + 1
          end do

          len = last - first + 1
          if (n_good > 0 .and. n_good < good_baseline_fraction * len) then
             pixels(first:last, idet) = dummy_pixel
             n_bad = n_bad + 1
             n_flagged = n_flagged + n_good
          end if
       end do
    end do

    call sum_mpi(n_flagged)
    call sum_mpi(n_bad)

    if (id == 0 .and. n_flagged > 0) write (*,'(a,i0,a,i0,a)') &
         'Flagged out additional ',n_flagged,' samples on ',n_bad,' poorly '&
         // 'sampled baselines'

  end subroutine flag_bad_baselines

  !---------------------------------------------------------------------------


  SUBROUTINE init_read_data

    return

  END SUBROUTINE init_read_data


  !-------------------------------------------------------------------------


  SUBROUTINE read_inmask(mask)

    integer,intent(out)  :: mask(nosubpix_cross,nosubmaps)
    integer              :: idr, i, m, n, isubmap, id_recv
    integer              :: nside_inmask,nosubpix_inmask
    integer(idp)         :: offset
    integer,allocatable  :: ibuffer(:)
    logical              :: ringflag
    character(len=80)    :: ordering
    type(fitshandle)     :: infile

    if (len_trim(file_inmask) == 0) return

    idr = 1
    if (ntasks == 1) idr = 0

    if (id == idr .and. info > 1) write(*,*) 'Reading mask...'

    if (id == idr) then
       call fits_open(infile, file_inmask, fits_readonly)
       call fits_get_key(infile, 'ORDERING', ordering)

       if (ordering == 'NESTED') then
          ringflag = .false.
       elseif (ordering == 'RING') then
          ringflag = .true.
          write(*,idf) ID,'Error in subroutine Read_inmask.'
          write(*,idf) ID,'Input mask should be in NESTED ordering.'
          call exit_with_status(1)
       else
          write(*,idf) ID,'Error in subroutine Read_inmask.'
          write(*,idf) ID,'Cannot determine the ordering (nested or ring).'
          write(*,idf) ID,ordering
          call exit_with_status(1)
       endif

       call fits_get_key(infile, 'NSIDE', nside_inmask)

       ! de/upgrade to resolution nside_cross
       if (info > 0) then
          if (nside_cross > nside_inmask) then
             write(*,'(a,i8)') ' Mask upgraded to resolution nside =', nside_cross
          elseif (nside_cross < nside_inmask) then
             write(*,'(a,i8)') ' Mask downgraded to resolution nside =', &
                  nside_cross
          end if
       endif

    endif

    call broadcast_mpi(nside_inmask, idr)

    ! up/downgrade
    if (nside_cross.ge.nside_inmask) then
       n = (nside_cross/nside_inmask)**2
    elseif (nside_cross.lt.nside_inmask) then
       n = (nside_inmask/nside_cross)**2
    endif

    nosubpix_inmask = (nside_inmask/nside_submap)**2
    allocate(ibuffer(nosubpix_inmask))

    offset = 0
    m = 0
    do isubmap = 0,nosubmaps_tot-1

       if (id == idr) then
          call fits_read_column(infile,1,ibuffer,offset)
          offset = offset+nosubpix_inmask
       endif

       id_recv = id_submap(isubmap)
       call send_mpi(ibuffer,nosubpix_inmask,idr,id_recv)

       if (id == id_recv) then

          m = m+1
          if (nside_cross.ge.nside_inmask) then

             do i = 1,nosubpix_inmask
                if (ibuffer(i).gt.0) then
                   mask((i-1)*n+1:i*n,m) = 1
                else
                   mask((i-1)*n+1:i*n,m) = 0
                endif
             enddo
          else

             do i = 1,nosubpix_cross
                if (any(ibuffer((i-1)*n+1:i*n).le.0)) then
                   mask(i,m) = 0
                else
                   mask(i,m) = 1
                endif
             enddo

          endif

       endif

    enddo

    if (id == idr) call fits_close(infile)

    deallocate(ibuffer)

    n = count(mask == 0)
    call sum_mpi(n)

    if (id == idr .and. info > 2) write(*,'(i9,a)') n,' pixels masked out'

  END SUBROUTINE read_inmask

  !-------------------------------------------------------------------------


  SUBROUTINE close_read_data()

  END SUBROUTINE close_read_data


  !--------------------------------------------------------------------------


  SUBROUTINE check_files( nperiod, periods, nsamp )

    integer(c_long), intent(in) :: nperiod
    integer(c_long), intent(in) :: periods(nperiod)
    integer(c_long), intent(in) :: nsamp

    integer :: ierr, isend
    integer(i8b) :: ninterval_max, n, len
    character(len=256) :: name

    integer :: iperiod, idet
    integer, parameter :: N_PERIOD_MAX=1e6
    integer(i8b) :: i, days, hours, minutes
    integer(i8b) :: ngood, nflag, nzeroweight
    real(dp) :: seconds
    real(c_double) :: overlap

    if (id == 0 .and. info > 1) write (*,'(/,a,/)') 'Examining periods'

    call reset_time(12)

    ninterval = nperiod
    allocate(intervals(ninterval), stat=ierr)
    if (ierr /= 0) call abort_mpi('Failed to allocate intervals')

    ! count the periods and their lengths

    if (nperiod > N_PERIOD_MAX) call abort_mpi('Please increase N_PERIOD_MAX')

    do iperiod = 1, nperiod
       if (iperiod < nperiod) then
          len = periods(iperiod+1) - periods(iperiod)
       else
          len = nsamp - periods(iperiod)
       end if
       if (basis_order > 0 .and. len < basis_order + 1) then
          ! The baseline hierarchy is degenerate for this period.  We could
          ! try to handle the issue through flagging but likely the inputs
          ! are wrong.
          print *, id, ' : ERROR : period # ', iperiod, ' is only ', len, &
               ' samples long but basis_order = ', basis_order
          call abort_mpi('Too short period')
       end if
       intervals(iperiod) = len
    end do

    ! Count the flagged and unflagged samples
    nzeroweight = 0
    do idet = 1,nodetectors
       do i = 1, nsamp
          if (all(weights(:, i, idet) == 0) .and. &
               .not. pixels(i, idet) < 0) then
             nzeroweight = nzeroweight + 1
             pixels(i,idet) = -1
          end if
       end do
    end do
    call sum_mpi(nzeroweight)

    ngood = 0
    nflag = 0
    do idet = 1,nodetectors
       do i = 1, nsamp
          if (pixels(i, idet) < 0) then
             nflag = nflag + 1
          else
             ngood = ngood + 1
          end if
       end do
    end do
    call sum_mpi(ngood)
    call sum_mpi(nflag)

    nosamples_proc = nsamp
    nosamples_proc_max = nosamples_proc
    call max_mpi(nosamples_proc_max)

    ninterval_max = ninterval
    call max_mpi(ninterval_max)

    ninterval_tot = ninterval
    call sum_mpi(ninterval_tot)

    nosamples_tot = nosamples_proc
    call sum_mpi(nosamples_tot)

    if (id == 0 .and. info > 0) then
       write (*,'(a,i0)') 'Total number of samples (single detector): ', &
            nosamples_tot
       write (*,'(a,f7.3," %")') 'Zero-weight fraction (all detectors): ', &
            nzeroweight*100. / (ngood+nflag)
       write (*,'(a,f7.3," %")') 'Flagged fraction (all detectors): ', &
            nflag*100. / (ngood+nflag)
       write (*,'(a,i0)') 'Total number of intervals: ', ninterval_tot
       write (*,'(a,i0)') 'Max number of samples per task: ', nosamples_proc_max
    end if

    cputime_read_pointing_periods = cputime_read_pointing_periods &
         + get_time_and_reset(12)

  END SUBROUTINE check_files


  !--------------------------------------------------------------------------


  SUBROUTINE get_pointing_periods

  END SUBROUTINE get_pointing_periods


  !-----------------------------------------------------------------


  SUBROUTINE check_todfile(ok,filename_in,key)

    logical :: ok
    character(len=*) :: filename_in(nofiles)
    character(len=*) :: key

  END SUBROUTINE check_todfile


  !-------------------------------------------------------------------------

END MODULE read_data
