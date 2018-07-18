MODULE noise_routines
  ! Build and apply a noise filter

  use, intrinsic :: iso_c_binding

  use commonparam
  use fourier
  use mpi_wrappers
  use timing

  use tod_storage, only : sampletime

  implicit none
  public

  include "fftw3.f03"

  integer, save :: nof = -1, nofilter = 0, notail = 0

  real(sp), save, public :: memory_filter = 0
  real(sp), save, public :: memory_precond = 0
  real(sp), save, public :: cputime_filter = 0
  real(sp), save, public :: cputime_precond = 0
  real(sp), save, public :: cputime_prec_construct = 0

  real(dp), allocatable :: xx(:)
  complex(dp), allocatable :: fx(:)
  complex(dp), allocatable, target :: fcov(:, :), fprecond(:, :)

  !Store interpolated spectra
  real(dp),allocatable :: freqtab(:), spectab(:, :)

  !Preconditioner

  integer :: nband
  logical :: use_diagonal = .false.
  real(dp), allocatable, target :: bandprec(:, :, :)
  real(dp), allocatable, target :: prec_diag(:, :)

  ! PSD information

  integer :: npsdtot
  integer, allocatable :: psddet(:), psdind(:)
  real(dp), allocatable :: psdstart(:), psdstop(:)

  public cinvmul, build_filter, close_filter, &
       construct_preconditioner, preconditioning_band, interpolate_psd

CONTAINS

  subroutine measure_noise_weights(radiometers)
    !
    ! Determine noise weights and split the PSDs into
    ! baselines and "white" (small scale) noise
    !
    logical, intent(in) :: radiometers
    real(dp) :: fstep, fbase, flower, fupper, sigma, rms, plateau
    integer :: nn, ipsd, nbin, ilower, iupper, ibase, idet, i
    real(dp), allocatable :: freqs(:), data(:)

    if (radiometers) then
       ! well-behaved PSD, just get the last PSD bin value
       nbin = 1
       allocate(freqs(nbin), data(nbin))
       freqs = fsample / 2

       do idet = 1, nodetectors
          do ipsd = 1, detectors(idet)%npsd
             call interpolate_psd(detectors(idet)%psdfreqs, &
                  detectors(idet)%psds(:,ipsd), freqs, data)
             plateau = data(1)
             rms = sqrt(plateau * fsample)
             detectors(idet)%sigmas(ipsd) = rms
             detectors(idet)%plateaus(ipsd) = plateau
          end do
       end do
       deallocate(freqs, data)
    else
       ! Measure the plateau value, VERY Planck Specific but will work
       ! for white noise filters
       nbin = 10000
       fstep = fsample * .5 / nbin
       allocate(freqs(nbin), data(nbin))
       freqs = (/ (fstep * i, i=1, nbin) /)
       fbase = 1 / (dnshort / fsample) ! baseline frequency
       ibase = fbase / fstep
       flower = 1 ! Plateau lower limit
       fupper = 10 ! Plateau upper limit
       ilower = flower / fstep
       iupper = fupper / fstep

       do idet = 1, nodetectors
          do ipsd = 1, detectors(idet)%npsd
             call interpolate_psd(detectors(idet)%psdfreqs, &
                  detectors(idet)%psds(:, ipsd), freqs, data)
             plateau = minval(data(ilower:iupper))
             ! integrate the RMS from the residual PSD:
             ! 1. constant PSD up to fbase
             rms = ibase * plateau
             ! 2. then integrate everything left
             do i = ibase + 1, nbin
                rms = rms + data(i)
             end do
             rms = sqrt(rms / nbin * fsample)
             detectors(idet)%sigmas(ipsd) = rms
             detectors(idet)%plateaus(ipsd) = plateau
          end do
       end do
       deallocate(freqs, data)
    end if

  end subroutine measure_noise_weights


  subroutine interpolate_psd(freq, psd, targetfreq, targetpsd)

    ! Use linear interpolation in the log-log plane

    real(dp), intent(in) :: freq(:), psd(:), targetfreq(:)
    real(dp), intent(out) :: targetpsd(:)

    real(dp), allocatable :: x(:), y(:)
    real(dp) :: flog, r

    integer :: n_in, n_out, i, ierr, j, startbin

    n_in = size(freq)
    n_out = size(targetfreq)

    startbin = 1
    if (freq(1) == 0) then
       n_in = n_in -1
       startbin = 2
    end if

    allocate(x(n_in), y(n_in), stat=ierr)
    if (ierr /= 0) stop 'No room to interpolate the PSD'

    x = log(freq(startbin:))
    y = log(psd(startbin:))

    j = 1
    do i = 1, n_out
       if (targetfreq(i) == 0) then
          targetpsd(i) = 0
          cycle
       end if
       flog = log(targetfreq(i))
       if (j < n_in - 1) then
          do while (x(j+1) < flog)
             j = j + 1
             if (j == n_in - 1) exit
          end do
       end if
       r = (flog - x(j)) / (x(j+1) - x(j))
       targetpsd(i) = exp(y(j) * (1._dp - r) + y(j + 1) * r)
    end do

    deallocate(x, y)

  end subroutine interpolate_psd


  function psd_index(idet, starttime) result(ipsd)
    !
    ! Returns the global PSD index (accross all detectors)
    !
    integer :: idet, ipsd
    real(dp) :: starttime

    do ipsd = 1, npsdtot
       if (psddet(ipsd) == idet .and. psdstart(ipsd) <= starttime &
            .and. psdstop(ipsd) >= starttime) return
    end do

    ipsd = -1
    return

  end function psd_index


  function psd_index_det(idet, starttime) result(ipsd)
    !
    ! Returns the PSD index for this detector
    !
    integer :: idet, ipsd
    real(dp) :: starttime
    integer :: ipsdtot

    if (.not. allocated(detectors(idet)%psds) &
         .and. detectors(idet)%npsd == 1) then
       ipsd = 1
       return
    end if

    ipsd = 1
    do
       if (ipsd == detectors(idet)%npsd) then
          if (detectors(idet)%psdstarts(ipsd) > starttime) ipsd = -1
          return
       end if
       if (detectors(idet)%psdstarts(ipsd+1) > starttime) exit
       ipsd = ipsd + 1
    end do

    return

  end function psd_index_det


  !-----------------------------------------------------------------------


  SUBROUTINE init_filter

    ! Stock version of Madam now has filter_time and tail_time hardcoded to their
    ! default values. We leave them free but retain matching default values

    integer :: nof_min, ierr

    if (.not. kfilter) return

    nofilter = (filter_time*fsample+.5) / nshort
    notail   = (tail_time*fsample+.5) / nshort

    !if (nofilter > noba_short_max) nofilter = noba_short_max
    if (notail > nofilter) notail = nofilter

    nof_min = nofilter + 2*notail
    ! Ensure that the filter is long enough for the preconditioner
    if (precond_width > nof_min-1) nof_min = precond_width+1

    nof = 1
    do
       if (nof >= nof_min) exit
       nof = 2*nof
    enddo
    if (ID == 0) write(*,*) 'FFT length =', nof

    allocate(fx(nof/2+1), xx(nof), stat=ierr) ! Work space
    if (ierr /= 0) stop 'No room for fx and xx'

    memory_filter = (nof/2+1)*16. + nof*8

    call init_fourier(nof)

  END SUBROUTINE init_filter


  !-----------------------------------------------------------------------


  SUBROUTINE close_filter

    if (.not. kfilter) return

    if (allocated(fcov)) deallocate(fcov)
    if (allocated(fprecond)) deallocate(fprecond)
    if (allocated(psddet)) deallocate(psddet)
    if (allocated(psdind)) deallocate(psdind)
    if (allocated(psdstart)) deallocate(psdstart)
    if (allocated(psdstop)) deallocate(psdstop)
    if (allocated(fx)) deallocate(fx)
    if (allocated(xx)) deallocate(xx)

    call close_fourier

    nof = -1

  END SUBROUTINE close_filter


  !-----------------------------------------------------------------------


  SUBROUTINE build_filter()
    ! Build the noise filter
    integer, parameter :: kmax = 2
    integer :: i, k, idet, nolines, nocols
    real(dp) :: fa, x
    real(dp) :: slope, fmin, fknee, sigma
    logical :: kread_file
    real(dp), allocatable :: aspec(:), f(:), g(:), spectrum(:), &
         ftable(:), spectrum_table(:, :)
    integer :: ipsd, ipsdtot, ierr

    if (.not.(kfilter)) return

    if (info == 3 .and. ID == 0) write(*,*) 'Building filter...'
    if (info > 4) write(*,*) 'Building filter...'

    if (nof < 0) call init_filter

    npsdtot = 0
    do idet = 1, nodetectors
       npsdtot = npsdtot + detectors(idet)%npsd
    end do

    ! Construct and store the filters
    if (allocated(fcov)) deallocate(fcov)
    allocate(fcov(nof/2+1, npsdtot), & ! Fourier inverse of the noise filter
         psddet(npsdtot), psdind(npsdtot), psdstart(npsdtot), psdstop(npsdtot), &
         stat=ierr)
    if (ierr /= 0) stop 'No room for fcov'
    memory_filter = memory_filter + (nof/2+1)*npsdtot*16. + npsdtot*24.

    allocate(aspec(nof/2+1), spectrum(nof), f(nof), g(nof), stat=ierr)
    if (ierr /= 0) stop 'No room for aspec'

    kread_file = (len_trim(file_spectrum) > 0)
    if (kread_file) call read_spectrum

    fa = fsample/nshort

    ipsdtot = 0
    do idet = 1, nodetectors
       if (.not. allocated(detectors(idet)%psds)) then
          stop 'All detectors must have at least one PSD to build the filter'
       end if
       do ipsd = 1, detectors(idet)%npsd
          ipsdtot = ipsdtot + 1

          psddet(ipsdtot) = idet
          psdind(ipsdtot) = ipsd
          psdstart(ipsdtot) = detectors(idet)%psdstarts(ipsd)
          if (ipsd < detectors(idet)%npsd) then
             psdstop(ipsdtot) = detectors(idet)%psdstarts(ipsd+1)
          else
             psdstop(ipsdtot) = 1e30
          end if

          aspec = 0
          do k = 0, kmax

             do i = 1, nof
                f(i) = k*fa +(i-1)*fa/nof
             enddo

             if (kread_file) then
                call get_spectrum_file
             else
                !call get_spectrum_powerlaw
                !call get_spectrum_toast(idet, ipsd)
                call get_spectrum_interp(idet, ipsd)
             endif

             do i = 1, nof
                x = pi*f(i)/fa
                if (x < 1.e-30) then
                   g(i) = 1.d0
                else
                   g(i) = sin(x)**2/x**2
                end if
             end do

             if (k == 0) then
                aspec(1) = spectrum(1)
             else
                aspec(1) = aspec(1) + 2*spectrum(1)*g(1)
             end if

             do i = 2, nof/2+1
                aspec(i) = aspec(i) + spectrum(i)*g(i) &
                     + spectrum(nof-i+2)*g(nof-i+2)
             end do

          end do

          aspec = 1.d0 / (fa*aspec)
          ! Madam/TOAST edit: do not constrain baseline mean
          if (.not. filter_mean) aspec(1) = 0
          fcov(:, ipsdtot) = aspec

          ! RK edit: check for NaNs
          do i = 1, nof/2+1
             if (isnan(abs(fcov(i, ipsdtot)))) &
                  call abort_mpi('NaN in fcov. Check input PSD and sigma.')
          end do

       end do ! ipsd
    end do ! idet

    deallocate(aspec, f, g, spectrum)
    if (kread_file) deallocate(ftable, spectrum_table)

    if (info > 4) write(*,*) 'Done'

  CONTAINS

    SUBROUTINE get_spectrum_powerlaw()

      do i = 1, nof
         if (f(i) > fmin) then
            spectrum(i) = sigma*sigma/fsample*(f(i)/fknee)**(slope)
         else
            spectrum(i) = sigma*sigma/fsample*(fmin/fknee)**(slope)
         endif
      enddo

    END SUBROUTINE get_spectrum_powerlaw


    SUBROUTINE read_spectrum()

      ! RK edit:
      ! To handle arbitrary M3 detector ordering,
      ! the spectrum file format is now this:
      !
      ! nobins nodetectors
      ! detname1 detname2 detname3 ...
      ! sigma1 sigma2 sigma3 ...
      ! freq psd1 psd2 psd3 psd4 ...
      !

      integer :: ierr, icol
      real(dp), allocatable :: dtemp(:, :), sigmas(:)
      character(len=SLEN), allocatable :: stemp(:)

      if (id == 0) then
         open(unit=17, file=file_spectrum, status='old')

         read (17, *, iostat=ierr) nolines, nocols
         if (ierr /= 0) then
            write (*,*) 'ERROR, ' // trim(file_spectrum) // &
                 ' is not formatted correctly'
            stop
         endif

         allocate(ftable(nolines), spectrum_table(nolines, nodetectors), &
              dtemp(nocols+1, nolines), stemp(nocols), sigmas(nocols), &
              stat=ierr)
         if (ierr /= 0) stop 'No room to read the noise spectra'

         ! get names of detectors
         read (17, *, iostat=ierr) stemp
         if (ierr /= 0) stop 'Failed to read detector names from spectrum file'

         ! get detector noise levels
         read (17, *, iostat=ierr) sigmas
         if (ierr /= 0) stop &
              'Failed to read detector noise levels from spectrum file'

         ! get the table of psds
         read (17, *, iostat=ierr) dtemp
         if (ierr /= 0) stop 'Failed to read the table from spectrum file'

         close(17)

         ! now insert the read columns for the correct detectors
         ftable = log(dtemp(1, :))
         do idet = 1, nodetectors
            ! find the column of the read table that corresponds to
            ! detector # idet
            do icol = 1, nocols+1
               if (icol > nocols) then
                  write (*,*) 'ERROR: ' // trim(detectors(idet)%name) // &
                       ' not in ' // trim(file_spectrum)
                  stop
               endif
               if (index(stemp(icol),trim(detectors(idet)%name)) > 0) exit
            end do

            ! insert the sigma
            detectors(idet)%sigmas = sigmas(icol)
            write (*,*) trim(detectors(idet)%name)//', SIGMA == ',sigmas(icol)

            ! prepare for logaritmic interpolation
            ! first column was the frequency
            spectrum_table(:, idet) = log(dtemp(icol+1, :))
         end do

         deallocate(stemp, dtemp)
      endif

      ! broadcast the results
      call broadcast_mpi(nolines, 0)
      do idet = 1, nodetectors
         call broadcast_mpi(detectors(idet)%sigmas, detectors(idet)%npsd, 0)
         detectors(idet)%weights = 1 / detectors(idet)%sigmas**2
      end do

      if (id /= 0) then
         allocate(ftable(nolines), spectrum_table(nolines, nodetectors), &
              stat=ierr)
         if (ierr /= 0) stop 'No room for ftable'
      endif

      call broadcast_mpi(ftable, nolines, 0)
      nocols = nodetectors
      do idet = 1, nodetectors
         call broadcast_mpi(spectrum_table(:, idet), nolines, 0)
      end do

    END SUBROUTINE read_spectrum


    SUBROUTINE get_spectrum_file()
      !
      ! Find the spectrum for detector idet by logarithmic interpolation
      !
      integer  :: i, k, icol
      real(dp) :: p, logf, logspec

      if (nocols > 1) then
         icol = idet
      else
         icol = 1
      endif

      k = 1
      do i = 1, nof
         logf = log(f(i))

         do
            if (logf <= ftable(k+1) .or. k == nolines-1) exit
            k = k + 1
         enddo

         p = (logf-ftable(k)) / (ftable(k+1)-ftable(k))

         if (p < 0) p = 0
         if (p > 1) p = 1

         logspec = (1.d0-p)*spectrum_table(k, icol) + p*spectrum_table(k+1, icol)
         spectrum(i) = exp(logspec) ! *sigma*sigma/fsample  -RK
      enddo

    END SUBROUTINE get_spectrum_file



    SUBROUTINE get_spectrum_interp(idet, ipsd)
      !
      ! Find the spectrum for detector idet by logarithmic interpolation
      !
      integer :: idet, ipsd

      integer  :: i, j, ierr, n
      real(dp) :: plateau
      real(dp), allocatable :: freqs(:), data(:)
      real(dp), pointer :: psd(:)

      psd => detectors(idet)%psds(:, ipsd)

      call interpolate_psd(detectors(idet)%psdfreqs, psd, f, spectrum)

      if (noise_weights_from_psd) then

         ! the sigmas have already been updated from the PSDs

         !plateau = detectors(idet)%sigmas(ipsd)**2 / fsample
         plateau = detectors(idet)%plateaus(ipsd)

      else

         ! Measure the white noise levels, update sigmas
         ! (noise weights remain constant)

         if (radiometers) then

            ! well-behaved PSD, just pick the last bin value

            n = 1
            allocate(freqs(n), data(n))
            freqs = fsample / 2

            call interpolate_psd(detectors(idet)%psdfreqs, psd, freqs, data)

            plateau = data(1)

         else

            ! Measure and subtract the plateau value, VERY Planck Specific

            n = 1000
            allocate(freqs(n), data(n))
            freqs = (/ (1 + dble(i*10)/n, i=0, n-1) /)

            call interpolate_psd(detectors(idet)%psdfreqs, psd, freqs, data)

            plateau = minval(data)

         end if

         deallocate(freqs, data)

         detectors(idet)%sigmas(ipsd) = sqrt(plateau * fsample)

      end if

      if (radiometers) then
         ! subtract with a small margin
         spectrum = spectrum - plateau * .999999
      else
         ! subtract with a larger margin
         spectrum = spectrum - plateau * .99
      end if

      if (any(spectrum <= 0)) then
         loop_bins : do i = 1, nof
            if (spectrum(i) <= 0) then
               ! First look for valid value in higher frequency
               loop_next_bins: do j = i+1, nof
                  if (spectrum(j) > 0) then
                     spectrum(i) = spectrum(j)
                     cycle loop_bins
                  end if
               end do loop_next_bins
            end if
            if (spectrum(i) <= 0) then
               ! If needed, get the value from lower frequency
               loop_previous_bins: do j = i-1, 1, -1
                  if (spectrum(j) > 0) then
                     spectrum(i) = spectrum(j)
                     cycle loop_bins
                  end if
               end do loop_previous_bins
            end if
         end do loop_bins
      end if

      if (f(1) == 0) spectrum(1) = 0 ! leave the mean unconstrained

    END SUBROUTINE get_spectrum_interp

  END SUBROUTINE build_filter


  !-----------------------------------------------------------------------


  SUBROUTINE cinvmul(ca, aa)

    real(dp), intent(out) :: ca(noba_short, nodetectors)
    real(dp), intent(in) :: aa(noba_short, nodetectors)

    if (kfilter) then
       call convolve_pp(ca, aa, fcov, 2)
    else
       ca = 0
    endif

  END SUBROUTINE cinvmul


  !-----------------------------------------------------------------------


  SUBROUTINE convolve_pp(y, x, fc, mode)
    ! Convolve a baseline vector with the noise prior, one pointing period
    ! at a time. FIXME: should have an option to apply the filter across the boundaries.

    real(dp), intent(out) :: y(noba_short, nodetectors)
    real(dp), intent(in) :: x(noba_short, nodetectors)
    complex(dp), intent(in) :: fc(nof/2+1, npsdtot)
    integer, intent(in) :: mode
    integer :: ichunk, idet, m, no, noba, kstart, i, ipsd, ierr
    real(dp) :: x0

    real(C_DOUBLE), pointer :: xx(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: fx(:) => NULL()
    type(C_PTR) :: pxx(0:nthreads-1), pfx(0:nthreads-1)
    integer :: id_threads

    call reset_time(14)

    y = 0

    ! Compiler bug workaround:
    ! GCC 4.7 up to at least 4.8.0 segfaults when compiling
    ! openMP pragmas that declare type(C_PTR) variables private
    ! We allocate enough pointers for every thread to use a separate
    ! (but shared) one

    !$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(nthreads) &
    !$OMP     PRIVATE(idet, ichunk, noba, kstart, x0, m, no, xx, fx, &
    !$OMP         id_thread, ipsd, ierr)

    id_thread = omp_get_thread_num()

    ! FFTW malloc is not thread safe but we want to allocate the workspace
    ! to be aligned to allow SIMD optimization in FFTW

!!$    !$OMP CRITICAL
!!$    pfx(id_thread) = fftw_alloc_complex(int(nof/2+1, C_SIZE_T))
!!$    call c_f_pointer(pfx(id_thread), fx, [nof/2+1])
!!$
!!$    pxx(id_thread) = fftw_alloc_real(int(nof, C_SIZE_T))
!!$    call c_f_pointer(pxx(id_thread), xx, [nof])
!!$    !$OMP END CRITICAL
    allocate(fx(nof/2+1), xx(nof), stat=ierr)
    if (ierr /= 0) stop 'No room for Fourier transform'

    !$OMP DO SCHEDULE(DYNAMIC,1)
    do idet = 1, nodetectors

       do ichunk = 1, ninterval
          noba = noba_short_pp(ichunk)
          kstart = sum(noba_short_pp(1:ichunk-1))
          ipsd = psd_index(idet, baselines_short_time(kstart+1))

          if (ipsd == -1) then
             ! no PSD
             y(kstart+1:kstart+noba, idet) = x(kstart+1:kstart+noba, idet)
             cycle
          end if

          ! Remove the PID average first.
          ! This ensures that solution is independent of a constant
          ! offset/pointing period
          x0 = 0
          if (noba > 0 .and. mode >= 1) then
             x0 = sum(x(kstart+1:kstart+noba, idet)) / noba
          end if

          ! overlap-save loop
          m = 0
          do
             if (m == noba) exit

             no = min(noba-m, nofilter)
             xx = 0.0

             ! Copy one data section+tails and convolve
             if (m > 0) then
                xx(1:notail) = x(kstart+m-notail+1:kstart+m, idet) - x0
             end if

             if (no+notail <= noba-m) then
                xx(notail+1:no+2*notail) = &
                     x(kstart+m+1:kstart+m+no+notail, idet) - x0
             else
                xx(notail+1:notail+noba-m) = &
                     x(kstart+m+1:kstart+noba, idet) - x0
             endif

             call dfft(fx, xx)
             fx = fx * fc(:, ipsd)
             call dfftinv(xx, fx)

             y(kstart+m+1:kstart+m+no, idet) = xx(notail+1:notail+no)
             m = m + no
          end do

          x0 = 0
          if (noba > 0 .and. mode == 2) then
             x0 = sum(y(kstart+1:kstart+noba, idet)) / noba
          end if
          y(kstart+1:kstart+noba, idet) = y(kstart+1:kstart+noba, idet) - x0

       end do
    end do
    !$OMP END DO

!!$    !$OMP CRITICAL
!!$    call fftw_free(pfx(id_thread))
!!$    call fftw_free(pxx(id_thread))
!!$    !$OMP END CRITICAL
    deallocate(fx)
    deallocate(xx)

    !$OMP END PARALLEL

    cputime_filter = cputime_filter + get_time(14)

  END SUBROUTINE convolve_pp



  !---------------------------------------------------------------------------


  SUBROUTINE construct_preconditioner(nna)

    real(dp), intent(in)  :: nna(noba_short_max, nodetectors)
    integer :: i, j, k, kstart, n, noba, idet, ichunk, ipsd, ierr, nbandmin, &
         try, ipsddet
    real(dp), allocatable :: invcov(:, :), blockm(:, :)
    real(dp) :: r
    integer, parameter :: trymax = 3
    integer :: ntries(0:trymax)
    logical :: neg

    if (precond_width < 1) return

    call reset_time(16)

    if (ID == 0 .and. info > 1) write(*,*) 'Constructing preconditioner'
    if (info > 4) write(*, idf) ID, 'Constructing preconditioner'

    if (.not. kfilter .or. precond_width == 1) then
       nband = 1
       use_diagonal = .true.
       if (.not. allocated(prec_diag)) then
          allocate(prec_diag(noba_short_max, nodetectors), stat=ierr)
          if (ierr /= 0) stop 'No room for prec_diag'
       end if
       memory_precond = noba_short_max*nodetectors*8.

       prec_diag = 0
       do idet = 1, nodetectors
          do ichunk = 1, ninterval
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(1:ichunk-1))
             ipsddet = psd_index_det(idet, baselines_short_time(kstart+1))
             if (ipsddet < 0) then
                ! No PSD available
                do k = kstart+1, kstart+noba
                   prec_diag(k, idet) = 1.
                end do
                cycle
             end if
             do k = kstart+1, kstart+noba
                if (nna(k, idet) == 0) then
                   prec_diag(k, idet) = 1. / detectors(idet)%weights(ipsddet)
                else
                   prec_diag(k, idet) = 1. / nna(k, idet)
                end if
             end do
          end do
       end do
       return
    end if

    use_diagonal = .false.

    if (use_fprecond) then
       if (allocated(fprecond)) deallocate(fprecond)
       allocate(fprecond(nof/2+1, npsdtot), stat=ierr)
       if (ierr /= 0) stop 'No room for fprecond'
       memory_precond = (nof/2+1)*npsdtot*16.
       fprecond = 0
       do idet = 1, nodetectors
          do ichunk = 1, ninterval
             noba = noba_short_pp(ichunk)
             if (noba == 0) cycle
             kstart = sum(noba_short_pp(1:ichunk-1))
             ipsd = psd_index(idet, baselines_short_time(kstart+1))
             ipsddet = psd_index_det(idet, baselines_short_time(kstart+1))
             if (ipsddet < 0) cycle
             fprecond(:, ipsd) = 1.d0 /( &
                  nshort*detectors(idet)%weights(ipsddet) + fcov(:, ipsd))
          end do
       end do
       return
    end if

    nband = min(precond_width, nof / 2 - 1)
    if (.not. allocated(bandprec)) then
       allocate(bandprec(nband+1, noba_short_max, nodetectors), stat=ierr)
       if (ierr /= 0) stop 'No room for bandprec'
    end if
    if (.not. allocated(prec_diag)) then
       allocate(prec_diag(noba_short_max, nodetectors), stat=ierr)
       if (ierr /= 0) stop 'No room for prec_diag'
    end if

    bandprec = 0
    prec_diag = 0
    memory_precond = noba_short_max*nodetectors*(nband*4.+8.)

    allocate(invcov(nband+1, npsdtot), stat=ierr)
    if (ierr /= 0) stop 'No room for invcov'

    do ipsd = 1, npsdtot
       call dfftinv(xx, fcov(:, ipsd)) ! C_a inverse into real domain
       invcov(:, ipsd) = xx(1:nband+1)
    end do

    ntries = 0

    !$OMP PARALLEL DO IF (nodetectors >= nthreads) &
    !$OMP     DEFAULT(NONE) &
    !$OMP     SHARED(nodetectors, ninterval, noba_short_pp, &
    !$OMP         baselines_short_time, nband, invcov, id, &
    !$OMP         sampletime, nof, baselines_short_start, &
    !$OMP         baselines_short_stop, bandprec, nna, prec_diag, &
    !$OMP         detectors, nthreads, ntries) &
    !$OMP     PRIVATE(idet, ichunk, noba, kstart, ipsd, blockm, &
    !$OMP         i, j, k, n, neg, ierr, nbandmin, try)
    do idet = 1,nodetectors

       !$OMP PARALLEL DO IF (nodetectors < nthreads) &
       !$OMP     DEFAULT(NONE) &
       !$OMP     SHARED(idet, nodetectors, ninterval, &
       !$OMP         noba_short_pp, baselines_short_time, nband, &
       !$OMP         invcov, id, sampletime, nof, baselines_short_start, &
       !$OMP         baselines_short_stop, bandprec, nna, prec_diag, &
       !$OMP         detectors, ntries) &
       !$OMP     PRIVATE(ichunk, noba, kstart, ipsd, blockm, &
       !$OMP         i, j, k, n, neg, ierr, nbandmin, try)
       do ichunk = 1, ninterval
          noba = noba_short_pp(ichunk)
          kstart = sum(noba_short_pp(1:ichunk-1))
          ipsd = psd_index(idet, baselines_short_time(kstart+1))

          if (noba < nband) then
             ! This chunk (pointing period) has fewer baselines than
             ! the width of the preconditioner.  This should be a very
             ! infrequent issue unless the user has chosen a too wide
             ! preconditioner.
             if (idet == 1 .and. omp_get_thread_num() == 0) then
                write(*,*) id, ' : WARNING : noba < nband : ', noba, ' < ', &
                     nband, ' : ichunk = ', ichunk, &
                     ', tstart = ', &
                     sampletime(baselines_short_start(kstart+1)), &
                     ', length = ', &
                     sampletime(baselines_short_stop(kstart+noba)) &
                     - sampletime(baselines_short_start(kstart+1))
             end if
             nbandmin = noba
          else
             nbandmin = nband
          end if

          if (ipsd == -1) then
             bandprec(1:nbandmin+1, kstart+1:kstart+noba, idet) = 0
             prec_diag(kstart+1:kstart+noba, idet) = 0
             cycle
          end if

          ! Rows from kstart+1 to kstart+noba for one band-diagonal
          ! submatrix.  Do the computation in double precision
          allocate(blockm(nbandmin+1, noba), stat=ierr)
          if (ierr /= 0) stop 'No room for blockm'

          try = 0
          loop_try : do
             blockm = 0
             blockm = spread(invcov(1:nbandmin+1, ipsd), 2, noba)
             blockm(1, :) = blockm(1, :) + nna(kstart+1:kstart+noba, idet)

             if (try /= 0) then
                ! Band diagonal failed to decompose.
                ! Regularize the matrix for decomposition by ensuring that the
                ! band power at nbandmin is negligible
                do i = 1, nbandmin+1
                   blockm(i, :) = blockm(i, :) * exp(-(dble(try*i) / nbandmin) ** 2)
                end do
             end if

             ! Cholesky decompose
             call DPBTRF('L', noba, nbandmin, blockm, nbandmin+1, ierr)
             if (ierr == 0 .or. try == trymax) then
                exit loop_try
             end if
             try = try + 1
          end do loop_try
          ntries(try) = ntries(try) + 1

          if (ierr /= 0) then
             if (ierr < 0) then
                print *, id, ' : failed to cholesky decompose. argument # ', &
                     -ierr, ' had an illegal value'
             end if
             print *, id, ' : Preconditioner was not positive definite. ', &
                  'Will use C_a preconditioner for ', &
                  'det = ', idet, ' = ', trim(detectors(idet)%name), &
                  ' : chunk = ', ichunk, ', ipsd = ', ipsd
             bandprec(1:nbandmin+1, kstart+1:kstart+noba, idet) = 0
             prec_diag(kstart+1:kstart+noba, idet) = 0
          else
             bandprec(1:nbandmin+1, kstart+1:kstart+noba, idet) = blockm(:, :)
             prec_diag(kstart+1:kstart+noba, idet) = 1 / blockm(1, :)
          end if

          deallocate(blockm)

       end do
       !$OMP END PARALLEL DO

    end do
    !$OMP END PARALLEL DO

    call sum_mpi(ntries, trymax+1)
    if (id == 0) then
       print *, 'Constructed ', sum(ntries), ' band preconditioners. ', &
            ntries(0), ' required no regularization'
       do try = 1, trymax
          if (ntries(try) > 0) print *, ntries(try), ' required level = ', &
               try, ' regularization'
       end do
    end if

    deallocate(invcov)
    cputime_prec_construct = cputime_prec_construct + get_time(16)

    if (info > 5) write(*, idf) ID, 'Done'

  END SUBROUTINE construct_preconditioner


  !------------------------------------------------------------------------------


  SUBROUTINE preconditioning_band(z, r)
    !Apply the preconditioner

    ! z and r may be 3-dimensional arrays in the calling code but
    ! here we collapse the first dimension.  The band preconditioner
    ! only works for basis_order == 0.
    real(dp), intent(out), target :: z(noba_short, nodetectors)
    real(dp), intent(in), target :: r(noba_short, nodetectors)
    integer :: i, j, k, idet, kstart, noba, ichunk, ierr, ipsd, ipsddet

    real(C_DOUBLE), pointer :: xx(:)
    complex(C_DOUBLE_COMPLEX), pointer :: fx(:)
    type(C_PTR) :: pxx(0:nthreads-1), pfx(0:nthreads-1)
    integer :: id_thread, m, no, nbandmin
    real(dp) :: x0

    if (precond_width < 1) then
       z = r
       return
    endif

    call reset_time(16)

    if (use_diagonal) then
       !$OMP PARALLEL DO IF (nodetectors >= nthreads) &
       !$OMP     DEFAULT(NONE) &
       !$OMP     SHARED(nodetectors, noba_short, z, r, prec_diag, nthreads) &
       !$OMP     PRIVATE(idet, k)
       do idet = 1, nodetectors
          !$OMP PARALLEL DO IF (nodetectors < nthreads) &
          !$OMP     DEFAULT(NONE) &
          !$OMP     SHARED(nodetectors, z, noba_short, r, prec_diag, idet) &
          !$OMP     PRIVATE(k)
          do k = 1, noba_short
             z(k, idet) = r(k, idet) * prec_diag(k, idet)
          end do
          !$OMP END PARALLEL DO
       end do
       !$OMP END PARALLEL DO
    else if (use_fprecond) then
       ! Apply the filter-based preconditioner
       call convolve_pp(z, r, fprecond, 0)
    else
       !$OMP PARALLEL DO IF (nodetectors >= nthreads) &
       !$OMP     DEFAULT(NONE) &
       !$OMP     SHARED(noba_short_pp, ninterval, nband, prec_diag, &
       !$OMP         bandprec, baselines_short_time, pfx, pxx, z, r, &
       !$OMP         nodetectors, id, nof, nofilter, notail, &
       !$OMP         fcov, nshort, detectors, nthreads) &
       !$OMP     PRIVATE(idet, ichunk, kstart, noba, j, k, &
       !$OMP         ierr, ipsd, x0, m, no, xx, fx, id_thread, ipsddet, &
       !$OMP         nbandmin)
       do idet = 1, nodetectors
          !$OMP PARALLEL DO IF (nodetectors < nthreads) &
          !$OMP     DEFAULT(NONE) &
          !$OMP     SHARED(idet, noba_short_pp, ninterval, nband, prec_diag, &
          !$OMP         bandprec, baselines_short_time, pfx, pxx, z, &
          !$OMP         r, nodetectors, id, nof, nofilter, notail, &
          !$OMP         fcov, nshort, detectors) &
          !$OMP     PRIVATE(ichunk, kstart, noba, j, k, ierr, &
          !$OMP         ipsd, x0, m, no, xx, fx, id_thread, ipsddet, nbandmin)
          do ichunk = 1, ninterval

             noba = noba_short_pp(ichunk)
             if (noba == 0) cycle
             kstart = sum(noba_short_pp(1:ichunk-1))

             if (noba < nband) then
                nbandmin = noba
             else
                nbandmin = nband
             end if

             if (any(prec_diag(kstart+1:kstart+noba, idet) /= 0)) then
                ! Use the precomputed Cholesky decomposition

                z(kstart+1:kstart+noba, idet) = r(kstart+1:kstart+noba, idet)
                call DPBTRS('L', noba, nbandmin, 1, &
                     bandprec(1:nbandmin+1, kstart+1:kstart+noba, idet), &
                     nbandmin+1, z(kstart+1:kstart+noba, idet), noba, ierr)
                if (ierr /= 0) then
                   print *, id, ' : failed to Cholesky decompose. argument # ', &
                        -ierr, ' had an illegal value'
                   call abort_mpi('bad preconditioner2')
                end if
             else
                ! The Cholesky decomposition is not available.
                ! Apply (F^T N_w^-1 F + C_a^-1) instead

                ipsd = psd_index(idet, baselines_short_time(kstart+1))
                ipsddet = psd_index_det(idet, baselines_short_time(kstart+1))
                if (ipsddet < 0) then
                   ! No PSD available
                   z(kstart+1:kstart+noba, idet) = r(kstart+1:kstart+noba, idet)
                   cycle
                end if

                id_thread = omp_get_thread_num()
!!$                !$OMP CRITICAL
!!$                pfx(id_thread) = fftw_alloc_complex(int(nof/2+1, C_SIZE_T))
!!$                call c_f_pointer(pfx(id_thread), fx, [nof/2+1])
!!$
!!$                pxx(id_thread) = fftw_alloc_real(int(nof, C_SIZE_T))
!!$                call c_f_pointer(pxx(id_thread), xx, [nof])
!!$                !$OMP END CRITICAL
                allocate(fx(nof/2+1), xx(nof), stat=ierr)
                if (ierr /= 0) stop 'No room for Fourier transform'

                !Remove the interval average first.
                !This ensures that solution is independent of a constant
                ! offset/pointing period
                x0 = sum(r(kstart+1:kstart+noba, idet)) / noba

                !overlap-save loop
                m = 0
                do
                   if (m == noba) exit

                   no = min(noba-m, nofilter)
                   xx = 0

                   !Copy one data section+tails and convolve
                   if (m > 0) then
                      xx(1:notail) = r(kstart+m-notail+1:kstart+m, idet) - x0
                   end if

                   if (no+notail <= noba-m) then
                      xx(notail+1:no+2*notail) = &
                           r(kstart+m+1:kstart+m+no+notail, idet) - x0
                   else
                      xx(notail+1:notail+noba-m) = &
                           r(kstart+m+1:kstart+noba, idet) - x0
                   endif

                   call dfft(fx, xx)
                   fx = fx / (fcov(:, ipsd) &
                        + nshort * detectors(idet)%weights(ipsddet))
                   call dfftinv(xx, fx)

                   z(kstart+m+1:kstart+m+no, idet) = xx(notail+1:notail+no)
                   m = m + no
                end do

                x0 = sum(z(kstart+1:kstart+noba, idet)) / noba
                z(kstart+1:kstart+noba, idet) = &
                     z(kstart+1:kstart+noba, idet) - x0

!!$                !$OMP CRITICAL
!!$                call fftw_free(pfx(id_thread))
!!$                call fftw_free(pxx(id_thread))
!!$                !$OMP END CRITICAL
                deallocate(fx)
                deallocate(xx)

             end if

          end do
          !$OMP END PARALLEL DO
       end do
       !$OMP END PARALLEL DO
    end if

    cputime_precond = cputime_precond + get_time(16)

  END SUBROUTINE preconditioning_band


  !---------------------------------------------------------------------------



END MODULE noise_routines
