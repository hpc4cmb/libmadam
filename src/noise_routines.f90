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

  integer, save :: nof = -1, nofilter = 0

  real(dp), save, public :: memory_filter = 0, memory_precond = 0, &
       cputime_filter = 0, cputime_precond = 0, &
       cputime_prec_construct = 0

  !real(dp), allocatable :: xx(:)
  !complex(dp), allocatable :: fx(:)
  complex(dp), allocatable, target :: fcov(:, :), fprecond(:, :)
  real(dp), allocatable :: invcov(:, :)

  !Store interpolated spectra
  real(dp), allocatable :: freqtab(:), spectab(:, :)

  !Preconditioner

  logical :: use_diagonal = .false.
  !real(dp), allocatable, target :: bandprec(:, :, :)
  type bandprec_pointer
     real(dp), allocatable :: data(:, :)
  end type bandprec_pointer
  ! Allow up to NSUBMAX sub-intervals for filtering and preconditioning
  integer, parameter :: NSUBMAX = 100
  type(bandprec_pointer), allocatable :: bandprec(:, :, :)
  real(dp), allocatable, target :: prec_diag(:, :)

  ! PSD information

  integer :: npsdtot
  integer, allocatable :: psddet(:), psdind(:)
  real(dp), allocatable :: psdstart(:), psdstop(:)

  public cinvmul, build_filter, close_filter, &
       construct_preconditioner, preconditioning_band, interpolate_psd

  ! LAPACK
  external DPBTRF, DPBTRS

CONTAINS

  subroutine measure_noise_weights(radiometers)
    !
    ! Determine noise weights and split the PSDs into
    ! baselines and "white" (small scale) noise
    !
    logical, intent(in) :: radiometers
    real(dp) :: fstep, fbase, flower, fupper, rms, plateau
    integer :: ipsd, nbin, ilower, iupper, ibase, i, idet
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
             detectors(idet)%weights(ipsd) = 1 / (plateau * fsample)
          end do
       end do
       deallocate(freqs, data)
    end if

  end subroutine measure_noise_weights


  subroutine interpolate_psd(freq, psd, targetfreq, targetpsd)

    ! Use linear interpolation in the log-log plane

    real(dp), intent(in) :: freq(:), psd(:), targetfreq(:)
    real(dp), intent(out) :: targetpsd(:)

    real(dp), allocatable :: logfreq(:), logpsd(:)
    real(dp) :: flog, r

    integer :: n_in, n_out, i, ierr, j, startbin

    if (size(freq) /= size(psd)) then
       call abort_mpi('freq and psd do not have the same shape')
    end if

    if (size(targetfreq) /= size(targetpsd)) then
       call abort_mpi('target freq and psd do not have the same shape')
    end if

    n_in = size(freq)
    n_out = size(targetfreq)
    allocate(logfreq(n_in), logpsd(n_in), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room to interpolate the PSD')

    if (all(freq == 0)) call abort_mpi('input freqs are all zero.')
    startbin = 1
    do while (freq(startbin) == 0)
       startbin = startbin + 1
    end do

    logfreq = log(freq)
    logpsd = log(psd)

    j = startbin
    do i = 1, n_out
       if (targetfreq(i) == 0) then
          targetpsd(i) = 0
          cycle
       end if
       flog = log(targetfreq(i))
       if (j < n_in - 1) then
          do while (logfreq(j+1) < flog)
             j = j + 1
             if (j == n_in - 1) exit
          end do
       end if
       if (psd(j) > 0 .and. psd(j+1) > 0) then
          ! Logarithmic interpolation
          r = (flog - logfreq(j)) / (logfreq(j+1) - logfreq(j))
          targetpsd(i) = exp(logpsd(j) * (1._dp - r) + logpsd(j + 1) * r)
       else
          ! Linear interpolation
          r = (targetfreq(i) - freq(j)) / (freq(j+1) - freq(j))
          targetpsd(i) = psd(j) * (1._dp - r) + psd(j + 1) * r
       end if
    end do

    deallocate(logfreq, logpsd)

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

    integer :: nof_min, ierr

    if (.not. kfilter) return

    nofilter = 2 * noba_short_pp_max

    nof_min = nofilter
    ! Ensure that the filter is long enough for the preconditioner
    if (precond_width_max > nof_min-1) nof_min = precond_width_max + 1

    nof = 1
    do
       if (nof >= nof_min) exit
       nof = 2*nof
    end do
    if (ID == 0) write(*,*) 'FFT length =', nof

    !allocate(fx(nof/2+1), xx(nof), stat=ierr) ! Work space
    !if (ierr /= 0) call abort_mpi('No room for fx and xx')

    memory_filter = ((nof/2+1)*16. + nof*8) * nthreads

    call init_fourier(nof)

  END SUBROUTINE init_filter


  !-----------------------------------------------------------------------


  subroutine allocate_bandprec
    integer :: ierr

    call free_bandprec
    allocate(bandprec(NSUBMAX, ninterval, nodetectors), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for bandprec')
  end subroutine allocate_bandprec



  subroutine free_bandprec
    integer :: ichunk, idet, isub

    if (allocated(bandprec)) then
       do idet = 1, nodetectors
          do ichunk = 1, ninterval
             do isub = 1, NSUBMAX
                if (allocated(bandprec(isub, ichunk, idet)%data)) then
                   deallocate(bandprec(isub, ichunk, idet)%data)
                end if
             end do
          end do
       end do
       deallocate(bandprec)
    end if
  end subroutine free_bandprec



  SUBROUTINE close_filter

    if (.not. kfilter) return

    if (allocated(fcov)) deallocate(fcov)
    if (allocated(invcov)) deallocate(invcov)
    if (allocated(fprecond)) deallocate(fprecond)
    if (allocated(psddet)) deallocate(psddet)
    if (allocated(psdind)) deallocate(psdind)
    if (allocated(psdstart)) deallocate(psdstart)
    if (allocated(psdstop)) deallocate(psdstop)
    !if (allocated(fx)) deallocate(fx)
    !if (allocated(xx)) deallocate(xx)
    if (allocated(prec_diag)) deallocate(prec_diag)
    call free_bandprec

    call close_fourier

    nof = -1

  END SUBROUTINE close_filter


  !-----------------------------------------------------------------------


  SUBROUTINE build_filter()
    ! Build the noise filter
    integer, parameter :: kmax = 2
    integer :: i, k, idet, nolines, nocols, num_threads, itask
    real(dp) :: fa, x
    logical :: kread_file
    real(dp), allocatable :: aspec(:), f(:), g(:), spectrum(:), &
         ftable(:), spectrum_table(:, :)
    integer :: ipsd, ipsdtot, ierr
    real(dp) :: memsum, mem_min, mem_max

    if (.not. kfilter) return

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
    if (ierr /= 0) call abort_mpi('No room for fcov')
    memory_filter = memory_filter + (nof/2+1)*npsdtot*16. + npsdtot*24.

    memsum = memory_filter / 2**20
    mem_min = memsum; mem_max = memsum
    call min_mpi(mem_min); call max_mpi(mem_max)
    call sum_mpi(memsum)
    if (ID == 0 .and. info > 0) then
       write(*,'(1x,a,t32,3(f9.1," MB"))') &
            'Allocated memory for filter:', memsum, mem_min, mem_max
    end if

    kread_file = (len_trim(file_spectrum) > 0)
    if (kread_file) call read_spectrum

    fa = fsample/nshort

    !$OMP PARALLEL NUM_THREADS(nthreads) &
    !$OMP     DEFAULT(NONE) &
    !$OMP     SHARED(nodetectors, detectors, psddet, psdind, psdstart, &
    !$OMP         psdstop, nof, fa, fcov, kread_file) &
    !$OMP     PRIVATE(id_thread, num_threads, ipsdtot, itask, idet, ipsd, &
    !$OMP         aspec, spectrum, f, g, ierr, k, x)

    id_thread = omp_get_thread_num()
    num_threads = omp_get_num_threads()

    allocate(aspec(nof/2+1), spectrum(nof), f(nof), g(nof), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for aspec')

    ipsdtot = 0
    itask = -1
    do idet = 1, nodetectors
       if (.not. allocated(detectors(idet)%psds)) then
          call abort_mpi( &
               'All detectors must have at least one PSD to build the filter')
       end if
       do ipsd = 1, detectors(idet)%npsd
          ipsdtot = ipsdtot + 1
          itask = itask + 1
          if (modulo(itask, num_threads) /= id_thread) cycle

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
                f(i) = k*fa + (i-1)*fa/nof
             end do

             if (kread_file) then
                call get_spectrum_file
             else
                call get_spectrum_interp(idet, ipsd, f, spectrum)
             end if

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

          where (abs(aspec) > 1e-30)
             aspec = 1.d0 / (fa * aspec)
          elsewhere
             aspec = 0
          end where
          aspec(1) = 0
          fcov(:, ipsdtot) = aspec

          ! RK edit: check for NaNs
          do i = 1, nof / 2 + 1
             if (isnan(abs(fcov(i, ipsdtot)))) &
                  call abort_mpi('NaN in fcov. Check input PSD and sigma.')
          end do

       end do ! ipsd
    end do ! idet

    deallocate(aspec, spectrum, f, g)

    !$OMP END PARALLEL

    if (kread_file) deallocate(ftable, spectrum_table)

    if (info > 4) write(*,*) 'Done'

  CONTAINS

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
            call abort_mpi('ERROR, ' // trim(file_spectrum) // &
                 ' is not formatted correctly')
         end if

         allocate(ftable(nolines), spectrum_table(nolines, nodetectors), &
              dtemp(nocols+1, nolines), stemp(nocols), sigmas(nocols), &
              stat=ierr)
         if (ierr /= 0) call abort_mpi('No room to read the noise spectra')

         ! get names of detectors
         read (17, *, iostat=ierr) stemp
         if (ierr /= 0) then
            call abort_mpi('Failed to read detector names from spectrum file')
         end if

         ! get detector noise levels
         read (17, *, iostat=ierr) sigmas
         if (ierr /= 0) then
            call abort_mpi( &
                 'Failed to read detector noise levels from spectrum file')
         end if

         ! get the table of psds
         read (17, *, iostat=ierr) dtemp
         if (ierr /= 0) then
            call abort_mpi('Failed to read the table from spectrum file')
         end if

         close(17)

         ! now insert the read columns for the correct detectors
         ftable = log(dtemp(1, :))
         do idet = 1, nodetectors
            ! find the column of the read table that corresponds to
            ! detector # idet
            do icol = 1, nocols+1
               if (icol > nocols) then
                  call abort_mpi('ERROR: ' // trim(detectors(idet)%name) // &
                       ' not in ' // trim(file_spectrum))
               end if
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
      end if

      ! broadcast the results
      call broadcast_mpi(nolines, 0)
      do idet = 1, nodetectors
         call broadcast_mpi(detectors(idet)%sigmas, detectors(idet)%npsd, 0)
         detectors(idet)%weights = 1 / detectors(idet)%sigmas**2
      end do

      if (id /= 0) then
         allocate(ftable(nolines), spectrum_table(nolines, nodetectors), &
              stat=ierr)
         if (ierr /= 0) call abort_mpi('No room for ftable')
      end if

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
      integer :: i, k, icol
      real(dp) :: p, logf, logspec

      if (nocols > 1) then
         icol = idet
      else
         icol = 1
      end if

      k = 1
      do i = 1, nof
         logf = log(f(i))

         do
            if (logf <= ftable(k+1) .or. k == nolines-1) exit
            k = k + 1
         end do

         p = (logf-ftable(k)) / (ftable(k+1)-ftable(k))

         if (p < 0) p = 0
         if (p > 1) p = 1

         logspec = (1.d0-p)*spectrum_table(k, icol) + p*spectrum_table(k+1, icol)
         spectrum(i) = exp(logspec) ! *sigma*sigma/fsample  -RK
      end do

    END SUBROUTINE get_spectrum_file



    SUBROUTINE get_spectrum_interp(idet, ipsd, f, spectrum)
      !
      ! Find the spectrum for detector idet by logarithmic interpolation
      !
      integer, intent(in) :: idet, ipsd
      real(dp), intent(inout) :: f(nof), spectrum(nof)

      integer :: i, j, n
      real(dp) :: plateau, margin
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
         margin = 1e-6
      else
         ! subtract with a larger margin
         margin = 1e-3
      end if

      where (spectrum > plateau * (1 + margin))
         spectrum = spectrum - plateau
      elsewhere
         spectrum = spectrum * margin
      end where

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


  SUBROUTINE cinvmul(ca, aa, nna)

    real(dp), intent(out) :: ca(0:basis_order, noba_short, nodetectors)
    real(dp), intent(in) :: aa(0:basis_order, noba_short, nodetectors)
    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)

    if (kfilter) then
       call convolve_pp(ca, aa, fcov, nna)
    else
       ca = 0
    end if

  END SUBROUTINE cinvmul


  !-----------------------------------------------------------------------


  subroutine trim_interval(kstart, kstop, noba, idet, nna)
    ! Adjust kstart and noba to exclude flagged baselines in
    ! either end of the interval.  If there are gaps larger than
    ! one minute, kstart and noba will reflect only the first segment
    ! up until the gap.
    integer, intent(inout) :: kstart, kstop
    integer, intent(out) :: noba
    integer, intent(in) :: idet
    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)
    integer :: gapmin, gaplen, i

    if (all(nna(:, :, kstart+1:kstop, idet) == 0)) then
       noba = 0
       return
    end if

    ! trim the beginning

    do while (kstart < kstop)
       if (any(nna(:, :, kstart+1, idet) /= 0)) exit
       kstart = kstart + 1
    end do

    ! trim the end

    do while (kstop > kstart)
       if (any(nna(:, :, kstop, idet) /= 0)) exit
       kstop = kstop - 1
    end do

    ! now determine if there are gaps

    noba = kstop - kstart
    gapmin = 60 * fsample / dnshort + 1
    gaplen = 0
    do i = kstart + 1, kstop
       if (any(nna(:, :, i, idet) /= 0)) then
          gaplen = 0
       else
          gaplen = gaplen + 1
          if (gaplen == gapmin) then
             ! Adjust noba to reach the beginning of the gap and return
             noba = i - gaplen - kstart
             return
          end if
       end if
    end do

  end subroutine trim_interval


  subroutine convolve_interval(ichunk, idet, x, y, xx, fx, fc, nna)
    integer, intent(in) :: ichunk, idet
    real(dp), intent(out) :: y(0:basis_order, noba_short, nodetectors)
    real(dp), intent(in) :: x(0:basis_order, noba_short, nodetectors)
    complex(dp), intent(in) :: fc(nof/2 + 1, npsdtot)
    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)

    real(C_DOUBLE) :: xx(nof)
    complex(C_DOUBLE_COMPLEX) :: fx(nof/2 + 1)
    integer :: noba, kstart, kstop, ipsd

    kstart = sum(noba_short_pp(1:ichunk-1))
    kstop = kstart + noba_short_pp(ichunk)
    ipsd = psd_index(idet, baselines_short_time(kstart+1))

    y(0, kstart+1:kstop, idet) = x(0, kstart+1:kstop, idet)

    if (ipsd == -1) then
       ! no PSD
       return
    end if

    ! Only apply the filter to the unflagged segments of the baseline vector

    do
       call trim_interval(kstart, kstop, noba, idet, nna)
       if (noba == 0) exit

       xx(1:noba) = x(0, kstart+1:kstart+noba, idet)
       xx(1:noba) = xx(1:noba) - sum(xx(1:noba)) / noba
       xx(noba+1:) = 0
       call dfft(fx, xx)
       fx = fx * fc(:, ipsd)
       call dfftinv(xx, fx)
       xx(1:noba) = xx(1:noba) - sum(xx(1:noba)) / noba
       y(0, kstart+1:kstart+noba, idet) = xx(1:noba)

       kstart = kstart + noba
    end do

  end subroutine convolve_interval



  SUBROUTINE convolve_pp(y, x, fc, nna)
    ! Convolve a baseline vector with the noise prior, one pointing period
    ! at a time. FIXME: should have an option to apply the filter across the
    ! boundaries.

    real(dp), intent(out) :: y(noba_short, nodetectors)
    real(dp), intent(in) :: x(noba_short, nodetectors)
    complex(dp), intent(in) :: fc(nof/2+1, npsdtot)
    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)

    integer :: ichunk, idet, m, ierr, id_thread, itask, num_threads

    real(C_DOUBLE), pointer :: xx(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: fx(:) => NULL()
    type(C_PTR) :: pxx, pfx

    call reset_time(14)

    y = x

    !$OMP PARALLEL NUM_THREADS(nthreads) &
    !$OMP     DEFAULT(NONE) &
    !$OMP     SHARED(nof, nodetectors, ninterval, noba_short_pp, &
    !$OMP         baselines_short_time, fc, x, y, nthreads, nna) &
    !$OMP     PRIVATE(idet, ichunk, m, xx, fx, ierr, id_thread, itask, &
    !$OMP         num_threads, pxx, pfx)

    id_thread = omp_get_thread_num()
    num_threads = omp_get_num_threads()

    !allocate(fx(nof/2+1), xx(nof), stat=ierr)
    !if (ierr /= 0) call abort_mpi('No room for Fourier transform')
    pxx = fftw_alloc_real(int(nof, C_SIZE_T))
    pfx = fftw_alloc_complex(int(nof/2 + 1, C_SIZE_T))
    call c_f_pointer(pxx, xx, [nof])
    call c_f_pointer(pfx, fx, [nof/2 + 1])

    itask = -1
    do idet = 1, nodetectors
       do ichunk = 1, ninterval
          itask = itask + 1
          if (modulo(itask, num_threads) /= id_thread) cycle
          call convolve_interval(ichunk, idet, x, y, xx, fx, fc, nna)
       end do
    end do

    !deallocate(fx, xx)
    call fftw_free(pxx)
    call fftw_free(pfx)

    !$OMP END PARALLEL

    cputime_filter = cputime_filter + get_time(14)

  END SUBROUTINE convolve_pp



  !---------------------------------------------------------------------------


  SUBROUTINE construct_preconditioner(nna)

    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)
    integer :: i, j, k, kstart, kstop, n, noba, idet, ichunk, ipsd, ierr, try, &
         ipsddet, nempty, nfail, nband, itask, id_thread, num_threads, isub
    integer, parameter :: trymax = 1000
    integer :: ntries(trymax)
    real(dp) :: memsum, mem_min, mem_max, p, p_tot
    real(dp), parameter :: plim = 1 - 1e-5

    real(C_DOUBLE), pointer :: xx(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: fx(:) => NULL()
    type(C_PTR) :: pxx, pfx

    if (precond_width_max < 1) return

    call reset_time(16)

    if (ID == 0 .and. info > 1) write(*,*) 'Constructing preconditioner'
    if (info > 4) write(*, idf) ID, 'Constructing preconditioner'

    memory_precond = 0

    if (.not. kfilter .or. precond_width_max == 1) then
       use_diagonal = .true.
       if (.not. allocated(prec_diag)) then
          allocate(prec_diag(noba_short, nodetectors), stat=ierr)
          if (ierr /= 0) call abort_mpi('No room for prec_diag')
       end if
       memory_precond = noba_short*nodetectors*8.

       prec_diag = 0
       do idet = 1, nodetectors
          do ichunk = 1, ninterval
             kstart = sum(noba_short_pp(1:ichunk-1))
             noba = noba_short_pp(ichunk)
             ipsddet = psd_index_det(idet, baselines_short_time(kstart+1))
             if (ipsddet < 0) then
                ! No PSD available
                do k = kstart+1, kstart+noba
                   prec_diag(k, idet) = 1.
                end do
                cycle
             end if
             do k = kstart+1, kstart+noba
                if (nna(0, 0, k, idet) == 0) then
                   prec_diag(k, idet) = 1. / detectors(idet)%weights(ipsddet)
                else
                   prec_diag(k, idet) = 1. / nna(0, 0, k, idet)
                end if
             end do
          end do
       end do
       return
    end if

    ! Build fprecond to act as a fall back option
    if (allocated(fprecond)) deallocate(fprecond)
    allocate(fprecond(nof/2+1, npsdtot), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for fprecond')
    fprecond = 0
    memory_precond = memory_precond + (nof/2+1)*npsdtot*16.
    ipsd = 0
    do idet = 1, nodetectors
       do ipsddet = 1, detectors(idet)%npsd
          ipsd = ipsd + 1
          fprecond(:, ipsd) = 1.d0 /( &
               nshort*detectors(idet)%weights(ipsddet) + fcov(:, ipsd))
       end do
    end do

    if (use_fprecond) then
       return
    end if

    use_diagonal = .false.

    if (allocated(invcov)) deallocate(invcov)
    allocate(invcov(nof, npsdtot), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for invcov')
    invcov = 0

    call allocate_bandprec

    !$OMP PARALLEL NUM_THREADS(nthreads) &
    !$OMP     DEFAULT(NONE) &
    !$OMP     SHARED(nof, npsdtot, fcov, invcov) &
    !$OMP     PRIVATE(id_thread, num_threads, ipsd, pxx, pfx, xx, fx)

    id_thread = omp_get_thread_num()
    num_threads = omp_get_num_threads()

    pxx = fftw_alloc_real(int(nof, C_SIZE_T))
    pfx = fftw_alloc_complex(int(nof/2 + 1, C_SIZE_T))
    call c_f_pointer(pxx, xx, [nof])
    call c_f_pointer(pfx, fx, [nof/2 + 1])

    do ipsd = 1, npsdtot
       if (modulo(ipsd-1, num_threads) /= id_thread) cycle
       fx = fcov(:, ipsd)
       call dfftinv(xx, fx) ! C_a inverse into real domain
       invcov(:, ipsd) = xx
    end do

    call fftw_free(pxx)
    call fftw_free(pfx)

    !$OMP END PARALLEL

    if (.not. use_cgprecond) then
       nempty = 0
       nfail = 0
       ntries = 0
       memory_precond = 0

       !$OMP PARALLEL NUM_THREADS(nthreads) &
       !$OMP     DEFAULT(NONE) &
       !$OMP     SHARED(nodetectors, ninterval, noba_short_pp, &
       !$OMP         baselines_short_time, invcov, id, nof, &
       !$OMP         baselines_short_start, baselines_short_stop, bandprec, &
       !$OMP         nna, detectors, nthreads, precond_width_min, &
       !$OMP         precond_width_max) &
       !$OMP     PRIVATE(idet, ichunk, noba, kstart, kstop, ipsd, i, j, k, n, &
       !$OMP         ierr, try, nband, id_thread, itask, num_threads, p, p_tot, &
       !$OMP         isub) &
       !$OMP     REDUCTION(+:memory_precond, nempty, ntries, nfail)

       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       itask = -1
       do idet = 1, nodetectors
          do ichunk = 1, ninterval
             itask = itask + 1
             if (modulo(itask, num_threads) /= id_thread) cycle

             kstart = sum(noba_short_pp(1:ichunk-1))
             kstop = kstart + noba_short_pp(ichunk)
             isub = 1
             loop_subchunk : do
                call trim_interval(kstart, kstop, noba, idet, nna)
                if (noba == 0) exit
                ipsd = psd_index(idet, baselines_short_time(kstart+1))

                if (ipsd == -1 .or. noba < 1 .or. &
                     all(nna(0, 0, kstart+1:kstart+noba, idet) == 0)) then
                   nempty = nempty + 1
                   cycle
                end if

                ! Rough guess for the width of the preconditioner based on
                ! total power in the band

                p_tot = sum(invcov(:nof/2, ipsd)**2)
                try = 1
                ierr = -1
                loop_try : do
                   ! Rows from kstart+1 to kstart+noba for one band-diagonal
                   ! submatrix.  Do the computation in double precision
                   nband = min(noba, try * precond_width_min)
                   nband = min(nband, precond_width_max)
                   p = sum(invcov(:nband, ipsd)**2)
                   if (p > plim * p_tot) then
                      ierr = 0
                      exit
                   end if
                   if (nband == noba  .or. nband == precond_width_max) exit
                   if (try == trymax) exit
                   try = try + 1
                end do loop_try

                if (ierr == 0) then
                   ! Compute the Cholesky decomposition
                   do
                      call get_cholesky_decomposition( &
                           bandprec(isub, ichunk, idet)%data, nband, noba, &
                           invcov(1, ipsd), nna(0, 0, kstart+1, idet), ierr)

                      if (ierr == 0) exit

                      ! Not positive definite

                      if (try == trymax) exit
                      if (nband == noba .or. nband == precond_width_max) exit

                      ! Increase preconditioner width

                      try = try + 1
                      nband = min(nband + precond_width_min, precond_width_max)
                   end do
                end if

                if (ierr /= 0) then
                   !write (*, '(1x,a,i0,a,es18.10,a,es18.10)') &
                   !     'Cholesky decomposition failed for ' // &
                   !     trim(detectors(idet)%name) // ', noba = ', noba, ', t = ', &
                   !     baselines_short_time(kstart+1), ' - ', &
                   !     baselines_short_time(kstart+noba)
                   nfail = nfail + 1
                else
                   ntries(try) = ntries(try) + 1
                   memory_precond = memory_precond + nband*noba*8
                end if

                kstart = kstart + noba
                isub = isub + 1
             end do loop_subchunk
          end do
       end do
       !$OMP END PARALLEL

       call sum_mpi(ntries)
       call sum_mpi(nempty)
       call sum_mpi(nfail)
       if (id == 0) then
          print *, 'Constructed ', sum(ntries), ' band preconditioners.'
          do try = 1, trymax
             if (ntries(try) > 0) print *, ntries(try), ' at width = ', &
                  min(try*precond_width_min, precond_width_max)
          end do
          if (nfail > 0) print *, '            ', nfail, ' failed (using CG)'
          if (nempty > 0) print *, '    Skipped ', nempty, ' empty intervals.'
       end if

    end if

    memsum = memory_precond / 2**20
    mem_min = memsum; mem_max = memsum
    call min_mpi(mem_min); call max_mpi(mem_max)
    call sum_mpi(memsum)
    if (ID == 0 .and. info > 0) then
       write(*,'(1x,a,t32,3(f9.1," MB"))') &
            'Allocated memory for precond:', memsum, mem_min, mem_max
    end if

    cputime_prec_construct = cputime_prec_construct + get_time(16)

    if (info > 5) write(*, idf) ID, 'Done'

  END SUBROUTINE construct_preconditioner


  !------------------------------------------------------------------------------


  SUBROUTINE preconditioning_band(z, r, nna)
    !Apply the preconditioner

    ! z and r may be 3-dimensional arrays in the calling code but
    ! here we collapse the first dimension.  The band preconditioner
    ! only works for basis_order == 0.
    real(dp), intent(out), target :: z(0:basis_order, noba_short, nodetectors)
    real(dp), intent(in), target :: r(0:basis_order, noba_short, nodetectors)
    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)

    integer :: j, k, idet, kstart, kstop, noba, ichunk, ierr, itask, &
         num_threads, m, nband, ipsd, isub
    real(dp) :: t1, t2, tf

    real(C_DOUBLE), pointer :: xx(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: fx(:) => NULL()
    type(C_PTR) :: pxx, pfx

    z = r

    if (precond_width_max < 1) then
       return
    end if

    call reset_time(16)

    if (use_diagonal) then
       z(0, :, :) = r(0, :, :) * prec_diag
    else if (use_fprecond) then
       ! Apply the filter-based preconditioner
       call convolve_pp(z, r, fprecond, nna)
    else
       !$OMP PARALLEL NUM_THREADS(nthreads) &
       !$OMP     DEFAULT(NONE) &
       !$OMP     SHARED(noba_short_pp, ninterval, bandprec, z, r, nodetectors, &
       !$OMP         id, nof, nshort, detectors, nthreads, fprecond, nna, &
       !$OMP         baselines_short_time, invcov, fcov, checknan) &
       !$OMP     PRIVATE(idet, ichunk, kstart, kstop, noba, j, k, ierr, m, &
       !$OMP         xx, fx, nband, id_thread, itask, num_threads, ipsd, isub, &
       !$OMP         t1, t2, tf, pxx, pfx)

       !t1 = get_time(55)
       !tf = 0
       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       !allocate(fx(nof/2+1), xx(nof), stat=ierr)
       !if (ierr /= 0) call abort_mpi('No room for Fourier transform')
       pxx = fftw_alloc_real(int(nof, C_SIZE_T))
       pfx = fftw_alloc_complex(int(nof/2 + 1, C_SIZE_T))
       call c_f_pointer(pxx, xx, [nof])
       call c_f_pointer(pfx, fx, [nof/2 + 1])

       itask = -1
       do idet = 1, nodetectors
          do ichunk = 1, ninterval
             kstart = sum(noba_short_pp(1:ichunk-1))
             kstop = kstart + noba_short_pp(ichunk)
             isub = 1
             loop_subchunk : do
                call trim_interval(kstart, kstop, noba, idet, nna)
                if (noba == 0) exit loop_subchunk
                itask = itask + 1
                if (modulo(itask, num_threads) == id_thread) then
                   ipsd = psd_index_det(idet, baselines_short_time(kstart+1))
                   if (ipsd < 0) cycle
                   if (detectors(idet)%weights(ipsd) == 0) cycle
                   ipsd = psd_index(idet, baselines_short_time(kstart+1))
                   if (allocated(bandprec(isub, ichunk, idet)%data)) then
                      ! Use the precomputed Cholesky decomposition
                      nband = size(bandprec(isub, ichunk, idet)%data, 1)
                      call apply_cholesky_decomposition( &
                           noba, z(0, kstart+1, idet), r(0, kstart+1, idet), &
                           nband, bandprec(isub, ichunk, idet)%data)
                   else
                      ! Use CG iteration to apply the preconditioner
                      call apply_CG_preconditioner( &
                           noba, nof, fcov(1, ipsd), nna(0, 0, kstart+1, idet), &
                           z(0, kstart+1, idet), r(0, kstart+1, idet), &
                           invcov(1, ipsd), fx, xx, tf)
                   end if
                end if
                kstart = kstart + noba
                isub = isub + 1
             end do loop_subchunk
          end do
       end do

       !deallocate(fx, xx)
       call fftw_free(pxx)
       call fftw_free(pfx)

       !t2 = get_time(55)
       !print *, id, ' : ', id_thread, ' : applied preconditioner in = ', &
       !     t2 - t1, ' s. Filter time = ', tf
       !$OMP END PARALLEL
    end if

    cputime_precond = cputime_precond + get_time(16)

  END SUBROUTINE preconditioning_band


  subroutine apply_cholesky_decomposition(noba, z, r, nband, blockm)
    integer, intent(in) :: noba
    real(dp), intent(out) :: z(noba)
    real(dp), intent(in) :: r(noba)
    integer, intent(in) :: nband
    real(dp), intent(in) :: blockm(nband, noba)
    integer :: ierr

    z = r
    call DPBTRS('L', noba, nband-1, 1, blockm, nband, z, noba, ierr)
    if (ierr /= 0) then
       print *, id, ' : failed to Cholesky solve. argument # ', &
            -ierr, ' had an illegal value'
       call abort_mpi('bad preconditioner2')
    end if

  end subroutine apply_cholesky_decomposition



  subroutine get_cholesky_decomposition(blockm, nband, noba, invcov, nna, ierr)
    !
    ! Calculate the Cholesky decomposition of a band-diagonal matrix
    !
    real(dp), allocatable, intent(out) :: blockm(:, :)
    integer, intent(in) :: nband, noba
    real(dp), intent(in) :: invcov(noba), nna(noba)
    integer, intent(out) :: ierr

    if (allocated(blockm)) deallocate(blockm)
    allocate(blockm(nband, noba), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for Cholesky decomposition')

    blockm = spread(invcov(1:nband), 2, noba)
    blockm(1, :) = blockm(1, :) + nna

    ! Cholesky decompose

    call DPBTRF('L', noba, nband-1, blockm, nband, ierr)

    if (ierr /= 0) deallocate(blockm)

  end subroutine get_cholesky_decomposition



  subroutine apply_CG_preconditioner( &
       noba, nof, fcov, nna, x, b, invcov, fx, xx, tf)
    !
    ! Apply the preconditioner using CG iteration
    !
    integer, intent(in) :: noba, nof
    complex(dp), intent(in) :: fcov(nof/2+1)
    real(dp), intent(in) :: nna(noba), invcov(noba)
    real(dp), intent(out) :: x(noba)
    real(dp), intent(in) :: b(noba)
    real(dp), intent(inout) :: tf
    real(C_DOUBLE), intent(inout) :: xx(nof)
    complex(C_DOUBLE_COMPLEX), intent(inout) :: fx(nof/2 + 1)

    real(dp), allocatable :: resid(:), prop(:), Aprop(:), zresid(:), precond(:)
    real(dp) :: alpha, beta, rr, rr0, rr_old, rz, rz_old, Anorm
    integer :: iter, ierr
    ! iteration limits for applying the preconditioner.  We do not need
    ! perfect convergence because the top level CG iteration uses the
    ! Polak-Ribiere formula to allow for a changing preconditioner
    real(dp), parameter :: cglimit = 1e-6
    integer, parameter :: itermax = 100

    allocate(resid(noba), prop(noba), Aprop(noba), zresid(noba), &
         precond(noba), stat=ierr)
    if (ierr /= 0) call abort_mpi('not enough memory to apply preconditioner')

    x = 0
    resid = b - x
    rr = dot_product(resid, resid)
    if (rr == 0) return
    rr0 = rr

    precond = nna + invcov(1)
    where (precond /= 0) precond = 1 / precond
    call apply_precond(resid, zresid)
    rz = dot_product(resid, zresid)

    prop = zresid
    iter = 1
    do
       call apply_A(prop, Aprop, tf)
       Anorm = dot_product(prop, Aprop)
       alpha = rz / Anorm
       x = x + alpha * prop
       resid = resid - alpha * Aprop
       rr_old = rr
       rr = dot_product(resid, resid)
       if (rr / rr0 < cglimit .or. iter == itermax) exit
       call apply_precond(resid, zresid)
       rz_old = rz
       rz = dot_product(resid, zresid)
       beta = rz / rz_old
       rz_old = rz
       prop = zresid + beta * prop
       iter = iter + 1
    end do

    deallocate(resid, prop, Aprop, zresid, precond)

  contains

    subroutine apply_precond(x, z)
      !
      ! Apply the preconditioner to `x` and return the result in `z`.
      !
      real(dp), intent(in) :: x(noba)
      real(dp), intent(out) :: z(noba)

      z = precond * x
    end subroutine apply_precond

    subroutine apply_A(x, y, tf)
      !
      ! Apply the square matrix `A` to `x` and place the result in `y`
      !
      real(dp), intent(in) :: x(noba)
      real(dp), intent(out) :: y(noba)
      real(dp), intent(inout) :: tf
      call convolve(x, fcov, y, tf)
      y = y + nna*x
    end subroutine apply_A

    subroutine convolve(x, fc, y, tf)
      !
      ! Convolve `x` with filter `fc` and place result in `y`.
      !
      real(dp), intent(in) :: x(noba)
      complex(dp) :: fc(nof/2 + 1)
      real(dp), intent(out) :: y(noba)
      real(dp), intent(inout) :: tf
      real(dp) :: t1, t2
      xx(:noba) = x - sum(x)/noba
      xx(noba+1:) = 0
      !t1 = get_time(56)
      call dfft(fx, xx)
      fx = fx*fc
      call dfftinv(xx, fx)
      !t2 = get_time(56)
      !tf = tf + t2 - t1
      y = xx(:noba)
      y = y - sum(y) / noba
    end subroutine convolve

  end subroutine apply_CG_preconditioner



  !---------------------------------------------------------------------------



END MODULE noise_routines
