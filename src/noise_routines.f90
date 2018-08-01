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

  real(dp), allocatable :: xx(:)
  complex(dp), allocatable :: fx(:)
  complex(dp), allocatable, target :: fcov(:, :), fprecond(:, :)

  !Store interpolated spectra
  real(dp), allocatable :: freqtab(:), spectab(:, :)

  !Preconditioner

  logical :: use_diagonal = .false.
  !real(dp), allocatable, target :: bandprec(:, :, :)
  type bandprec_pointer
     real(dp), allocatable :: data(:, :)
  end type bandprec_pointer
  type(bandprec_pointer), allocatable :: bandprec(:, :)
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
    if (ierr /= 0) call abort_mpi('No room to interpolate the PSD')

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

    nofilter = noba_short_pp_max

    nof_min = nofilter
    ! Ensure that the filter is long enough for the preconditioner
    if (precond_width_max > nof_min-1) nof_min = precond_width_max + 1

    nof = 1
    do
       if (nof >= nof_min) exit
       nof = 2*nof
    end do
    if (ID == 0) write(*,*) 'FFT length =', nof

    allocate(fx(nof/2+1), xx(nof), stat=ierr) ! Work space
    if (ierr /= 0) call abort_mpi('No room for fx and xx')

    memory_filter = (nof/2+1)*16. + nof*8

    call init_fourier(nof)

  END SUBROUTINE init_filter


  !-----------------------------------------------------------------------


  subroutine allocate_bandprec
    integer :: ierr

    call free_bandprec
    allocate(bandprec(ninterval, nodetectors), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for bandprec')
  end subroutine allocate_bandprec



  subroutine free_bandprec
    integer :: ichunk, idet

    if (allocated(bandprec)) then
       do ichunk = 1, ninterval
          do idet = 1, nodetectors
             if (allocated(bandprec(ichunk, idet)%data)) then
                deallocate(bandprec(ichunk, idet)%data)
             end if
          end do
       end do
       deallocate(bandprec)
    end if
  end subroutine free_bandprec



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
       write(*,'(x,a,t32,3(f9.1," MB"))') &
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

          aspec = 1.d0 / (fa*aspec)
          aspec(1) = 0
          fcov(:, ipsdtot) = aspec

          ! RK edit: check for NaNs
          do i = 1, nof/2+1
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


  subroutine trim_interval(kstart, noba, idet, nna)
    ! Adjust kstart and noba to exclude flagged baselines in
    ! either end of the interval
    integer, intent(inout) :: kstart, noba
    integer, intent(in) :: idet
    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)

    if (all(nna(:, :, kstart+1:kstart+noba, idet) == 0)) then
       noba = 0
       return
    end if

    do while (noba > 0)
       if (any(nna(:, :, kstart+1, idet) /= 0)) exit
       kstart = kstart + 1
       noba = noba - 1
    end do

    do while (noba > 0)
       if (any(nna(:, :, kstart+noba, idet) /= 0)) exit
       noba = noba - 1
    end do

    return
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
    integer :: noba, kstart, ipsd

    noba = noba_short_pp(ichunk)
    kstart = sum(noba_short_pp(1:ichunk-1))
    ipsd = psd_index(idet, baselines_short_time(kstart+1))

    y(0, kstart+1:kstart+noba, idet) = x(0, kstart+1:kstart+noba, idet)

    if (ipsd == -1) then
       ! no PSD
       return
    end if

    ! Only apply the filter to the unflagged center of the baseline vector

    call trim_interval(kstart, noba, idet, nna)

    xx(1:noba) = x(0, kstart+1:kstart+noba, idet)
    xx(1:noba) = xx(1:noba) - sum(xx(1:noba)) / noba
    xx(noba+1:) = 0
    call dfft(fx, xx)
    fx = fx * fc(:, ipsd)
    call dfftinv(xx, fx)
    xx(1:noba) = xx(1:noba) - sum(xx(1:noba)) / noba
    y(0, kstart+1:kstart+noba, idet) = xx(1:noba)

    return
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

    integer :: ichunk, idet, m, no, noba, kstart, ipsd, ierr, id_thread, &
         itask, num_threads
    real(dp) :: x0

    real(C_DOUBLE), pointer :: xx(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: fx(:) => NULL()

    call reset_time(14)

    y = 0

    !$OMP PARALLEL NUM_THREADS(nthreads) &
    !$OMP     DEFAULT(NONE) &
    !$OMP     SHARED(nof, nodetectors, ninterval, noba_short_pp, &
    !$OMP         baselines_short_time, fc, x, y, nthreads, nna) &
    !$OMP     PRIVATE(idet, ichunk, noba, kstart, x0, m, no, xx, fx, &
    !$OMP         ipsd, ierr, id_thread, itask, num_threads)

    id_thread = omp_get_thread_num()
    num_threads = omp_get_num_threads()

    allocate(fx(nof/2+1), xx(nof), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for Fourier transform')

    itask = -1
    do idet = 1, nodetectors
       do ichunk = 1, ninterval
          itask = itask + 1
          if (modulo(itask, num_threads) /= id_thread) cycle
          call convolve_interval(ichunk, idet, x, y, xx, fx, fc, nna)
       end do
    end do

    deallocate(fx)
    deallocate(xx)

    !$OMP END PARALLEL

    cputime_filter = cputime_filter + get_time(14)

  END SUBROUTINE convolve_pp



  !---------------------------------------------------------------------------


  SUBROUTINE construct_preconditioner(nna)

    real(dp), intent(in)  :: &
         nna(0:basis_order, 0:basis_order, noba_short, nodetectors)
    integer :: i, j, k, kstart, n, noba, idet, ichunk, ipsd, ierr, try, ipsddet
    real(dp), allocatable :: invcov(:, :)
    integer :: nempty, nfail, nband, itask, id_thread, num_threads
    integer, parameter :: trymax = 1000
    integer :: ntries(trymax)
    real(dp) :: memsum, mem_min, mem_max, p, p_tot

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

    allocate(invcov(nof, npsdtot), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for invcov')
    invcov = 0

    call allocate_bandprec

    do ipsd = 1, npsdtot
       call dfftinv(xx, fcov(:, ipsd)) ! C_a inverse into real domain
       ! DEBUG begin
       !do i = 1, nof/2
       !   write (1000+100*id+ipsd, *) i, &
       !        real(fcov(i, ipsd)), imag(fcov(i, ipsd)), xx(i)
       !end do
       ! DEBUG end
       invcov(:, ipsd) = xx
    end do

    nempty = 0
    nfail = 0
    ntries = 0
    memory_precond = 0

    !$OMP PARALLEL NUM_THREADS(nthreads) &
    !$OMP     DEFAULT(NONE) &
    !$OMP     SHARED(nodetectors, ninterval, noba_short_pp, &
    !$OMP         baselines_short_time, invcov, id, sampletime, nof, &
    !$OMP         baselines_short_start, baselines_short_stop, bandprec, nna, &
    !$OMP         detectors, nthreads, precond_width_min, precond_width_max) &
    !$OMP     PRIVATE(idet, ichunk, noba, kstart, ipsd, i, j, k, n, ierr, try, &
    !$OMP         nband, id_thread, itask, num_threads, p, p_tot) &
    !$OMP     REDUCTION(+:memory_precond, nempty, ntries, nfail)

    id_thread = omp_get_thread_num()
    num_threads = omp_get_num_threads()

    itask = -1
    do idet = 1, nodetectors
       do ichunk = 1, ninterval
          itask = itask + 1
          if (modulo(itask, num_threads) /= id_thread) cycle

          noba = noba_short_pp(ichunk)
          kstart = sum(noba_short_pp(1:ichunk-1))
          call trim_interval(kstart, noba, idet, nna)
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
          loop_try : do
             ! Rows from kstart+1 to kstart+noba for one band-diagonal
             ! submatrix.  Do the computation in double precision
             nband = min(noba, try * precond_width_min)
             if (nband == noba  .or. nband == precond_width_max) exit
             p = sum(invcov(:nband, ipsd)**2)
             if (1 - p / p_tot < 1e-5) exit
             if (try == trymax) exit
             try = try + 1
          end do loop_try

          ! Compute the Cholesky decomposition

          do
             allocate(bandprec(ichunk, idet)%data(nband, noba), stat=ierr)
             if (ierr /= 0) call abort_mpi('No room for bandprec block')

             bandprec(ichunk, idet)%data = spread( &
                  invcov(1:nband, ipsd), 2, noba)
             bandprec(ichunk, idet)%data(1, :) = &
                  bandprec(ichunk, idet)%data(1, :) &
                  + nna(0, 0, kstart+1:kstart+noba, idet)

             ! Cholesky decompose

             call DPBTRF( &
                  'L', noba, nband-1, bandprec(ichunk, idet)%data, nband, ierr)

             if (ierr == 0) exit

             ! Not positive definite

             deallocate(bandprec(ichunk, idet)%data)

             if (try == trymax) exit
             if (nband == noba .or. nband == precond_width_max) exit

             ! Increase preconditioner width

             try = try + 1
             nband = min(nband + precond_width_min, precond_width_max)
          end do

          if (ierr /= 0) then
             print *,'Cholesky decomposition failed for ', &
                  trim(detectors(idet)%name), ',  noba = ', noba
             nfail = nfail + 1
!!$             ! DEBUG begin
!!$             write (1000+id, *) trim(detectors(idet)%name), ierr, noba
!!$             do i = 1, noba
!!$                write (1000+id, *) &
!!$                     i, invcov(i, ipsd), nna(0, 0, kstart+i, idet), &
!!$                     baselines_short_time(kstart+i)
!!$             end do
!!$             ! DEBUG end
          else
             ntries(try) = ntries(try) + 1
             memory_precond = memory_precond + nband*noba*8
          end if
       end do
    end do
    !$OMP END PARALLEL

    memsum = memory_precond / 2**20
    mem_min = memsum; mem_max = memsum
    call min_mpi(mem_min); call max_mpi(mem_max)
    call sum_mpi(memsum)
    if (ID == 0 .and. info > 0) then
       write(*,'(x,a,t32,3(f9.1," MB"))') &
            'Allocated memory for precond:', memsum, mem_min, mem_max
    end if

    call sum_mpi(ntries, trymax)
    call sum_mpi(nempty)
    call sum_mpi(nfail)
    if (id == 0) then
       print *, 'Constructed ', sum(ntries), ' band preconditioners.'
       do try = 1, trymax
          if (ntries(try) > 0) print *, ntries(try), ' at width = ', &
               try*precond_width_min
       end do
       print *, '            ', nfail, ' failed (using C_a)'
       print *, '    Skipped ', nempty, ' empty intervals.'
    end if

    deallocate(invcov)
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

    integer :: j, k, idet, kstart, noba, ichunk, ierr, itask, num_threads
    integer :: m, no, nband
    real(dp) :: x0

    real(C_DOUBLE), pointer :: xx(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: fx(:) => NULL()

    if (precond_width_max < 1) then
       z = r
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
       !$OMP         id, nof, nshort, detectors, nthreads, fprecond, nna) &
       !$OMP     PRIVATE(idet, ichunk, kstart, noba, j, k, ierr, x0, m, no, &
       !$OMP         xx, fx, nband, id_thread, itask, num_threads)

       id_thread = omp_get_thread_num()
       num_threads = omp_get_num_threads()

       allocate(fx(nof/2+1), xx(nof), stat=ierr)
       if (ierr /= 0) call abort_mpi('No room for Fourier transform')

       itask = -1
       do idet = 1, nodetectors
          do ichunk = 1, ninterval
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(1:ichunk-1))
             call trim_interval(kstart, noba, idet, nna)
             if (noba == 0) cycle

             itask = itask + 1
             if (modulo(itask, num_threads) /= id_thread) cycle

             if (allocated(bandprec(ichunk, idet)%data)) then
                ! Use the precomputed Cholesky decomposition
                nband = size(bandprec(ichunk, idet)%data, 1)
                z(0, kstart+1:kstart+noba, idet) = &
                     r(0, kstart+1:kstart+noba, idet)
                call DPBTRS( &
                     'L', noba, nband-1, 1, bandprec(ichunk, idet)%data, &
                     nband, z(0, kstart+1, idet), noba, ierr)
                if (ierr /= 0) then
                   print *, id, ' : failed to Cholesky solve. argument # ', &
                        -ierr, ' had an illegal value'
                   call abort_mpi('bad preconditioner2')
                end if
             else
                call convolve_interval( &
                     ichunk, idet, r, z, xx, fx, fprecond, nna)
             end if

          end do
       end do

       deallocate(fx)
       deallocate(xx)

       !$OMP END PARALLEL
    end if

    cputime_precond = cputime_precond + get_time(16)

  END SUBROUTINE preconditioning_band


  !---------------------------------------------------------------------------



END MODULE noise_routines
