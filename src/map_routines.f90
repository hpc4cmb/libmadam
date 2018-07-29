MODULE map_routines
  !
  ! Computation in matrix domain
  !
  use commonparam
  use matrix
  use mpi_wrappers
  use timing

  implicit none
  private

  real(dp), save, public :: cputime_inv=0
  character(len=15) :: stdstr(99,99)=''

  public stdstr, ccmultiply, invert_pixelmatrix_cross, &
       invert_pixelmatrix_map, makemaps, map_analysis

CONTAINS

  !-----------------------------------------------------------------------


  SUBROUTINE ccmultiply(cc, map, nopix)
    ! Multiply a map by the inverse of the pixel matrix

    real(dp),intent(in) :: cc(nmap, nmap, nopix)
    real(dp),intent(inout) :: map(nmap, nopix)
    integer, intent(in) :: nopix
    integer :: ip
    real(dp) :: tm, qm, um

    if (nmap == 1) then
       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(map, cc, nopix) PRIVATE(ip)
       do ip = 1, nopix
          if (cc(1, 1, ip) == 0) then
             map(:, ip) = 0
          else
             map(1, ip) = cc(1, 1, ip) * map(1, ip)
          end if
       end do
       !$OMP END PARALLEL DO
    else if (nmap == 3) then
       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(map, cc, nopix) &
       !$OMP     PRIVATE(ip, tm, qm, um)
       do ip = 1, nopix
          if (cc(1, 1, ip) == 0) then
             map(:, ip) = 0
          else
             tm = cc(1, 1, ip)*map(1, ip) + &
                  cc(2, 1, ip)*map(2, ip) + &
                  cc(3, 1, ip)*map(3, ip)
             qm = cc(2, 1, ip)*map(1, ip) + &
                  cc(2, 2, ip)*map(2, ip) + &
                  cc(3, 2, ip)*map(3, ip)
             um = cc(3, 1, ip)*map(1, ip) + &
                  cc(3, 2, ip)*map(2, ip) + &
                  cc(3, 3, ip)*map(3, ip)
             map(1, ip) = tm
             map(2, ip) = qm
             map(3, ip) = um
          end if
       end do
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO DEFAULT(NONE) SHARED(map, cc, nopix) PRIVATE(ip)
       do ip = 1,nopix
          if (cc(1, 1, ip) == 0) then
             map(:, ip) = 0
          else
             map(:, ip) = matmul(cc(:, :, ip), map(:, ip))
          end if
       end do
       !$OMP END PARALLEL DO
    end if

  END SUBROUTINE ccmultiply


  !-----------------------------------------------------------------------


  SUBROUTINE invert_pixelmatrix_cross(cc, inmask)
    !
    ! Invert the pixel matrices
    !
    real(dp), intent(inout) :: cc(nmap, nmap, 0:nopix_cross-1)
    integer, intent(in) :: inmask(0:nopix_cross-1)

    integer  :: ip, noeig, imap
    integer  :: cov0, cov1, covn
    real(dp) :: cdet, trace, tracen, rcnd
    logical  :: sing

    if (info == 3 .and. ID == 0) write(*,*) 'Inverting pixel matrices...'
    if (info > 4) write(*,idf) ID, 'Inverting pixel matrices...'

    call reset_time(10)

    if (use_inmask) then
       do ip = 0, nopix_cross-1
          if (inmask(ip) == 0) cc(:, :, ip) = 0
       end do
    end if

    cov0 = 0
    cov1 = 0
    covn = 0

    if (nmap == 1) then
       where (cc(1, 1, :) /= 0)
          cc(1, 1, :) = 1._dp / cc(1, 1, :)
       elsewhere
          cc(1, 1, :) = 0
       end where
       cov1 = count(cc(1, 1, :) /= 0)
    else
       do ip = 0, nopix_cross-1

          if (all(cc(:, :, ip) == 0)) cycle

          select case (pixmode_cross)
          case (0) ! absolute determinant
             call invert_LU(cc(:, :, ip), cdet, sing, nmap)

             if (cdet < pixlim_cross .or. sing) then
                cc(:, :, ip) = 0
             else
                covn = covn + 1
             end if
          case (1) ! scaled determinant
             trace = 0
             do imap = 1, nmap
                trace = trace + cc(imap, imap, ip)
             end do

             if (abs(trace) < 1.e-30) then
                tracen = 0
             else
                tracen = (nmap/trace)**nmap
             end if

             call invert_LU(cc(:, :, ip), cdet, sing, nmap)

             if (cdet*tracen < pixlim_cross .or. sing) then
                cc(:, :, ip) = 0
             else
                covn = covn + 1
             end if
          case (2) ! rcond
             call invert_eig(cc(:, :, ip), cdet, rcnd, sing, nmap)

             if (rcnd < pixlim_cross .or. sing) then
                cc(:, :, ip) = 0
             else
                covn = covn + 1
             end if
          case (3) ! use only T
             if (cc(1, 1, ip) > 1.e-20) then
                cc(1, 1, ip) = 1.d0 / cc(1, 1, ip)
                cc(:, 2:nmap, ip) = 0
                cc(2:nmap, 1, ip) = 0
                covn = covn + 1
             else
                cc(:, :, ip) = 0
             end if
          case (4) ! invert the non-singular part
             call invert_eig_s( &
                  cc(:, :, ip), cdet, rcnd, noeig, pixlim_cross, nmap)

             if (noeig == nmap) then
                covn = covn + 1
             else
                if (noeig > 0) cov1 = cov1 + 1
             end if
          end select
       end do
    end if

    cov0 = nopix_cross - cov1 - covn

    call sum_mpi(covn)
    call sum_mpi(cov1)
    call sum_mpi(cov0)

    if (ID == 0 .and. info > 1) then
       if (pixmode_cross == 4 .and. nmap > 1) then
          write(*,*)
          write(*,'(i9,a)') covn,' pixels fully solved'
          write(*,'(i9,a)') cov1,' pixels partly solved'
          write(*,'(i9,a)') cov0,' pixels unsolved'
       else
          write(*,*)
          write(*,'(i9,a)') cov1+covn,' pixels solved'
          write(*,'(i9,a)') cov0,     ' pixels unsolved'
       end if
    end if

    cputime_inv = cputime_inv + get_time(10)

    if (info > 4) write(*,idf) ID,'Done'

  END SUBROUTINE invert_pixelmatrix_cross


  !-----------------------------------------------------------------------


  SUBROUTINE invert_pixelmatrix_map(cc, mask, crit)
    ! Invert the pixel matrices

    real(dp),intent(inout) :: cc(nmap,nmap, 0:nopix_map-1)
    integer, intent(out) :: mask(0:nopix_map-1)
    real(sp),intent(out) :: crit(0:nopix_map-1)
    integer :: nstatic ! -RK

    integer  :: ip, cover, imap
    real(dp) :: cdet, trace, tracen, rcnd, pcrit
    logical  :: sing

    if (info == 3 .and. ID == 0) write(*,*) 'Inverting pixel matrices...'
    if (info > 4) write(*,idf) ID, 'Inverting pixel matrices...'

    call reset_time(10)

    cover = 0

    nstatic = max(1, ceiling(dble(nopix_map)/nthreads)) ! -RK
!!$OMP PARALLEL DEFAULT(PRIVATE) &
!!$OMP      SHARED(nopix_map,nmap,pixmode_map,cc,mask,pixlim_map,crit)&
!!$OMP      REDUCTION(+:cover)
    if (nmap == 1) then
!!$OMP DO SCHEDULE(STATIC,nstatic)
       do ip = 0,nopix_map-1
          if (cc(1, 1, ip) /= 0) then
             cc(1, 1, ip) = 1._dp / cc(1, 1, ip)
             mask(ip) = 1
             cover = cover + 1
          else
             cc(1, 1, ip) = 0
             mask(ip) = 0
          end if
       end do
!!$OMP END DO
    else
!!$OMP DO SCHEDULE(STATIC,nstatic)
       do ip = 0, nopix_map-1

          if (all(cc(:, :, ip) == 0)) cycle

          select case(pixmode_map)
          case (0) ! determinant
             call invert_LU(cc(:, :, ip), cdet, sing, nmap)
             pcrit = cdet
          case(1) ! scaled determinant
             trace = 0
             do imap = 1,nmap
                trace = trace + cc(imap, imap, ip)
             end do

             if (abs(trace) < 1.e-30) then
                tracen = 0
             else
                tracen = (nmap/trace)**nmap
             end if

             call invert_LU(cc(:, :, ip), cdet, sing, nmap)

             pcrit = cdet*tracen

          case(2) ! rcond
             call invert_eig(cc(:, :, ip), cdet, rcnd, sing, nmap)

             pcrit = rcnd
          case default
             stop 'Unknown pixmode_map'
          end select

          if (pcrit < pixlim_map .or. sing) then
             cc(:, :, ip) = 0
             mask(ip) = 0
          else
             mask(ip) = 1
             cover = cover + 1
          end if

          if (do_mask) crit(ip) = pcrit

       end do
!!$OMP END DO
    end if
!!$OMP END PARALLEL

    call sum_mpi(cover)

    if (ID == 0 .and. info >= 2) then
       write(*,*)
       write(*,'(i9,a)') cover,                ' pixels solved'
       write(*,'(i9,a)') 12*nside_map**2-cover,' pixels unsolved'
    end if

    cputime_inv = cputime_inv + get_time(10)

    if (info > 4) write(*,idf) ID,'Done'

  END SUBROUTINE invert_pixelmatrix_map


  !------------------------------------------------------------------------


  SUBROUTINE makemaps(map, cc, mask)
    !
    ! Construct the final output map
    real(dp),intent(inout) :: map(nmap, 0:nopix_map-1)
    real(dp),intent(in) :: cc(nmap,nmap, 0:nopix_map-1)
    integer, intent(in) :: mask(0:nopix_map-1)
    real(dp) :: monopole
    integer :: cover

    if (info == 3 .and. ID == 0) write(*,*) 'Constructing output map...'
    if (info > 4) write(*,idf) ID, 'Constructing output map...'

    call ccmultiply(cc, map, nopix_map)

    if (.not. rm_monopole) return

    cover = count(mask > 0)
    call sum_mpi(cover)

    if (cover == 0) return

    ! Remove the monopole from temperature map only
    monopole = sum(map(1, :))
    call sum_mpi(monopole)

    monopole = monopole/cover
    where(mask > 0) map(1, :) = map(1, :) - monopole

    if (info <= 5) write(*,idf) ID,'Done'

  END SUBROUTINE makemaps


  !-----------------------------------------------------------------------


  SUBROUTINE map_analysis(map, mask)
    ! Analyze the output map: Compute min, max, std etc,

    real(dp),intent(in) :: map(nmap, 0:nosubpix_map-1, nosubmaps)
    integer, intent(in) :: mask(0:nosubpix_map-1, nosubmaps)
    character(len=16) :: stokes(nmap)
    integer :: cover, imap

    cover = count(mask > 0)
    call sum_mpi(cover)

    do imap = 1, nmap
       write (stokes(imap), '(i9)') imap
    end do

    stdstr = ''

    if (id == 0 .and. info > 0) then
       write (*,*)
       write (*,*) 'Destriped map:'
    end if

    call write_std(stdstr(:, 1), map)

  CONTAINS

    SUBROUTINE write_std(stdstr, map)

      character(len=15) :: stdstr(nmap)
      real(dp),intent(in) :: map(nmap, 0:nosubpix_map-1, nosubmaps)
      integer :: i, j, k
      logical :: kcover(0:nosubpix_map-1)
      real(dp) :: summ(nmap), summ2(nmap), std(nmap)
      real(dp) :: mean(nmap), minv(nmap), maxv(nmap), dmap(nmap)
      character(len=15) :: str(nmap)

      summ = 0
      summ2 = 0
      dmap = 0

      do i = 1, nosubmaps
         do j = 0, nosubpix_map-1
            do k = 1, nmap
               dmap(k) = map(k, j, i)
            end do
            summ = summ + dmap
            summ2 = summ2 + dmap*dmap
         end do
      end do

      minv = 1.e30
      maxv =-1.e30
      do i = 1,nosubmaps
         kcover = (mask(:, i) > 0)
         do k = 1, nmap
            minv(k) = min(minv(k), minval(map(k, :, i), kcover))
            maxv(k) = max(maxv(k), maxval(map(k, :, i), kcover))
         end do
      end do

      do k = 1, nmap
         call sum_mpi(summ(k))
         call sum_mpi(summ2(k))
         call min_mpi(minv(k))
         call max_mpi(maxv(k))
      end do

      std = 0
      mean = 0
      if (cover > 1) then
         do k = 1, nmap
            mean(k) = summ(k) / cover
            std(k) = summ2(k) - summ(k)*summ(k)/cover
            std(k) = sqrt(std(k) / (cover-1))
         end do
      end if

      if (ID == 0 .and. info > 0) then
         write(*,'(a9,6x,99(2x,a))') 'Map     ',(stokes(k),k=1,nmap)

         call tempstr(stdstr, std)
         write(*,'(a9,2x,99(4x,a))') 'Std     ',(stdstr(k),k=1,nmap)

         if (info > 0) then
            call tempstr(str, mean)
            write(*,'(a9,2x,99(4x,a))') 'Mean    ',(str(k),k=1,nmap)

            call tempstr(str, minv)
            write(*,'(a9,2x,99(4x,a))') 'Min     ',(str(k),k=1,nmap)

            call tempstr(str, maxv)
            write(*,'(a9,2x,99(4x,a))') 'Max     ',(str(k),k=1,nmap)
         end if
         write(*,*)
      end if

    END SUBROUTINE write_std

    SUBROUTINE tempstr(str, t)

      real(dp) :: t(nmap), abst
      character(len=15) :: str(nmap)
      integer :: k

      str = ''
      do k = 1, nmap
         abst = abs(t(k))

         if (abst > 1.e-10 .and. abst <= 1.e-6) then
            write(str(k),'(f12.8," uK")') t(k)*1.e6

         elseif (abst > 1e-6 .and. abst < 0.001) then
            write(str(k),'(f12.6," uK")') t(k)*1.e6

         elseif (abst >= 0.001 .and. abst < 1) then
            write(str(k),'(f12.3," uK")') t(k)*1.e6

         elseif (abst >= 1 .and. abst < 1000) then
            write(str(k),'(f12.6," K")') t(k)

         elseif (abst >= 1000 .and. abst < 1000000) then
            write(str(k),'(f12.3," K")') t(k)

         else
            write(str(k),'(es12.5," K")') t(k)
         end if
      end do

    END SUBROUTINE tempstr

  END SUBROUTINE map_analysis


  !-----------------------------------------------------------------------

END MODULE map_routines
