
MODULE madam_routines

  use commonparam
  use mpi_wrappers
  use maptod_transfer
  use noise_routines
  use map_routines
  use pointing
  use timing
  use maps_and_baselines, only : nna_inv
  use output, only : write_leakmatrix

  implicit none
  private

  real(sp), save, public :: memory_cg = 0.0

  real(sp), save, public :: cputime_build_matrix = 0.0,  &
       cputime_send_matrix  = 0.0,  &
       cputime_bin_maps     = 0.0,  &
       cputime_send_maps    = 0.0,  &
       cputime_count_hits   = 0.0,  &
       cputime_cga          = 0.0,  &
       cputime_cga_init     = 0.0,  &
       cputime_cga_mpi      = 0.0,  &
       cputime_cga_1        = 0.0,  &
       cputime_cga_2        = 0.0,  &
       cputime_cga_cc       = 0.0,  &
       cputime_cgb          = 0.0,  &
       cputime_cgb_init     = 0.0,  &
       cputime_cgb_mpi      = 0.0,  &
       cputime_cgb_1        = 0.0,  &
       cputime_cgb_2        = 0.0,  &
       cputime_clean_tod    = 0.0,  &
       cputime_unclean_tod  = 0.0, &
       cputime_leakmatrix   = 0.0

  public pixel_matrix, bin_tod, &
       count_hits, initialize_a, iterate_a, &
       subtract_baselines_a, clean_tod, unclean_tod, &
       leakmatrix

CONTAINS

  !-------------------------------------------------------------------------


  SUBROUTINE leakmatrix(idet, cc, mask)
    !
    ! Build and write the detector-specific leakage matrix
    !
    integer, intent(in) :: idet
    real(dp),intent(in) :: cc(nmap,nmap,0:nopix_map-1)
    integer, intent(in), optional :: mask(nosubpix_map,nosubmaps)
    real(dp) :: cca_dummy(1,1,1)
    logical :: detflags_save(NDETMAX), kfirst_save
    real(dp), allocatable :: cc_det(:,:,:), loccc_save(:,:,:)
    integer :: ierr, ip

    if (info == 3 .and. id == 0) &
         write(*,*) 'Building leakage matrices...'
    if (info > 4) write(*,idf) ID, 'Building leakage matrices...'

    detflags_save = detflags
    detflags = .false.
    detflags(idet) = .true.
    kfirst_save = kfirst
    kfirst = .false.

    allocate(cc_det(nmap,nmap,0:nopix_map-1), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for cc_det')
    allocate(loccc_save(nmap,nmap,0:nsize_locmap), stat=ierr)
    if (ierr /= 0) call abort_mpi('No room for loccc_save')
    loccc_save = loccc

    call pixel_matrix(cca_dummy, cc_det)

    call reset_time(10)

    do ip = 0,nopix_map-1
       cc_det(:,:,ip) = matmul(cc(:,:,ip), cc_det(:,:,ip))
    end do
    
    cputime_leakmatrix = cputime_leakmatrix + get_time_and_reset(10)

    call write_leakmatrix(cc_det, mask, detectors(idet)%name)

    deallocate(cc_det)
    loccc = loccc_save
    deallocate(loccc_save)
    detflags = detflags_save
    kfirst = kfirst_save

    if (info > 4) write(*,idf) id,'Done'

  END SUBROUTINE leakmatrix


  !----------------------------------------------------------------------------


  SUBROUTINE pixel_matrix(cca, cc)
    !
    ! Build pixel matrices
    !
    real(dp),intent(out) :: cca(nmap,nmap,0:nopix_cross-1)
    real(dp),intent(inout) :: cc(nmap,nmap,0:nopix_map-1)
    integer :: i, ip, n, ierr, idet, num_threads, ichunk, noba, kstart, ipsd
    integer :: k, firstpix, lastpix
    real(dp) :: detweight, sqrtweight
    real(dp) :: w(nmap)

    if (info == 3 .and. id == 0) &
         write(*,*) 'Building pixel matrices...'
    if (info > 4) write(*,idf) ID, 'Building pixel matrices...'

    call reset_time(10)

    loccc = 0

    !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
    !$OMP   SHARED(nodetectors,nmap,first_chunk,last_chunk,noba_short_pp,&
    !$OMP          baselines_short_time,detectors,baselines_short_start,&
    !$OMP          baselines_short_stop,isubchunk,subchunkpp,pixels,id,npsdtot,&
    !$OMP          dummy_pixel,loccc,weights,nthreads,nsize_locmap,detflags,&
    !$OMP          surveyflags) &
    !$OMP   PRIVATE(firstpix,lastpix,i,ip,id_thread,num_threads,idet,detweight,&
    !$OMP           sqrtweight,w,ichunk,noba,kstart,ipsd,k)

    if ( id == 0 .and. nthreads /= omp_get_num_threads() ) &
         print *, 'WARNING : building pixel matrices with only ', &
         omp_get_num_threads(),' / ',nthreads,' threads'

    id_thread = omp_get_thread_num()
    firstpix = id_thread * nsize_locmap / omp_get_num_threads()
    lastpix = (id_thread+1) * nsize_locmap / omp_get_num_threads() - 1

    do idet = 1,nodetectors
       if ( .not. detflags(idet) ) cycle ! not included in this subset
       if (nmap == 1) then
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) then
                print *, id,' : WARNING: there is no PSD for det # ', idet, &
                     ' at ', baselines_short_time(kstart+1), '. npsd = ', &
                     detectors(idet)%npsd, ', start times = ', &
                     detectors(idet)%psdstarts
                cycle
             end if
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle       
             do k = kstart+1,kstart+noba
                do i = baselines_short_start(k), baselines_short_stop(k)
                   if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
                   if ( .not. surveyflags(i) ) cycle
                   ip = pixels(i,idet)
                   !if (ip == dummy_pixel) cycle
                   if ( ip < firstpix .or. ip > lastpix ) cycle
                   loccc(1,1,ip) = loccc(1,1,ip) + detweight
                end do
             end do
          end do
       else if ( nmap == 3 ) then
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) then
                print *, id,' : WARNING: there is no PSD for det # ', idet, &
                     ' at ', baselines_short_time(kstart+1), '. npsd = ', &
                     detectors(idet)%npsd, ', start times = ', &
                     detectors(idet)%psdstarts
                cycle
             end if
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle
             sqrtweight = sqrt(detweight)
             do k = kstart+1,kstart+noba
                do i = baselines_short_start(k), baselines_short_stop(k)
                   if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
                   if ( .not. surveyflags(i) ) cycle
                   ip = pixels(i,idet)
                   !if (ip == dummy_pixel) cycle
                   if ( ip < firstpix .or. ip > lastpix ) cycle
                   w = weights(1:nmap,i,idet) * sqrtweight
                   loccc(1,1,ip) = loccc(1,1,ip) + w(1)**2
                   loccc(1,2,ip) = loccc(1,2,ip) + w(1) * w(2)
                   loccc(1,3,ip) = loccc(1,3,ip) + w(1) * w(3)
                   loccc(2,2,ip) = loccc(2,2,ip) + w(2)**2
                   loccc(2,3,ip) = loccc(2,3,ip) + w(2) * w(3)
                   loccc(3,3,ip) = loccc(3,3,ip) + w(3)**2
                end do
             end do
          end do
       else
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) then
                print *, id,' : WARNING: there is no PSD for det # ', idet, &
                     ' at ', baselines_short_time(kstart+1), '. npsd = ', &
                     detectors(idet)%npsd, ', start times = ', &
                     detectors(idet)%psdstarts
                cycle
             end if
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle       
             do k = kstart+1,kstart+noba
                do i = baselines_short_start(k), baselines_short_stop(k)
                   if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
                   if ( .not. surveyflags(i) ) cycle
                   ip = pixels(i,idet)
                   !if (ip == dummy_pixel) cycle
                   if ( ip < firstpix .or. ip > lastpix ) cycle
                   call dsyr('U', nmap, detweight, weights(:,i,idet), 1, &
                        loccc(:,:,ip), nmap)
                end do
             end do
          end do
       end if
    end do
    !$OMP END PARALLEL

    cputime_build_matrix = cputime_build_matrix + get_time_and_reset(10)

    call collect_cc(cc, nosubpix_map)

    if (kfirst) call collect_cc(cca, nosubpix_cross)

    cputime_send_matrix = cputime_send_matrix + get_time(10)

    if (info > 4) write(*,idf) id,'Done'

  END SUBROUTINE pixel_matrix


  !----------------------------------------------------------------------------


  SUBROUTINE bin_tod(map, binmap, wamap, tod_stored)
    !
    ! Bin TOD onto map
    !
    real(dp),intent(inout) :: map(nmap,0:nopix_map-1)
    real(dp),intent(inout) :: binmap(nmap,0:nopix_map-1)
    real(dp),intent(inout) :: wamap(nmap,0:nopix_cross-1)
    real(dp), intent(in) :: tod_stored(nosamples_proc,nodetectors)
    integer :: i, n, ip, m, ierr, firstpix, lastpix, idet, ichunk, noba, kstart
    integer :: ipsd, k
    real(dp) :: pw, detweight

    if (info == 3 .and. ID == 0) write(*,*) 'Binning TOD...'

    call reset_time(10)

    locmap = 0

    !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
    !$OMP     PRIVATE(i,ip,id_thread,idet,detweight,ierr,&
    !$OMP             ichunk,noba,kstart,ipsd,k,firstpix,lastpix) &
    !$OMP     SHARED(nsize_locmap,nodetectors,detflags,nmap,first_chunk,&
    !$OMP            last_chunk,noba_short_pp,detectors,baselines_short_start,&
    !$OMP            baselines_short_stop,isubchunk,subchunkpp,pixels,&
    !$OMP            dummy_pixel,baselines_short_time,surveyflags,&
    !$OMP            tod_stored,weights,nthreads,locmap,id,nosamples_proc)

    if ( id == 0 .and. nthreads /= omp_get_num_threads() ) &
         print *,'WARNING : binning TOD with only ', omp_get_num_threads(), &
         ' / ', nthreads, ' threads'

    id_thread = omp_get_thread_num()
    firstpix = id_thread * nsize_locmap / omp_get_num_threads()
    lastpix = (id_thread+1) * nsize_locmap / omp_get_num_threads() - 1

    do idet = 1,nodetectors
       if ( .not. detflags(idet) ) cycle ! not included in this subset
       if (nmap == 1) then
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) cycle
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle       
             do k = kstart+1,kstart+noba
                do i = baselines_short_start(k), baselines_short_stop(k)
                   if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
                   if ( .not. surveyflags(i) ) cycle
                   ip = pixels(i,idet)
                   !if (ip == dummy_pixel) cycle
                   if ( ip < firstpix .or. ip > lastpix ) cycle
                   locmap(1,ip) = locmap(1,ip) + tod_stored(i,idet) * detweight
                end do
             end do
          end do
       else
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) cycle
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle       
             do k = kstart+1,kstart+noba
                do i = baselines_short_start(k), baselines_short_stop(k)
                   if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
                   if ( .not. surveyflags(i) ) cycle
                   ip = pixels(i,idet)
                   !if (ip == dummy_pixel) cycle
                   if ( ip < firstpix .or. ip > lastpix ) cycle
                   locmap(:,ip) = locmap(:,ip) &
                        + weights(:,i,idet) * tod_stored(i,idet) * detweight
                end do
             end do
          end do
       end if
    end do
    !$OMP END PARALLEL

    cputime_bin_maps = cputime_bin_maps + get_time_and_reset(10)

    if (kfirst) call collect_map(map, nosubpix_map)

    if (do_binmap) call collect_map(binmap, nosubpix_map)
    if (kfirst)    call collect_map(wamap, nosubpix_cross)

    cputime_send_maps = cputime_send_maps + get_time(10)

    if (info > 4) write(*,idf) ID,'Done'

  END SUBROUTINE bin_tod


  !---------------------------------------------------------------------------


  SUBROUTINE count_hits(nohits)

    ! Count the number of hits/pixel

    integer,intent(out) :: nohits(0:nopix_map-1,*)
    integer :: i, n, ip, idet, ichunk, noba, kstart, ipsd, k, firstpix, lastpix
    real(dp) :: detweight

    if (.not. do_hits) return

    if (info == 3 .and. id == 0) write(*,*) 'Counting hits...'
    if (info > 4) write(*,idf) id,'Counting hits'

    call reset_time(10)

    lochits = 0

    do idet = 1,nodetectors
       if ( .not. detflags(idet) ) cycle ! not included in this subset

       !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
       !$OMP   SHARED(nodetectors,first_chunk,last_chunk,noba_short_pp,&
       !$OMP          baselines_short_time,baselines_short_start,&
       !$OMP          baselines_short_stop,isubchunk,subchunkpp,surveyflags,&
       !$OMP          lochits,dummy_pixel,do_dethits,nohits,nosubpix_map,idet,&
       !$OMP          nsize_locmap,nthreads,detectors,pixels,id,nopix_map) &
       !$OMP   PRIVATE(ichunk,noba,kstart,ipsd,detweight,k,i,ip,firstpix,&
       !$OMP          lastpix,id_thread)
       
       if ( id == 0 .and. nthreads /= omp_get_num_threads() ) &
            print *, 'WARNING : counting hits with only ', &
            omp_get_num_threads(), ' / ', nthreads, ' threads'

       id_thread = omp_get_thread_num()
       firstpix = id_thread * nsize_locmap / omp_get_num_threads()
       lastpix = (id_thread+1) * nsize_locmap / omp_get_num_threads() - 1

       do ichunk = first_chunk, last_chunk
          noba = noba_short_pp(ichunk)
          kstart = sum(noba_short_pp(first_chunk:ichunk-1))
          ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
          if ( ipsd < 0 ) then
             print *, id,' : WARNING: there is no PSD for det # ', idet, &
                  ' at ', baselines_short_time(kstart+1), '. npsd = ', &
                  detectors(idet)%npsd, ', start times = ', &
                  detectors(idet)%psdstarts
             cycle
          end if
          detweight = detectors(idet)%weights(ipsd)
          if (detweight == 0) cycle
          do k = kstart+1,kstart+noba
             do i = baselines_short_start(k), baselines_short_stop(k)
                if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
                if ( .not. surveyflags(i) ) cycle
                ip = pixels(i,idet)
                if (ip == dummy_pixel) cycle
                if ( ip < firstpix .or. ip > lastpix ) cycle
                lochits(ip) = lochits(ip) + 1
             end do
          end do
       end do

       !$OMP END PARALLEL 

       if (do_dethits) then
          call collect_hits(nohits(0,idet), nosubpix_map)
          lochits = 0
       end if
    end do

    if (.not. do_dethits) call collect_hits(nohits(0,1), nosubpix_map)

    cputime_count_hits = cputime_count_hits + get_time(10)

    if (info > 4) write(*,idf) id,'Done'

  END SUBROUTINE count_hits


  !---------------------------------------------------------------------------


  SUBROUTINE initialize_a(yba, nna, wamap, cca, tod_stored)

    real(dp),intent(inout) :: yba(0:basis_order, noba_short_max, nodetectors)
    real(dp), intent(out) :: nna(0:basis_order, 0:basis_order, noba_short_max, nodetectors)
    real(dp),intent(inout) :: wamap(nmap,0:nopix_cross-1)
    real(dp),intent(in) :: cca(nmap,nmap,0:nopix_cross-1)
    real(dp), intent(in) :: tod_stored(nosamples_proc,nodetectors)
    real(dp) :: detweight, bf
    integer(i8b) :: i, k, ip, j, order, order2, i0, imap, n, workspace_length
    integer(i8b) :: ntot, nbad, ngood
    integer(i4b) :: ichunk, noba, kstart, ipsd, idet
    real(dp), pointer :: basis_function(:,:)
    REAL(dp), allocatable :: workspace(:), eigenvectors(:,:)
    REAL(dp), allocatable :: eigenvectorsT(:,:), eigenvalues(:)
    integer :: ierr
    logical :: bad_baseline

    if (info == 3 .and. id == 0) write(*,*) 'Building RHS...'
    if (info > 4) write(*,idf) id,'Building RHS...'

    call reset_time(10)

    ! globally, apply pixel matrices to the binned map
    call ccmultiply(cca, wamap, nopix_cross)
    ! global -> local, update locmap
    call scatter_map(wamap, nosubpix_cross)

    if (checknan) then
       loop_idet : do idet = 1,nodetectors
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) then
                print *, id, ' : WARNING: there is no PSD for det # ', idet, &
                     ' at ', baselines_short_time(kstart+1), '. npsd = ', &
                     detectors(idet)%npsd, ', start times = ', &
                     detectors(idet)%psdstarts
                cycle
             end if
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle
             do k = kstart+1,kstart+noba
                do i = baselines_short_start(k), baselines_short_stop(k)
                   if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                   do order = 0,basis_order
                      if (isnan(yba(order, k, idet))) then
                         print *,id,' : NaN in uninitialized yba'
                         exit loop_idet
                      end if
                   end do
                   ip = pixels(i,idet)
                   if (ip == dummy_pixel) cycle
                   !do ip = 0,nsize_locmap
                   do j = 1,nmap
                      if (isnan(locmap(j,ip))) then
                         print *,id,' : NaN in locmap, stokes = ',j
                         exit loop_idet
                      end if
                   end do
                   !end do
                   if (nmap /= 1) then
                      do imap = 1,nmap
                         if ( isnan(weights(imap,i,idet)) ) then
                            print *,id,' : NaN in weights'
                            exit loop_idet
                         end if
                      end do
                   end if
                end do
             end do
          end do
       end do loop_idet
    end if
    
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP    PRIVATE(idet,ichunk,noba,kstart,ipsd,detweight,&
    !$OMP            k,i0,basis_function,i,ip,order,bf) &
    !$OMP    SHARED(nodetectors,first_chunk,last_chunk,noba_short_pp,detectors,&
    !$OMP           yba,nna,baselines_short_start,basis_functions,&
    !$OMP           baselines_short_stop,isubchunk,subchunk,pixels,dummy_pixel,&
    !$OMP           locmap,basis_order,tod_stored,nmap,baselines_short_time,&
    !$OMP           weights,order2,checknan,id,noba_short_max,nosamples_proc)
    !$OMP DO SCHEDULE(DYNAMIC,1)
    do idet = 1,nodetectors

       do ichunk = first_chunk, last_chunk
          noba = noba_short_pp(ichunk)
          kstart = sum(noba_short_pp(first_chunk:ichunk-1))
          ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
          if ( ipsd < 0 ) cycle
          detweight = detectors(idet)%weights(ipsd)
          if (detweight == 0) cycle       
          do k = kstart+1,kstart+noba

             yba(:,k,idet) = 0
             nna(:,:,k,idet) = 0

             i0 = baselines_short_start(k)
             basis_function => basis_functions(k)%arr
             do i = baselines_short_start(k), baselines_short_stop(k)
                if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                ip = pixels(i,idet)
                if (ip == dummy_pixel) cycle
                
                ! Use of locmask rather than locmap==0 agrees with the old
                ! versions of Madam and does not exclude poorly conditioned
                ! pixels.
                if ( all(locmap(:,ip) == 0) ) then
                   cycle ! masked by inmask or pixel rejection
                end if
                !if ( use_inmask .and. locmask(ip) == 0 ) cycle
                
                do order = 0,basis_order
                   bf = basis_function( order, i-i0 )
                   yba(order, k, idet) = yba(order, k, idet) &
                        + bf * tod_stored(i,idet)
                   nna(:, order, k, idet) = nna(:, order, k, idet) &
                        + bf * basis_function(:,i-i0)
                end do

                ! This loop implicitly applies the locmap==0 (loccc==0) criterion
                do order = 0,basis_order
                   bf = basis_function( order, i-i0 )
                   if (nmap == 1) then
                      yba(order,k,idet) = yba(order,k,idet) - bf * locmap(1,ip)
                   else
                      yba(order,k,idet) = yba(order,k,idet) &
                           - bf * dot_product(locmap(:,ip), weights(:,i,idet))
                   end if
                end do
             end do

             yba(:,k,idet) = yba(:,k,idet) * detweight
             ! Stock version applied detweight later
             nna(:,:,k,idet) = nna(:,:,k,idet) * detweight

          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! invert the F^T F blocks for the diagonal preconditioner
    ! FIXME add timer for this
    ! FIXME handle singular matrices gracefully
    ntot = nodetectors * noba_short
    nbad = 0
    ngood = 0
    do idet = 1,nodetectors
       n = basis_order + 1
       !if (n == 1) then
       !   do k = 1, noba_short
       !      if (nna(1,1,k,idet) /= 0) nna_inv(1,1,k,idet) = 1 / nna(1,1,k,idet)
       !   end do
       !else
       workspace_length = max(5*n, 1000)
       allocate(eigenvalues(n), eigenvectors(n,n), eigenvectorsT(n,n), &
            workspace(workspace_length), stat=ierr)
       if (ierr /= 0) stop 'No room to invert F^T F'
       do ichunk = first_chunk, last_chunk
          noba = noba_short_pp(ichunk)
          kstart = sum(noba_short_pp(first_chunk:ichunk-1))
          ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
          if ( ipsd < 0 ) cycle
          detweight = detectors(idet)%weights(ipsd)
          if (detweight == 0) cycle       
          do k = kstart+1,kstart+noba
             nna_inv(:,:,k,idet) = 0
             if ( abs(nna(0,0,k,idet)) < 1e-30 ) then
                eigenvalues = 0
             else
                eigenvectors = nna(:, :, k, idet)
                call dsyev('V', 'U', n, eigenvectors, n, eigenvalues, &
                     workspace, workspace_length, ierr)
             end if
             bad_baseline = .true.
             if (ierr /= 0) then
                print *,id,' : detector: ',trim(detectors(idet)%name)
                print *,id,' : nna'
                do i = 1,n
                   print *,id," : ",nna(:,i,k,idet)
                end do
                !stop 'DSYEV failed!'
                nna(:,:,k,idet) = 0
                yba(:,k,idet) = 0
                nbad = nbad + 1
             else if (minval(eigenvalues) == 0 &
                  .and. maxval(eigenvalues) == 0) then
                nna(:,:,k,idet) = 0
                yba(:,k,idet) = 0
                cycle
             else if (minval(eigenvalues) / maxval(eigenvalues) < 1e-6) then
                ! Degenerate basis functions
                print *,id,' : rcond too small: ',eigenvalues
                print *,id,' : hits : ',nna(0,0,k,idet),' / ',baselines_short(k)
                print *,id,' : ',nna(:,:,k,idet)
                nna(:,:,k,idet) = 0
                yba(:,k,idet) = 0
                nbad = nbad + 1
             else
                bad_baseline = .false.
             end if

             if (bad_baseline) then
                if ( .not. kfilter ) then
                   call abort_mpi('Masking bad baselines not implemented. Increase good_baseline_fraction')
                   ! There will be no baseline for this segment, subtract the corresponding
                   ! binned signal
                   do i = baselines_short_start(k), baselines_short_stop(k)
                      ip = pixels(i, idet)
                      if (ip == dummy_pixel) cycle
                      if (nmap == 1) then
                         locmap(1,ip) = locmap(1,ip) &
                              - tod_stored(i,idet) * detweight
                         loccc(1,1,ip) = loccc(1,1,ip) - detweight
                      else
                         locmap(:,ip) = locmap(:,ip) &
                              - weights(:,i,idet) * tod_stored(i,idet) * detweight
                         call dsyr('U', nmap, detweight, weights(:,i,idet), &
                              -1, loccc(:,:,ip), nmap)
                      end if
                      if (do_hits) lochits(ip) = lochits(ip) - 1
                   end do
                end if
                cycle
             end if

             eigenvectorsT = transpose(eigenvectors)
             do i = 1,n
                eigenvectors(:,i) = eigenvectors(:,i) / eigenvalues(i)          
             end do
             nna_inv(:, :, k, idet) = matmul( eigenvectors, eigenvectorsT)
             
             ngood = ngood + 1
          end do
       end do
       
       deallocate(eigenvalues, eigenvectors, eigenvectorsT, workspace)
       !end if
    end do

    call sum_mpi(nbad)
    call sum_mpi(ngood)
    call sum_mpi(ntot)
    if (id == 0 .and. info > 0) then
       print *,'Succesfully inverted ', 100.0*ngood/dble(ntot), &
            ' % of the baseline matrices'
       print *,' Inverse rejected on ', 100.0*nbad/dble(ntot), &
            ' % of the baseline matrices'
       print *,'                     ', 100.0*(ntot-nbad-ngood) / dble(ntot), &
            ' % were empty'
    end if
    
    cputime_cga_init = cputime_cga_init + get_time(10)

    if (info > 4) write(*,idf) id,'Done'

  END SUBROUTINE initialize_a


  !-------------------------------------------------------------------------


  SUBROUTINE iterate_a(aa, yba, nna, wamap, cca, tod_stored)
    !
    ! Solution of the destriping equation by conjugate gradient algorithm
    ! This routine uses tables pixel,weights directly to speed up the computation.
    !
    real(dp), intent(out)   :: aa(0:basis_order, noba_short_max, nodetectors)
    real(dp), intent(in)    :: yba(0:basis_order, noba_short_max, nodetectors)
    real(dp), intent(in)    :: nna(0:basis_order, 0:basis_order, noba_short_max, nodetectors)
    real(dp), intent(inout) :: wamap(nmap,0:nopix_cross-1)
    real(dp), intent(in)    :: cca(nmap,nmap,0:nopix_cross-1)
    real(dp), intent(in) :: tod_stored(nosamples_proc,nodetectors)

    real(dp),allocatable   :: r(:,:,:), p(:,:,:), z(:,:,:), ap(:,:,:), proj(:)
    real(dp)               :: rz, rzinit, rzo, pap
    real(dp)               :: alpha, beta, pw, apn, bf, detweight
    integer :: i, j, k, n, m, ip, istep, idet, first, last, order, i0, m0
    integer :: order2, ichunk, noba, kstart, ipsd, itask
    real(dp), pointer :: basis_function(:,:)

    ! for openmp -RK
    integer :: ierr, imap
    real(dp), allocatable, target :: locmap_all_threads(:,:,:)
    real(dp), allocatable, target :: ap_all_threads(:,:,:,:)
    real(dp), pointer :: locmap_thread(:,:), ap_thread(:,:,:)
    real(dp) :: apnv(0:basis_order)    


    if (info > 4) write(*,idf) ID,'Begin iteration...'

    call reset_time(10)

    allocate( &
         r(0:basis_order, noba_short, nodetectors), &
         p(0:basis_order, noba_short, nodetectors), &
         ap(0:basis_order, noba_short, nodetectors), &
         z(0:basis_order, noba_short, nodetectors), &
         proj(0:basis_order), stat=ierr)
    if (ierr /= 0) stop 'iterate_a: no room for CG iteration'

    memory_cg = max(noba_short*nodetectors*32.,memory_cg)

    r = 0.0
    p = 0.0
    ap = 0.0
    z = 0.0

    aa = 0.0

    if (initialize_by_regression) then
       ! Set first guess aa from a template fit of the basis functions
       
       do idet = 1, nodetectors
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) cycle
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle       
             do k = kstart+1,kstart+noba
                m0 = baselines_short_start(k)
                proj = 0
                if (basis_order == 0) then
                   ! Shortcut
                   do m = baselines_short_start(k), baselines_short_stop(k)
                      if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                      ip = pixels(m, idet)
                      if (ip == dummy_pixel) cycle
                      
                      proj = proj + tod_stored(m,idet)
                   end do
                else
                   ! basis_order > 0
                   basis_function => basis_functions(k)%arr
                   do m = baselines_short_start(k), baselines_short_stop(k)
                      if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                      ip = pixels(m, idet)
                      if (ip == dummy_pixel) cycle
                      
                      proj = proj + basis_function(:,m-m0) * tod_stored(m,idet)
                   end do
                end if ! basis_order == 0 if
                
                proj = proj * detweight ! nna_inv includes detector weights
             
                proj = matmul( nna_inv(:,:,k,idet), proj )
                aa(:,k,idet) = proj
             end do  !loop over baselines
          end do
       end do
       
       ! Evaluate the first residual accordingly
       !r = yba - A * aa
       locmap = 0
       do idet = 1, nodetectors
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) cycle
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle       
             do k = kstart+1,kstart+noba
                m0 = baselines_short_start(k)
                basis_function => basis_functions(k)%arr
                do m = baselines_short_start(k), baselines_short_stop(k)
                   if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                   ip = pixels(m, idet)
                   if (ip == dummy_pixel) cycle
                   
                   pw = dot_product(basis_function(:,m-m0), &
                        aa(:,k,idet)) * detweight
                   if (nmap == 1) then
                      locmap(:,ip) = locmap(:,ip) + pw
                   else
                      locmap(:,ip) = locmap(:,ip) + weights(:,m,idet)*pw
                   end if
                end do
             end do
          end do
       end do

       wamap = 0.0
       call collect_map(wamap, nosubpix_cross)
       call ccmultiply(cca, wamap,nopix_cross)
       call scatter_map(wamap, nosubpix_cross)
       
       if (kfilter) then
          call cinvmul(ap,aa)
       else
          ap = 0
       endif
       
       do idet = 1, nodetectors
          do ichunk = first_chunk, last_chunk
             noba = noba_short_pp(ichunk)
             kstart = sum(noba_short_pp(first_chunk:ichunk-1))
             ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
             if ( ipsd < 0 ) cycle
             detweight = detectors(idet)%weights(ipsd)
             if (detweight == 0) cycle       
             do k = kstart+1,kstart+noba
                do order = 0,basis_order
                   apn = dot_product( nna(:,order,k,idet), aa(:,k,idet) )
                
                   i0 = baselines_short_start(k)
                   basis_function => basis_functions(k)%arr
                   do i = baselines_short_start(k), baselines_short_stop(k)
                      if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                      ip = pixels(i, idet)
                      if (ip == dummy_pixel) cycle
                      
                      if (nmap == 1) then
                         apn = apn - locmap(1,ip) &
                              * basis_function(order,i-i0) * detweight
                      else
                         apn = apn &
                              - dot_product(weights(:,i,idet), locmap(:,ip)) &
                              * basis_function(order,i-i0) * detweight
                      end if
                   end do
                   ap(order, k, idet) = ap(order, k, idet) + apn
                end do
             end do
          end do
       end do

       r = yba(:,1:noba_short,1:nodetectors) - ap(:,1:noba_short,1:nodetectors)
       ap = 0
    else
       ! Standard aa=0 first guess
       r = yba(:,1:noba_short,1:nodetectors)
    endif
    
    if (basis_order == 0) then
       call preconditioning_band(z(0,:,:), r(0,:,:))
    else
       do idet = 1,nodetectors
          do i = 1,noba_short
             z(0:basis_order,i,idet) = matmul( &
                  nna_inv(0:basis_order,0:basis_order,i,idet), &
                  r(0:basis_order,i,idet))
          end do
       end do
    end if
    
    p = z
    rz = sum(r*z)

    if (checknan) then
       if (isnan(rz)) print *,id,' : ERROR: rz is nan'          
          
       ybaloop : do idet = 1,nodetectors
          do i = 1,noba_short
             do order = 0, basis_order
                if (isnan(yba(order, i, idet))) then
                   print *,id,' : yba has NaN(s) : idet, base, order : ', &
                        idet, i, order
                   exit ybaloop
                end if
             end do
          end do
       end do ybaloop ! yba

       rloop : do idet = 1,nodetectors
          do i = 1,noba_short
             do order = 0, basis_order
                if (isnan(r(order,i,idet))) then
                   print *,id,' : r has NaN(s) : idet, base, order : ', &
                        idet, i, order
                   exit rloop
                end if
             end do
          end do
       end do rloop ! r

       zloop : do idet = 1,nodetectors
          do i = 1,noba_short
             do order = 0, basis_order
                if (isnan(z(order,i,idet))) then
                   print *,id,' : z has NaN(s) : idet, base, order : ', &
                        idet, i, order
                   print *,id,' : nna_inv = ', &
                        nna_inv(0:basis_order,0:basis_order,i,idet)
                   print *,id,' : nna = ',nna(0:basis_order,0:basis_order,i,idet)
                   print *,id,' : r = ',r(0:basis_order,i,idet)
                   exit zloop
                end if
             end do
          end do
       end do zloop ! z
          
       ploop : do idet = 1,nodetectors
          do i = 1,noba_short
             do order = 0, basis_order
                if (isnan(p(order,i,idet))) then
                   print *,id,' : p has NaN(s) : idet, base, order : ', &
                        idet, i, order
                   exit ploop
                end if
             end do
          end do
       end do ploop ! p          
    end if ! checknan

    call sum_mpi(rz)

    rzinit = rz
    if (rzinit.lt.1.e-30) then
       if (id == 0) then
          write(*,*)
          write(*,'(x,a,es25.15,a)') 'Anomalously low rzinit =', &
               rzinit,' NO CG ITERATION'
       end if
       return
    end if

    if (ID==0.and.info.ge.2) then
       write(*,*)
       write(*,*) 'CG iteration begins'
       write(*,'(x,a,es25.15)') 'rzinit =',rzinit
       write(*,*) 'istep, rz/rzinit'
    endif

    if (isnan(rz)) then
       if (id == 0) write(*,'(a)') 'WARNING: residual is NaN, no iterations'
       return
    endif


    allocate(locmap_all_threads(nmap,0:nsize_locmap,0:nthreads-1), stat=ierr)
    if (ierr /= 0) stop 'iterate_a: no room for locmap_all_threads'

    allocate(ap_all_threads(0:basis_order,noba_short,nodetectors,0:nthreads-1), &
         stat=ierr)
    if (ierr /= 0) stop 'iterate_a: no room for ap_all_threads'

    call reset_time(99) ! -RK

    istep = 0
    do
       if (istep >= iter_max) exit 
       istep = istep+1

       call reset_time(12)

       ! From baseline to map

       locmap = 0
       locmap_all_threads = 0

       if ( basis_order == 0 ) then
          if ( nmap == 1 ) then
             call baseline_to_map_order0_nopol()
          else if ( nmap == 3 ) then
             call baseline_to_map_order0_pol()
          else
             call baseline_to_map_order0_general()
          end if
       else
          if ( nmap == 1 ) then
             call baseline_to_map_nopol()
          else
             call baseline_to_map_general()
          end if
       end if

       locmap = sum( locmap_all_threads, dim=3 )

       ! Communicate maps between processes
       call wait_mpi
       cputime_cga_1 = cputime_cga_1 +get_time_and_reset(12)

       wamap = 0.0
       call collect_map(wamap,nosubpix_cross) ! locmap -> wamap
       call wait_mpi
       cputime_cga_mpi = cputime_cga_mpi + get_time_and_reset(12)
       call ccmultiply(cca,wamap,nopix_cross) ! apply cca, rejects masked pixels
       call wait_mpi
       cputime_cga_cc = cputime_cga_cc + get_time_and_reset(12)
       call scatter_map(wamap,nosubpix_cross) ! wamap -> locmap
       call wait_mpi
       cputime_cga_mpi = cputime_cga_mpi +get_time_and_reset(12)

       if (kfilter) then
          call cinvmul(ap,p)
       else
          ap = 0
       endif

       ! From map to baseline

       call reset_time(12)

       ap_all_threads = 0

       if ( basis_order == 0 ) then
          if ( nmap == 1 ) then
             call map_to_baseline_order0_nopol()
          else if ( nmap == 3 ) then
             call map_to_baseline_order0_pol()
          else
             call map_to_baseline_order0_general()
          end if
       else
          if ( nmap == 1 ) then
             call map_to_baseline_nopol()
          else
             call map_to_baseline_general()
          end if
       end if

       ap = ap + sum( ap_all_threads, dim=4 )

       Cputime_cga_2 = cputime_cga_2 + get_time(12)

       pap = sum(p*ap)
       call sum_mpi(pap)

       alpha = rz/pap
       r = r - alpha*ap

       aa(:, 1:noba_short, 1:nodetectors) = &
            aa(:, 1:noba_short, 1:nodetectors) + alpha*p

       if (basis_order == 0) then
          call preconditioning_band(z(0,:,:), r(0,:,:))
       else
          do idet = 1,nodetectors
             do i = 1,noba_short
                z(0:basis_order,i,idet) = matmul( &
                     nna_inv(0:basis_order,0:basis_order,i,idet), &
                     r(0:basis_order,i,idet))
             end do
          end do
       end if

       rzo = rz

       rz = sum(r*z)
       call sum_mpi(rz)

       if (ID==0.and.info.ge.2) write(*,'(i4,es16.6," (",f6.3,"s)")') &
            istep, rz/rzinit, get_time_and_reset(99)

       if (rz/rzinit < cglimit .and. istep > iter_min) exit
       if (rz == 0) exit

       beta = rz/rzo
       p = z + beta*p

    end do

    deallocate(locmap_all_threads, ap_all_threads)

    deallocate(r, p, ap, z, proj)

    noiter = istep
    if (id == 0 .and. info > 0) then
       write(*,*) 'Iteration done'
       write(*,'(i6,a)') noiter,' iteration steps'
    endif
    call wait_mpi

    cputime_cga = cputime_cga +get_time(10)

    if (info >= 5) write(*,idf) ID,'Done'

  contains
    
    ! Specific implementations of the CG iteration loop

    subroutine baseline_to_map_order0_nopol()
      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(locmap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,p,&
      !$OMP          baselines_short_start,baselines_short_stop,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel) &
      !$OMP   PRIVATE(itask,id_thread,locmap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,pw,m,ip)
      itask = -1
      id_thread = omp_get_thread_num()
      locmap_thread => locmap_all_threads(:,:,id_thread)
      loop_detector : do idet = 1, nodetectors
         loop_chunk : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk

            loop_baseline : do k = kstart+1,kstart+noba
               pw = p(0,k,idet) * detweight
               do m = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                  ip = pixels(m,idet)
                  if (ip == dummy_pixel) cycle
                  ip = ip + 1 ! Start locmap indexing from 1
                  
                  locmap_thread(1,ip) = locmap_thread(1,ip) + pw
               end do
            end do loop_baseline
         end do loop_chunk
      end do loop_detector
      !$OMP END PARALLEL      
      
    end subroutine baseline_to_map_order0_nopol


    subroutine baseline_to_map_order0_pol()

      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(locmap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,p,&
      !$OMP          baselines_short_start,baselines_short_stop,locmap,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel,weights) &
      !$OMP   PRIVATE(itask,id_thread,locmap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,pw,m,ip)

      itask = -1
      id_thread = omp_get_thread_num()
      locmap_thread => locmap_all_threads(:,:,id_thread)

      loop_detector : do idet = 1, nodetectors
         loop_chunk : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk

            loop_baseline : do k = kstart+1,kstart+noba
               pw = p(0,k,idet) * detweight
               do m = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                  ip = pixels(m,idet)
                  if (ip == dummy_pixel) cycle
                  ip = ip + 1 ! Start locmap indexing from 1
                  
                  locmap_thread(1,ip) = locmap_thread(1,ip) &
                       + weights(1,m,idet) * pw
                  locmap_thread(2,ip) = locmap_thread(2,ip) &
                       + weights(2,m,idet) * pw
                  locmap_thread(3,ip) = locmap_thread(3,ip) &
                       + weights(3,m,idet) * pw                         
               end do
            end do loop_baseline
         end do loop_chunk
      end do loop_detector

      !$OMP END PARALLEL      
      
    end subroutine baseline_to_map_order0_pol


    subroutine baseline_to_map_order0_general()
      
      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(locmap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,p,&
      !$OMP          baselines_short_start,baselines_short_stop,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel,weights) &
      !$OMP   PRIVATE(itask,id_thread,locmap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,pw,m,ip)
      itask = -1
      id_thread = omp_get_thread_num()
      locmap_thread => locmap_all_threads(:,:,id_thread)
      loop_detector : do idet = 1, nodetectors
         loop_chunk : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk

            loop_baseline : do k = kstart+1,kstart+noba
               pw = p(0, k, idet) * detweight
               do m = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                  ip = pixels(m, idet)
                  if (ip == dummy_pixel) cycle
                  ip = ip + 1 ! Start locmap indexing from 1
                  
                  locmap_thread(:,ip) = locmap_thread(:,ip) &
                       + weights(:,m,idet) * pw
               end do
            end do loop_baseline
         end do loop_chunk
      end do loop_detector
      !$OMP END PARALLEL      
      
    end subroutine baseline_to_map_order0_general


    subroutine baseline_to_map_nopol()
      
      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(locmap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,p,&
      !$OMP          baselines_short_start,baselines_short_stop,&
      !$OMP          basis_functions,basis_order,isubchunk,subchunk,pixels,&
      !$OMP          dummy_pixel) &
      !$OMP   PRIVATE(itask,id_thread,locmap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,m0,basis_function,m,ip)
      itask = -1
      id_thread = omp_get_thread_num()
      locmap_thread => locmap_all_threads(:,:,id_thread)
      loop_detector : do idet = 1, nodetectors
         loop_chunk : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk

            loop_baseline : do k = kstart+1,kstart+noba
               m0 = baselines_short_start(k)
               basis_function => basis_functions(k)%arr
               do m = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                  ip = pixels(m, idet)
                  if (ip == dummy_pixel) cycle
                  ip = ip + 1 ! Start locmap indexing from 1
                  
                  locmap_thread(1,ip) = locmap_thread(1,ip) + detweight &
                       * dot_product(basis_function(0:basis_order,m-m0), &
                       p(0:basis_order,k,idet))
               end do
            end do loop_baseline
         end do loop_chunk
      end do loop_detector
      !$OMP END PARALLEL      
      
    end subroutine baseline_to_map_nopol

    
    subroutine baseline_to_map_general()
      
      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(locmap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,p,&
      !$OMP          baselines_short_start,baselines_short_stop,&
      !$OMP          basis_functions,basis_order,isubchunk,subchunk,pixels,&
      !$OMP          dummy_pixel,nmap,weights) &
      !$OMP   PRIVATE(itask,id_thread,locmap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,m0,basis_function,m,ip,pw)
      itask = -1
      id_thread = omp_get_thread_num()
      locmap_thread => locmap_all_threads(:,:,id_thread)
      loop_detector : do idet = 1, nodetectors
         loop_chunk : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk

            loop_baseline : do k = kstart+1,kstart+noba
               m0 = baselines_short_start(k)
               basis_function => basis_functions(k)%arr
               do m = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(m) /= isubchunk ) cycle
                  ip = pixels(m, idet)
                  if (ip == dummy_pixel) cycle
                  ip = ip + 1 ! Start locmap indexing from 1
                  
                  pw = dot_product(basis_function(0:basis_order,m-m0), &
                       p(0:basis_order,k,idet)) * detweight
                  locmap_thread(1:nmap,ip) = locmap_thread(1:nmap,ip) &
                       + weights(1:nmap,m,idet)*pw
               end do
            end do loop_baseline
         end do loop_chunk
      end do loop_detector
      !$OMP END PARALLEL      
      
    end subroutine baseline_to_map_general


    subroutine map_to_baseline_order0_nopol()

      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(ap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,&
      !$OMP          nna,p,baselines_short_start,baselines_short_stop,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel,locmap) &
      !$OMP   PRIVATE(itask,id_thread,ap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,apn,i,ip)
      itask = -1
      id_thread = omp_get_thread_num()
      ap_thread => ap_all_threads(:,:,:,id_thread)

      loop_detector_ap : do idet = 1, nodetectors
         loop_chunk_ap : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))            
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk_ap
            
            loop_baseline_ap : do k = kstart+1,kstart+noba
               ! / detweight is canceled later
               apn = nna(0,0,k,idet) * p(0,k,idet) / detweight
               
               do i = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                  ip = pixels(i, idet)
                  if (ip == dummy_pixel) cycle
                  
                  apn = apn - locmap(1,ip)
               end do
                  
               ap_thread(1, k, idet) = ap_thread(1, k, idet) + apn * detweight
            end do loop_baseline_ap
         end do loop_chunk_ap
      end do loop_detector_ap
      !$OMP END PARALLEL

    end subroutine map_to_baseline_order0_nopol


    subroutine map_to_baseline_order0_pol()

      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(ap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,&
      !$OMP          nna,p,baselines_short_start,baselines_short_stop,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel,weights,locmap) &
      !$OMP   PRIVATE(itask,id_thread,ap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,apn,i,ip)

      itask = -1
      id_thread = omp_get_thread_num()
      ap_thread => ap_all_threads(:,:,:,id_thread)

      loop_detector_ap : do idet = 1, nodetectors
         loop_chunk_ap : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle loop_chunk_ap
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk_ap
            
            loop_baseline_ap : do k = kstart+1,kstart+noba
               ! / detweight is canceled later
               apn = nna(0,0,k,idet) * p(0,k,idet) / detweight
               
               do i = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                  ip = pixels(i,idet)
                  if (ip == dummy_pixel) cycle
                  
                  apn = apn - weights(1,i,idet) * locmap(1,ip) &
                       - weights(2,i,idet) * locmap(2,ip) &
                       - weights(3,i,idet) * locmap(3,ip)
               end do
               
               ap_thread(1,k,idet) = ap_thread(1,k,idet) + apn * detweight
            end do loop_baseline_ap
         end do loop_chunk_ap
      end do loop_detector_ap

      !$OMP END PARALLEL

    end subroutine map_to_baseline_order0_pol


    subroutine map_to_baseline_order0_general()

      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(ap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,&
      !$OMP          nna,p,baselines_short_start,baselines_short_stop,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel,weights,locmap,nmap) &
      !$OMP   PRIVATE(itask,id_thread,ap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,apn,i,ip)
      itask = -1
      id_thread = omp_get_thread_num()
      ap_thread => ap_all_threads(:,:,:,id_thread)

      loop_detector_ap : do idet = 1, nodetectors
         loop_chunk_ap : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))            
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk_ap
            
            loop_baseline_ap : do k = kstart+1,kstart+noba
               ! / detweight is canceled later
               apn = nna(0,0,k,idet) * p(0,k,idet) / detweight
               
               do i = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                  ip = pixels(i,idet)
                  if (ip == dummy_pixel) cycle
                  
                  apn = apn &
                       - dot_product(weights(1:nmap,i,idet), locmap(1:nmap,ip))
               end do

               ap_thread(1, k, idet) = ap_thread(1, k, idet) + apn * detweight
            end do loop_baseline_ap
         end do loop_chunk_ap
      end do loop_detector_ap
      !$OMP END PARALLEL

    end subroutine map_to_baseline_order0_general


    subroutine map_to_baseline_nopol()
      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(ap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,&
      !$OMP          basis_order,basis_functions,&
      !$OMP          nna,p,baselines_short_start,baselines_short_stop,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel,weights,locmap) &
      !$OMP   PRIVATE(itask,id_thread,ap_thread,idet,ichunk,noba,kstart,ipsd,&
      !$OMP          detweight,k,apnv,i0,basis_function,i,ip)
      itask = -1
      id_thread = omp_get_thread_num()
      ap_thread => ap_all_threads(:,:,:,id_thread)

      loop_detector_ap : do idet = 1, nodetectors
         loop_chunk_ap : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))            
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk_ap
            
            loop_baseline_ap : do k = kstart+1,kstart+noba
               ! / detweight is canceled later
               apnv = matmul( nna(0:basis_order,0:basis_order,k,idet), &
                    p(0:basis_order,k,idet) ) / detweight
                  
               i0 = baselines_short_start(k)
               basis_function => basis_functions(k)%arr
               do i = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                  ip = pixels(i,idet)
                  if (ip == dummy_pixel) cycle
                  
                  apnv = apnv - locmap(1,ip) * basis_function(0:basis_order,i-i0)
               end do
                  
               ap_thread(:,k,idet) = ap_thread(:,k,idet) + apnv * detweight
            end do loop_baseline_ap
         end do loop_chunk_ap
      end do loop_detector_ap
      !$OMP END PARALLEL

    end subroutine map_to_baseline_nopol


    subroutine map_to_baseline_general()
      
      !$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
      !$OMP   SHARED(ap_all_threads,nodetectors,first_chunk,last_chunk,&
      !$OMP          noba_short_pp,baselines_short_time,detectors,nthreads,&
      !$OMP          basis_order,basis_functions,&
      !$OMP          nna,p,baselines_short_start,baselines_short_stop,&
      !$OMP          isubchunk,subchunk,pixels,dummy_pixel,weights,locmap,nmap) &
      !$OMP   PRIVATE(itask,id_thread,ap_thread,idet,ichunk,noba,kstart,&
      !$OMP          ipsd,detweight,k,apnv,i0,basis_function,i,ip)
      itask = -1
      id_thread = omp_get_thread_num()
      ap_thread => ap_all_threads(:,:,:,id_thread)

      loop_detector_ap : do idet = 1, nodetectors
         loop_chunk_ap : do ichunk = first_chunk, last_chunk
            noba = noba_short_pp(ichunk)
            kstart = sum(noba_short_pp(first_chunk:ichunk-1))            
            ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
            if ( ipsd < 0 ) cycle
            detweight = detectors(idet)%weights(ipsd)
            if (detweight == 0) cycle
            
            itask = itask + 1
            if ( modulo( itask, omp_get_num_threads() ) /= id_thread ) &
                 cycle loop_chunk_ap

            loop_baseline_ap : do k = kstart+1,kstart+noba
               ! / detweight is canceled later
               apnv = matmul(nna(0:basis_order,0:basis_order,k,idet), &
                    p(0:basis_order,k,idet)) / detweight
                  
               i0 = baselines_short_start(k)
               basis_function => basis_functions(k)%arr
               do i = baselines_short_start(k), baselines_short_stop(k)
                  if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                  ip = pixels(i,idet)
                  if (ip == dummy_pixel) cycle
                  
                  apnv = apnv - dot_product(weights(1:nmap,i,idet), &
                       locmap(1:nmap,ip)) * basis_function(0:basis_order,i-i0)
               end do
                  
               ap_thread(:,k,idet) = ap_thread(:,k,idet) + apnv * detweight
            end do loop_baseline_ap
         end do loop_chunk_ap
      end do loop_detector_ap
      !$OMP END PARALLEL

    end subroutine map_to_baseline_general


  END SUBROUTINE iterate_a


  !---------------------------------------------------------------------------


  SUBROUTINE subtract_baselines_a(map, aa)
    !
    ! Subtract baselines and compute the final map
    real(dp), intent(inout) :: map(nmap,0:nopix_map-1)
    real(dp), intent(in) :: aa(0:basis_order, noba_short_max, nodetectors)
    integer :: i, k, ip, idet, i0, order, ichunk, noba, kstart, ipsd
    real(dp) :: aw, detweight
    real(dp), pointer :: basis_function(:,:)

    if (info == 3 .and. id == 0) write(*,*) 'Subtracting baselines...'
    if (info > 4) write(*,idf) id,'Subtracting baselines...'

    locmap = 0

    do idet = 1,nodetectors
       do ichunk = first_chunk, last_chunk
          noba = noba_short_pp(ichunk)
          kstart = sum(noba_short_pp(first_chunk:ichunk-1))
          ipsd = psd_index_det( idet, baselines_short_time(kstart+1) )
          if ( ipsd < 0 ) cycle
          detweight = detectors(idet)%weights(ipsd)
          if (detweight == 0) cycle       
          do k = kstart+1,kstart+noba
             i0 = baselines_short_start(k)
             basis_function => basis_functions(k)%arr
             do i = baselines_short_start(k), baselines_short_stop(k)
                if ( isubchunk /= 0 .and. subchunk(i) /= isubchunk ) cycle
                
                ip = pixels(i,idet)
                if (ip == dummy_pixel) cycle
                
                aw = dot_product(basis_function(:,i-i0), aa(:,k,idet)) * detweight
                
                if (nmap == 1) then
                   locmap(1,ip) = locmap(1,ip) - aw
                else
                   locmap(:,ip) = locmap(:,ip) - aw*weights(:,i,idet)
                end if

             end do
          end do
       end do
    end do

    call collect_map(map, nosubpix_map)

    if (info > 4) write(*,idf) id,'Done'

  END SUBROUTINE subtract_baselines_a


  !---------------------------------------------------------------------------


  SUBROUTINE clean_tod(tod,aa)
    !
    ! Subtract baselines from the TOD

    real(dp),intent(inout) :: tod(nosamples_proc,nodetectors)
    real(dp),intent(in) :: aa(0:basis_order,noba_short_max,nodetectors)
    integer :: i, k, idet, order, i0
    real(dp), pointer :: basis_function(:,:)

    if (.not. write_tod .and. ndetset == 0 .and. nsurvey == 0) return

    if (.not. kfirst) return

    if ( tod_is_clean ) return

    if (info == 3 .and. id == 0) write(*,*) 'Subtracting baselines from TOD...'
    if (info > 4) write(*,idf) ID,'Subtracting baselines from TOD...'

    call reset_time(10)

    do idet = 1,nodetectors
       do k = 1,noba_short
          i0 = baselines_short_start(k)
          basis_function => basis_functions(k)%arr
          do i = baselines_short_start(k), baselines_short_stop(k)
             if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
             do order = 0, basis_order
                tod(i,idet) = tod(i,idet) &
                     - basis_function(order,i-i0) * aa(order,k,idet)
             end do
          end do
       end do
    end do

    tod_is_clean = .true.

    cputime_clean_tod = cputime_clean_tod  + get_time(10)

    if (info > 4) write(*,idf) ID,'Done'

  END SUBROUTINE clean_tod


  !--------------------------------------------------------------------------


  SUBROUTINE unclean_tod(tod, aa)
    !
    ! Subtract baselines from the TOD

    real(dp),intent(inout) :: tod(nosamples_proc,nodetectors)
    real(dp),intent(in) :: aa(0:basis_order,noba_short_max,nodetectors)
    integer :: i, k, idet, order, i0
    real(dp), pointer :: basis_function(:,:)

    if (.not. write_tod .and. ndetset == 0 .and. nsurvey == 0) return

    if ( .not. tod_is_clean ) return

    if (info == 3 .and. id == 0) write(*,*) 'Adding baselines to TOD...'
    if (info > 4) write(*,idf) ID,'Adding baselines to TOD...'

    call reset_time(10)

    do idet = 1,nodetectors
       do k = 1,noba_short
          i0 = baselines_short_start(k)
          basis_function => basis_functions(k)%arr
          do i = baselines_short_start(k), baselines_short_stop(k)
             if ( isubchunk /= 0 .and. subchunkpp(i) /= isubchunk ) cycle
             do order = 0, basis_order
                tod(i,idet) = tod(i,idet) &
                     + basis_function(order,i-i0) * aa(order,k,idet)
             end do
          end do
       end do
    end do

    tod_is_clean = .false.

    cputime_clean_tod = cputime_clean_tod + get_time(10)

    if (info > 4) write(*,idf) ID,'Done'

  END SUBROUTINE unclean_tod


  !--------------------------------------------------------------------------

END MODULE madam_routines
