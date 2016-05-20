PROGRAM hole_filler

   use fits_routines
   use matrix
   use inout_holefiller

   implicit none
   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: isp = selected_int_kind(9)
   integer, parameter :: idp = selected_int_kind(18)

   integer  :: iargc, nopar, i, k, ipix, nside, ibuff
   integer  :: nbuff, nosteps, degrade, level
   real(dp) :: cdet, rcond
   real(dp),allocatable :: cc(:), cdegraded(:), bmap(:)

   integer  :: nopix(12)   = 0, &
               nopix_good  = 0, &
               nopix_bad   = 0

   integer  :: nside_out = 0,  &
               nside_min = 16
   real(dp) :: rcond_min = 0.0

   character(len=200) :: file_map_in     = '', &
                         file_map_out    = '', &
                         file_matrix_in  = '', &
                         file_matrix_out = '', &
                         file_nside      = ''

   real(dp),allocatable :: cc_in(:,:),  cc_out(:,:)
   real(dp),allocatable :: map_in(:,:), map_out(:,:)
   integer, allocatable :: nside_map(:)

   call read_parameters
   call open_infiles(file_map_in,file_matrix_in)
   call create_outfiles(file_map_out,file_matrix_out,file_nside,  &
                        nside_out,rcond_min)

   nbuff = nside_out*nside_out
   nosteps = 12*nside_out*nside_out/nbuff

   allocate(map_in(nostokes,nbuff))
   allocate(map_out(nostokes,nbuff))
   allocate(cc_in(ncc,nbuff))
   allocate(cc_out(ncc,nbuff))
   allocate(nside_map(nbuff))
   allocate(cc(ncc),cdegraded(ncc),bmap(nostokes))


   do ibuff = 1,nosteps

      write(*,*) 'ibuff =',ibuff

      call read_infiles(map_in,cc_in,nbuff,nside_out)

      do i = 1,nbuff

         if (nostokes==3) then
            cc = cc_in(:,i)
            call invert3_eig(cc,cdet,rcond)
         else
            cc = cc_in(1,i)
            if (cc(1).lt.1.e-20) then
               cc(1) = 0.0
               rcond = 0.0
            else
               cc(1) = 1.d0/cc(1)
               rcond = 1.d0
            endif
         endif

         if (rcond.gt.rcond_min) then

            nopix_good  = nopix_good+1
            cc_out(:,i) = cc_in(:,i)
            bmap = map_in(:,i)
            nside_map(i) = nside_out
         else
! Fill the pixel from lower resolution
            level   = 0
            degrade = 1
            nside   = nside_out
            do
               level   = level+1
               degrade = degrade*4
               nside   = nside/2
               write(*,*) 'i,nside =',i,nside

               if (nside.lt.nside_min) then
                   nopix_bad = nopix_bad+1
                   cc   = 0.0
                   bmap = 0.0
                   nside_map(i) = 0
                   cc_out(:,i)  = 0.0
                   exit
               endif

               ipix = (i-1)/degrade
               ipix = ipix*degrade

               do k = 1,nostokes
                  bmap(k) = sum(map_in(k,ipix+1:ipix+degrade))
               enddo
               do k = 1,ncc
                  cdegraded(k) = sum(cc_in(k,ipix+1:ipix+degrade))
               enddo

              if (nostokes==3) then

                  cc = cdegraded
                  call invert3_eig(cc,cdet,rcond)
               else

                  if (cdegraded(1).gt.1.e-20) then
                     cc = 1.d0/cdegraded
                     rcond = 1.d0
                  else
                     cc = 0.0
                     rcond = 0.0
                  endif
               endif

               if (rcond.gt.rcond_min) then
                  nopix(level) = nopix(level)+1
                  cc_out(:,i)  = cdegraded
                  nside_map(i) = nside
                  exit
               endif

            enddo
         endif

         if (nostokes==3) then
            map_out(1,i) = cc(1)*bmap(1)+cc(2)*bmap(2)+cc(3)*bmap(3)
            map_out(2,i) = cc(2)*bmap(1)+cc(4)*bmap(2)+cc(5)*bmap(3)
            map_out(3,i) = cc(3)*bmap(1)+cc(5)*bmap(2)+cc(6)*bmap(3)
         else
            map_out(1,i) = cc(1)*bmap(1)
         endif
      enddo

      call write_outfiles(map_out,cc_out,nside_map,nbuff)

   enddo

   call close_files(sum(nopix))

   write(*,'(i9,x,a)') nopix_good,'pixels accepted'
   nside = nside_out
   do k = 1,12
      nside = nside/2
      if (nopix(k).gt.0)   &
         write(*,'(i9,x,a,i5)') nopix(k),'pixels filled from nside',nside
   enddo
   write(*,'(i9,x,a)') nopix_bad,'pixels rejected'
   write(*,'(i9,x,a)') 12*nside_out*nside_out,'pixels total'
   write(*,*)

   write(*,*) 'Output written in files'
   write(*,*) trim(file_map_out)
   if (do_matrix) write(*,*) trim(file_matrix_out)
   if (do_nside)  write(*,*) trim(file_nside)
   write(*,*)

!END

CONTAINS

!------------------------------------------------------------------------------


   SUBROUTINE read_parameters

      integer            :: i, k, n, iend, nolines, ipar, nopar, iargc
      character(len=200) :: lines(100), key, value
      character(len=80)  :: par

      nopar = iargc()

      k = 1
      do ipar = 1,nopar
         call getarg(ipar,par)

         n = scan(par,'=')

         if (n==0) then
            open(unit=15,file=par,status='old')
            do
              read(15,'(a200)',iostat=iend) lines(k)
              if (iend.lt.0) exit

              k = k+1
           enddo
           close(15)
        else
           lines(k) = par
           k = k+1
        endif
      enddo

      nolines = k-1

      do i = 1,nolines
         k = scan(lines(i),'=')
         if (k==0.or.lines(i)(1:1)=='#') cycle

         key = trim(adjustl(lines(i)(1:k-1)))
         value = trim(adjustl(lines(i)(k+1:200)))

         if (key=='file_matrix_in')  file_matrix_in = trim(value)
         if (key=='file_matrix_out') file_matrix_out = trim(value)
         if (key=='file_map_in')     file_map_in = trim(value)
         if (key=='file_map_out')    file_map_out = trim(value)
         if (key=='file_nside')      file_nside = trim(value)

         if (key=='nside_out')       read(value,*) nside_out
         if (key=='nside_min')       read(value,*) nside_min
         if (key=='rcond_min')       read(value,*) rcond_min
      enddo

      write(*,*)
      write(*,*) 'Input files:'
      write(*,*) 'file_map_in     = ',trim(file_map_in)
      write(*,*) 'file_matrix_in  = ',trim(file_matrix_in)
      write(*,*)
      write(*,*) 'Output files:'
      write(*,*) 'file_map_out    = ',trim(file_map_out)
      write(*,*) 'file_matrix_out = ',trim(file_matrix_out)
      write(*,*) 'file_nside      = ',trim(file_nside)
      write(*,*)
      write(*,'(x,a,i5)')     'nside_out =',nside_out
      write(*,'(x,a,i5)')     'nside_min =',nside_min
      write(*,'(x,a,es14.5)') 'rcond_min =',rcond_min
      write(*,*)

   END SUBROUTINE

!------------------------------------------------------------------------------

END PROGRAM