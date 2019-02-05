MODULE fits_routines

! Routines for handling of fitsio files
! Partly copied from LevelS

   implicit none
   private

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: isp = selected_int_kind(9)
   integer, parameter :: idp = selected_int_kind(18)

   character, parameter :: datatype_short  = 'I'
   character, parameter :: datatype_int    = 'J'
   character, parameter :: datatype_long   = 'K'
   character, parameter :: datatype_real   = 'E'
   character, parameter :: datatype_double = 'D'
   character, parameter :: datatype_char   = 'A'

! FITS file access types
   integer, parameter, public :: FITS_READONLY  = 0
   integer, parameter, public :: FITS_READWRITE = 1

! FITS HDU types
   integer, parameter, public :: FITS_IMAGE  = 0
   integer, parameter, public :: FITS_ASCTAB = 1
   integer, parameter, public :: FITS_BINTAB = 2

! constant for yet undefined values
   integer, parameter, public :: FITS_INVALID = -123456

!! data type containing information about the currently connected HDU
   type handle_hdu
      integer                     :: unit    = -1
      integer(idp)                :: nrows   = -1
      integer                     :: ncols   = -1
      integer                     :: hdutype = -1
      type(handle_column),pointer :: columns(:) =>NULL()
   end type

!! data type containing information about a single column in a FITS table
   type handle_column
      character(len=80) :: name = ''
      character(len=80) :: unit = ''
      integer           :: repcount = 1
      character         :: datatype = ''
      integer(idp)      :: length = 0
   end type

   public fits_open,              &
          fits_close,             &
          fits_create,            &
          fits_open_for_write,    &
          fits_goto_hdu

   public fits_init_bintab,       &
          fits_set_name,          &
          fits_set_unit

   public fits_find_column,       &
          fits_read_bincolumn,    &
          fits_write_bincolumn

   public fits_put_key,           &
          fits_get_key,           &
          fits_put_comment,       &
          addi

   public handle_hdu, handle_column

   interface fits_set_name
      module procedure set_name_s,    &
                       set_name_vec
   end interface

   interface fits_set_unit
      module procedure set_unit_s,    &
                       set_unit_vec
   end interface

   interface fits_read_bincolumn
      module procedure read_column_dbin,  &
                       read_column_sbin,  &
                       read_column_ibin
   end interface

   interface fits_write_bincolumn
      module procedure write_column_dbin,  &
                       write_column_sbin,  &
                       write_column_ibin
   end interface

   interface fits_get_key
      module procedure get_key_int,     &
                       get_key_logical, &
                       get_key_string,  &
                       get_key_real,    &
                       get_key_double
   end interface

   interface fits_put_key
      module procedure put_key_int,     &
                       put_key_logical, &
                       put_key_string,  &
                       put_key_real,    &
                       put_key_double
   end interface

CONTAINS

!------------------------------------------------------------------------------


   SUBROUTINE fits_open(handle,file)
! Open an existing file for read

      type(handle_hdu),intent(out) :: handle
      character(len=*),intent(in)  :: file
      integer                      :: stat, block
      logical                      :: there

      inquire(file=trim(file),exist=there)

      if (.not.there) then
         write(*,*) 'ERROR in fits_open: File not found'
         write(*,*) trim(file)
         stop
      endif

      stat = 0
      call FTGIOU(handle%unit,stat)
      call FTOPEN(handle%unit,trim(file),fits_readonly,block,stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE fits_open_for_write(handle,file)
! Open an existing file for write

      type(handle_hdu),intent(out) :: handle
      character(len=*),intent(in)  :: file
      integer                      :: stat, block
      logical                      :: there

      inquire(file=trim(file),exist=there)

      if (.not.there) then
         write(*,*) 'ERROR in fits_open: File not found'
         write(*,*) trim(file)
         stop
      endif

      stat = 0
      call FTGIOU(handle%unit,stat)
      call FTOPEN(handle%unit,trim(file),fits_readwrite,block,stat)
      call check_errors(stat)

   END SUBROUTINE



!------------------------------------------------------------------------------


   SUBROUTINE fits_goto_hdu(handle,nhdu)

      type(handle_hdu),intent(inout) :: handle
      integer,         intent(in)    :: nhdu
      integer                        :: stat, i, repcount, tnull
      real(dp)                       :: tscal,tzero
      character(len=80)              :: tdisp, ttype, tunit, comment
      character                      :: datatype

      stat = 0
      call FTMAHD(handle%unit,nhdu,handle%hdutype,stat)
      call check_errors(stat)

      call FTGNCL(handle%unit,handle%ncols,stat)
!      call FTGNRW(handle%unit,handle%nrows,stat)
      call FTGKYK(handle%unit,'NAXIS2',handle%nrows,comment,stat)
! patch 29.5.08 to handle file sizes exceeding 2^31 rows

      call check_errors(stat)

      if (associated(handle%columns)) deallocate(handle%columns)

      allocate(handle%columns(handle%ncols))

      stat = 0
      do i = 1,handle%ncols

         call FTGBCL(handle%unit,i,ttype,tunit,datatype,repcount,   &
                     tscal,tzero,tnull,tdisp,stat)

         handle%columns(i)%datatype = uppercase(datatype)
         handle%columns(i)%repcount = repcount
         handle%columns(i)%length   = repcount*handle%nrows
         handle%columns(i)%name     = trim(ttype)
         handle%columns(i)%unit     = trim(tunit)

         call check_errors(stat)
      enddo

  CONTAINS

      FUNCTION uppercase(c) RESULT(cu)
         character :: c, cu
         integer   :: ic

         ic = ichar(c)
         if (ic.ge.ichar('a').and.ic.le.ichar('z')) &
            ic = ic +ichar('A')-ichar('a')
         cu = char(ic)

      END FUNCTION

   END SUBROUTINE


!------------------------------------------------------------------------------


   FUNCTION fits_find_column(handle,colname) RESULT(colnum)
! Find column with given name
! Return -1 if file or column not found.

     integer          :: colnum
     type(handle_hdu) :: handle
     character(len=*) :: colname
     integer          :: i

     colnum = -1
     do i = 1,handle%ncols
        if (handle%columns(i)%name==colname) then
           colnum = i
           exit
        endif
     enddo

   END FUNCTION


!------------------------------------------------------------------------------


   SUBROUTINE check_column_read(handle,icol,frow,felem,n,offset)

      type(handle_hdu),intent(in)   :: handle
      integer(idp),    intent(out)  :: frow
      integer,         intent(out)  :: felem
      integer,         intent(in)   :: n, icol
      integer(idp),    intent(in)   :: offset
      integer(idp)                  :: idrep, ilast

      if (icol.gt.handle%ncols) then
         write(*,*) 'ERROR: Index exceeds number of columns.'
         write(*,'(a,2i5)') ': icol, ncols = ',icol,handle%ncols
         stop
      endif

      idrep = handle%columns(icol)%repcount
      frow  = offset/idrep+1
      felem = mod(offset,idrep)+1

      ilast = offset+n

      if (ilast.gt.handle%columns(icol)%length) then
         write(*,*) 'ERROR: Reading past end of file'
         write(*,'(a,2i9)') 'offset, n =',offset,n
         write(*,'(a,i9)')  'length    =',handle%columns(icol)%length
         stop
      endif

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE read_column_dbin(handle,icol,dbuffer,n,offset)
! Read n samples from column icol, skipping the first offset samples

      type(handle_hdu),intent(in)          :: handle
      integer,         intent(in)          :: n, icol
      real(dp),        intent(out)         :: dbuffer(n)
      integer(idp),    intent(in),optional :: offset
      integer                              :: stat, felem
      integer(idp)                         :: ioff, frow
      logical                              :: anynul

      if (present(offset)) then
         ioff = offset
      else
         ioff = 0
      endif

      call check_column_read(handle,icol,frow,felem,n,ioff)

      stat = 0
      call FTGCVD(handle%unit,icol,int(frow),felem,n,0._dp,dbuffer,anynul,stat)

      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE read_column_sbin(handle,icol,sbuffer,n,offset)
! Read n samples from column icol, skipping the first offset samples

      type(handle_hdu),intent(in)          :: handle
      integer,         intent(in)          :: n, icol
      real(sp),        intent(out)         :: sbuffer(n)
      integer(idp),    intent(in),optional :: offset
      integer                              :: stat, felem
      integer(idp)                         :: ioff, frow
      logical                              :: anynul

      if (present(offset)) then
         ioff = offset
      else
         ioff = 0
      endif

      call check_column_read(handle,icol,frow,felem,n,ioff)

      stat = 0
      call FTGCVE(handle%unit,icol,int(frow),felem,n,0.,sbuffer,anynul,stat)

      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE read_column_ibin(handle,icol,ibuffer,n,offset)
! Read n samples from column icol, skipping the first offset samples

      type(handle_hdu),intent(in)          :: handle
      integer,         intent(in)          :: n, icol
      integer,         intent(out)         :: ibuffer(n)
      integer(idp),    intent(in),optional :: offset
      integer                              :: stat, felem
      integer(idp)                         :: ioff, frow
      logical                              :: anynul

      if (present(offset)) then
         ioff = offset
      else
         ioff = 0
      endif

      call check_column_read(handle,icol,frow,felem,n,ioff)

      stat = 0
      call FTGCVJ(handle%unit,icol,int(frow),felem,n,0,ibuffer,anynul,stat)

      call check_errors(stat)

   END SUBROUTINE


!-------------------------------------------------------------------------------


   SUBROUTINE fits_close(handle)

     type(handle_hdu),intent(inout) :: handle
     integer                        :: stat

     stat = 0
     call FTCLOS(handle%unit,stat)
     call FTFIOU(handle%unit,stat)
     call check_errors(stat)

     handle%unit = -1
     if (associated(handle%columns)) deallocate(handle%columns)

     handle%columns =>NULL()
     handle%hdutype = -1
     handle%nrows = -1
     handle%ncols = -1

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE fits_create(handle,file)

      type(handle_hdu), intent(out) :: handle
      character(len=*), intent(in)  :: file
      integer :: stat, zeroff(2)

      stat = 0
      call ftgiou(handle%unit,stat)
      call ftinit(handle%unit,trim(file),1,stat)
      zeroff = 0
      call ftphps(handle%unit,8,0,zeroff,stat)
      call ftpdat(handle%unit,stat)
      call ftplsw(handle%unit,stat)

      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE fits_init_bintab(handle,datatype,repcount)

      type(handle_hdu), intent(inout)       :: handle
      character(len=*), intent(in)          :: datatype
      integer,          intent(in),optional :: repcount
      character(len=80),allocatable         :: tform(:),tunit(:),ttype(:)
      integer                               :: stat, i, n, rep

      n = len(datatype)

      handle%ncols = n
      handle%nrows = 0

      allocate(tform(n),tunit(n),ttype(n))

      if (associated(handle%columns)) deallocate(handle%columns)
      allocate (handle%columns(n))

      ttype = ''
      tunit = ''

      rep = 1
      if (present(repcount)) rep=repcount

      do i = 1, n
         write (tform(i),*) rep
         tform(i) = trim(adjustl(tform(i)))//datatype(i:i)

         handle%columns(i)%repcount = rep
      enddo

      do i = 1,n
         handle%columns(i)%datatype = datatype(i:i)
      enddo

      stat = 0
      call FTIBIN(handle%unit,0,n,ttype,tform,tunit,'xtension',0,stat)
      call check_errors(stat)

      handle%hdutype = FITS_BINTAB

      deallocate(tform,tunit,ttype)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE set_name_vec(handle,name)

      type(handle_hdu),intent(inout) :: handle
      character(len=*),intent(in)    :: name(:)
      integer                        :: stat, k, n
      character(len=8)               :: key

      n = size(name)

      if (n.ne.handle%ncols) then
         write(*,*) 'ERROR in fits_set_name:'
         write(*,*) 'Number of entries doe not match number of columns.'
         stop
      endif

      do k = 1,handle%ncols

         key = addi('TTYPE',k)

         call FTUKYS(handle%unit,key,name(k),'&',stat)

         handle%columns(k)%name = name(k)
      enddo

   END SUBROUTINE


   SUBROUTINE set_name_s(handle,name,icol)

      type(handle_hdu),intent(inout)       :: handle
      character(len=*),intent(in)          :: name
      integer,         intent(in),optional :: icol
      integer                              :: stat, k
      character(len=8)                     :: key

      if (present(icol)) then
         if (icol.gt.handle%ncols) then
            write(*,*) 'ERROR in fits_set_name:'
            write(*,*) 'Column index icol =',icol
            write(*,*) 'exceeds number of columns =',handle%ncols
            stop
         endif
      endif

      do k = 1,handle%ncols

         if (present(icol)) then
             if (k.ne.icol) cycle
         endif

         key = addi('TTYPE',k)

         call FTUKYS(handle%unit,key,name,'&',stat)

         handle%columns(k)%name = name
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE set_unit_vec(handle,name)

      type(handle_hdu),intent(inout) :: handle
      character(len=*),intent(in)    :: name(:)
      integer                        :: stat, k, n
      character(len=8)               :: key

      n = size(name)

      if (n.ne.handle%ncols) then
         write(*,*) 'ERROR in fits_set_name:'
         write(*,*) 'Number of entries doe not match number of columns.'
         stop
      endif

      do k = 1,handle%ncols

         key = addi('TUNIT',k)

         call FTUKYS(handle%unit,key,name(k),'&',stat)

         handle%columns(k)%name = name(k)
      enddo

   END SUBROUTINE


   SUBROUTINE set_unit_s(handle,name,icol)

      type(handle_hdu),intent(inout)       :: handle
      character(len=*),intent(in)          :: name
      integer,         intent(in),optional :: icol
      integer                              :: stat, k
      character(len=8)                     :: key

      if (present(icol)) then
         if (icol.gt.handle%ncols) then
            write(*,*) 'ERROR in fits_set_unit:'
            write(*,*) 'Column index icol =',icol
            write(*,*) 'exceeds number of columns =',handle%ncols
            stop
         endif
      endif

      do k = 1,handle%ncols

         if (present(icol)) then
             if (k.ne.icol) cycle
         endif

         key = addi('TUNIT',k)

         call FTUKYS(handle%unit,key,name,'&',stat)

         handle%columns(k)%unit = name
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   FUNCTION addi(key_in,i) result(key)

      character(len=*) :: key_in
      character(len=8) :: key
      integer          :: i

      if (i.lt.10) then
         write(key,'(a,i1)') key_in,i
      elseif (i.lt.100) then
         write(key,'(a,i2)') key_in,i
      else
         write(key,'(a,i3)') key_in,i
      endif

   END FUNCTION


!------------------------------------------------------------------------------


   SUBROUTINE check_column_write(handle,icol,frow,felem,offset)

      type(handle_hdu),intent(in)  :: handle
      integer(idp),    intent(out) :: frow
      integer,         intent(out) :: felem
      integer,         intent(in)  :: icol
      integer(idp),    intent(in)  :: offset
      integer(idp)                 :: idrep

      if (icol.gt.handle%ncols) then
         write(*,*) 'ERROR: Index exceeds number of columns.'
         write(*,'(a,2i5)') ': icol, ncols = ',icol,handle%ncols
         stop
      endif

      idrep = handle%columns(icol)%repcount
      frow  = offset/idrep+1
      felem = mod(offset,idrep)+1

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE write_column_dbin(handle,colnum,dbuffer,n,offset)

      type(handle_hdu),intent(inout)   :: handle
      integer,         intent(in)      :: colnum, n
      real(dp),        intent(in)      :: dbuffer(n)
      integer(idp),intent(in),optional :: offset
      integer                          :: stat, felem, i
      integer(idp)                     :: frow

      if (present(offset)) then
         call check_column_write(handle,colnum,frow,felem,offset)
      else
         frow = 1
         felem = 1
      endif

      stat = 0
      call FTPCLD(handle%unit,colnum,int(frow),felem,n,dbuffer,stat)
      call FTGNRW(handle%unit,handle%nrows,stat)

      call check_errors(stat)

      do i = 1,handle%ncols
         handle%columns(i)%length = handle%columns(i)%repcount*handle%nrows
      enddo

   END SUBROUTINE


   SUBROUTINE write_column_sbin(handle,colnum,sbuffer,n,offset)

      type(handle_hdu),intent(inout)   :: handle
      integer,         intent(in)      :: colnum, n
      real(sp),        intent(in)      :: sbuffer(n)
      integer(idp),intent(in),optional :: offset
      integer                          :: stat, felem, i
      integer(idp)                     :: frow

      if (present(offset)) then
         call check_column_write(handle,colnum,frow,felem,offset)
      else
         frow = 1
         felem = 1
      endif

      stat = 0
      call FTPCLE(handle%unit,colnum,int(frow),felem,n,sbuffer,stat)
      call FTGNRW(handle%unit,handle%nrows,stat)

      call check_errors(stat)

      do i = 1,handle%ncols
         handle%columns(i)%length = handle%columns(i)%repcount*handle%nrows
      enddo

   END SUBROUTINE


   SUBROUTINE write_column_ibin(handle,colnum,ibuffer,n,offset)

      type(handle_hdu),intent(inout)   :: handle
      integer,         intent(in)      :: colnum, n
      integer,         intent(in)      :: ibuffer(n)
      integer(idp),intent(in),optional :: offset
      integer                          :: stat, felem, i
      integer(idp)                     :: frow

      if (present(offset)) then
         call check_column_write(handle,colnum,frow,felem,offset)
      else
         frow = 1
         felem = 1
      endif

      stat = 0
      call FTPCLJ(handle%unit,colnum,int(frow),felem,n,ibuffer(1:n),stat)
      call FTGNRW(handle%unit,handle%nrows,stat)

      call check_errors(stat)

      do i = 1,handle%ncols
         handle%columns(i)%length = handle%columns(i)%repcount*handle%nrows
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE get_key_int(handle,keyname,value)

      integer           :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYJ(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE get_key_logical(handle,keyname,value)

      logical           :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYL(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE get_key_string(handle,keyname,value)

      character(len=*),intent(out) :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYS(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE get_key_real(handle,keyname,value)

      real              :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYE(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE get_key_double(handle,keyname,value)

      double precision  :: value
      type(handle_hdu)  :: handle
      character(len=*)  :: keyname
      character(len=80) :: comment
      integer           :: stat

      stat = 0
      call FTGKYD(handle%unit,trim(keyname),value,comment,stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE put_key_int(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      integer          :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8) keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYJ(handle%unit,trim(keyname),value,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE put_key_logical(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      logical          :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8)  keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYL(handle%unit,trim(keyname),value,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE put_key_string(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      character(len=*) :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8)  keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
!      call FTUKYS(handle%unit,trim(keyname),value,trim(comment),stat)
      call FTUKLS(handle%unit,trim(keyname),trim(adjustl(value)),  &
                  trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE put_key_real(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      real             :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8)  keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYE(handle%unit,trim(keyname),value,-8,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


   SUBROUTINE put_key_double(handle,keyname,value,comment)

      type(handle_hdu) :: handle
      double precision :: value
      character(len=*) :: keyname,comment
      integer          :: stat

      if (len_trim(keyname).gt.8) keyname=keyname(1:8)
      if (len_trim(comment).gt.47) comment=comment(1:47)

      stat = 0
      call FTUKYD(handle%unit,trim(keyname),value,-15,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE fits_put_comment(handle,comment)

      type(handle_hdu) :: handle
      character(len=*) :: comment
      integer          :: stat

      stat = 0
      call FTPCOM(handle%unit,trim(comment),stat)
      call check_errors(stat)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE check_errors(stat)
! Print error messages and stop

      integer,intent(inout) :: stat

      if (stat<=0) return

      call write_errors(stat)
      stop

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE write_errors(stat)
! Print error messages and proceed

      integer,intent(inout) :: stat
      character(len=30)     :: errtext
      character(len=80)     :: errmessage

      if (stat<=0) return

      call ftgerr(stat,errtext)
      write(*,*) 'FITSIO Error Status =',stat
      write(*,*) errtext

      call ftgmsg(errmessage)
      do while (errmessage /= ' ')
         write(*,*) errmessage
         call ftgmsg(errmessage)
      enddo

   END SUBROUTINE


!------------------------------------------------------------------------------

END MODULE
