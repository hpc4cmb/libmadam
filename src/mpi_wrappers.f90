MODULE mpi_wrappers

  use mpi

  ! This module defines wrapper routines for some common mpi routines
  implicit none
  !private
  !include 'mpif.h'

  integer, parameter, public :: null_mpi = mpi_proc_null
  logical, parameter, public :: parallel_version = .true.

  integer, public :: comm

  integer, parameter, private :: sp = kind(1.0)
  integer, parameter, private :: dp = kind(1.0d0)
  integer, parameter, private :: isp = selected_int_kind(9)
  integer, parameter, private :: idp = selected_int_kind(18)

  integer, save, private :: rc, status(mpi_status_size)
  integer, save, private :: id_task, ntasks

   ! RK edit begins
   interface abort_mpi
      module procedure abort_mpi
   end interface
   ! RK edit ends

  interface sum_mpi
     module procedure sum_mpi_dp,   &
          sum_mpi_sp,   &
          sum_mpi_int,  &
          sum_mpi_long, &
          sum_mpi_cmplx, &
          sum_mpi_vec_int, &
          sum_mpi_vec_long, &
          sum_mpi_vec_double, &
          sum_mpi_map_double, &
          sum_mpi_matrix_double, &
          sum_mpi_vec_logical
  end interface

  interface min_mpi
     module procedure min_mpi_dp,   &
          min_mpi_sp,   &
          min_mpi_int
  end interface

  interface max_mpi
     module procedure max_mpi_dp,   &
          max_mpi_sp,   &
          max_mpi_int, &
          max_mpi_long ! -RK
  end interface

  interface send_mpi
     module procedure send_mpi_dp,       &
          send_mpi_sp,       &
          send_mpi_int,      &
          send_mpi_long,     &
          send_mpi_log,      &
          send_mpi_str,      &
          send_mpi_cmplx,    &
          send_mpi_vec_dp,   &
          send_mpi_vec_sp,   &
          send_mpi_vec_int,  &
          send_mpi_vec_long, &
          send_mpi_vec_log,  &
          send_mpi_vec_str,  &
          send_mpi_vec_cmplx
  end interface

  interface broadcast_mpi
     module procedure broadcast_mpi_dp,        &
          broadcast_mpi_sp,        &
          broadcast_mpi_int,       &
          broadcast_mpi_long,      &
          broadcast_mpi_log,       &
          broadcast_mpi_str,       &
          broadcast_mpi_cmplx,     &
          broadcast_mpi_vec_dp,    &
          broadcast_mpi_vec_sp,    &
          broadcast_mpi_vec_int,   &
          broadcast_mpi_vec_long,  &
          broadcast_mpi_vec_log,   &
          broadcast_mpi_vec_str,   &
          broadcast_mpi_vec_cmplx
  end interface

  interface collect_mpi
     module procedure collect_mpi_dp,    &
          collect_mpi_sp,    &
          collect_mpi_int,   &
          collect_mpi_cmplx
  end interface

  interface sendrecv_mpi
     module procedure sendrecv_mpi_dp,   &
          sendrecv_mpi_sp,   &
          sendrecv_mpi_int
  end interface

  integer,allocatable :: sendtag(:), recvtag(:)

  public sum_mpi, min_mpi, max_mpi, broadcast_mpi, send_mpi, sendrecv_mpi, &
       collect_mpi, init_mpi, wait_mpi, close_mpi, &
       broadcast_mpi_vec_dp, broadcast_mpi_vec_sp, broadcast_mpi_vec_int, &
       broadcast_mpi_vec_long, broadcast_mpi_vec_log, &
       broadcast_mpi_vec_str, broadcast_mpi_vec_cmplx, &
       send_mpi_vec_dp, send_mpi_vec_sp, send_mpi_vec_int, &
       send_mpi_vec_long, send_mpi_vec_log, send_mpi_vec_str, &
       send_mpi_vec_cmplx, &
       collect_mpi_dp, collect_mpi_sp, collect_mpi_int, collect_mpi_cmplx, &
       abort_mpi

CONTAINS


  ! RK edit begins
  subroutine abort_mpi(mess)
    character(len=*), optional :: mess

    write (*,*) ''
    if (present(mess)) then
       write (*,*) ' Madam fatal error: '// trim(mess)
    else
       write (*,*) ' Madam fatal error'
    end if
    write (*,*) ''

    call mpi_abort(comm, -1, rc)
    if (rc /= 0) stop 'Even MPI_abort failed ...'

  end subroutine abort_mpi
  ! RK edit ends


  !---------------------------------------------------------------------------


  SUBROUTINE init_mpi(comm_in, ntasks_out, id_out)

    integer, intent(in) :: comm_in
    integer, intent(out) :: ntasks_out, id_out
    !integer :: required, provided

    comm = comm_in

    !call mpi_init(rc)

    !required = MPI_THREAD_FUNNELED
    !call mpi_init_thread(required, provided, rc)
    !if ( rc /= 0 ) stop 'MPI_init failed'

    call mpi_comm_size(comm, ntasks, rc)
    call mpi_comm_rank(comm, id_task, rc)

    !if ( provided < required .and. id_task == 0 ) &
    !     print *,'WARNING: MPI environment did not provide thread support'

    ntasks_out = ntasks
    id_out = id_task

    allocate(sendtag(0:ntasks-1))
    allocate(recvtag(0:ntasks-1))
    sendtag = 0
    recvtag = 0

  END SUBROUTINE init_mpi


  !---------------------------------------------------------------------------


  SUBROUTINE close_mpi()

    !call mpi_finalize(rc)

    deallocate(sendtag,recvtag)

  END SUBROUTINE close_mpi


  !---------------------------------------------------------------------------


  SUBROUTINE wait_mpi()

    call mpi_barrier(comm,rc)

  END SUBROUTINE wait_mpi


  !---------------------------------------------------------------------------
  !
  !  interface sum_mpi
  !

  SUBROUTINE sum_mpi_dp(dsum)

    double precision,intent(inout) :: dsum
    double precision               :: dtemp

    call mpi_allreduce(dsum,dtemp,1,mpi_double_precision,mpi_sum,  &
         comm,rc)
    dsum = dtemp

  END SUBROUTINE sum_mpi_dp


  SUBROUTINE sum_mpi_sp(ssum)

    real,intent(inout) :: ssum
    real               :: stemp

    call mpi_allreduce(ssum,stemp,1,mpi_real,mpi_sum,comm,rc)
    ssum = stemp

  END SUBROUTINE sum_mpi_sp


  SUBROUTINE sum_mpi_int(isum)

    integer,intent(inout) :: isum
    integer               :: itemp

    call mpi_allreduce(isum,itemp,1,mpi_integer,mpi_sum,comm,rc)
    isum = itemp

  END SUBROUTINE sum_mpi_int


  SUBROUTINE sum_mpi_long(isum)

    integer(idp),intent(inout) :: isum
    integer(idp)               :: itemp

    call mpi_allreduce(isum,itemp,1,mpi_integer8,mpi_sum,comm,rc)
    isum = itemp

  END SUBROUTINE sum_mpi_long


  SUBROUTINE sum_mpi_cmplx(csum)

    complex(dp),intent(inout) :: csum
    complex(dp)               :: ctemp

    call mpi_allreduce(csum,ctemp,1,mpi_double_complex,mpi_sum,  &
         comm,rc)
    csum = ctemp

  END SUBROUTINE sum_mpi_cmplx

  ! RK additions

  SUBROUTINE sum_mpi_vec_int(ivec)

    integer, intent(inout) :: ivec(:)
    integer, allocatable :: itemp(:)
    integer :: ierr, n

    n = size(ivec)

    allocate(itemp(n), stat=ierr)
    if (ierr /= 0) call abort_mpi('sum_mpi_vec_int failed to allocate.')

    call mpi_allreduce(ivec, itemp, n, mpi_integer, mpi_sum, comm, rc)

    ivec = itemp

    deallocate(itemp)

  END SUBROUTINE sum_mpi_vec_int


  SUBROUTINE sum_mpi_vec_long(lvec)

    integer(idp), intent(inout) :: lvec(:)
    integer(idp), allocatable :: ltemp(:)
    integer :: ierr, n

    n = size(lvec)

    allocate(ltemp(n), stat=ierr)
    if (ierr /= 0) call abort_mpi('sum_mpi_vec_long failed to allocate.')

    call mpi_allreduce(lvec, ltemp, n, mpi_integer8, mpi_sum, comm, rc)

    lvec = ltemp

    deallocate(ltemp)

  END SUBROUTINE sum_mpi_vec_long


  SUBROUTINE sum_mpi_vec_double(dvec)

    real(dp), intent(inout) :: dvec(:)
    real(dp), allocatable :: dtemp(:)
    integer :: ierr, n

    n = size(dvec)

    allocate(dtemp(n), stat=ierr)
    if (ierr /= 0) call abort_mpi('sum_mpi_vec_double failed to allocate.')

    call mpi_allreduce(dvec, dtemp, n, MPI_DOUBLE_PRECISION, mpi_sum, comm, rc)

    dvec = dtemp

    deallocate(dtemp)

  END SUBROUTINE sum_mpi_vec_double


  SUBROUTINE sum_mpi_map_double(dmap)

    real(dp), intent(inout) :: dmap(:, :)
    real(dp), allocatable :: dtemp(:, :)
    integer :: ierr, nmap, npix

    nmap = size(dmap, 1)
    npix = size(dmap, 2)

    allocate(dtemp(nmap, npix), stat=ierr)
    if (ierr /= 0) call abort_mpi('sum_mpi_map_double failed to allocate.')

    call mpi_allreduce(dmap, dtemp, nmap * npix, MPI_DOUBLE_PRECISION, &
         mpi_sum, comm, rc)

    dmap = dtemp

    deallocate(dtemp)

  END SUBROUTINE sum_mpi_map_double


  SUBROUTINE sum_mpi_matrix_double(dmatrix)

    real(dp), intent(inout) :: dmatrix(:, :, :)
    real(dp), allocatable :: dtemp(:, :, :)
    integer :: ierr, nmap1, nmap2, npix

    nmap1 = size(dmatrix, 1)
    nmap2 = size(dmatrix, 2)
    npix = size(dmatrix, 3)

    allocate(dtemp(nmap1, nmap2, npix), stat=ierr)
    if (ierr /= 0) call abort_mpi('sum_mpi_matrix_double failed to allocate.')

    call mpi_allreduce(dmatrix, dtemp, nmap1 * nmap2 * npix, &
         MPI_DOUBLE_PRECISION, mpi_sum, comm, rc)

    dmatrix = dtemp

    deallocate(dtemp)

  END SUBROUTINE sum_mpi_matrix_double


  SUBROUTINE sum_mpi_vec_logical(bvec)

    logical, intent(inout) :: bvec(:)
    logical, allocatable :: btemp(:)
    integer :: ierr, n

    n = size(bvec)

    allocate(btemp(n), stat=ierr)
    if (ierr /= 0) call abort_mpi('sum_mpi_vec_logical failed to allocate.')

    call mpi_allreduce(bvec, btemp, n, MPI_LOGICAL, MPI_LOR, comm, rc)

    bvec = btemp

    deallocate(btemp)

  END SUBROUTINE sum_mpi_vec_logical


  !-------------------------------------------------------------------------
  !
  !  interface max_mpi
  !

  SUBROUTINE max_mpi_dp(dmax)

    double precision,intent(inout) :: dmax
    double precision               :: dtemp

    call mpi_allreduce(dmax,dtemp,1,mpi_double_precision,mpi_max,  &
         comm,rc)
    dmax = dtemp

  END SUBROUTINE max_mpi_dp


  SUBROUTINE max_mpi_sp(smax)

    real,intent(inout) :: smax
    real               :: stemp

    call mpi_allreduce(smax,stemp,1,mpi_real,mpi_max,comm,rc)
    smax = stemp

  END SUBROUTINE max_mpi_sp


  SUBROUTINE max_mpi_int(imax)

    integer,intent(inout) :: imax
    integer               :: itemp

    call mpi_allreduce(imax,itemp,1,mpi_integer,mpi_max,comm,rc)
    imax = itemp

  END SUBROUTINE max_mpi_int


  SUBROUTINE max_mpi_long(imax)

    integer(idp),intent(inout) :: imax
    integer(idp)               :: itemp

    call mpi_allreduce(imax,itemp,1,mpi_integer8,mpi_max,comm,rc)
    imax = itemp

  END SUBROUTINE max_mpi_long


  !---------------------------------------------------------------------------
  !
  !  interface min_mpi
  !

  SUBROUTINE min_mpi_dp(dmin)

    double precision,intent(inout) :: dmin
    double precision               :: dtemp

    call mpi_allreduce(dmin,dtemp,1,mpi_double_precision,mpi_min,  &
         comm,rc)
    dmin = dtemp

  END SUBROUTINE min_mpi_dp


  SUBROUTINE min_mpi_sp(smin)

    real,intent(inout) :: smin
    real               :: stemp

    call mpi_allreduce(smin,stemp,1,mpi_real,mpi_min,comm,rc)
    smin = stemp

  END SUBROUTINE min_mpi_sp


  SUBROUTINE min_mpi_int(imin)

    integer,intent(inout) :: imin
    integer               :: itemp

    call mpi_allreduce(imin,itemp,1,mpi_integer,mpi_min,comm,rc)
    imin = itemp

  END SUBROUTINE min_mpi_int


  !---------------------------------------------------------------------------
  !
  !  interface send_mpi
  !

  SUBROUTINE send_mpi_dp(dx,id_send,id_recv)

    integer,         intent(in)    :: id_send, id_recv
    double precision,intent(inout) :: dx
    integer                        :: tag

    if (id_recv==id_send) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(dx,1,mpi_double_precision,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(dx,1,mpi_double_precision,id_send,tag,comm, &
            status,rc)
    endif

  END SUBROUTINE send_mpi_dp


  SUBROUTINE send_mpi_sp(sx,id_send,id_recv)

    integer,intent(in)    :: id_send, id_recv
    real,   intent(inout) :: sx
    integer               :: tag

    if (id_recv==id_send) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(sx,1,mpi_real,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(sx,1,mpi_real,id_send,tag,comm,status,rc)

    endif

  END SUBROUTINE send_mpi_sp


  SUBROUTINE send_mpi_int(n,id_send,id_recv)

    integer,intent(in)    :: id_send, id_recv
    integer,intent(inout) :: n
    integer               :: tag

    if (id_recv==id_send) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(n,1,mpi_integer,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(n,1,mpi_integer,id_send,tag,comm,status,rc)

    endif

  END SUBROUTINE send_mpi_int


  SUBROUTINE send_mpi_long(n,id_send,id_recv)

    integer,     intent(in)    :: id_send, id_recv
    integer(idp),intent(inout) :: n
    integer                    :: tag

    if (id_recv==id_send) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(n,1,mpi_integer8,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(n,1,mpi_integer8,id_send,tag,comm, &
            status,rc)

    endif

  END SUBROUTINE send_mpi_long


  SUBROUTINE send_mpi_log(n,id_send,id_recv)

    integer,intent(in)    :: id_send, id_recv
    logical,intent(inout) :: n
    integer               :: tag

    if (id_recv==id_send) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(n,1,mpi_logical,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(n,1,mpi_logical,id_send,tag,comm,status,rc)

    endif

  END SUBROUTINE send_mpi_log


  SUBROUTINE send_mpi_str(s,id_send,id_recv)

    integer,         intent(in)    :: id_send, id_recv
    character(len=*),intent(inout) :: s
    integer                        :: tag, n

    if (id_recv==id_send) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       n = len_trim(s)
       call mpi_send(n,1,mpi_integer,id_recv,tag,comm,rc) !fix 23.06.08
       call mpi_send(s,n,mpi_character,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       s = ''
       call mpi_recv(n,1,mpi_integer,id_send,tag,comm,status,rc)
       call mpi_recv(s,n,mpi_character,id_send,tag,comm,status,rc)

    endif

  END SUBROUTINE send_mpi_str


  SUBROUTINE send_mpi_cmplx(cx,id_send,id_recv)

    integer,    intent(in)    :: id_send, id_recv
    complex(dp),intent(inout) :: cx
    integer                   :: tag

    if (id_recv==id_send) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(cx,1,mpi_double_complex,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(cx,1,mpi_double_complex,id_send,tag,    &
            comm,status,rc)
    endif

  END SUBROUTINE send_mpi_cmplx


  !------------------------------------------------------------------------------
  !

  SUBROUTINE send_mpi_vec_dp(dbuffer,n,id_send,id_recv)

    integer,         intent(in)    :: n,id_send, id_recv
    double precision,intent(inout) :: dbuffer(n)
    integer                        :: tag

    if (id_recv==id_send) return
    if (n.le.0) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(dbuffer,n,mpi_double_precision,id_recv,tag,  &
            comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(dbuffer,n,mpi_double_precision,id_send,tag,   &
            comm,status,rc)
    endif

  END SUBROUTINE send_mpi_vec_dp


  SUBROUTINE send_mpi_vec_sp(sbuffer,n,id_send,id_recv)

    integer,intent(in)    :: n,id_send, id_recv
    real,   intent(inout) :: sbuffer(n)
    integer               :: tag

    if (id_recv==id_send) return
    if (n.le.0) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(sbuffer,n,mpi_real,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(sbuffer,n,mpi_real,id_send,tag,comm,status,rc)
    endif

  END SUBROUTINE send_mpi_vec_sp


  SUBROUTINE send_mpi_vec_int(ibuffer,n,id_send,id_recv)

    integer,intent(in)    :: n,id_send, id_recv
    integer,intent(inout) :: ibuffer(n)
    integer               :: tag

    if (id_recv==id_send) return
    if (n.le.0) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(ibuffer,n,mpi_integer,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(ibuffer,n,mpi_integer,id_send,tag,comm,  &
            status,rc)
    endif

  END SUBROUTINE send_mpi_vec_int


  SUBROUTINE send_mpi_vec_long(ibuffer,n,id_send,id_recv)

    integer,     intent(in)    :: n, id_send, id_recv
    integer(idp),intent(inout) :: ibuffer(n)
    integer                    :: tag

    if (id_recv==id_send) return
    if (n.le.0) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(ibuffer,n,mpi_integer8,id_recv,tag,  &
            comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(ibuffer,n,mpi_integer8,id_send,tag,  &
            comm,status,rc)
    endif

  END SUBROUTINE send_mpi_vec_long


  SUBROUTINE send_mpi_vec_log(ibuffer,n,id_send,id_recv)

    integer, intent(in)    :: n, id_send, id_recv
    logical, intent(inout) :: ibuffer(n)
    integer                :: tag

    if (id_recv==id_send) return
    if (n.le.0) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(ibuffer,n,mpi_logical,id_recv,tag,comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(ibuffer,n,mpi_logical,id_send,tag,comm,  &
            status,rc)
    endif

  END SUBROUTINE send_mpi_vec_log


  SUBROUTINE send_mpi_vec_str(sbuffer,n,id_send,id_recv)

    integer,         intent(in)    :: n, id_send, id_recv
    character(len=*),intent(inout) :: sbuffer(n)
    integer                        :: tag, i, m

    if (id_recv==id_send) return
    if (n.le.0) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       do i = 1,n
          m = len_trim(sbuffer(i))
          call mpi_send(m,1,mpi_integer,id_recv,tag,comm,rc)
          call mpi_send(sbuffer(i),m,mpi_character,id_recv,tag,  &
               comm,rc)
       enddo

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       do i = 1,n
          sbuffer(i) = ''
          call mpi_recv(m,1,mpi_integer,id_send,tag,comm,  &
               status,rc)
          call mpi_recv(sbuffer(i),m,mpi_character,id_send,tag,      &
               comm,status,rc)
       enddo

    endif

  END SUBROUTINE send_mpi_vec_str


  SUBROUTINE send_mpi_vec_cmplx(cbuffer,n,id_send,id_recv)

    integer,    intent(in)    :: n,id_send, id_recv
    complex(dp),intent(inout) :: cbuffer(n)
    integer                   :: tag

    if (id_recv==id_send) return
    if (n.le.0) return

    if (id_task==id_send) then

       tag = sendtag(id_recv)
       tag = tag+1
       sendtag(id_recv) = tag

       call mpi_send(cbuffer,n,mpi_double_complex,id_recv,tag,  &
            comm,rc)

    elseif (id_task==id_recv) then

       tag = recvtag(id_send)
       tag = tag+1
       recvtag(id_send) = tag

       call mpi_recv(cbuffer,n,mpi_double_complex,id_send,tag,   &
            comm,status,rc)
    endif

  END SUBROUTINE send_mpi_vec_cmplx


  !------------------------------------------------------------------------------

  ! interface broadcast:
  ! Send a parameter to all processes

  SUBROUTINE broadcast_mpi_dp(dx,id)

    double precision :: dx
    integer          :: id

    call mpi_bcast(dx,1,mpi_double_precision,id,comm,rc)

  END SUBROUTINE broadcast_mpi_dp


  SUBROUTINE broadcast_mpi_sp(sx,id)

    real    :: sx
    integer :: id

    call mpi_bcast(sx,1,mpi_real,id,comm,rc)

  END SUBROUTINE broadcast_mpi_sp


  SUBROUTINE broadcast_mpi_int(n,id)

    integer :: n, id

    call mpi_bcast(n,1,mpi_integer,id,comm,rc)

  END SUBROUTINE broadcast_mpi_int


  SUBROUTINE broadcast_mpi_long(n,id)

    integer(idp) :: n
    integer      :: id

    call mpi_bcast(n,1,mpi_integer8,id,comm,rc)

  END SUBROUTINE broadcast_mpi_long


  SUBROUTINE broadcast_mpi_log(flag,id)

    logical :: flag
    integer :: id

    call mpi_bcast(flag,1,mpi_logical,id,comm,rc)

  END SUBROUTINE broadcast_mpi_log


  SUBROUTINE broadcast_mpi_str(s,id)

    character(len=*) :: s
    integer          :: n, id

    n = len(s)
    call mpi_bcast(s,n,mpi_character,id,comm,rc)

  END SUBROUTINE broadcast_mpi_str


  SUBROUTINE broadcast_mpi_cmplx(c,id)

    complex(dp) :: c
    integer     :: id

    call mpi_bcast(c,1,mpi_double_complex,id,comm,rc)

  END SUBROUTINE broadcast_mpi_cmplx


  !------------------------------------------------------------------------------


  SUBROUTINE broadcast_mpi_vec_dp(dbuffer,n,id)

    double precision :: dbuffer(n)
    integer          :: n, id

    call mpi_bcast(dbuffer,n,mpi_double_precision,id,comm,rc)

  END SUBROUTINE broadcast_mpi_vec_dp


  SUBROUTINE broadcast_mpi_vec_sp(sbuffer,n,id)

    real    :: sbuffer(n)
    integer :: n, id

    call mpi_bcast(sbuffer,n,mpi_real,id,comm,rc)

  END SUBROUTINE broadcast_mpi_vec_sp


  SUBROUTINE broadcast_mpi_vec_int(ibuffer,n,id)

    integer :: ibuffer(n)
    integer :: n, id

    call mpi_bcast(ibuffer,n,mpi_integer,id,comm,rc)

  END SUBROUTINE broadcast_mpi_vec_int


  SUBROUTINE broadcast_mpi_vec_long(ibuffer,n,id)

    integer(idp) :: ibuffer(n)
    integer      :: n, id

    call mpi_bcast(ibuffer,n,mpi_integer8,id,comm,rc)

  END SUBROUTINE broadcast_mpi_vec_long


  SUBROUTINE broadcast_mpi_vec_log(lbuffer,n,id)

    logical :: lbuffer(n)
    integer :: n, id

    call mpi_bcast(lbuffer,n,mpi_logical,id,comm,rc)

  END SUBROUTINE broadcast_mpi_vec_log


  SUBROUTINE broadcast_mpi_vec_str(sbuffer,n,id)

    character(len=*) :: sbuffer(n)
    integer          :: n, id, i, m

    m = len(sbuffer)
    do i = 1,n
       call mpi_bcast(sbuffer(i),m,mpi_character,id,comm,rc)
    enddo

  END SUBROUTINE broadcast_mpi_vec_str


  SUBROUTINE broadcast_mpi_vec_cmplx(cbuffer,n,id)

    complex(dp) :: cbuffer(n)
    integer     :: n, id

    call mpi_bcast(cbuffer,n,mpi_double_complex,id,comm,rc)

  END SUBROUTINE broadcast_mpi_vec_cmplx


  !------------------------------------------------------------------------------
  !
  !  Interface collect_mpi
  !
  !  Sum the input array over all processes
  !  and store the result into another array.

  SUBROUTINE collect_mpi_int(ibuffer,isum,n,id)

    integer             :: n, id
    integer,intent(in)  :: ibuffer(n)
    integer,intent(out) :: isum(n)

    call mpi_reduce(ibuffer,isum,n,mpi_integer,mpi_sum,id,comm,rc)

  END SUBROUTINE collect_mpi_int


  SUBROUTINE collect_mpi_sp(sbuffer,ssum,n,id)

    integer          :: n, id
    real,intent(in)  :: sbuffer(n)
    real,intent(out) :: ssum(n)

    call mpi_reduce(sbuffer,ssum,n,mpi_real,mpi_sum,id,comm,rc)

  END SUBROUTINE collect_mpi_sp


  SUBROUTINE collect_mpi_dp(dbuffer,dsum,n,id)

    integer                      :: n, id
    double precision,intent(in)  :: dbuffer(n)
    double precision,intent(out) :: dsum(n)

    call mpi_reduce(dbuffer,dsum,n,mpi_double_precision,mpi_sum,id,  &
         comm,rc)

  END SUBROUTINE collect_mpi_dp


  SUBROUTINE collect_mpi_cmplx(cbuffer,csum,n,id)

    integer                 :: n, id
    complex(dp),intent(in)  :: cbuffer(n)
    complex(dp),intent(out) :: csum(n)

    call mpi_reduce(cbuffer,csum,n,mpi_double_complex,mpi_sum,id,  &
         comm,rc)

  END SUBROUTINE collect_mpi_cmplx


  !------------------------------------------------------------------------------


  SUBROUTINE sendrecv_mpi_dp(sendbuffer,recvbuffer,n,id_recv,id_send)

    integer          :: n, id_send, id_recv
    double precision :: sendbuffer(n), recvbuffer(n)

    recvbuffer = 0.0

    call mpi_sendrecv(sendbuffer,n,mpi_double_precision,id_recv,18,   &
         recvbuffer,n,mpi_double_precision,id_send,18,   &
         comm,status,rc)

  END SUBROUTINE sendrecv_mpi_dp


  SUBROUTINE sendrecv_mpi_sp(sendbuffer,recvbuffer,n,id_recv,id_send)

    integer :: n, id_send, id_recv
    real    :: sendbuffer(n), recvbuffer(n)

    recvbuffer = 0.0

    call mpi_sendrecv(sendbuffer,n,mpi_real,id_recv,18,   &
         recvbuffer,n,mpi_real,id_send,18,   &
         comm,status,rc)

  END SUBROUTINE sendrecv_mpi_sp


  SUBROUTINE sendrecv_mpi_int(sendbuffer,recvbuffer,n,id_recv,id_send)

    integer :: n, id_send, id_recv
    integer :: sendbuffer(n), recvbuffer(n)

    recvbuffer = 0.0

    call mpi_sendrecv(sendbuffer,n,mpi_integer,id_recv,18,   &
         recvbuffer,n,mpi_integer,id_send,18,   &
         comm,status,rc)

  END SUBROUTINE sendrecv_mpi_int


  !------------------------------------------------------------------------------

END MODULE mpi_wrappers
