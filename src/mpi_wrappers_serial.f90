MODULE mpi_wrappers

! This module defines wrapper routines for some common mpi routines
! Serial version

   implicit none
   private

   integer, parameter, public :: null_mpi = -1
   logical, parameter, public :: parallel_version = .false.

   integer, parameter :: sp = kind(1.0)
   integer, parameter :: dp = kind(1.0d0)
   integer, parameter :: isp = selected_int_kind(9)
   integer, parameter :: idp = selected_int_kind(18)

   integer, save :: id_task, ntasks

   interface sum_mpi
      module procedure sum_mpi_dp,   &
                       sum_mpi_sp,   &
                       sum_mpi_int,  &
                       sum_mpi_long, &
                       sum_mpi_cmplx
   end interface

   interface min_mpi
      module procedure min_mpi_dp,   &
                       min_mpi_sp,   &
                       min_mpi_int
   end interface

   interface max_mpi
      module procedure max_mpi_dp,   &
                       max_mpi_sp,   &
                       max_mpi_int
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
          collect_mpi_dp, collect_mpi_sp, collect_mpi_int, collect_mpi_cmplx

 CONTAINS

!------------------------------------------------------------------------------


   SUBROUTINE init_mpi(ntasks_out,id_out)

      integer,intent(out) :: ntasks_out, id_out

      ntasks = 1
      id_task = 0

      ntasks_out = ntasks
      id_out = id_task

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE close_mpi()

      return

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE wait_mpi()

      return

   END SUBROUTINE


!------------------------------------------------------------------------------
!
!  interface sum_mpi
!

   SUBROUTINE sum_mpi_dp(dsum)

      double precision,intent(inout) :: dsum
      return

   END SUBROUTINE


   SUBROUTINE sum_mpi_sp(ssum)

      real,intent(inout) :: ssum
      return

   END SUBROUTINE


   SUBROUTINE sum_mpi_int(isum)

      integer,intent(inout) :: isum
      return

   END SUBROUTINE


   SUBROUTINE sum_mpi_long(isum)

      integer(idp),intent(inout) :: isum
      return

   END SUBROUTINE


   SUBROUTINE sum_mpi_cmplx(csum)

      complex(dp),intent(inout) :: csum
      return

   END SUBROUTINE


!------------------------------------------------------------------------------
!
!  interface max_mpi
!

   SUBROUTINE max_mpi_dp(dmax)

      double precision,intent(inout) :: dmax
      return

   END SUBROUTINE


   SUBROUTINE max_mpi_sp(smax)

      real,intent(inout) :: smax
      return

   END SUBROUTINE


   SUBROUTINE max_mpi_int(imax)

      integer,intent(inout) :: imax
      return

   END SUBROUTINE


!------------------------------------------------------------------------------
!
!  interface min_mpi
!

   SUBROUTINE min_mpi_dp(dmin)

      double precision,intent(inout) :: dmin
      return

   END SUBROUTINE


   SUBROUTINE min_mpi_sp(smin)

      real,intent(inout) :: smin
      return

   END SUBROUTINE


   SUBROUTINE min_mpi_int(imin)

      integer,intent(inout) :: imin
      return

   END SUBROUTINE


!------------------------------------------------------------------------------
!
!  interface send_mpi
!

   SUBROUTINE send_mpi_dp(dx,id_send,id_recv)

      integer,         intent(in)    :: id_send, id_recv
      double precision,intent(inout) :: dx
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_sp(sx,id_send,id_recv)

      integer,intent(in)    :: id_send, id_recv
      real,   intent(inout) :: sx
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_int(n,id_send,id_recv)

      integer,intent(in)    :: id_send, id_recv
      integer,intent(inout) :: n
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_long(n,id_send,id_recv)

      integer,     intent(in)    :: id_send, id_recv
      integer(idp),intent(inout) :: n
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_log(n,id_send,id_recv)

      integer,intent(in)    :: id_send, id_recv
      logical,intent(inout) :: n
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_str(s,id_send,id_recv)

      integer,         intent(in)    :: id_send, id_recv
      character(len=*),intent(inout) :: s
      return

   END SUBROUTINE



   SUBROUTINE send_mpi_cmplx(cx,id_send,id_recv)

      integer,    intent(in)    :: id_send, id_recv
      complex(dp),intent(inout) :: cx
      return

   END SUBROUTINE


!------------------------------------------------------------------------------
!

   SUBROUTINE send_mpi_vec_dp(dbuffer,n,id_send,id_recv)

      integer,         intent(in)    :: n,id_send, id_recv
      double precision,intent(inout) :: dbuffer(n)
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_vec_sp(sbuffer,n,id_send,id_recv)

      integer,intent(in)    :: n,id_send, id_recv
      real,   intent(inout) :: sbuffer(n)
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_vec_int(ibuffer,n,id_send,id_recv)

      integer,intent(in)    :: n,id_send, id_recv
      integer,intent(inout) :: ibuffer(n)
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_vec_long(ibuffer,n,id_send,id_recv)

      integer,     intent(in)    :: n, id_send, id_recv
      integer(idp),intent(inout) :: ibuffer(n)
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_vec_log(ibuffer,n,id_send,id_recv)

      integer, intent(in)    :: n, id_send, id_recv
      logical, intent(inout) :: ibuffer(n)
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_vec_str(sbuffer,n,id_send,id_recv)

      integer,         intent(in)    :: n, id_send, id_recv
      character(len=*),intent(inout) :: sbuffer(n)
      return

   END SUBROUTINE


   SUBROUTINE send_mpi_vec_cmplx(cbuffer,n,id_send,id_recv)

      integer,    intent(in)    :: n, id_send, id_recv
      complex(dp),intent(inout) :: cbuffer(n)
      return

   END SUBROUTINE


!------------------------------------------------------------------------------

! interface broadcast:
! Send a parameter to all processes

   SUBROUTINE broadcast_mpi_dp(dx,id)

      double precision :: dx
      integer          :: id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_sp(sx,id)

      real    :: sx
      integer :: id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_int(n,id)

      integer :: n, id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_long(n,id)

      integer(idp) :: n
      integer      :: id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_log(flag,id)

      logical :: flag
      integer :: id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_str(s,id)

      character(len=*) :: s
      integer          :: id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_cmplx(c,id)

      complex(dp) :: c
      integer     :: id
      return

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE broadcast_mpi_vec_dp(dbuffer,n,id)

      double precision :: dbuffer(n)
      integer          :: n, id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_vec_sp(sbuffer,n,id)

      real    :: sbuffer(n)
      integer :: n, id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_vec_int(ibuffer,n,id)

      integer :: ibuffer(n)
      integer :: n, id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_vec_long(ibuffer,n,id)

      integer(idp) :: ibuffer(n)
      integer      :: n, id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_vec_log(lbuffer,n,id)

      logical :: lbuffer(n)
      integer :: n, id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_vec_str(sbuffer,n,id)

      character(len=*) :: sbuffer(n)
      integer          :: n, id
      return

   END SUBROUTINE


   SUBROUTINE broadcast_mpi_vec_cmplx(cbuffer,n,id)

      complex(dp) :: cbuffer(n)
      integer     :: n, id
      return

   END SUBROUTINE


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

      isum(1:n) = ibuffer(1:n)

   END SUBROUTINE


   SUBROUTINE collect_mpi_sp(sbuffer,ssum,n,id)

      integer          :: n, id
      real,intent(in)  :: sbuffer(n)
      real,intent(out) :: ssum(n)

      ssum(1:n) = sbuffer(1:n)

   END SUBROUTINE


   SUBROUTINE collect_mpi_dp(dbuffer,dsum,n,id)

      integer                      :: n, id
      double precision,intent(in)  :: dbuffer(n)
      double precision,intent(out) :: dsum(n)

      dsum(1:n) = dbuffer(1:n)

   END SUBROUTINE


   SUBROUTINE collect_mpi_cmplx(cbuffer,csum,n,id)

      integer                 :: n, id
      complex(dp),intent(in)  :: cbuffer(n)
      complex(dp),intent(out) :: csum(n)

      csum(1:n) = cbuffer(1:n)

   END SUBROUTINE


!------------------------------------------------------------------------------


   SUBROUTINE sendrecv_mpi_dp(sendbuffer,recvbuffer,n,id_recv,id_send)

      integer          :: n, id_send, id_recv
      double precision :: sendbuffer(n), recvbuffer(n)

      recvbuffer(1:n) = 0.0

   END SUBROUTINE


   SUBROUTINE sendrecv_mpi_sp(sendbuffer,recvbuffer,n,id_recv,id_send)

      integer :: n, id_send, id_recv
      real    :: sendbuffer(n), recvbuffer(n)

      recvbuffer(1:n) = 0.0

   END SUBROUTINE


   SUBROUTINE sendrecv_mpi_int(sendbuffer,recvbuffer,n,id_recv,id_send)

      integer :: n, id_send, id_recv
      integer :: sendbuffer(n), recvbuffer(n)

      recvbuffer(1:n) = 0

   END SUBROUTINE



!------------------------------------------------------------------------------

END MODULE
