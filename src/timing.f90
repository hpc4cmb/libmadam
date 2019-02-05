MODULE timing

   implicit none
   private

   !double precision,save :: tagtime(0:19) = 0.0 -RK
   double precision, save :: tagtime(0:99) = 0.0 ! -RK

   public get_time, get_time_and_reset, reset_time

 CONTAINS

!------------------------------------------------------------------------------


   double precision FUNCTION get_time(tag)

      integer, intent(in) :: tag
      double precision :: time
      integer :: t(8)

      call date_and_time(values=t)

      time = t(3)*86400.d0 +t(5)*3600.d0 +t(6)*60.d0 +t(7) +t(8)*1.d-3

      get_time = time-tagtime(tag)

      do
         if (get_time.ge.0) exit
         get_time = get_time+86400.d0
      enddo

   END FUNCTION


!------------------------------------------------------------------------------


   double precision FUNCTION get_time_and_reset(tag)

      integer, intent(in) :: tag
      double precision :: time
      integer :: t(8)

      call date_and_time(values=t)

      time = t(3)*86400.d0 +t(5)*3600.d0 +t(6)*60.d0 +t(7) +t(8)*1.d-3

      get_time_and_reset = time-tagtime(tag)

      do
         if (get_time_and_reset.ge.0) exit
         get_time_and_reset = get_time_and_reset+86400.d0
      end do

      tagtime(tag) = time

   END FUNCTION


!------------------------------------------------------------------------------


   SUBROUTINE reset_time(tag)

      integer, intent(in), optional :: tag
      double precision :: time
      integer :: t(8)

      call date_and_time(values=t)

      time = t(3)*86400.d0 +t(5)*3600.d0 +t(6)*60.d0 +t(7) +t(8)*1.d-3

      if (present(tag)) then
         tagtime(tag) = time
      else
         tagtime = time
      end if

   END SUBROUTINE


!------------------------------------------------------------------------------

END MODULE
