MODULE submap_transfer

  use commonparam
  use mpi_wrappers
  implicit none
  private

  interface get_submap
     module procedure get_submap_dp,   get_submap_sp,    &
          get_submap_dp2sp,get_submap_int,   &
          get_submap_dp_s, get_submap_sp_s,  &
          get_submap_int_s, get_submap_dp_2d
  end interface

  public get_submap, get_submap_int

CONTAINS

  !---------------------------------------------------------------------------

  ! Identify the index of a given submap in the local array
  FUNCTION get_m (isubmap) result(m)
    integer, intent(in) :: isubmap
    integer m

    m = count (id_submap(0:isubmap)==ID)
  END FUNCTION get_m

  SUBROUTINE get_submap_dp(dbuffer,nosubpix,map,nomaps,imap,isubmap,id_recv)
    integer, intent(in)  :: nomaps, nosubpix
    real(dp),intent(out) :: dbuffer(nosubpix)
    real(dp),intent(in)  :: map(nomaps,nosubpix,nosubmaps)
    integer, intent(in)  :: imap, isubmap, id_recv

    if (ID==id_submap(isubmap)) dbuffer = map(imap,:,get_m(isubmap))
    call send_mpi(dbuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_dp


  SUBROUTINE get_submap_dp_2d(dbuffer, nosubpix, map, nomaps, imap, jmap, isubmap, id_recv)
    integer, intent(in)  :: nomaps, nosubpix
    real(dp),intent(out) :: dbuffer(nosubpix)
    real(dp),intent(in)  :: map(nomaps,nomaps,nosubpix,nosubmaps)
    integer, intent(in)  :: imap, jmap, isubmap, id_recv

    !print *,id,' : get_submap_dp_2d : ',isubmap,imap,id_submap(isubmap),id_recv ! debug

    if (ID==id_submap(isubmap)) dbuffer = map(imap,jmap,:,get_m(isubmap))
    call send_mpi(dbuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_dp_2d


  SUBROUTINE get_submap_dp2sp(sbuffer,nosubpix,map,nomaps,imap,isubmap,id_recv)
    integer, intent(in)  :: nomaps, nosubpix
    real(sp),intent(out) :: sbuffer(nosubpix)
    real(dp),intent(in)  :: map(nomaps,nosubpix,nosubmaps)
    integer, intent(in)  :: imap, isubmap, id_recv

    if (ID==id_submap(isubmap)) sbuffer = map(imap,:,get_m(isubmap))
    call send_mpi(sbuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_dp2sp


  SUBROUTINE get_submap_sp(sbuffer,nosubpix,map,nomaps,imap,isubmap,id_recv)
    integer, intent(in)  :: nomaps, nosubpix
    real(sp),intent(out) :: sbuffer(nosubpix)
    real(sp),intent(in)  :: map(nomaps,nosubpix,nosubmaps)
    integer, intent(in)  :: imap, isubmap, id_recv

    !print *,id,' : get_submap_sp : ',isubmap,id_submap(isubmap) ! debug

    if (ID==id_submap(isubmap)) sbuffer = map(imap,:,get_m(isubmap))
    call send_mpi(sbuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_sp


  SUBROUTINE get_submap_int(ibuffer,nosubpix,map,nomaps,imap,isubmap,id_recv)
    integer,intent(in)  :: nomaps, nosubpix
    integer,intent(out) :: ibuffer(nosubpix)
    integer,intent(in)  :: map(nomaps,nosubpix,nosubmaps)
    integer,intent(in)  :: imap, isubmap, id_recv

    if (ID==id_submap(isubmap)) ibuffer = map(imap,:,get_m(isubmap))
    call send_mpi(ibuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_int


  SUBROUTINE get_submap_dp_s(dbuffer,nosubpix,map,isubmap,id_recv)
    integer, intent(in)  :: nosubpix
    real(dp),intent(out) :: dbuffer(nosubpix)
    real(dp),intent(in)  :: map(nosubpix,nosubmaps)
    integer, intent(in)  :: isubmap, id_recv

    if (ID==id_submap(isubmap)) dbuffer = map(:,get_m(isubmap))
    call send_mpi(dbuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_dp_s


  SUBROUTINE get_submap_sp_s(sbuffer,nosubpix,map,isubmap,id_recv)
    integer, intent(in)  :: nosubpix
    real(sp),intent(out) :: sbuffer(nosubpix)
    real(sp),intent(in)  :: map(nosubpix,nosubmaps)
    integer, intent(in)  :: isubmap, id_recv

    if (ID==id_submap(isubmap)) sbuffer = map(:,get_m(isubmap))
    call send_mpi(sbuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_sp_s


  SUBROUTINE get_submap_int_s(ibuffer,nosubpix,map,isubmap,id_recv)
    integer, intent(in)  :: nosubpix
    integer, intent(out) :: ibuffer(nosubpix)
    integer, intent(in)  :: map(nosubpix,nosubmaps)
    integer, intent(in)  :: isubmap, id_recv

    if (ID==id_submap(isubmap)) ibuffer = map(:,get_m(isubmap))
    call send_mpi(ibuffer,nosubpix,id_submap(isubmap),id_recv)
  END SUBROUTINE get_submap_int_s

  !------------------------------------------------------------------------------

END MODULE submap_transfer
