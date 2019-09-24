MODULE maps_and_baselines

  ! This module defines some allocatable arrays needed by Madam
  !  and provides routines for the (de)allocation of them.

  use planck_config, only : sp, dp, i4b, i8b
  use commonparam, only : basis_functions, basis_order, id, info, kfirst, &
       noba_short, nodetectors, do_binmap, do_hits, do_mask, &
       ndetset, nmap, nopix_cross, nopix_map, nsurvey, use_inmask
  use mpi_wrappers, only : min_mpi, max_mpi, sum_mpi
  use memory_and_time, only : check_stat, write_memory

  implicit none
  private

  real(dp), allocatable, public :: nna(:, :, :, :), nna_inv(:, :, :, :)
  real(dp), allocatable, public :: aa(:, :, :), ab(:, :)
  real(dp), allocatable, public :: yba(:, :, :)

  integer(i4b), allocatable, public :: inmask(:), outmask(:)
  integer(i4b), allocatable, public :: nohits(:, :)
  real(dp), allocatable, public :: cca(:, :, :), cc(:, :, :)
  real(dp), allocatable, public :: wamap(:, :)
  real(dp), allocatable, public :: map(:, :), binmap(:, :)
  real(sp), allocatable, public :: crit(:)

  real(dp), save, public :: memory_maps = 0, memory_baselines = 0, &
        memory_basis_functions = 0

  character(len=40), parameter :: mstr='(x,a,t32,f9.1," MB")'
  character(len=40), parameter :: mstr3='(x,a,t32,3(f9.1," MB"))'

  public allocate_maps, free_maps, allocate_baselines, free_baselines, free_mask

CONTAINS

  !---------------------------------------------------------------------------


  SUBROUTINE allocate_maps

    integer :: allocstat
    real(sp) :: memsum, mem_min, mem_max

    memory_maps = 0

    allocate(map(nmap, 0:nopix_map-1), cc(nmap, nmap, 0:nopix_map-1), &
         outmask(0:nopix_map-1), stat=allocstat)
    call check_stat(allocstat, 'map, cc and outmask')

    map = 0
    cc = 0
    outmask = 0
    memory_maps = memory_maps + (nmap ** 2. * 8. + nmap * 8. + 4.) * nopix_map

    if (do_binmap .or. nsurvey > 0 .or. ndetset > 0) then
       allocate(binmap(nmap, 0:nopix_map-1), stat=allocstat)
       call check_stat(allocstat, 'binmap')
       binmap = 0
       memory_maps = memory_maps + nmap * nopix_map * 8
    else
       allocate(binmap(nmap, 0:0))
    end if

    if (do_hits) then
       allocate(nohits(0:nopix_map-1, 1), stat=allocstat)
       memory_maps = memory_maps + nopix_map * 4.
       call check_stat(allocstat, 'nohits')
       nohits = 0
    else
       allocate(nohits(0:0, 1))
    end if

    if (do_mask) then
       allocate(crit(0:nopix_map - 1), stat=allocstat)
       call check_stat(allocstat, 'crit')
       crit = 0
       memory_maps = memory_maps + nopix_map * 4.
    else
       allocate(crit(0:0))
    end if

    if (use_inmask) then
       if (.not. allocated(inmask)) then
          allocate(inmask(0:nopix_cross - 1), stat=allocstat)
          call check_stat(allocstat, 'inmask')
          inmask = 0
          memory_maps = memory_maps + nopix_cross * 4.
       end if
    else
       if (.not. allocated(inmask)) allocate(inmask(0:0))
    end if

    if (kfirst) then
       allocate(cca(nmap, nmap, 0:nopix_cross - 1),  &
            wamap(nmap, 0:nopix_cross - 1), stat=allocstat)
       call check_stat(allocstat, 'cca and wamap')
       cca = 0
       wamap = 0
       memory_maps = memory_maps + (nmap ** 2. * 8. + nmap * 8.) * nopix_cross
    else
       allocate(cca(nmap, nmap, 0:0), wamap(nmap, 0:0))
    endif

    call write_memory("Map memory", memory_maps)

  END SUBROUTINE allocate_maps


  !--------------------------------------------------------------------------


  SUBROUTINE free_maps

    if (allocated(map)) deallocate(map)
    if (allocated(nohits)) deallocate(nohits)
    if (allocated(cc)) deallocate(cc)
    if (allocated(cca)) deallocate(cca)
    if (allocated(wamap)) deallocate(wamap)
    if (allocated(outmask)) deallocate(outmask)
    if (allocated(crit)) deallocate(crit)
    if (allocated(binmap)) deallocate(binmap)

  END SUBROUTINE free_maps


  SUBROUTINE free_mask

    if (allocated(inmask))  deallocate(inmask)

  END SUBROUTINE free_mask


  !---------------------------------------------------------------------------


  SUBROUTINE allocate_baselines

    integer :: allocstat
    real(sp) :: memsum, mem_min, mem_max

    memory_baselines = 0

    if (kfirst) then
       allocate(aa(0:basis_order, noba_short, nodetectors), &
            yba(0:basis_order, noba_short, nodetectors), &
            nna(0:basis_order, 0:basis_order, noba_short, nodetectors), &
            nna_inv(0:basis_order, 0:basis_order, noba_short, nodetectors), &
            stat=allocstat)
       call check_stat(allocstat, 'aa, yba, nna and nna_niv')

       aa = 0
       yba = 0
       nna = 0
       nna_inv = 0

       memory_baselines = memory_baselines &
            + nodetectors*noba_short*(basis_order+1)*16
       memory_baselines = memory_baselines &
            + nodetectors*noba_short*(basis_order+1)**2*16
    else
       allocate(yba(1, 1, 1), nna(1, 1, 1, 1))
    endif

    call write_memory("Baseline memory", memory_baselines)

  END SUBROUTINE allocate_baselines


  !------------------------------------------------------------------------------


  SUBROUTINE free_baselines

    integer :: k

    if (allocated(yba)) deallocate(yba)
    if (allocated(aa)) deallocate(aa)
    if (allocated(nna)) deallocate(nna)
    if (allocated(nna_inv)) deallocate(nna_inv)
    if (allocated(basis_functions)) then
       do k = 1, noba_short
          if (.not. basis_functions(k)%copy) then
             deallocate(basis_functions(k)%arr)
          end if
          basis_functions(k)%arr => NULL()
       end do
       deallocate(basis_functions)
       memory_basis_functions = 0
    end if

  END SUBROUTINE free_baselines


  !------------------------------------------------------------------------------

END MODULE maps_and_baselines
