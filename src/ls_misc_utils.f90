module ls_misc_utils
use planck_config
implicit none
private

public fatal_error, assert, assert_present, assert_not_present, &
       assert_alloc, file_present, toString

interface toString
  module procedure toString_i4b, toString_i8b, toString_sp, toString_dp, &
                   toString_l
end interface

contains

subroutine fatal_error (msg)
  character(len=*), intent(in), optional :: msg

  if (present(msg)) then
    print *,'Fatal error: ', trim(msg)
  else
    print *,'Fatal error'
  endif
  call exit_with_status(1)
end subroutine fatal_error

function file_present (filename)
  character(len=*), intent(in) :: filename
  logical :: file_present

  inquire(file=trim(filename),exist=file_present)
end function

subroutine assert_present (filename)
  character(len=*), intent(in) :: filename

  if (.not. file_present(trim(filename))) then
    print *, "Error:  file '" // trim(filename) // "' does not exist!"
    call exit_with_status(1)
  end if
end subroutine assert_present

subroutine assert_not_present (filename)
  character(len=*), intent(in) :: filename

  if (file_present(trim(filename))) then
    print *, "Error:  file '" // trim(filename) // "' already exists!"
    call exit_with_status(1)
  end if
end subroutine assert_not_present

subroutine assert_alloc (stat,code,arr)
  integer, intent(in) :: stat
  character(len=*), intent(in) :: code, arr

  if (stat==0) return

  print *, trim(code)//'> cannot allocate memory for array: '//trim(arr)
  call exit_with_status(1)
end subroutine assert_alloc

subroutine assert (testval,msg,errcode)
  logical, intent(in) :: testval
  character(len=*), intent(in), optional :: msg
  integer, intent(in), optional :: errcode

  if (testval) return

  print *,"Assertion failed: "
  if (present(msg)) print *, trim(msg)
  if (present(errcode)) call exit_with_status (errcode)
  call exit_with_status(1)
end subroutine

function toString_i4b (num)
  integer(i4b), intent(in) :: num
  character (len=30) :: toString_i4b

  write(toString_i4b,*) num
  toString_i4b = trim(adjustl(toString_i4b))
end function

function toString_i8b (num)
  integer(i8b), intent(in) :: num
  character (len=30) :: toString_i8b

  write(toString_i8b,*) num
  toString_i8b = trim(adjustl(toString_i8b))
end function

function toString_sp (num)
  real(sp), intent(in) :: num
  character (len=30) :: toString_sp

  write(toString_sp,*) num
  toString_sp = trim(adjustl(toString_sp))
end function

function toString_dp (num)
  real(dp), intent(in) :: num
  character (len=30) :: toString_dp

  write(toString_dp,*) num
  toString_dp = trim(adjustl(toString_dp))
end function

function toString_l (val)
  logical, intent(in) :: val
  character (len=10) :: toString_l

  write(toString_l,*) val
  toString_l = trim(adjustl(toString_l))
end function

end module ls_misc_utils
