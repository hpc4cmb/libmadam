module linebufmod
implicit none
private

public linebuf, linebuf_add, linebuf_clear, linebuf_set

type linebuf
  character(len=80), dimension(:), pointer :: data=>NULL()
  integer :: size=0
end type linebuf

contains

subroutine linebuf_add (buf, line)
  type(linebuf), intent(inout) :: buf
  character(len=*), intent(in) :: line
  character(len=80), dimension(:), pointer :: data2
  data2=>NULL()

  if (.not. associated(buf%data)) allocate(buf%data(10))

  if (size(buf%data)==buf%size) then
    allocate (data2(2*size(buf%data)))
    data2(1:size(buf%data)) = buf%data
    deallocate (buf%data)
    buf%data => data2
  endif

  buf%size=buf%size+1
  buf%data(buf%size) = trim (line)
end subroutine linebuf_add

subroutine linebuf_clear (buf)
  type(linebuf), intent(inout) :: buf

  if (associated(buf%data)) deallocate(buf%data)
  buf%size=0
end subroutine linebuf_clear

subroutine linebuf_set (buf, orig)
  type(linebuf), intent(inout) :: buf
  character(len=80), dimension(:), intent(in) :: orig

  if (associated(buf%data)) deallocate(buf%data)

  allocate(buf%data(size(orig)))
  buf%data=orig
  buf%size=size(orig)
end subroutine linebuf_set

end module linebufmod
