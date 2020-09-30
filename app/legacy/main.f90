program main
  use mod1, only: dp
  implicit none

  interface
    function legacy_rocket(input)
      import dp
      character(len=*), intent(in) :: input
      real(dp), allocatable :: legacy_rocket(:,:)
    end function
  end interface

  integer :: i
  real(dp), allocatable :: output(:,:)

  output = legacy_rocket("rocket.inp")
  do i = 1, size(output, 1)
      print'(11e15.6,1x)', output(i,:)
  end do
end program
