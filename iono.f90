subroutine nOut(n0, T, z, H, n)

  implicit none
  real, intent(in) :: n0
  real, dimension(1000), intent(in) :: T, z, H
  real, intent(out), dimension(1000) :: n

  integer :: i
  !----------------------------------------------------------------
  n(1) = n0
  do i=2, 1000
     n(i) = (n(i-1)*T(i-1))/(T(i))*exp(-(z(i)-z(i-1))/H(i))
  end do

end subroutine nOut

subroutine col_nOut(n, H, z, col_n)

  implicit none
  real, dimension(1000), intent(in) :: n, H, z
  real, intent(out), dimension(1000) :: col_n

  integer :: i
  !----------------------------------------------------------------
  col_n(1000) = n(1000)*H(1000)
  !print*,'z 1000: ', z(1000)
  do i=999, 1, -1
     col_n(i) = n(i)*(z(i+1) - z(i)) + col_n(i+1)
  end do

end subroutine col_nOut
