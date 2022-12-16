    program hermites
    implicit NONE

    integer, parameter :: wp = kind(1.0d0)
    integer, parameter :: nr = 50
    integer :: i, j, k, n, max_n, max_l, fact
    real(wp), parameter :: pi = acos(-1.d0)
    real(wp) :: b, rmax, dr, beta
    real(wp) :: X(nr), r(nr), Y(nr), C(nr), F(nr)
    real(wp) :: rho(nr)

! User input
    print *, "#Enter values for l and n:"
    read *, max_l, n

! Initialize
    b = 1.010_wp*(16.0_wp**(1.0_wp/6.0_wp))
    rmax = 10.0_wp
    dr = rmax/nr
    n = 0

    r = 0.0_wp

    do i = 1, nr   ! Populate r(nr) array
    r(i) = i*dr
    end do

    C = 0.0_wp


    do j = 0, max_l

     do i = 1, nr
      call chi(r, n, nr, i, j, b, X, C)
      C(i) = X(i)
      call harmonic(n, nr, i, j, r, Y, F)
      F(i) = Y(i)
     end do
     !rho =
    end do


end program hermites


subroutine chi(r, n, nr, i, j, b, X, C)
implicit none

integer, parameter :: wp = kind(1.0d0)
integer, intent(in) :: nr
integer, intent(in) :: n
integer, intent(in) :: i, j
integer :: Fact
real(wp), intent(in) :: r(nr)
real(wp), intent(in) :: b
real(wp), intent(out) :: X(nr)
real(wp), intent(inout) :: C(nr)
real(wp) :: beta

    beta = gamma(n + j + 3.0_wp/2.0_wp)
    X(i) = C(i) + (b**(-3/2)) * (((2.0_wp*fact(n))/beta)**(1.0_wp/2.0_wp)) * ((r(i)/b)**j) * exp((-1.0_wp/2.0_wp) * ((r(i)/b)**2))
end subroutine chi

subroutine harmonic(n, nr, i, j, r, Y, F)
implicit none

integer, parameter :: wp = kind(1.0d0)
integer, intent(in) :: nr
integer, intent(in) :: n
integer, intent(in) :: i, j
real(wp), parameter :: pi = acos(-1.d0)
real(wp), intent(in) :: r(nr)
real(wp), intent(out) :: Y(nr)
real(wp), intent(inout) :: F(nr)

  Y(i) = F(i) + (2.0_wp * ((2.0_wp * j) + 1.0_wp)) / (4.0_wp * pi * r(i)*r(i))

end subroutine harmonic

FUNCTION Fact(n)
IMPLICIT NONE
INTEGER :: Fact
INTEGER :: i
INTEGER, INTENT(IN) :: n

IF (n == 0) THEN
   Fact = 1
ELSE
   Fact = 1
   Do i = 1, n
      fact = i * fact
   End Do
END IF

END FUNCTION Fact
