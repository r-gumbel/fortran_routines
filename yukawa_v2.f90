program yukawa
implicit none

  integer, parameter :: wp = kind(1.0d0)
  integer, parameter :: nc = 100
  real(wp), parameter :: pi = acos(-1.0d0)
  integer :: i
  real(wp) :: rmax, dr, mu
  real(wp) :: phi_yc(nc)      ! Calculated of Yukawa potential
  real(wp) :: phi_ye(nc)      ! Exact value of Yukawa potential
  real(wp) :: diff(nc), r(nc)

    rmax = 10.0_wp
    dr = rmax/nc
    mu = 2.0_wp

! Loop for radial position r
    do i = 1, nc
      r(i) = i * dr
    end do

  call yukawa_potential(nc, pi, r, rmax, phi_yc, dr, mu)
  call yukawa_exact(nc, pi, r, dr, mu, phi_ye)

! Calculate percentage error between exact potential and approximated
    do i = 1, nc
      diff(i) = abs((phi_ye(i) - phi_yc(i))/phi_ye(i)) * 100
    end do

! Display results in a table:

  write( *, '(a,/,a,/,a,/a)')&
&"=================================================================================================",&
&" Position Value:         Exact Result:            Approximate Result:          Percent Error:    ",&
&"================================================================================================="

    do i = 1, nc
      write(*, '(f10.5,16x,f10.5,16x,f10.5,16x,f10.5)')r(i), phi_ye(i), phi_yc(i), diff(i)
    end do


end program yukawa

!===============================================================================
subroutine yukawa_potential(nc, pi, r, rmax, phi_yc, dr, mu)
implicit none
!-------------------------------------------------------------------------------
! Calculation of the Yukawa potential, in a spherically symmetrical system,
! accomplished by a modified solution to the Helmholtz equation:
! (Laplacian - mu**2)*Phi(r) = -4*pi*Vy*rho(r). In this case we'll set Vy = 1.
!-------------------------------------------------------------------------------

! Internal parameters:
      integer , parameter :: wp = kind(1.0d0)
      integer :: i
      real(wp) :: sum1
      real(wp), intent(in) :: mu
      real(wp) :: bi

! Global parameters:
      integer, intent(in) :: nc
      real(wp), intent(in) :: pi
      real(wp), intent(in) :: rmax, dr


! Internal arrays:
      real(wp), dimension(:), allocatable :: Ui, di, rho
      real(wp), dimension(:), allocatable :: Ei, Fi

! Global arrays:
      real(wp), intent(out) :: phi_yc(nc)
      real(wp), intent(in)  :: r(nc)

! Allocate dimensions to internal arrays
      allocate(Ui(0:nc), di(nc), rho(nc))
      allocate(Ei(0:nc+1), Fi(0:nc+1))

! Initialize
      Ei = 0.0_wp
      Ei(nc+1) = 0.0_wp
      Fi = 0.0_wp
      Fi(nc+1) = 0.0_wp
      bi = -(mu*mu)*(dr*dr) - 2.0_wp
      Ui(0) = 0.0_wp
      sum1 = 0.0_wp

      open(unit = 11, file = "yukawa_calc.dat")


! Loop for density rho, and sum:
      do i = 1 , nc
        rho(i) = 1.0_wp/(8.0_wp * pi) * exp(-r(i))
        sum1 = sum1 + 4.0_wp * pi * r(i)**2 * rho(i) * dr
      end do
      rho = rho/sum1




! Loop for array di:
      do i = 1, nc
        di(i) = -4.0 * pi * rho(i) * r(i) * (dr*dr)
      end do

! Backwards loop for arrays Ei , Fi:
       do i = nc, 1, -1
        Ei(i) = -1.0_wp/(bi + Ei(i+1))
        Fi(i) = (di(i) - Fi(i+1))/(bi + Ei(i+1))
      end do

! Loop to calculate Ui:
      do i = 1, nc
        Ui(i) = (Ei(i) * Ui(i-1)) + Fi(i)
      end do


! Calculate contribution from Yukawa potential, phi_yc, using Ui and r:
      do i = 1, nc
        phi_yc(i) = Ui(i)/r(i)
        !write(11, *)r(i), phi_yc(i)
      end do

      close(11)

end subroutine yukawa_potential

!-------------------------------------------------------------------------------

subroutine yukawa_exact(nc, pi, r, dr, mu, phi_ye)
implicit none

! Subroutine to calculate exact value of the Yukawa potential, using the
! solution to the Helmholtz equation.

! Local parameters/variables:
      integer, parameter :: wp = kind(1.0d0)
      integer :: i

! Global parameters/variables:
      integer, intent(in) :: nc
      real(wp), intent(in) :: pi
      real(wp), intent(in) :: dr, mu
      real(wp), intent(out) :: phi_ye(nc)
      real(wp) :: p(nc)
      real(wp), intent(in) :: r(nc)


      p = 0.0_wp
      phi_ye = 0.0_wp
      open(unit = 12, file = "yukawa_exact.dat")

! Loop to calculate exact Yukawa potential
      do i = 1, nc
        p(i) = ((1.0_wp/(1.0_wp - mu**2))**2)*(exp(-mu*r(i)) - (exp(-r(i))*(1.0_wp + ((1.0_wp/2.0_wp)*(1.0_wp - mu**2))*r(i))))
        phi_ye(i) = p(i)/(r(i))
        !write(12,*)r(i), phi_ye(i)
end do

    close(12)

end subroutine yukawa_exact
