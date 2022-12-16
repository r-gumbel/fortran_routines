  program poisson
  implicit none

! Define parameters:
    integer, parameter  :: wp = kind(1.0d0)
    integer, parameter  :: nc = 100
    integer :: i
    real(wp), parameter :: pi = acos(-1.0d0)

! Arrays:
    real(wp) :: F(0:nc+1), E(0:nc+1), Wc(0:nc+1)
    real(wp) :: di(nc), r(nc), rho(nc)
    real(wp) :: phi_a(nc), phi(nc), r_exact(nc), exact_phi(nc), diff(nc)

  call loops(nc, Wc, F, E, di, r, rho, phi_a)

  call phi_exact(nc, r_exact, phi, exact_phi)

  do i = 1, nc
     diff(i) =  abs((exact_phi(i) - phi_a(i)) / exact_phi(i)) *100
  end do


! Display results in a table:

write(*,'(a,/,a,/,a,/a)')&
   &"=================================================================================================",&
   &" Position Value:         Exact Result:            Approximate Result:          Percent Error:    ",&
   &"================================================================================================="

      do i = 1, nc
   write (*,'(f10.5,16x,f10.5,16x,f10.5,16x,f10.5)') r(i), exact_phi(i), phi_a(i), diff(i)
      end do

  end program poisson

!======================================================================================================
  subroutine loops(nc, Wc, F, E, di, r, rho, phi_a)
  implicit none

! Declare parameters:
    integer, parameter :: wp = kind(1.d0)
    integer, intent(in) :: nc
    integer :: i
    real(wp) :: pi = acos(-1.0d0), dr, rmax, sum

! Declare Arrays:
    real(wp) :: F(0:nc+1), E(0:nc+1), Wc(0:nc+1)
    real(wp), intent(inout) :: di(nc), r(nc), rho(nc), phi_a(nc)

! Initialize
    E = 0.0_wp
    r = 0.0_wp
    F = 0.0_wp
    Wc = 0.0_wp
    Wc(nc+1) = 1.0_wp
    rmax = 10.0_wp
    dr = rmax/nc
    phi_a = 0.0_wp

  open(unit = 11, file = 'Poisson.dat')

! Loop for position r, density rho, and sum:
    do i = 1,nc
      r(i) = dr * i
      rho(i) = 1.0_wp/(8.0_wp*pi) * exp(-r(i))
      sum = sum + 4.0_wp * pi * r(i)**2 * rho(i) * dr
    end do
    rho = rho/sum

! Loop for arrays, E, F, and di:
    do i = 1, nc
      di(i) = -4.0_wp * pi * rho(i) * r(i) * dr**2
      E(i) = -1.0_wp / (E(i-1) -2.0_wp)
      F(i) = (di(i) - F(i-1)) / (E(i-1) -2.0_wp)
    end do

! Backwards loop to fill array Wc using arrays E and F:
    do i = nc, 1, -1
      Wc(i) = E(i) * Wc(i+1) + F(i)
    end do

! Loop to calculate phi_c using Wc and r:
    do i = 1, nc
      phi_a(i) = Wc(i) / r(i)
      write(11,*)r(i),phi_a(i)
    end do

  close(11)

end subroutine loops

subroutine phi_exact(nc, r_exact, phi, exact_phi)
implicit none

! Declare parameters and arrays
    integer, parameter :: wp = kind(1.0d0)
    integer, intent(in) :: nc
    real(wp) :: rmax, dr
    integer :: i
    real(wp), intent(inout) :: r_exact(nc), phi(nc), exact_phi(nc)

! Initialize:
    rmax = 10.0_wp
    dr = rmax / nc
    r_exact = 0.0_wp
    exact_phi = 0.0_wp
    phi = 0.0

  open(unit = 12, file = "Exact_potential.dat")

    do i = 1, nc
      r_exact(i) = dr * i
      phi(i) = 1.0_wp - 1.0_wp/2.0_wp * (r_exact(i) + 2.0_wp) * exp(-1.0_wp*r_exact(i))
      exact_phi(i) = phi(i) / r_exact(i)
      write(12,*)r_exact(i), exact_phi(i)
    end do

  close(12)

end subroutine phi_exact
