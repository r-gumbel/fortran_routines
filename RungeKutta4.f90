program runge_kutta4
    implicit none
!=========================================================================
! Routine implements Fourth-Order Runge Kutta Method to solve an initial
! value problem: dy/dx = xy, y(0) = 4.0, on the interval of 0 < x < 1
! Including procedures for the exact solution, y(x) = 4.0*exp(x**2/2)
! and to calculate relative error.
! Results will write to file 'RK4.dat'.
! ========================================================================
    integer, parameter :: wp = kind(1.d0)
    integer  :: i, N
    real(wp) :: k1, k2, k3, k4
    real(wp) :: xmin, xmax, dx
    real(wp), allocatable :: x(:), y(:), exact(:)

    xmin = 0.0_wp; xmax = 1.0_wp        ! interval of interest
    dx   = 0.01_wp                      ! step size
    N    = (xmax - xmin) / dx + 1       ! no. grid points

    allocate(x(N), y(N), exact(N))
    x = xmin + [(i*dx, i = 0, N-1)]
    y = 0.0_wp; y(1) = 4.0_wp
    exact = [(sol(x(i)), i = 1, N)]

    open(unit = 11, file = "RK4.dat")
    write(11,'(f12.9, 2x, f12.9, 2x, f12.9, 2x, f12.9)')x(1), y(1), exact(1), 0.0_wp

! Iterate 4th-order Runge-Kutta
    do i = 1, N - 1
        k1 = f(y(i), x(i))
        k2 = f(y(i) + 1.0_wp/2.0_wp * k1 * dx, x(i) + 1.0_wp/2.0_wp * dx)
        k3 = f(y(i) + 1.0_wp/2.0_wp * k2 * dx, x(i) + 1.0_wp/2.0_wp * dx)
        k4 = f(y(i) + k3 * dx, x(i) + dx)
        y(i+1) = y(i) + dx * (k1 + 2.0_wp * k2 + 2.0_wp * k3 + k4) / 6.0_wp
        write(11,'(f12.9, 2x, f12.9, 2x, f12.9, 2x, f12.9)')x(i+1), y(i+1), exact(i+1), error(exact(i+1),y(i+1))
    end do

    close(11)

    contains

!==========================================================================
!   ODE TO BE SOLVED, f(y,x) = xy
!==========================================================================
    real(wp) function f(y,x)

    real(wp), intent(in) :: x, y

    f = x * y

    end function f

!==========================================================================
! EXACT SOLUTION OF ODE, i.e. y = 4 * exp(1/2 * x^2)
!==========================================================================
    real(wp) function sol(x)

    real(wp), intent(in) :: x

    sol = 4.0_wp * exp(0.5_wp * x**2)

    end function sol

!==========================================================================
! FUNCTION TO CALCULATE RELATIVE ERROR
!==========================================================================
    real(wp) function error(a,b)

    real(wp), intent(in) :: a, b

    error = abs(a - b) / abs(a) * 100.0_wp

    end function error

end program runge_kutta4

