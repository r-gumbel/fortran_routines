program midpoint
    implicit none
!=====================================================================================
! Integration of a single variable function using composite midpoint method
! Test function in this case will be a simple quadratic (x-2)**2
! User supplies inputs. Four different values of step dx will print to the stdout
! dx = 1, dx = 0.1, dx = 0.01, dx = 0.001-- Along with relative error associated
! each
!=====================================================================================
    integer, parameter :: wp = kind(1.d0)
    integer  :: j, N
    real(wp) :: xmin, xmax, dx, area, error

    xmin = 0.0_wp; xmax = 10.0_wp

    print *, "Enter lower bound xmin:"
    read  *, xmin
    print *, "Enter upper bound xmax:"
    read  *, xmax
    print *, "=============================================================================="
    print *, "         dx              N         area          %error"
    print *, "=============================================================================="

    do j = 0, 3
        dx = 10.0_wp**(-j)
        call mid(xmin, xmax, dx, area, N)
        call exact(xmin, xmax, area, error)
        write(*,'(5x, f9.4, 5x, i9, 5x, f9.4, 5x, f9.4)') dx, N, area, error
    end do
    print *, "==============================================================================="

    contains
!----------------------------------------------------------------------------------------
!   FUNCTION TO BE INTEGRATED,  i.e. f(x) = (x-2)^2
!----------------------------------------------------------------------------------------
    real(wp) function f(x)
    real(wp), intent(in) :: x
    f = (x - 2.0_wp)**2
    end function f

!----------------------------------------------------------------------------------------
!   EXACT VALUE OF ANTI-DERITATIVE,  i.e. F(x) = 2x - 4
!----------------------------------------------------------------------------------------
    real(wp) function antiF(x)
    real(wp), intent(in) :: x
    antiF = 1.0_wp/3.0_wp * x**3 - 2.0_wp * x**2 + 4.0_wp * x
    end function antiF

!----------------------------------------------------------------------------------------
!   SUBROUTINE TO CALCULATE INTEGRAL USING COMPOSITE MIDPOINT METHOD
!----------------------------------------------------------------------------------------
    subroutine mid(xmin, xmax, dx, area, N)

    real(wp), intent(in)  :: xmin, xmax, dx
    real(wp), intent(out) :: area
    integer, intent(out)  :: N
    integer :: i

    N = (xmax - xmin) / dx
    area = 0.0_wp

!   Iterate over midpoint sum
    do i = 0, N - 1
        area = area + dx * f((xmin + (i-1) *dx + i*dx)/2.0_wp)
    end do

    end subroutine mid

!----------------------------------------------------------------------------------------
!   Subroutine to calculate error at given step-size dx
!----------------------------------------------------------------------------------------
    subroutine exact(xmin, xmax, area, error)

    real(wp), intent(in)  :: xmin, xmax, area
    real(wp), intent(out) :: error

    error = antiF(xmax) - antiF(xmin)        ! Calculate relative error of trapezoidal
    error = abs(error - area)/ abs(error)    ! method for a given step-size dx
    error = error * 100.0_wp
    end subroutine exact

end program midpoint

