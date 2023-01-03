program simpson
    implicit none
!=====================================================================================
! Integration of a single variable function using composite simpson's rule
! Test function in this case will be a simple quadratic (x-2)**2
! User supplies inputs. Four different values of step dx will print to the stdout
! dx = 1, dx = 0.1, dx = 0.01, dx = 0.001-- Along with relative error associated
! each
!=====================================================================================
    integer, parameter :: wp = kind(1.d0)
    integer  :: j, N
    real(wp) :: xmin, xmax, area, dx, error

    print *, "Enter lower bound xmin:"
    read  *,  xmin
    print *, "Enter upper bound xmax:"
    read *,   xmax

    print *, "=============================================================================="
    print *, "         dx              N         area          %error"
    print *, "=============================================================================="

    do j = 0, 3
        dx = 10.0_wp**(-j)
        N  = (xmax - xmin) / dx
        call simp(xmin, xmax, dx, area, N)
        call exact(xmin, xmax, area, error)
        write(*, '(5x, f9.4, 5x, i9, 5x, f9.4, 5x, f9.4)') dx, N, area, error
    end do
    print *, "==============================================================================="

    contains

!-------------------------------------------------------------------------------------
!   FUNCTION TO BE INTEGRATED, i.e. f(x) = (x-2)^2
!-------------------------------------------------------------------------------------
    real(wp) function f(x)
    real(wp), intent(in) :: x
    f = (x - 2.0_wp)**2
    end function f

!-------------------------------------------------------------------------------------
!   ANTI-DERIVATIVE OF f(x), TO BE USED IN CALCULATING RELATIVE ERROR
!-------------------------------------------------------------------------------------
    real(wp) function antiF(x)
    real(wp), intent(in) :: x
    antiF = 1.0_wp/3.0_wp * x**3 - 2.0_wp * x**2 + 4.0 * x
    end function antiF

!-------------------------------------------------------------------------------------
!   SUBROUTINE THAT USES SIMPSONS RULE TO CALCULATE DEFINITE INTEGRAL
!   FOR A SPECIFIED STEP SIZE (dx), NUMBER OF POINTS (N), UPPER AND
!   LOWER BOUNDS (xmin and xmax)
!-------------------------------------------------------------------------------------
    subroutine simp(xmin, xmax, dx, area, N)

    integer, intent(in)   :: N
    real(wp), intent(in)  :: xmin, xmax, dx
    real(wp), intent(out) :: area
    integer :: i

    area = f(xmin) + f(xmax)
    do i = 2, N-2
        if(mod(i,2) == 0) then
            area = area + 4.0_wp * (f(xmin + i * dx))
        else
            area = area + 2.0_wp * (f(xmin + i * dx))
        end if
    end do
    area = dx/3.0_wp * area

    end subroutine simp

!-------------------------------------------------------------------------------------
!   Subroutine to calculate error at specified step-size dx
!-------------------------------------------------------------------------------------
    subroutine exact(xmin, xmax, area, error)

    real(wp), intent(in)  :: xmin, xmax, area
    real(wp), intent(out) :: error

    error = antiF(xmax) - antiF(xmin)           ! Calculate relative error of simpsons
    error = abs(error - area) / abs(error)      ! rule method at step size dx

    end subroutine exact

end program simpson


