program integration
    implicit none
!========================================================================================
! Routine that implements and compares 4 methods of integration: Composite Trapezoidal,
! Composite Midpoint, Simpson's, and Monte Carlo.
! Function f(x) = 4.0_wp * exp(-x)
! xmin = -4.0, xmax = 4.0, step size dx = 0.01
!========================================================================================
    integer, parameter :: wp = kind(1.d0)
    integer  :: N
    real(wp) :: xmin, xmax, dx
    real(wp) :: area1, area2, area3, area4, area5
    real(wp) :: err1, err2, err3, err4

    print *, "Enter value for lower bound xmin: "
    read  *, xmin
    print *, "Enter value for upper bound xmax:"
    read  *, xmax
    dx = 0.01_wp
    N = (xmax - xmin)/dx + 1

    call trap(xmin,dx,area1,N)              ! Calculate integral using trapezoidal method
    call midpoint(xmin,dx,area2,N)          ! Calculate integral using midpoint method
    call simpsons(xmin,xmax,dx,area3,N)     ! Calculate integral using Simpson's rule
    call monte(xmin,xmax,area4)             ! Calculate integral using monte carlo method
    area5 = exact(xmax) - exact(xmin)       ! Calculate integral using expression for its anti-derivative

    call error(area5,area1,err1)            ! Calculate
    call error(area5,area2,err2)            ! relative error
    call error(area5,area3,err3)            ! associated with
    call error(area5,area4,err4)            ! each integration method

    write(*,*)" "
    write(*,'(a, f6.4)')"Step size dx = ", dx
    write(*,'(a, i6)') "Number of points N =", N
    write(*,*)" "
    write(*,'(a)')" The relative error associated with each method is:"
    write(*,*)"==============================================================================="
    write(*,*)"   Trapezoidal     Midpoint        Simpsons        Monte Carlo (10^5 points) "
    write(*,*)"==============================================================================="
    write(*,'(3x,f7.4, x, a, 7x, f7.4, x, a, 7x, f7.4, x, a, 7x, f7.4, x, a)') err1,"%", err2,"%", err3,"%", err4, "%"

contains

!-----------------------------------------------------------------------------
    real(wp) function f(x)
    real(wp), intent(in) :: x
    f = 4.0_wp * exp(-x)
    end function f
!-----------------------------------------------------------------------------
    real(wp) function exact(x)
    real(wp), intent(in) :: x
    exact = -4.0_wp * exp(-x)
    end function exact
!-----------------------------------------------------------------------------
    subroutine trap(xmin,dx,area,N)

    real(wp), intent(in) :: xmin, dx
    real(wp), intent(inout) :: area
    integer, intent(in) :: N
    integer :: i

    area = 0.0_wp
    do i = 0, N-1
        area = area + dx/2.0_wp * (f(xmin + (i+1)*dx) + f(xmin + i*dx))
    end do

    end subroutine trap
!-----------------------------------------------------------------------------
    subroutine midpoint(xmin,dx,area,N)

    real(wp), intent(in) :: xmin, dx
    real(wp), intent(inout) :: area
    integer, intent(in) :: N
    integer :: i

    area = 0.0_wp
    do i = 0, N
        area = area + dx * (f((xmin + i*dx )/2.0_wp + (xmin + (i+1)*dx)/2.0_wp))
    end do

    end subroutine midpoint
!------------------------------------------------------------------------------
    subroutine simpsons(xmin,xmax,dx,area,N)

    real(wp), intent(in) :: xmin, xmax, dx
    real(wp), intent(inout) :: area
    integer, intent(in) :: N
    integer :: i

    area = 0.0_wp
    area = f(xmin) + f(xmax)
    do i = 2, N-1
        if(mod(i,2) == 0) then
            area = area + 4.0_wp * f(xmin + i*dx)
        else
            area = area + 2.0_wp * f(xmin + i*dx)
        end if
    end do
    area = dx/3.0_wp * area

    end subroutine simpsons
!--------------------------------------------------------------------------------
    subroutine monte(xmin,xmax,area)

    real(wp), intent(in) :: xmin, xmax
    real(wp), intent(inout) :: area
    real(wp), allocatable :: x(:), y(:)
    integer  :: i, M = 100000
    real(wp) :: avg

    allocate(x(M), y(M))
    call random_number(x)
    x = xmin + (xmax - xmin) * x
    y = [(f(x(i)), i = 1, M)]
    avg = sum(y)/M
    area = (xmax - xmin) * avg

    end subroutine monte
!--------------------------------------------------------------------------------
    subroutine error(a, b, err)

    real(wp), intent(in)  :: a, b
    real(wp), intent(out) :: err

    err = abs(a - b) / abs(a) * 100.0_wp

    end subroutine error
!--------------------------------------------------------------------------------
end program integration
