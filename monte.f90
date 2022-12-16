program monte
!==============================================================================================================================
! Monte Carlo integration routine for a single variable function
!==============================================================================================================================
    implicit none

    integer, parameter :: wp = kind(1.d0)
    integer  :: i, M
    real(wp) :: xmin, xmax
    real(wp), allocatable :: x(:), y(:)
    real(wp) :: avg, area

    M = 10000
    xmin = -2.0_wp; xmax = 2.0_wp
    allocate(x(M), y(M))
    call random_number(x)
    x = xmin + (xmax - xmin) * x
    y = [(f(x(i)), i = 1, M)]
    avg = sum(y)/M
    area = (xmax - xmin)*avg

    print *, area

    contains

    real(wp) function f(x)

    real(wp), intent(in) :: x

    f = exp(-x**2) * x**2

    end function f

end program monte
