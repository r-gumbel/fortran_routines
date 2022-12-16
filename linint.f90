program neville
    implicit none

    integer, parameter :: wp = kind(1.d0)
    integer  :: i, N
    real(wp) :: xmin = -1.5_wp*acos(-1.d0), xmax = 1.5*acos(-1.d0), xin
    real(wp) :: a(6), b(6)
    real(wp) :: dx = 0.01
    real(wp), allocatable :: x(:), y(:), z(:)
    real(wp) :: f123456

    a = [0.0_wp, 0.20_wp, 0.40_wp, 0.60_wp, 0.80_wp, 1.0_wp]
    a = xmin + (xmax - xmin) * a

    print *, a

    b = [(f(a(i)), i = 1, 6)]
    N = 100
    allocate(x(N), y(N), z(0:N))


    x = [(i*dx, i = 1, N)]
    x(1:N) =  xmin + (xmax - xmin) * x

!     do i = 1, N
!         print *, x(i)
!     end do

    open(unit = 11, file = 'test.dat')
    do i = 1, N-1
        xin = x(i)
        call nelville(a, b, xin, f123456)
        write(11,*)xin, f123456
    end do
    close(11)

    open(unit = 12, file = "exact.dat")
    do i = 0, N - 1
        z(i) = f(x(i))
        write(12,*) x(i), z(i)
    end do
    close(12)




    contains

    real(wp) function f(x)
    real(wp), intent(in) :: x

    f = sin(x)

    end function f

    subroutine nelville(a, b, x, f123456)

    real(wp), intent(in) :: a(6), b(6)
    real(wp), intent(in) :: x
    real(wp) :: f12, f23, f34, f45, f56
    real(wp) :: f123, f234, f345, f456
    real(wp) :: f1234, f2345, f3456
    real(wp) :: f12345, f23456
    real(wp), intent(out) :: f123456

    f12 = 1._wp / (a(1) - a(2)) * ( (x - a(2))*b(1) - (x - a(1)) * b(2))

    f23 = 1._wp / (a(2) - a(3)) * ( (x - a(3))*b(2) - (x - a(2)) * b(3))

    f34 = 1._wp / (a(3) - a(4)) * ( (x - a(4))*b(3) - (x - a(3)) * b(4))

    f45 = 1._wp / (a(4) - a(5)) * ( (x - a(5))*b(4) - (x - a(4)) * b(5))

    f56 = 1._wp / (a(5) - a(6)) * ( (x - a(6))*b(5) - (x - a(5)) * b(6))



    f123 = 1._wp / (a(1) - a(3)) * ( (x - a(3)) * f12 - (x - a(1)) * f23)

    f234 = 1._wp / (a(2) - a(4)) * ( (x - a(4)) * f23 - (x - a(2)) * f34)

    f345 = 1._wp / (a(3) - a(5)) * ( (x - a(5)) * f34 - (x - a(3)) * f45)

    f456 = 1._wp / (a(4) - a(6)) * ( (x - a(6)) * f45 - (x - a(4)) * f56)


    f1234 = 1._wp / (a(1) - a(4)) * ( (x - a(4)) * f123 - (x - a(1)) * f234)

    f2345 = 1._wp / (a(2) - a(5)) * ( (x - a(5)) * f234 - (x - a(2)) * f345)

    f3456 = 1._wp / (a(3) - a(6)) * ( (x - a(6)) * f345 - (x - a(3)) * f456)


    f12345 = 1._wp / (a(1) - a(5)) * ( (x - a(5)) * f1234 - (x - a(1)) * f2345)

    f23456 = 1._wp / (a(2) - a(6)) * ( (x - a(6)) * f2345 - (x - a(2)) * f3456)


    f123456 = 1._wp / (a(1) - a(6)) * ( (x - a(6)) * f12345 - (x - a(1)) * f23456)


    end subroutine nelville
end program neville
