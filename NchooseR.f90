program combination
    implicit none

    integer :: N, R
    integer :: combo

    write(*,*) "Enter total number of objects to choose from:"
    read (*,*) N
    write(*,*) "Enter total number of choices:"
    read (*,*) R
    call NchooseR(N,R,combo)
    write(*,*) "Total number of possible choices = ", combo


contains

    subroutine NchooseR(n,r,combo)

    integer, intent(in)  :: n, r
    integer, intent(out) :: combo

    combo = fact(n) / ((fact(r)) * (fact(n - r)))

    end subroutine NchooseR

    integer function fact(a)
    integer, intent(in) :: a
    integer :: b

    fact = 1
    b = a

    do while (b .ge. 1)
        fact = b * fact
        b = b - 1
    end do

    end function
end program combination
