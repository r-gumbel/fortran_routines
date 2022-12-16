program primes
    implicit none
!===================================================================================================
! Prime number generation routine.  Written from general algorithm by Samuel Wong
!---------------------------------------------------------------------------------------------------
!    G E N E R A L  A L G O R I T H M
!---------------------------------------------------------------------------------------------------
!   Initialization:
!       (a) Set K as the total number of prime numbers required
!       (b) Set N_BEG as the number prime numbers to start with
!       (c) Input the first N_BEG prime numbers explicitly
!   1. Store the last prime number in the list as NEW
!   2. Construct a new integer by adding 2 to NEW
!   3. Test if NEW is prime by dividing it by all the existing primes <= NEW/2:
!       (a) If divisible, try a new integer by adding 2 to NEW and repeat the test
!       (b) If not divisible, it is a new prime number
!   4. Repeat steps 2 and 3 to find the next prime number in the list.
!   5. Stop if the total number is K
!===================================================================================================
!
    integer :: i,j, K, M, N
    integer, allocatable :: N_BEG(:) ! Initial list of prime numbers
    integer, allocatable :: P(:)     ! Array used to store found prime values
    integer :: NEW
!
    N_BEG = [2,3,5]
    NEW = maxval(N_BEG)
!
    write(*,*) "Enter total number of primes desired, K:"
    read(*,*) K
    if( K .le. 3) then
        print *, "List of primes: 2, 3, 5"; stop
    else
        allocate(P(K))
    end if
    P(1:size(N_BEG)) = N_BEG; P(size(N_BEG)+1:K) = 0

! test NEW
    do j = 4, K
        do
            NEW = NEW + 2
            i = 2
            do while (P(i) .le. NEW/2 .and. P(i) .gt. 0)
                N = P(i)
                if (mod(NEW,N) == 0) exit
                i = i + 1
            end do
            if(mod(NEW,N) /= 0) exit
        end do
        P(j) = NEW
    end do


    do i = 1, K
        print *, P(i)
    end do

end program primes
