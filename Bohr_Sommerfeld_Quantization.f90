!===============================================================================
! program N_search
! Program to find number of bound states N using the Bohr-Sommerfeld classical
! approximation.  User will be prompted to enter values for an initial guess
! for E_last, and a value for gamma. This program will contain and call
! subroutines E_search & int_roots.
!===============================================================================
    program N_search
    implicit none
!-------------------------------------------------------------------------------
! Parameters:
    integer , parameter :: wp = kind(1.d0) !machine double default
    real , parameter  ::  PI = 3.1415927
!-------------------------------------------------------------------------------
! Input/Output variables:
    real(8)  :: Nmax
    integer  :: MaxN
    integer  :: N             !current level (input)
    real(wp) :: E             !trial energies (input)
    real(wp) :: gamma         !=sqrt(2*m*l**2*potential/hbar**2) (input)
    real(wp) :: answer_s      !action integral S(E)
    real(wp) :: xmin, xmax
!-------------------------------------------------------------------------------
! Prompt use for the following values: E, N, and gamma
!-------------------------------------------------------------------------------
    print *, '#Enter initial guess for E_last :'
    read *, E
    print *, '#Enter value for gamma:' !(note: you've been using 21.7)
    read *, gamma
!-------------------------------------------------------------------------------
!  Use E_last to find number of bound states
!-------------------------------------------------------------------------------
   call int_roots(E,xmin,xmax,answer_s,gamma)
   Nmax = answer_s/PI - 1.0/2.0
   MaxN = NINT(Nmax)
   print *, '#For gamma=', gamma, 'Nmax =', MaxN + 1
   print *,'==================================================================='
!-------------------------------------------------------------------------------
    E = -0.9999
    do N = 0, MaxN
      call E_search(E,N,gamma)
    end do
    end program N_search
!-------------------------------------------------------------------------------
    subroutine E_search(E,N,gamma)
    implicit none
!===============================================================================
! Subroutine Energy_Search(E,N,gamma)
! Finds the Bohr-Sommerfeld approximation quantization for a given energy bound
! state n.  Uses a step search of the action integral S(E), to find where
! S(E) = (n + 1/2)*pi, using an initial guess for E.
! This entails nesting the subroutine int_roots within Energy_Search.
! int_roots will calculate both turning points x_in and x_out, as well as
! output a new value for S(E + dE), as Energy_Search incrementally adjusts
! the value of E (until S(E +dE) = (n +1/2)*pi
!===============================================================================
! Establish parameters
!-------------------------------------------------------------------------------
    integer , parameter  ::  wp = kind(1.d0)   !machine default double
    real , parameter  ::  PI = 3.1415927
    real(wp)          ::  Etol           !convergence criteria for energy search
    integer           ::  Itmax          !max number of iterations
!-------------------------------------------------------------------------------
! Input/Output variables
    integer , intent(in)  :: N           !current level (input)
    real(wp) , intent(inout) :: E        !trial energies
    real(wp) , intent(in) :: gamma       !varies according to molecule
    real(wp) :: energy_n                 !result of energy search
    real(wp) :: func                     !potential function v(x)
    real(wp) :: answer_s                 !action (output)
    real(wp) :: xmin, xmax               !turning points (outputs)
    integer  :: iter                     !iteration of search
!-------------------------------------------------------------------------------
! Local variables:
    real(wp) :: dE                       !increment in energy search
!-------------------------------------------------------------------------------
! Prompt user for values E and N
!-------------------------------------------------------------------------------
    Etol  = 1.0d-5
    iter  = 0
    Itmax = 100
    iter  = 0
    dE    = abs(E)/4.0_wp

  do while (abs(dE) >= Etol .and. iter < Itmax)

    iter  = iter + 1
    E     = E + dE
    call int_roots(E,xmin,xmax,answer_s,gamma)
    if(answer_S .ge. (N + 1.0/2.0)*PI .or. E > 0.0) then
      E  = E - dE
      dE = dE/2.0
    end if

  end do
!
  energy_n = E
  print *, '#Bound state n =', N
  print *, '#energy at n =', energy_n
  print *, '#x in =', xmin, 'x out =', xmax

!
    end subroutine E_search
!===============================================================================
    subroutine int_roots(E,xmin,xmax,answer_s,gamma)
    implicit none
!-------------------------------------------------------------------------------
!   Subroutine to calculate the integral S(E) and the classical
!   turning points xmin, xmax for a given energy E
!-------------------------------------------------------------------------------
!   Establish the parameters we'll use
!-------------------------------------------------------------------------------
!   Shared parameters:
  integer, parameter :: wp = kind(1.d0) !machine defaut double
!
!   For calculating xmin and xmax:
  real(wp)  :: x_initial                !x value at equilibrium
  real(wp)  :: Fold                     !E - func(x_initial)
  real(wp)  :: xtol                     !tolerance for turning points search
  real(wp)  :: iter                     !iterations of step search
  real(wp)  :: Itmax                    !max number of iterations allowed
!
!    For calculating the integral S:
  real(wp) , intent(in) :: gamma                    !=sqrt(2*m*l**2*potential/hbar**2)
  integer   :: N                        !Number of integration points
!-------------------------------------------------------------------------------
!   Establish the variables we'll use:
!-------------------------------------------------------------------------------
!   Input/Output variables:
  real(wp) , intent(in)  :: E            !energy (input)
  real(wp) , intent(out) :: answer_S     !action (output)
  real(wp) , intent(out) :: xmin, xmax   !turning points (output)
  real(wp)               :: func         !potential as function of x
!
!   Shared variables:
  real(wp)  :: x                         !general independent variable
  real(wp)  :: dx                        !incremental change of x
!
!   Integration variables:
  integer   :: NP
  integer   :: i
!
!   Assign parameter values:
!   Roots:
  x_initial = 2.0_wp**(1.0_wp/6.0_wp)
  xtol      = 1.0d-9
  Itmax     = 100
!   Integral:
  !gamma     = 21.7_wp
!
!===============================================================================
! find x1 (inner turning point), begin step search at bottom of potential well.
!===============================================================================
! Setup initializations
!-------------------------------------------------------------------------------
  x = x_initial                           !equilibrium value for x
  dx = 0.5_wp                             !initial step size
  Fold = func(x_initial)                  !value of pot(x) at equilibrium
  iter = 0                                !starting point for iterations
!-------------------------------------------------------------------------------
! Start the search and come back if not converged
!-------------------------------------------------------------------------------
!
  do while (abs(dx) >= xtol .and. iter < Itmax)
!
    iter = iter + 1
    x    = x - dx
    if(func(x) .ge. E) then
      x  = x + dx
      dx = dx/2.0
    end if
!
  end do
!
  xmin = x
!
!===============================================================================
! find x2 (outer turning point), begin step search a bottom of potential well.
!===============================================================================
! Setup initializations
!-------------------------------------------------------------------------------
  x = x_initial                         !equilibrium value for x
  dx = 0.5_wp                           !initial step size
  Fold = func(x_initial)                !value of pot(x) at equilibrium
  iter = 0                              !starting point for iterations
!-------------------------------------------------------------------------------
! Start the search and come back if not converged
!-------------------------------------------------------------------------------
!
do while (abs(dx) >= xtol .and. iter < Itmax)
!
  iter = iter + 1
  x    = x + dx
  if(func(x) .ge. E) then
    x  = x - dx
    dx = dx/2.0
  end if
!
end do
!
  xmax = x
!===============================================================================
! Use Simpson's rule to calculate integral
! S(E) = integral of (gamma*(E -v(x)) using bounds xmin and xout
! from roots search of function E - v(x) for a provided energy level E
!   Note: current code prompts user to imput number of integration
!   points, but this is a feature I've kept from Dr. Umar's example
!   code, to experiment with different values for 'N'
!===============================================================================
!-------------------------------------------------------------------------------
! Read in the number of integration points and perform checking of input
!-------------------------------------------------------------------------------
N = 200
!
!-------------------------------------------------------------------------------
! If number of points even make it odd
!-------------------------------------------------------------------------------
    if(mod(N,2) == 0) then
      NP = N + 1
    else
      NP = N
    end if
!
    dx  = (xmax-xmin)/(NP-1)
!
    answer_s = 0.0
    do i = 2, NP-1, 2
      x = xmin + (i-1)*dx
      answer_s = answer_s + 4.0*sqrt(E - func(x))
    end do
    do i = 3, NP-2, 2
      x = xmin + (i-1)*dx
      answer_s = answer_s + 2.0*sqrt(E - func(x))
    end do
!
answer_s = dx/3.0*(answer_s + (sqrt(E - func(xmin)) + sqrt( E-func(xmax))))
!
answer_s = gamma*answer_s         !multiply integrated function by gamma
!
end subroutine int_roots
!-------------------------------------------------------------------------------
!       Function func(x)
!  evaluates the Lennard-Jones potential at x
!-------------------------------------------------------------------------------
function func(x)
implicit none
!-------------------------------------------------------------------------------
  integer, parameter  :: wp = kind(1.d0)  ! machine default double
  real(wp),intent(in) :: x
  real(wp)            :: func
!
  func = 4.0_wp*(x**(-12.0_wp) - x**(-6.0_wp))
!
end function func
