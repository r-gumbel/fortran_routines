program test
    implicit none
!-------------------------------------------------------------------------
! Driver for tri-linear interpolation procedure,
! tplina(nx, ny, nz, xar, yar, zar, f, xcu, ycu, zcu).
! Using a three variable test function, f(x,y,z) = x**2 + y**2 + z**2.
!-------------------------------------------------------------------------

    integer, parameter :: wp = kind(1.d0)
    integer  :: i, j, k
    integer  :: nx, ny, nz, nr
    real(wp) :: xcu, ycu, zcu
    real(wp) :: dx, dy, dz, dr
    real(wp) :: answer, error
    real(wp), allocatable :: xar(:), yar(:), zar(:)
    real(wp), allocatable :: f(:,:,:)

    nx = 100; ny = 100; nz = 100
    nr = 60
    dx = 0.1_wp; dy = 0.1_wp; dz = 0.1_wp  ! step size for x, y, z meshes

    allocate(xar(nx))                      ! x-mesh      Inputs to
    allocate(yar(ny))                      ! y-mesh      generates points
    allocate(zar(nz))                      ! z-mesh      for interpolation

    allocate(f(nx,ny,nz))                  ! Known function values array for interpolation


    xar = [(i*dx, i = 1, nx)]
    yar = [(i*dy, i = 1, ny)]
    zar = [(i*dz, i = 1, nz)]

    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
            f(i,j,k) = g(xar(i), yar(j), zar(k))
        end do
      end do
    end do

    xcu = 4.5_wp
    ycu = 4.33_wp
    zcu = 5.32_wp

   !====================================================================================!
   ! T E S T : Using a do-loop to range through 10 desired values of x, while keeping   !
   ! values for y and z fixed.  Comparing to exact value of a function, and calculating !
   ! the absolute error.                                                                !
   !====================================================================================!

    write(*,*)'=============================================================================='
    write(*,'(t10, a1, t20, a17, t37, a16, t68, a5)') 'x', 'f(x,4.33,5.32)', 'interp', 'error'
    write(*,*)'=============================================================================='
    do i = 1, 10
      xcu = i*0.33_wp
      answer = tplina(nx, ny, nz, xar, yar, zar, f, xcu, ycu, zcu)
      error  = abs(g(xcu,ycu,zcu) - answer) / abs(g(xcu,ycu,zcu)) * 100._wp
      write(*,'(1x, f14.8, 6x, f14.8, 6x, f14.8, 6x, e15.4)')xcu, g(xcu,ycu,zcu), answer, error
    end do

contains

Function tplina(nx, ny, nz, xar, yar, zar, f, xcu, ycu, zcu)
!---------------------------------------------------------------------------------
!       to do a trilinear inter. of f(nx,ny,nz), on xar(nx), yar(ny), and zar(nz)
!---------------------------------------------------------------------------------
      Implicit None
      Integer, Parameter :: wp = Kind(1.0D0)
!-------------------------------------------------------
!   A r g u m e n t s
!-------------------------------------------------------
      Integer, Intent(In)  :: nx
      Integer, Intent(In)  :: ny
      Integer, Intent(In)  :: nz
      Real(wp), Intent(In) :: xcu
      Real(wp), Intent(In) :: ycu
      Real(wp), Intent(In) :: zcu
      Real(wp), Intent(In) :: xar(nx)
      Real(wp), Intent(In) :: yar(ny)
      Real(wp), Intent(In) :: zar(nz)
      Real(wp), Intent(In) :: f(nx, ny, nz)
!-------------------------------------------------------
!   L o c a l   V a r i a b l e s
!-------------------------------------------------------
      Integer :: nxm, nym, nzm, i, icu, jcu, kcu
      Real(wp) :: tplina, dxf, dyf, dzf, dxc, dyc, dzc
      Real(wp) :: f111, f112, f121, f122
      Real(wp) :: f211, f212, f221, f222
!-------------------------------------------------------
      tplina = 0.0_wp
      nxm = nx - 1
      nym = ny - 1
      nzm = nz - 1
      Do i = 1, nxm
         If(xar(i) <= xcu .And. xcu <= xar(i+1)) Go To 20
      End Do
      Return
20    Continue
      icu = i
      Do i = 1, nym
         If(yar(i) <= ycu .And. ycu <= yar(i+1)) Go To 40
      End Do
      Return
40    Continue
      jcu = i
      Do i = 1, nzm
         If(zar(i) <= zcu .And. zcu <= zar(i+1)) Go To 60
      End Do
      Return
60    Continue
      kcu = i

      dxf = (xcu-xar(icu)) / (xar(icu+1)-xar(icu))
      dxf = (ycu-yar(jcu)) / (yar(jcu+1)-yar(jcu))
      dzf = (zcu-zar(kcu)) / (zar(kcu+1)-zar(kcu))

      dxc = 1.0_wp - dxf
      dyc = 1.0_wp - dyf
      dzc = 1.0_wp - dzf

      f111 = f(icu, jcu, kcu)
      f112 = f(icu, jcu, kcu + 1)
      f121 = f(icu, jcu + 1, kcu)
      f122 = f(icu, jcu + 1, kcu + 1)
      f211 = f(icu + 1, jcu, kcu)
      f212 = f(icu + 1, jcu, kcu + 1)
      f221 = f(icu + 1, jcu + 1, kcu)
      f222 = f(icu + 1, jcu + 1, kcu + 1)

      tplina = f111*dxc*dyc*dzc + f112*dxc*dyc*dzf + f121*dxc*dyf*dzc + f122*dxc*dyf*dzf &
      + f211*dxf*dyc*dzc + f212*dxf*dyc*dzf + f221*dxf*dyf*dzc + f222*dxf*dyf*dzf
!
      Return
End Function tplina

      real(wp) function g(x,y,z)

      real(wp), intent(in) :: x, y, z

      g  = x**2 + y**2 + z**2

      end function g

end program test
