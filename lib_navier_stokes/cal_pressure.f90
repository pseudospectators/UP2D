subroutine cal_pressure (time,u,uk,pk)
  use share_vars
  use rhs
  implicit none
  real(kind=pr), intent(in) :: time
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: uk, u
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (out) :: pk
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2) :: nlk
  real(kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2
  integer :: iy
  
  call cal_nlk (time, u, uk, work1, nlk, .true.)
  
  ! divergence
  call cofdx ( nlk(:,:,1), work1 )
  call cofdy ( nlk(:,:,2), work2 )
  ! solve poisson eqn
  call poisson ( work1+work2, pk)
     
end subroutine cal_pressure

