

subroutine cal_pressure (time,u,uk,pk)
  use share_vars
  use rhs
  implicit none
  real(kind=pr),intent(in) :: time
  real(kind=pr),dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: uk, u
  real(kind=pr),dimension(0:nx-1,0:ny-1), intent (out) :: pk
  real(kind=pr),dimension(:,:,:), allocatable :: nlk
  real(kind=pr),dimension(:,:), allocatable :: work1, work2
  integer :: iy

  allocate( nlk(0:nx-1,0:ny-1,1:2) )
  allocate( work1(0:nx-1, 0:ny-1), work2(0:nx-1, 0:ny-1) )

  !-- compute RHS with penalization
  call cal_nlk (time, u, uk, work1, nlk)

  !-- divergence
  call cofdx ( nlk(:,:,1), work1 )
  call cofdy ( nlk(:,:,2), work2 )

  !-- solve poisson eqn
  call poisson ( work1+work2, pk)

  deallocate ( nlk, work1, work2 )
end subroutine cal_pressure
