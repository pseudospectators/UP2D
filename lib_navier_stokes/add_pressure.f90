subroutine add_pressure (nlk)
  use share_vars
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: nlk
  real(kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2, work3
  integer :: iy
  
  ! divergence
  call cofdx ( nlk(:,:,1), work1 )
  call cofdy ( nlk(:,:,2), work2 )
  ! solve poisson eqn
  call poisson ( work1+work2, work3)
  ! gradient
  call cofdx (work3, work1)
  call cofdy (work3, work2)
  ! add gradient
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1) - work1(:,iy)    
    nlk(:,iy,2) = nlk(:,iy,2) - work2(:,iy)
  enddo
  !$omp end parallel do     
end subroutine add_pressure

