subroutine add_diffusion (nlk, uk)
  ! adds the explicit diffusion term nu*laplace(u) to nlk
  use share_vars
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (in) :: uk
  real(kind=pr), dimension(:,:), allocatable :: work1, work2
  integer :: iy

  allocate( work1(0:nx-1, 0:ny-1), work2(0:nx-1, 0:ny-1) )

  ! first component
  call cofdxdx( uk(:,:,1), work1 )
  call cofdydy( uk(:,:,1), work2 )
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1) + nu*( work1(:,iy) + work2(:,iy) )
  enddo
  !$omp end parallel do

  ! second component
  call cofdxdx( uk(:,:,2), work1 )
  call cofdydy( uk(:,:,2), work2 )
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,2) = nlk(:,iy,2) + nu*( work1(:,iy) + work2(:,iy) )
  enddo
  !$omp end parallel do

  deallocate( work1, work2 )
end subroutine add_diffusion
