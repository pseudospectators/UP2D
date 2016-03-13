subroutine add_pressure_grad( nlk, pk )
  use share_vars
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (in) :: pk
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: work1, work2
  integer :: iy
  
  !-- compute gradient
  call cofdx (pk, work1)
  call cofdy (pk, work2)
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1) + work1(:,iy)    
    nlk(:,iy,2) = nlk(:,iy,2) + work2(:,iy)
  enddo
  !$omp end parallel do   
end subroutine add_pressure_grad
