module rhs

  implicit none
  contains

subroutine cal_nlk (time, u, uk, vor, nlk, penalization)
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent (in) :: 					time
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: 		vor
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	uk
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: 	nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2
  logical, intent (in) :: penalization
  integer :: iy  

  !-- compute vorticity
  call curl (uk, work1)  
  call cofitxy(work1, vor)
    
  !-- NL term + penal
  if ( penalization ) then
      !$omp parallel do private(iy)
      do iy=0,ny-1
        work1(:,iy) = +vor(:,iy)*u(:,iy,2) -mask(:,iy)*(u(:,iy,1)-u_BC(:,iy,1))
        work2(:,iy) = -vor(:,iy)*u(:,iy,1) -mask(:,iy)*(u(:,iy,2)-u_BC(:,iy,2))
      enddo
      !$omp end parallel do
  else
      !$omp parallel do private(iy)
      do iy=0,ny-1
        work1(:,iy) = +vor(:,iy)*u(:,iy,2)
        work2(:,iy) = -vor(:,iy)*u(:,iy,1)
      enddo
      !$omp end parallel do  
  endif
  
  !-- NLK to fourier space
  call coftxy ( work1, nlk(:,:,1) )
  call coftxy ( work2, nlk(:,:,2) )  
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    nlk(:,iy,1) = nlk(:,iy,1)*dealiase(:,iy)
    nlk(:,iy,2) = nlk(:,iy,2)*dealiase(:,iy)
  enddo
  !$omp end parallel do  
end subroutine cal_nlk

end module rhs