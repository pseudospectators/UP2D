subroutine RK2_implicit (time, dt,it, u, uk, pk, vort, nlk)
  use share_vars
  use rhs
  use FieldExport
  use masks
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (inout) :: vort, pk
  real(kind=pr), intent (out) :: dt
  real(kind=pr), intent (in) :: time
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1,work2, workvis,div,divk
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2) :: nlk2, uk_tmp, u_tmp
  integer :: iy
  integer, intent(in) :: it
  real(kind=pr) :: timestep

  !-- modify: no dt<eps
  dt = timestep(time,it, u )
  !-----------------------------------------------------------------------------
  !-- 1st strang step: half time step for penalization
  !-- equation is solved exactly
  !-----------------------------------------------------------------------------
!  call create_mask(time)
  !$omp parallel do private (iy)
  do iy=0,ny-1
    u(:,iy,1) = (u(:,iy,1)-mask(:,iy)*eps*us(:,iy,1))*exp(-0.5d0*dt*mask(:,iy)) &
                 + mask(:,iy)*eps*us(:,iy,1)
    u(:,iy,2) = (u(:,iy,2)-mask(:,iy)*eps*us(:,iy,2))*exp(-0.5d0*dt*mask(:,iy)) &
                 + mask(:,iy)*eps*us(:,iy,2) 
  enddo
  !$omp end parallel do   
  
  call coftxy( u(:,:,1), uk(:,:,1) )
  call coftxy( u(:,:,2), uk(:,:,2) )
  
  !$omp parallel do private (iy)
  do iy=0,ny-1
    uk(:,iy,1) = dealiase(:,iy)*uk(:,iy,1)
    uk(:,iy,2) = dealiase(:,iy)*uk(:,iy,2)
  enddo
  !$omp end parallel do   
  
  !-----------------------------------------------------------------------------
  !-- 2nd Strang step: a full time step with RK2 and 
  !-- integrating factor, but without the penalization and
  !-- the pressure. So just diffusion and non-linear transport
  !-----------------------------------------------------------------------------
  call cofitxy ( uk(:,:,1),u(:,:,1) )
  call cofitxy ( uk(:,:,2),u(:,:,2) )
  
  !-- RHS, without penalization
  call cal_nlk(time, u, uk, vort, nlk, .false.)
  !-- add old pressure term
  call add_pressure_grad( nlk, pk )
  !-- compute integrating factor
  call cal_vis( dt, workvis )
  
  !$omp parallel do private (iy)
  do iy=0,ny-1
    uk_tmp(:,iy,1) = (uk(:,iy,1) + dt*nlk(:,iy,1) )*dealiase(:,iy)*workvis(:,iy)
    uk_tmp(:,iy,2) = (uk(:,iy,2) + dt*nlk(:,iy,2) )*dealiase(:,iy)*workvis(:,iy)
  enddo
  !$omp end parallel do
  
  !-- mean flow forcing
  call mean_flow(uk_tmp)
  !-- to phys space again
  call cofitxy ( uk_tmp(:,:,1),u_tmp(:,:,1) )
  call cofitxy ( uk_tmp(:,:,2),u_tmp(:,:,2) ) 
   
  !--RHS, without penalization
  call cal_nlk(time, u_tmp, uk_tmp, vort, nlk2, .false.)
  !--add old pressure term
  call add_pressure_grad( nlk2, pk )   
  
  !-- sum up all the terms.
  !$omp parallel do private (iy)
  do iy=0,ny-1
      uk(:,iy,1) = (uk(:,iy,1)*workvis(:,iy) &
          + 0.5d0*dt*(nlk(:,iy,1)*workvis(:,iy) + nlk2(:,iy,1)) )*dealiase(:,iy)      
      uk(:,iy,2) = (uk(:,iy,2)*workvis(:,iy) &
          + 0.5d0*dt*(nlk(:,iy,2)*workvis(:,iy) + nlk2(:,iy,2)) )*dealiase(:,iy)
  enddo
  !$omp end parallel do 
  
  !-- mean flow forcing
  call mean_flow(uk)    
  !-- velocity in phys. space
  call cofitxy (uk(:,:,1), u(:,:,1))
  call cofitxy (uk(:,:,2), u(:,:,2))    
  
  !-----------------------------------------------------------------------------
  !-- 3rd strang step: half time step for penalization
  !-----------------------------------------------------------------------------
!  call create_mask (time+dt)
  !$omp parallel do private (iy)
  do iy=0,ny-1
    u(:,iy,1) = (u(:,iy,1)-mask(:,iy)*eps*us(:,iy,1))*exp(-0.5*dt*mask(:,iy)) &
                 + mask(:,iy)*eps*us(:,iy,1)
    u(:,iy,2) = (u(:,iy,2)-mask(:,iy)*eps*us(:,iy,2))*exp(-0.5*dt*mask(:,iy)) &
                 + mask(:,iy)*eps*us(:,iy,2) 
  enddo
  !$omp end parallel do   
  
  call coftxy( u(:,:,1), uk(:,:,1) )
  call coftxy( u(:,:,2), uk(:,:,2) )
  
  !$omp parallel do private (iy)
  do iy=0,ny-1
    uk(:,iy,1) = dealiase(:,iy)*uk(:,iy,1)
    uk(:,iy,2) = dealiase(:,iy)*uk(:,iy,2)
  enddo
  !$omp end parallel do 
    
  !-----------------------------------------------------------------------------
  !-- Projection step
  !-----------------------------------------------------------------------------
  !-- compute divergence (of predicted velocity field)
  call cofdx(uk(:,:,1),work1)
  call cofdy(uk(:,:,2),work2)
  divk = work1 + work2  
  
  !-----------------------
!   call cofitxy (divk,div) 
!   !$omp parallel do private (iy)  
!   do iy=0,ny-1
!     div(:,iy) = div(:,iy)*(1.d0-mask(:,iy)*eps)  
!   enddo
!   !$omp end parallel do
!   call coftxy  (div,divk)
  !-----------------------  
  
  call poisson ( divk, workvis )
  
  !-- add pressure gradient
  call cofdx(workvis,work1)
  call cofdy(workvis,work2)  
  !$omp parallel do private (iy)
  do iy=0,ny-1
    uk(:,iy,1) = uk(:,iy,1) + work1(:,iy)
    uk(:,iy,2) = uk(:,iy,2) + work2(:,iy)
    ! add pressure increment to old pressure
    pk(:,iy) = pk(:,iy) + workvis(:,iy) / dt
  enddo
  !$omp end parallel do   
  
  !-- uk should now be divergence free
  
  !-- velocity in phys. space (for consistent output)
  call cofitxy (uk(:,:,1), u(:,:,1))
  call cofitxy (uk(:,:,2), u(:,:,2))      
end subroutine RK2_implicit

