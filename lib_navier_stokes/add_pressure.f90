subroutine add_pressure (nlk, uk, u, vor)
  use share_vars
  use fieldexport
  implicit none
  real(kind=pr), dimension (0:nx-1,0:ny-1,1:2), intent (inout) :: nlk
  real(kind=pr), dimension (0:nx-1,0:ny-1,1:2), intent (in) :: u, uk
  real(kind=pr), dimension (0:nx-1,0:ny-1), intent (in) :: vor
  real(kind=pr), dimension (0:nx-1,0:ny-1) :: work1, work2, work3, work4, work5
  integer :: iy
  real(kind=pr)::max_divergence
  
  if (ipressure == "classic") then
    !--------------------------------------------
    ! classic pressure projects all source terms
    !-------------------------------------------
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
      nlk(:,iy,1) = nlk(:,iy,1) + work1(:,iy)    
      nlk(:,iy,2) = nlk(:,iy,2) + work2(:,iy)
    enddo
    !$omp end parallel do   
    
  elseif (ipressure == "modified" ) then
  
    !---------------------------------------------------------------------------
    ! modified pressure relaxes incompressibility constraint
    ! div(u) is not zero to machine precision
    ! The divergence then obeys a penalized diffusion law:
    !     d/dt DIV = nu*laplace(DIV) - chi/eta * DIV
    !---------------------------------------------------------------------------    
    !-- we need to recompute the non-linear transport term
    !$omp parallel do private(iy)
    do iy=0,ny-1
      work1(:,iy) = +vor(:,iy)*u(:,iy,2)
      work2(:,iy) = -vor(:,iy)*u(:,iy,1)
    enddo
    !$omp end parallel do  
    !-- NL term to fourier space
    call coftxy ( work1, work3 )
    call coftxy ( work2, work4 )
    
    !-- divergence of the NL term:
    call cofdx ( work3, work1 )
    call cofdy ( work4, work2 )
    work5 = work1 + work2
    
    !---------------------------------------------------------------------------
    !-- penalized divergence of velocity field (1-chi)*(div(u))
    call cofdx(uk(:,:,1),work1)
    call cofdy(uk(:,:,2),work2)
    !-- div to phys space
    call cofitxy( work1+work2,work3 )    
    !-- penalize it
    !$omp parallel do private(iy)
    do iy=0,ny-1
      work4(:,iy) = work3(:,iy)*(1.d0-mask(:,iy)*eps)/eps
    enddo    
    !$omp end parallel do     
    
    !back to Fourier space
    call coftxy(work4, work3)    
    
    ! this is the RHS of poisson:
    !$omp parallel do private(iy)
    do iy=0,ny-1
      work3(:,iy) = (work5(:,iy) + work3(:,iy) )*dealiase(:,iy)
    enddo    
    !$omp end parallel do      
    
    !---------------------------------------------------------------------------
    !-- solve poisson eqn
    call poisson (work3, work1)
    !-- gradient
    call cofdx (work1, work3)
    call cofdy (work1, work2)
    !-- add gradient
    !$omp parallel do private(iy)
    do iy=0,ny-1
      nlk(:,iy,1) = nlk(:,iy,1) + work3(:,iy)    
      nlk(:,iy,2) = nlk(:,iy,2) + work2(:,iy)
    enddo
    !$omp end parallel do  
  else
    write(*,*) "ipressure undefined..."
    stop
  endif

end subroutine add_pressure