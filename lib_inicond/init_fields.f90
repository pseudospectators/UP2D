subroutine init_fields (u, uk, pk, vor, nlk)
  use share_vars
  use FieldExport
  implicit none
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2), intent (inout) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1), intent (inout) :: vor, pk
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: vortk, work1, work2
  real(kind=pr) :: r0,we,d,r1,r2
  integer :: ix,iy
  
  u = 0.d0
  uk = 0.d0
  vor = 0.d0
  pk = 0.d0
  
  select case (inicond)
  case ('quiescent')
    u = 0.d0
    uk= 0.d0
    vor = 0.d0
    pk = 0.d0
  case ('lamballais')
    call lamballais(u,uk,pk,vor,nlk)
  case ('couette')
  case ('turbulent')
    call random_seed()
    do ix=0,nx-1
    do iy=0,ny-1
       call RANDOM_NUMBER(d)
       vor(ix,iy) = 200.d0*(2.0d0*d - 1.d0)
    enddo
    enddo    
    
    call coftxy (vor, vortk)    
    call vorticity2velocity ( vortk, u )    
    call coftxy( u(:,:,1), uk(:,:,1))
    call coftxy( u(:,:,2), uk(:,:,2))
  case ('dipole')  
    x0 = 0.5d0*xl
    y0 = 0.5d0*yl
    r0 = 0.1d0
    we = 299.528385375226d0
    d = 0.1d0
    
    do ix=0,nx-1
    do iy=0,ny-1
       r1 = sqrt( (real(ix)*dx-x0)**2 + (real(iy)*dy-y0-d)**2 ) / r0
       r2 = sqrt( (real(ix)*dx-x0)**2 + (real(iy)*dy-y0+d)**2 ) / r0
       vor(ix,iy) = we * (1.d0-r1**2)*exp(-r1**2) - we * (1.d0-r2**2)*exp(-r2**2)
    enddo
    enddo
    
    call coftxy (vor, vortk)    
    call vorticity2velocity ( vortk, u )    
    call coftxy( u(:,:,1), uk(:,:,1))
    call coftxy( u(:,:,2), uk(:,:,2))
     
  end select

  !-- compute initial pressure (for implicit penalization)
  call cal_pressure ( 0.d0, u, uk, pk )
  
  call cofdx(uk(:,:,1),work1)
  call cofdy(uk(:,:,2),work2)
  call cofitxy(work1+work2,vortk)
  write(*,*) "ini2", maxval(vortk)
  write(*,*) kind(maxval(vortk))
  write(*,*) inicond
  write(*,'("inicond field divergence ",es12.4)') maxval(vortk)
  
end subroutine 