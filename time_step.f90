subroutine time_step 
  use share_vars
  use FieldExport
  use spectral_solver
  use spectral_essentials
  use PerformanceMeasurement

  implicit none
  real(kind=pr) :: time, dt1
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: p, vort
  real(kind=pr) :: dt, T_lastdrag, T_lastsave, t1, time_left
  integer :: it=0, iy
  character(len=17) :: timestring
  character (len=11) :: name

  time = 0.0
  dt1 = 0.0
  it = 0
  T_lastdrag = 0.0
  T_lastsave = 0.0
  
  !----------------------------------------------------------------  
  ! Initialize vorticity or read values from a backup file
  !----------------------------------------------------------------
  call init_fields (u, uk, p, vort, nlk)

  !----------------------------------------------------------------
  ! create startup mask
  !---------------------------------------------------------------- 
  call create_mask (time)
  call SaveGIF(mask, trim(dir_name)//'/'//trim(simulation_name)//"startup_mask", 13)  

  
  !----------------------------------------------------------------
  ! loop over time steps
  !---------------------------------------------------------------- 
  do while ((time<Tmax) .and. (it<=nt))
      t1 = Performance("start",1)
      
      call Runge (time, dt1,it, u, uk, p, vort, nlk) 
      
      time = time + dt1  ! Advance in time
      it = it + 1
      
      if (time-T_lastdrag>tdrag) then 
        !----------------------------------------------------------------
        !-- video snapshots
        !----------------------------------------------------------------
        write (timestring,'(i5.5)') nint(time*100.d0)         
        colorscale = 0.10*max(maxval(vort),abs(minval(vort))) 
        call SaveGIF(vort, trim(dir_name)//"/vor/"//trim(timestring)//".vor", 1, -colorscale, colorscale)
                        
        write (*,*) "time=", time
        T_lastdrag=time
      endif
!   
!       !--------------------------------------------------------------------------------------------
!       !--   save fields
!       !--------------------------------------------------------------------------------------------
!       if (time-T_lastsave>tsave) then
! 	  write (name, '(es10.4)') time
! 	  call cofdx ( uk(:,:,2), vort )
! 	  call cofdy ( uk(:,:,1), p )
! 	  !$omp parallel do private(iy)
! 	  do iy=0,ny-1
! 	      p(:,iy) = vort(:,iy) - p(:,iy)
! 	  enddo
! 	  !$omp end parallel do  
! 	  
! 	  call cofitxy( p, vort)	  
! 	  call SaveField (trim(dir_name)//'/fields/'//trim(simulation_name)//'vor_'//name, vort, 1, xl,yl, 'precision')
!   
! 	  T_lastsave=time
!       endif 
      t1 = Performance("stop",1)
      !----------------------------------------------------------------
      !-- remaining time
      !----------------------------------------------------------------      
      if (modulo(it,200)==0) then
        time_left = t1*(Tmax-time)/dt1
        write (*,'("time left=",i3,"d ",i2,"h ",i2,"m ",i2,"s (",i2,"%)")') &
          floor(time_left/(24.*3600.)),&
          floor(mod(time_left,24.*3600.)/3600.),&
          floor(mod(time_left,3600.)/60.),&
          floor(mod(mod(time_left,3600.),60.)),&
          nint(100.*time/Tmax)
      endif
  enddo
end subroutine time_step











