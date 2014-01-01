subroutine time_step 
  use share_vars
  use masks
  use FieldExport  
  use RK2_module
  use timing
  implicit none
  real(kind=pr) :: time=0.0d0, dt1=0.0d0, max_divergence
  real(kind=pr), dimension(0:nx-1,0:ny-1,1:2) :: u, uk, nlk
  real(kind=pr), dimension(0:nx-1,0:ny-1) :: pk, vort
  real(kind=pr) :: T_lastdrag=0.0d0, T_lastsave=0.0d0, t1, time_left
  integer :: it=0
  character(len=5) :: timestring
  character(len=11) :: name
  
  !-- Initialize vorticity or read values from a backup file
  call init_fields (u, uk, pk, vort, nlk)
  
  !-- create startup mask
  call create_mask (time)
  call SaveGIF(mask, trim(simulation_name)//"startup_mask", 13)   
  call SaveGif(us(:,:,1),trim(simulation_name)//'usx')
  call SaveGif(us(:,:,2),trim(simulation_name)//'usy')
  
  !----------------------------------------------------------------
  ! loop over time steps
  !---------------------------------------------------------------- 
  do while ((time<Tmax) .and. (it<=nt))
      t1 = Performance("start",1)
      
      !----------------------------------------------------------------
      !-- Actual time step
      !----------------------------------------------------------------
      select case (iMethod)
      case ('RK2')
        call RK2 (time, dt1,it, u, uk, pk, vort, nlk) 
      case ('RK2_implicit')
        call RK2_implicit (time, dt1,it, u, uk, pk, vort, nlk) 
      case default
        write (*,*) "Error: iMethod undefined"
      end select
      
      
      time = time + dt1  ! Advance in time
      it = it + 1      
      
      if (time-T_lastdrag>tdrag) then 
        !----------------------------------------------------------------
        !-- video snapshots
        !----------------------------------------------------------------
        write (timestring,'(i5.5)') nint(time*100.d0)         
        colorscale = 0.25d0*max(maxval(vort),dabs(minval(vort))) 
        call SaveGIF(vort, "vor/"//trim(timestring)//".vor", 1, -colorscale, colorscale)
                        
        write (*,'("Snapshot. time=",es12.4," vormax=",es12.4)') time, max(maxval(vort),dabs(minval(vort))) 
        T_lastdrag=time
      endif
 
      t1 = Performance("stop",1)
      !----------------------------------------------------------------
      !-- remaining time
      !----------------------------------------------------------------      
      if (modulo(it,200)==0) then
        time_left = t1*(Tmax-time)/dt1
        write (*,'("time left=",i3,"d ",i2,"h ",i2,"m ",i2,"s (",es8.2," s/dt) [",i2,"%] divu=",es12.4)') &
          floor(time_left/(24.d0*3600.d0)),&
          floor(mod(time_left,24.d0*3600.d0)/3600.d0),&
          floor(mod(time_left,3600.d0)/60.d0),&
          floor(mod(mod(time_left,3600.d0),60.d0)),&
          t1, &
          nint(100.d0*time/Tmax), max_divergence(uk)
      endif
  enddo
end subroutine time_step