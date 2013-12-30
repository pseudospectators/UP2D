! +bei jedem aufruf der function soll gestoppt werden. 
! +schreibe dazu, wenn called mit "start", die aktuelle zeit in eine save-variable.
! +dannach berechne bei "stop" die differenz in sekunden. 
! +schreibe den wert in ein array von 100 einträgen. wenn voll, fange wieder von vorne an (save-variable)
! +gebe den mittelwert der nicht-null einträge zurück (in sekunden). klappt also auch für wenige (<100) aufrufe. 
! +die zurückegegebene zeitmessung stellt also immer den durchschnittlichen zeitbedarf von 100 calls dar. 
! ----------------------------------------
! do ...
! 	time_NST = AvgPerformance("start","avg",1)
! 	call DoNavierStokesStep()
! 	time_NST = AvgPerformance("stop","avg",1)
! 	if (SaveNow) then write(33,*) time_NST
! enddo
! ---------------------------------------
! verzichte auf globale variable. das modul bekommt ein array 10x100, für 10 kanäle mit 100 werten
! ----------
! modul soll auch umrechnung in d:h:m:s enthalten. 
! ---------
! globale messungen (laufzeit, gesamtperformance) können ersetzt werden. für die laufzeit braucht man kein avg
! auch timestep performance kann gut in 100er schritte umgemünzt werden, muss nicht global gemittelt sein.
! ------------
! speichern: time/dt; summe; time_NSt; time_press; time_solid; time_mask



module PerformanceMeasurement
  !-------------------------------------------------------
  implicit none
  save 
  integer, parameter                        :: N_avg =20, N_chans = 20
  real, dimension (1:N_chans, 1:N_avg)      :: timetable = 0.0 !contains 20 values for each channel (for averaging)
  integer*8, dimension (1:N_chans)          :: times    = 0  !contains the time when the measurement started
  integer, dimension (1:N_chans)            :: counters = 1  !counts the measurements

  !-------------------------------------------------------
  
  contains

  real function Performance(Start, channel)
    implicit none
    integer, intent (in) :: channel
    integer*8 :: rate, dummy, time_temp
    character (len=*), intent(in) :: Start
    
    if (Start == "start") then
      call system_clock(times(channel), rate, dummy) !save start time to the "times"-array  
      performance = 0.0
    elseif (Start == "stop") then
      call system_clock(time_temp, rate, dummy) !save start time to the "times"-array       
      timetable(channel,counters(channel)) =  real(time_temp-times(channel)) / real(rate) !difference between stop and start in seconds
      ! now return the value. note that before the array is full, the value is too small. I don't mind, do you?
      performance = sum(timetable(channel,:)) / real(N_avg)
      counters(channel) = counters(channel) + 1 !increase counter
      if (counters(channel)>N_avg) counters(channel) = 1
    else
      write (33,*) "Performance encountered a problem (invalid command!)"
      stop
    endif
    
    return    
  end function Performance

  !----------------------------------------------------------------------------

  
  real function GetRuntime(Start)
    character (len=*), intent(in) :: Start
    integer*8, save :: rate, dummy,starttime, currenttime

    if (Start == "start") then
      call system_clock(starttime, rate, dummy) !save start time to the "times"-array
      GetRuntime = 0.0
    elseif (Start == "now") then
      call system_clock(currenttime, rate, dummy) !save start time to the "times"-array
      GetRuntime = real(currenttime-starttime) / real(rate) !difference between stop and start in seconds
    endif
    
    return
  end function GetRuntime
  
  !----------------------------------------------------------------------------

  subroutine SavePerformance(time,time_dt, time_nst, time_pressure, time_mask, time_solid, steps_left,dt1)
    use share_vars
    implicit none
    
    real, intent(in) :: time, time_dt, time_nst, time_pressure, time_mask, time_solid,dt1
    integer, intent(in) :: steps_left
    real :: totaltime, time_left, time_sum
    integer :: days,mins,hours,secs,days1,mins1,hours1,secs1

    totaltime=GetRuntime('now')

    days   =floor(totaltime/(24.*3600.))
    hours  =floor(mod(totaltime,24.*3600.)/3600.)
    mins   =floor(mod(totaltime,3600.)/60.)
    secs   =floor(mod(mod(totaltime,3600.),60.))

    time_left = steps_left*time_dt
    
    days1   =floor(time_left/(24.*3600.))
    hours1  =floor(mod(time_left,24.*3600.)/3600.)
    mins1   =floor(mod(time_left,3600.)/60.)
    secs1   =floor(mod(mod(time_left,3600.),60.))
    
    time_sum=time_nst+time_pressure+time_mask+time_solid
    
    open  (10, file = trim(dir_name)//'/'//trim(simulation_name)//"performance_details", status = 'unknown', access = 'append')
    write (10,'( 2(i3,"d ",i2,"h ",i2,"m ",i2,"s "), 2(es11.4,1x)," - ",2(es11.4,1x),4(es11.4,1x,i3,"% "))') &
        days1,hours1,mins1,secs1,days,hours,mins,secs,&
        time,dt1,time_dt,time_sum,&
        time_nst, nint(100.*time_nst/time_sum), time_pressure, &
        nint(100.*time_pressure/time_sum), time_mask, nint(100.*time_mask/time_sum),&
        time_solid, nint(100.*time_solid/time_sum)
    close (10)
   
  end subroutine SavePerformance

end module PerformanceMeasurement

