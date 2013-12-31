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
  end module PerformanceMeasurement

