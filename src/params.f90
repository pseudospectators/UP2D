! Wrapper to read parameters from an ini file for fsi.  Also reads
! parameters which are set in the vars module.
subroutine get_params(paramsfile)
  use share_vars
  use ini_files_parser

  type(inifile) :: PARAMS
  character(len=*) :: paramsfile

  ! Read the paramsfile and put the length i and the text in PARAMS
  call read_ini_file(PARAMS, paramsfile, .true.)

  ! Resolution section
  call read_param(PARAMS,"Resolution","nx",nx, 4)
  call read_param(PARAMS,"Resolution","ny",ny, 4)
  call read_param(PARAMS,"Time","nt",nt, 9999999)
  call read_param(PARAMS,"Time","Tmax",Tmax,1.d9)
  call read_param(PARAMS,"Time","CFL",cfl,0.1d0)
  call read_param(PARAMS,"Time","iMethod",iMethod,"RK2")
  call read_param(PARAMS,"Time","dt_fixed",dt_fixed,0.0d0)
  call read_param(PARAMS,"Time","dt_max",dt_max,0.0d0)
  call read_param(PARAMS,"ReynoldsNumber","nu",nu,1.d-2)

  ! Initial conditions section
  call read_param(PARAMS,"InitialCondition","inicond",inicond, "none")

  ! Dealasing section
  call read_param(PARAMS,"Dealiasing","iDealias",iDealias, 1)

  ! Penalization section
  call read_param(PARAMS,"Penalization","iMask",iMask, "none")
  call read_param(PARAMS,"Penalization","eps",eps, 1.d-2)

  ! Geometry section
  call read_param(PARAMS,"Geometry","xl",xl, 1.d0)
  call read_param(PARAMS,"Geometry","yl",yl, 1.d0)

  ! saving section
  call read_param(PARAMS,"Saving","tsave",tsave, 1.d0)
  call read_param(PARAMS,"Saving","itsave",itsave, 999999)
  call read_param(PARAMS,"Saving","tdrag",tdrag, 1.d0)
  call read_param(PARAMS,"Saving","iSaveVorticity",iSaveVorticity, 1)
  call read_param(PARAMS,"Saving","iSaveVelocity",iSaveVelocity, 1)
  call read_param(PARAMS,"Saving","iSavePressure",iSavePressure, 1)
  call read_param(PARAMS,"Saving","iSaveMask",iSaveMask, 1)

  ! sponge
  call read_param(PARAMS,"Sponge","iSpongeType",iSpongeType, "none")
  call read_param(PARAMS,"Sponge","use_sponge",use_sponge, 0)
  call read_param(PARAMS,"Sponge","eps_sponge",eps_sponge, 1.0d0)

  ! mean flow
  call read_param(PARAMS,"MeanFlow","iMeanFlow",iMeanFlow, "none")
  call read_param(PARAMS,"MeanFlow","ux_mean",ux_mean, 0.d0)
  call read_param(PARAMS,"MeanFlow","uy_mean",uy_mean, 0.d0)

  ! Set other parameters (all procs)
  pi = 4.d0 * datan(1.d0)
  ! lattice spacing is global
  dx = xl/dble(nx)
  dy = yl/dble(ny)


  ! clean ini file
  call clean_ini_file(PARAMS)

end subroutine get_params
