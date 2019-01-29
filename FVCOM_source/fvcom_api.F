 
!                                                                      !
!  FVCOM Main Driver:                                                  !
!                                                                      !
!  This driver executes FVCOM model running.  It                       !
!  controls the initialization, time-stepping, and finalization        !
!  of the model execution following ESMF conventions:        !
!                                                                      !
!     FVCOM_initialize                                                  !
!     FVCOM_run                                                         !
!     FVCOM_finalize                                                    !
!                                                                      !
!=======================================================================
!
MODULE  mod_ocean_control


    !====================================================================!
    !                        USE MODULES                                 !
    !=====================================================================!
    USE MOD_UTILS
    USE CONTROL
    USE MOD_PAR   
    USE MOD_STARTUP
    USE MOD_TIME
    USE MOD_CLOCK
    USE MOD_INPUT
    USE MOD_NCDIO
    USE MOD_NCLL
    USE MOD_SETUP
    USE MOD_SET_TIME
    USE MOD_FORCE
    USE MOD_OBCS
    USE MOD_NESTING
    USE MOD_REPORT
# if defined (DATA_ASSIM)
    USE MOD_ASSIM
# endif  
  
# if defined (LAG_PARTICLE)
    USE MOD_LAG
# endif
    USE PROBES
    USE MOD_BOUNDSCHK !bounds checking
# if defined (SEMI_IMPLICIT)
    USE MOD_SEMI_IMPLICIT
# endif
# if defined (BALANCE_2D)
    USE MOD_BALANCE_2D
# endif  
# if defined (ONE_D_MODEL)
    USE MOD_ONEDTIDE
# endif
# if defined (WATER_QUALITY)
    USE MOD_WQM
# endif    
# if defined (BioGen)
    USE MOD_BIO_3D
# endif    
# if defined (FABM)
    USE MOD_FABM_3D
# endif
# if defined (DYE_RELEASE)
    USE MOD_DYE
# endif
# if defined (RRKF)
   USE RRKVAL
   USE MOD_RRK
   USE MOD_RRKA
   USE MOD_RRKF_OBS
# endif
# if defined (ENKF)  
   USE ENKFVAL 
   USE MOD_ENKF
   USE MOD_ENKFA
   USE MOD_ENKF_OBS
# endif


# if defined (GOTM)
  USE MOD_GOTM 
# endif  
  
# if defined (SEDIMENT)
# if defined (ORIG_SED)
  USE MOD_SED 
# elif defined (CSTMS_SED)
  USE MOD_SED_CSTMS 
# endif
# endif

# if defined (MEAN_FLOW)
  USE MOD_MEANFLOW
  USE MOD_OBCS2
  USE MOD_OBCS3
# endif

# if defined (ICE)
  USE MOD_ICE
# endif

! Added by researchers at Akvaplan-niva 2018, idealized tests give promising results.
# if defined (TVD)
    USE MOD_TVD
# endif

! Added by Qin at Akvaplan-niva 2019
# if defined (ICESHELF)
    USE MOD_ISF
# endif
  
# if defined (NH)
   USE NON_HYDRO
# endif

# if defined (SEMI_IMPLICIT) || defined (NH) || (defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT))
  USE MOD_PETSc, ONLY : PETSc_SET, PETSc_CLEANUP
# endif

! for periodic lateral boundary conditions
# if defined (PLBC)
  USE MOD_PERIODIC_LBC
# endif

# if defined (WAVE_CURRENT_INTERACTION)
  USE TIMECOMM
  USE SWCOMM3
  USE VARS_WAVE
  USE MOD_WAVE_CURRENT_INTERACTION
# endif 

# if defined (WAVE_CURRENT_INTERACTION) && (SEDIMENT)
  USE MOD_BBL
# endif

# if defined (THIN_DAM)
  USE MOD_DAM
# endif

# if defined (PWP)
  USE MOD_PWP
# endif

  USE MOD_STATION_TIMESERIES 
  USE MOD_SPARSE_TIMESERIES
  
 !==============================================================================!
 !==============================================================================! 
	
      IMPLICIT NONE
      PRIVATE
      PUBLIC  :: FVCOM_initialize
      PUBLIC  :: FVCOM_run
      PUBLIC  :: FVCOM_finalize

      CONTAINS
!#ifdef FISOC
!      SUBROUTINE ROMS_initialize (first, mpiCOMM, inputFile)
!#else
!      SUBROUTINE ROMS_initialize (first, mpiCOMM)
!#endif
        SUBROUTINE FVCOM_initialize
!=======================================================================
!                                                                      !
!  This routine allocates and initializes FVCOM state variables    !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
!------------------------------------------------------------------------------|
  IMPLICIT NONE

  character(len=*),parameter::CVS_Id="$Id$" ! [sng] CVS Identification
  character(len=*),parameter::CVS_Date="$Date$" ! [sng] Date string
  character(len=*),parameter::CVS_Name="$Name$" ! [sng] File name string
  character(len=*),parameter::CVS_Revision="$Revision$" ! [sng] File revision string

# if defined (RRKF)
  CHARACTER(LEN=120):: RRKFILES
  CHARACTER(LEN=8)  :: RRKINP1
  CHARACTER(LEN=4)  :: RRKINP2
  CHARACTER(LEN=100)::  MKDIR
# endif
  INTEGER :: IERR

  type(watch) Timer
  
  type(TIME) ::GET_BEGIN
  integer status

  !==============================================================================!
  ! INITIALIZE ALL CONTROL VARIABLES
  !==============================================================================!
  CALL INITIALIZE_CONTROL("FVCOM")

#  if defined (MULTIPROCESSOR)
  ! INTIALIZE MPI CONTROL VARIABLES
  CALL INIT_MPI_ENV(MYID,NPROCS,SERIAL,PAR,MSR,MSRID)
  MPI_FVCOM_GROUP = MPI_COMM_WORLD ! FOR NOW MAKE THEM EQUAL
#  endif

  !==============================================================================!
  !   INITIALIZE A STOP WATCH TIMER FOR TESTING SUBROUTINE EFFICENCY             !
  !==============================================================================!
  CALL WATCH_INIT(TIMER)

  !==============================================================================!
  !   IMPORT CASENAME AND COMMAND LINE ARGUMENTS AND START LOG FILE              !
  !==============================================================================!
  CALL COMMANDLINEIO(CVS_ID,CVS_Date,CVS_Name,CVS_Revision)       
  if(DBG_SET(dbg_log)) Call WRITE_BANNER(PAR,NPROCS,MYID)

  !==============================================================================!
  ! SET DEFAULT VALUES AND READ NAME LISTS                                            
  !==============================================================================!

  CALL NAMELIST

  !==============================================================================!
  !   SET MODEL CONTROL PARAMTERS BASED ON NAME LIST HERE                        !
  !==============================================================================!
  CALL CNTRL_PRMTRS

  !==============================================================================!
  !   SET THE STARTUP TYPE TO BE USED!                                           !
  !==============================================================================!
  CALL SET_STARTUP_TYPE ! see: startup_type.F

  !==============================================================================!
  !   OPEN ALL FILES NEEDED BASED ON THE RUN PARAMETERS                          !
  !==============================================================================!
  CALL OPEN_ALL

  !==============================================================================!
  !   SET MODEL TIME BASED ON THE NAMELIST TIME STRINGS OR RESTART FILE          !
  !==============================================================================!
  CALL SETUP_TIME

  !==============================================================================!
  !   LOAD GRID CONNECTIVITY AND OBC LIST FOR METIS DECOMPOSITION                !
  !==============================================================================!
  CALL LOAD_GRID

  !==============================================================================!
  !   SETUP THE DOMAIN FOR PARALLEL OR SERIAL RUNNING                            !
  !==============================================================================!
  CALL SETUP_DOMAIN

  !==============================================================================!
  !  GET THE PARAMETERS OF FABM                                      !
  !==============================================================================!
#  if defined (FABM)
  IF (FABM_MODEL) CALL FABM_FVCOM_INITIALIZE
#  endif

  !==============================================================================!
  !   ALLOCATE ALL DOMAIN SIZE VARIABLES HERE                                    !
  !==============================================================================!
  CALL ALLOCATE_ALL
# if defined (DYE_RELEASE)
  CALL ALLOC_VARS_DYE
# endif
# if defined (BALANCE_2D)
  CALL ALLOC_BALANCE_VARS
# endif    
  !==============================================================================!
  !   LOAD/SETUP PHYSICAL QUANTITIES (CORIOLIS, GRAVITY, SPONGE LAYER, XY/LATLON)!
  !==============================================================================!
  CALL COORDS_N_CONST

  !==============================================================================!
  ! CALCULATE GRID METRICS - NEIGHBORS, GRADIENTS, CELL AREA, INTERP COEFF'S     !
  !==============================================================================!
  CALL GRID_METRICS

  !==============================================================================!
  ! SETUP THE SEDIMENT MODEL (MUST COME BEFORE SETUP_FORCING)                    ! 
  !==============================================================================!
#  if defined (SEDIMENT)
   IF(SEDIMENT_MODEL) &
  &  CALL SETUP_SED(SEDIMENT_MODEL_FILE,IOBCN,IOBCN_GL,N_SED,N_SED_MAX,SED_NAMES)
#  endif

   !SETUP TVD ADVECTION
# if defined (TVD)
   CALL SETUP_TVD
# endif

!JQI  !==============================================================================!
!JQI  !  SETUP THE MODEL FORCING                                                     !
!JQI  !==============================================================================!
!JQI  CALL SETUP_FORCING

  !==============================================================================!
  !  GET THE PARAMETERS OF BIOLOGICAL MODEL                                      !
  !==============================================================================!
#  if defined (BioGen)
   IF(BIOLOGICAL_MODEL)THEN
     KBV=KB
!     CALL GET_PARAMETER("./"//trim(input_dir)//"/"//trim(BIOLOGICAL_MODEL_FILE))
     CALL GET_PARAMETER_NEW("./"//trim(input_dir)//"/"//trim(BIOLOGICAL_MODEL_FILE))
     CALL BIO_PARAMETER_REPORT
   END IF
!   RIVER_TS_SETTING = 'NONE' !No bio river yet
#  endif

  !==============================================================================!
  !  SETUP THE MODEL FORCING                                                     !
  !==============================================================================!
  CALL SETUP_FORCING

  !==============================================================================!
  !  SETUP OTHER TOOLS, MODELS AND DATA ASSIMILATION                             !
  !==============================================================================!

# if defined (ICE)
  IF(ICE_MODEL) CALL ICE_INIT_0
# endif

# if defined (ICESHELF)
    CALL ISF_INIT_0
# endif

  !  SETUP PETSc FOR SEMI_IMPLICIT AND NON-HYDROSTATIC MODULE
# if defined (SEMI_IMPLICIT) || (NH) || (defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT))
  CALL PETSc_SET
# endif

#  if defined (WAVE_CURRENT_INTERACTION) 
#  if defined (WAVE_OFFLINE)
   CALL WAVE_CURRENT_SETUP
#  else
   CALL SWMAIN_SETUP
   CALL WAVE_CURRENT_SETUP
#  endif
#  endif  

#  if defined (PLBC)
   CALL FIND_NODE_CELL
#  endif

# if defined (WAVE_CURRENT_INTERACTION) && (SEDIMENT)
   CALL INIT_BBL
# endif


  ! SETUP DATA ASSIMILATION MODE
# if defined (DATA_ASSIM) 
  CALL SETUP_DATA_ASSIMILATION
# if defined (PWP)
  CALL SETUP_PWP
# endif
# endif
#  if defined (ENKF)
   if(msr)  print *, 'before enkf_set_time'
     CALL ENKF_SET_TIME
     CALL SET_ASSIM_ENKF_EVE
        IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   if(msr)  print *, 'finish enkf_set_time'
     !check point done
#  endif
! New Open Boundary Condition ----2
#  if defined (MEAN_FLOW)
     CALL FIND_OBSIDE
     IF(MSR)WRITE(IPT,*)'FIND_OBSIDE COMPLETED....'
     CALL ALLOC_OBC3_DATA
     IF(MSR)WRITE(IPT,*)'OBC3_DATA COMPLETED....'
     CALL SETUP_OBC3
     IF(MSR)WRITE(IPT,*)'SETUP_OBC3 COMPLETED....'
     CALL ALLOC_OBC2_DATA
     IF(MSR)WRITE(IPT,*)'OBC2_DATA COMPLETED....'
     MF_RST_STCNT = 0
#    endif
#  if defined (RRKF)

     CALL RRK_SET_TIME
     print *, 'finish rrkf_set_time'
!=============================================================
! THIS IS NOT STANDARD FORTRAN - REMOVE BEFORE PUBLIC RELEASE!
     IF(MSR) THEN 
        MKDIR = "mkdir -p "//TRIM(OUTPUT_DIR)//"/rrktemp"
        CALL SYSTEM( TRIM(MKDIR) )         
     ENDIF
# if defined(MULTIPROCESSOR)     
     CALL MPI_BARRIER(MPI_FVCOM_GROUP,ierr)
# endif
!=============================================================
! IF (RRK_ON) then
    print *, 'are you do here?' ! should be marked
     CALL RRK_SET_STARTUP
    print *, 'finish rrk_set_startup'
     CALL READ_SSH
    print *, 'finish read_ssh' 
     IF(WETTING_DRYING_ON) CALL READ_WETDRY
    print *, 'finish read_read_wetdry' 
     CALL READ_UV
    print *, 'finish read_uv' 
     CALL READ_TURB
    print *, 'finish read_turb' 
     CALL READ_TS
    print *, 'finish read_TS'
!  END IF
# else

  !==============================================================================!
  !  SET THE INITIAL CONDITIONS FOR THE MODEL                                    !
  !==============================================================================!
  CALL STARTUP
# endif

!qxu{ test function READ_DATETIME
!  GET_BEGIN = READ_DATETIME('2009-01-01T00:00:00.0','ymd','UTC',status)
!  CALL PRINT_TIME(GET_BEGIN,IPT,'2009-01-01T00:00:00.0')
!qxu}

  !==============================================================================!
  !  CALL ARCHIVE TO SETUP THE OUTPUT AND DUMP CONSTANT VALUES                   !
  !==============================================================================!
  CALL BCOND_GCN(8,0)
  
  CALL ARCHIVE
  ! ORDER MATTERS - ARCHIVE_NEST MUST GO AFTER ARCHIVE DURING SETUP!
  CALL ARCHIVE_NEST
# if defined(ENKF)
   enkf_out_Time=NC_DAT%FTIME%NEXT_IO
   DO ENKF_memberCONTR =1, ENKF_NENS
   write(enkf_num,'(i2.2)')  ENKF_memberCONTR
    IF(NC_ON)THEN
             call    ENKF_SETUP_NC(NC_DAT)
             IF (ENKF_memberCONTR == ENKF_NENS) ARCHIVE_TIMES=ARCHIVE_TIMES+1

    END IF
   END DO
# endif
# if defined (WAVE_CURRENT_INTERACTION)
  CALL ARCHIVE_NEST_WAVE
# endif

  CALL SET_PROBES(PROBES_ON,PROBES_NUMBER,PROBES_FILE)
  IF(OUT_STATION_TIMESERIES_ON)CALL READ_STATION_FILE
# if defined (WAVE_CURRENT_INTERACTION)
  IF(OUT_WAVE_SPARSE_TIMESERIES_ON)CALL SPARSE_STATION
# endif

  ! Setup Bounds checking (shutdown if variables exceed threshold)
  CALL SETUP_BOUNDSCHK !bounds checking

# if defined (LAG_PARTICLE)
  CALL SET_LAG
# endif

# if defined (GOTM)
  CALL INIT_GOTM
# endif
  
  IF(OUT_STATION_TIMESERIES_ON)THEN
    CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL),INTERVAL_TIME_SERIES)
    CALL OUT_STATION_TIMESERIES
    TIME_SERIES = STARTTIME + INTERVAL_TIME_SERIES
  END IF
# if defined (WAVE_CURRENT_INTERACTION)
  IF(OUT_WAVE_SPARSE_TIMESERIES_ON)THEN
    CALL GET_OUTPUT_FILE_INTERVAL(TRIM(OUT_INTERVAL_SPARSE),WAVE_INTERVAL_TIME_SERIES)
    CALL OUT_WAVE_SPARSE_TIMESERIES
    WAVE_TIME_SERIES = STARTTIME + WAVE_INTERVAL_TIME_SERIES
  END IF
# endif
  
      RETURN
   END SUBROUTINE FVCOM_initialize
   

   !=========================================================================!
   !                        RUN                                              !
   !=========================================================================!
   SUBROUTINE    FVCOM_run
     RETURN
   END SUBROUTINE FVCOM_run
      
   !=========================================================================!
   !                        FINALIZE                                          !
   !=========================================================================!

   SUBROUTINE    FVCOM_finalize
     RETURN
   END SUBROUTINE FVCOM_finalize

   END MODULE mod_ocean_control
