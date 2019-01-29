
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
MODULE mod_ocean_control


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

  SUBROUTINE FVCOM_initialize(casename_opt)
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
    character(len=*),INTENT(IN),OPTIONAL :: casename_opt

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

    IF (PRESENT(casename_opt)) THEN
       casename = casename_opt
    END IF

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
  SUBROUTINE    FVCOM_run(START_DATE_opt,END_DATE_opt)
    INTEGER STATUS
    character(len=*),INTENT(IN),OPTIONAL :: START_DATE_opt, END_DATE_opt

    TYPE(TIME) :: StartTime_sub
    TYPE(TIME) :: EndTime_sub

    IF (PRESENT(START_DATE_opt)) THEN ! thinf: compute coupled partial intergration period

          ! GET THE FISOC SUBSTEP START TIME                                                                                       
          StartTime_sub = READ_DATETIME(START_DATE_opt,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the FISOC substep date string START_DATE: "//trim(START_DATE_opt))

          ! GET THE FISOC SUBSTEP END TIME                                                                                         
          EndTime_sub = READ_DATETIME(END_DATE_opt,DATE_FORMAT,TIMEZONE,status)
          if (status == 0) &
               & Call Fatal_Error("Could not read the FISOC substep date string END_DATE: "&
                      &//trim(END_DATE_opt))
          ! SANITY CHECK                                                                                             
          if(StartTime_sub .GT. EndTime_sub) &
               & Call Fatal_Error("Runfile Start_Date exceeds or equal to End_Date")

          
          !! FIND THE START TIME IN THE RESTART FILE AND SET ISTART                                                           !     
          !CALL SET_STARTUP_FILE_STACK(StartTime,IINT)
          !
          ! ADVANCE THE TO THE FIRST TIME STEP OF THE MODEL FROM THE                                                 
          ! INITIAL CONDITION
          ! thinf: reuse previous timestep
          ISTART = IINT +1

          !CALCULATE THE NUMBER OF STEPS AND IEND                                                                    
          NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime_sub,EndTime_sub)
          IEND = ISTART + NSTEPS


    END IF

    !==============================================================================!
    !  SELECT THE RUN MODE AND EXECUTE THE MAIN LOOP
    !==============================================================================!
    SELECT CASE(FVCOM_RUN_MODE)
       ! RUN MODE SET IN mod_assim.F(set_assim_param)

       ! =============================================================================!
       ! == PURE SIMULATION MODE - Instantanious data assimilation only ==============!
    CASE(FVCOM_PURE_SIM)
       ! =============================================================================!


       !==============================================================================!
       !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
       !==============================================================================!
       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "!===================================================="
          write(IPT,*) "!===================================================="
          write(IPT,*) "!============== STARTING MAIN LOOP AT:==============="
          if(DBG_SET(dbg_log)) &
               & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
          write(IPT,*) "!===================================================="
       end if

       CALL REPORT('INITIAL CONDITIONS')

       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "!===================================================="
          write(IPT,*) "!===================================================="
          write(IPT,*) "!===================================================="
       end if

       !////////////////////////// MAIN LOOP //////////////////////////////////////////
       DO IINT=ISTART,IEND

          IntTime=IntTime + IMDTI

          CALL INTERNAL_STEP

          !==============================================================================!
          !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
          !==============================================================================!
          if(DBG_SET(dbg_log)) &
               & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)

          IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

          !==============================================================================!
          !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
          !==============================================================================!

          !for Eulerian velocity output
# if defined (WAVE_CURRENT_INTERACTION)
          !U = U-U_STOKES_3D
          !V = V-V_STOKES_3D
          !UA=UA-U_STOKES_2D
          !VA=VA-V_STOKES_2D
# endif

          CALL ARCHIVE

#  if defined (WAVE_CURRENT_INTERACTION)
          !Convert Eulerian velocity back to Lagrangian velocity to maintain consistent with Mellor's equation
          !U = U+U_STOKES_3D
          !V = V+V_STOKES_3D
          !UA=UA+U_STOKES_2D
          !VA=VA+V_STOKES_2D
#  endif

          CALL DUMP_PROBE_DATA
          IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
# if defined (WAVE_CURRENT_INTERACTION)
          IF(OUT_WAVE_SPARSE_TIMESERIES_ON) CALL OUT_WAVE_SPARSE_TIMESERIES
# endif
          !==============================================================================!
          !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
          !==============================================================================!
          CALL SHUTDOWN_CHECK(D1)

          !==============================================================================!
          !  CALL BOUNDS CHECK TO SEE IF VARIABLES EXCEED USER-DEFINED THRESHOLDS 
          !==============================================================================!
          CALL BOUNDSCHK  !bounds checking

          !==============================================================================!
          !    LAGRANGIAN PARTICLE TRACKING                                              |
          !==============================================================================!
# if defined (LAG_PARTICLE)
          CALL LAG_UPDATE
# endif

          !==============================================================================!
          !    NESTING OUTPUT                                                            |
          !==============================================================================!
          IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
          IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif

       END DO
       !////////////////////////// END MAIN LOOP //////////////////////////////////////

       ! ================================================================================!
       ! == THIS MAIN LOOP USES NUDGING OR OI METHODS TO ASSIMILATE =====================!
    CASE(FVCOM_NUDGE_OI_ASSIM)
       ! ================================================================================!
# if defined (DATA_ASSIM) 

       CALL ALLOC_BUFFER

       ! SET THE ASSIMILATION/SIMULATION RESET TIME
       ASSIM_RESTART_TIME = StartTime
       IF(TSGRD_ASSIM) THEN
          CALL Now_2_month_days(StartTime,Pyear,Pmonth,Pmdays)
          TSC_SAVE_TIME2  = IntTime + days2time(real(Pmdays,kind=DP))
       ENDIF

       IF(SSHGRD_ASSIM) SSH_SAVED(1:M,0) = EL(1:M)

       ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
       INT_START = ISTART
       INT_COUNT = ISTART
       IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
          INT_END = ISTART + 2*(IEND-ISTART)
       ELSEIF(TS_NGASSIM .OR. CUR_NGASSIM .OR.      &
            TS_OIASSIM .OR. CUR_OIASSIM)THEN
          INT_END = IEND
       END IF

       !==============================================================================!
       !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
       !==============================================================================!
       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
          Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
          write(IPT,*) "===================================================="
       end if

       CALL REPORT('INITIAL CONDITIONS')

       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
       end if

       IF(TSGRD_ASSIM) THEN
          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "SETUP ASSIMILATION RSTFILE"
          CALL SETUP_RSTFILE_ASSIM
       ENDIF

       DO WHILE(IntTime < EndTime)
          IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
             ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY
          ELSEIF(TS_NGASSIM.OR.CUR_NGASSIM .OR.       &
               TS_OIASSIM.OR.CUR_OIASSIM)THEN
             ASSIM_RESTART_TIME = EndTime
          END IF

          IF(SST_ASSIM.OR.SSTGRD_ASSIM) THEN
             SST_SAVE_INDEX = 0
             SST_SAVE_TIME  = IntTime + sst_save_interval
          ENDIF

          IF(TSGRD_ASSIM) THEN
             TSC_SAVE_INDEX = 0
             TSC_SAVE_TIME = IntTime + TSC_SAVE_INTERVAL 

             CALL Now_2_month_days(IntTime+IMDTI,Pyear,Pmonth,Pmdays)

             IF(ASSIM_RESTART_TIME > TSC_SAVE_TIME2) THEN
                TSC_SAVE_INDEX2 = 1
                TSC_SAVE_TIME2 = TSC_SAVE_TIME2 + days2time(real(Pmdays,kind=DP))
             ELSE
                TSC_SAVE_INDEX2 = 0
             ENDIF
          ENDIF

          IF(SSHGRD_ASSIM) THEN
             SSH_SAVE_INDEX = 0
             SSH_SAVE_TIME  = IntTime + ssh_save_interval
          ENDIF

          CALL SAVE_STATE

          if(DBG_SET(dbg_log)) THEN
             Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
             write(IPT,*) "======= Start 1 Day Simulation  ===================="
          end if

          !==============================================================================!
          !    RUN PURE SIMULATION MODE:
          !==============================================================================!

          ASSIM_FLAG = 0
          DO WHILE(IntTime < ASSIM_RESTART_TIME)

             IINT = IINT + 1
             IntTime=IntTime + IMDTI

             CALL INTERNAL_STEP

             !==============================================================================!
             !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
             !==============================================================================!
             if(DBG_SET(dbg_log)) &
                  & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

             IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

             IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
                ! do nothing
             ELSE
                CALL ARCHIVE

                CALL DUMP_PROBE_DATA
                IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
#            if defined (WAVE_CURRENT_INTERACTION)
                IF(OUT_WAVE_SPARSE_TIMESERIES_ON) CALL OUT_WAVE_SPARSE_TIMESERIES
#            endif


#            if defined (LAG_PARTICLE)
                CALL LAG_UPDATE
#            endif

                IF(NCNEST_ON)      CALL ARCHIVE_NEST     
#            if defined (WAVE_CURRENT_INTERACTION)
                IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#            endif
             END IF

             !==============================================================================!
             !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
             !==============================================================================!
             CALL SHUTDOWN_CHECK(D1)

          END DO

          !==============================================================================!
          !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
          !==============================================================================!
          IF(SST_ASSIM.OR.SSTGRD_ASSIM.OR.SSHGRD_ASSIM.OR.TSGRD_ASSIM)THEN
             if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
             CALL RESTORE_STATE
             CALL RHO_PMEAN

             IF(SST_ASSIM.OR.SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE
             IF(TSGRD_ASSIM) CALL TSGRD_OBSERVATION_UPDATE(IntTime+(days2time(1.0_DP))/2)
             IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

             if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="

             if(DBG_SET(dbg_log)) THEN
                Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
                write(IPT,*) "======= Start 1 Day Assimilation  ===================="
             end if

          ELSEIF(TS_NGASSIM .OR. CUR_NGASSIM .OR. TS_OIASSIM .OR. CUR_OIASSIM)THEN
             IntTime = ASSIM_RESTART_TIME
          END IF
          !==============================================================================!
          !    RUN DATA ASSIMILATION MODE:
          !==============================================================================!
          ASSIM_FLAG = 1
          DO WHILE(IntTime < ASSIM_RESTART_TIME)

             IINT = IINT + 1
             IntTime=IntTime + IMDTI

             CALL INTERNAL_STEP

             !==============================================================================!
             !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
             !==============================================================================!
             if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
             if(DBG_SET(dbg_log)) &
                  & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

             IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

             !==============================================================================!
             !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
             !==============================================================================!
             CALL ARCHIVE

             CALL DUMP_PROBE_DATA
             IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
#          if defined (WAVE_CURRENT_INTERACTION)
             IF(OUT_WAVE_SPARSE_TIMESERIES_ON)CALL OUT_WAVE_SPARSE_TIMESERIES
#          endif
             !==============================================================================!
             !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
             !==============================================================================!
             CALL SHUTDOWN_CHECK(D1)

             !==============================================================================!
             !    LAGRANGIAN PARTICLE TRACKING                                              |
             !==============================================================================!
# if defined (LAG_PARTICLE)
             CALL LAG_UPDATE
# endif

             !==============================================================================!
             !    NESTING OUTPUT                                                            |
             !==============================================================================!
             IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
             IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif

          END DO

          ! RUN THE NEXT DAY

          IF(SSHGRD_ASSIM) SSH_SAVED(1:M,0) = SSH_SAVED(1:M,SSHGRD_N_PER_INTERVAL)

       END DO
# endif


    CASE(FVCOM_RRKF_WITHOUT_SSA)
# if defined (RRKF)     

       IF(MSR) WRITE(IPT,*) 
       IF(MSR) WRITE(IPT,*) 'Starting Kalman filter input data processing......'

       !    IF(MSR) CALL RRK_REF
       CALL RRK_REF
       print *, 'finish rrk_ref'
       IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       print *, 'call setup_rrk_restart'   
       !    IF(LOCAL_DISK)THEN
       !       NCF_RRKRE => SETUP_RESTART(LOCAL_FILE)
       !    ELSE
       NCF_RRKRE => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE))
       !        NCF_RRKRE=>new_file()
       !        NCF_RRKRE%FNAME=TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE)
       !        call NC_LOAD(NCF_RRKRE)
       !        NCF_rrkre%FTIME=>NEW_FTIME()
       !        print *, 'finish nc_load'
       !    END IF

       !    IF(LOCAL_DISK)THEN
       !       NCF_RRKFCT => SETUP_RESTART(LOCAL_FILE1)
       !    ELSE
       NCF_RRKFCT => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE1))
       !    END IF
       print *, 'call rrk_rrk' 
       CALL RRK_RRK(1)    !  perturbation
       !        NCF_RRKRE=>new_file()
       !        NCF_RRKRE%FNAME=TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE)
       !        call NC_LOAD(NCF_RRKRE)
       !        NCF_rrkre%FTIME=>NEW_FTIME()
       print *, 'finish nc_load'
       print *, 'finish rrk_rrki'
       print *, 'before nc open'
       CALL NC_OPEN(NCF_RRKRE)
       !    CALL NC_LOAD(NCF_RRKRE)
       print *, 'finish nc_open'
       DO RRK_EOFCONTR =1, RRK_NEOF
          print *, 'before readrestart' 
          IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          print *, NCF_RRKRE%FNAME
          CALL READRESTART(NCF_RRKRE,RRK_EOFCONTR)
          IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          print *, 'finish readrestart'
          IntTime = REF_START_TIME
          ExtTime = REF_START_TIME

          iint = rkint_start         

          !==============================================================================!
          !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
          !==============================================================================!
          if(DBG_SET(dbg_log)) THEN
             write(IPT,*) "===================================================="
             write(IPT,*) "============   ", 'RRK_NEOF =', RRK_EOFCONTR,"=============="
             write(IPT,*) "===================================================="
             if(DBG_SET(dbg_log)) &
                  & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
             write(IPT,*) "===================================================="
          end if

          CALL REPORT('INITIAL CONDITIONS')

          if(DBG_SET(dbg_log)) THEN
             write(IPT,*) "===================================================="
             write(IPT,*) "===================================================="
             write(IPT,*) "===================================================="
          end if



          !////////////////////////// MAIN LOOP //////////////////////////////////////////
          DO WHILE(IntTime < (REF_START_TIME + RRK_INTERVAL) )

             IntTime=IntTime + IMDTI
             iint = iint +1

             CALL INTERNAL_STEP

             !==============================================================================!
             !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
             !==============================================================================!
             if(DBG_SET(dbg_log)) &
                  & Call REPORT_TIME(IINT,rkInt_start,rkInt_END,IntTime)

             !        IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS: Perturbation Mode')

             !==============================================================================!
             !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
             !==============================================================================!
             CALL SHUTDOWN_CHECK(D1)

          END DO
          !////////////////////////// END MAIN LOOP /////////////////////
          !/////////////////


          CALL WRITERESTART(NCF_RRKFCT,RRK_EOFCONTR)  

       END DO  ! RRK_NEOF
       print *, 'before rrk_rrk(2)'  
       CALL RRK_RRK(2)   ! calculate linearized model matrix
       print *, 'before rrk_rrk(4)'

       IF(MSR) CALL RRK_RRK(4) ! doubling algorithm


       !JQI4000 CONTINUE   

       IF(MSR) WRITE(IPT,*)
       IF(MSR) WRITE(IPT,*) 'Finish Kalman filter preparation! '
       IF(MSR) WRITE(IPT,*) 'Start data assimilation! '
       IF(MSR) WRITE(IPT,*)

       CALL NC_OPEN(NC_START)  


       CALL READRESTART(NC_START,REF_START_TIME)


       IntTime = RRK_START_TIME
       ExtTime = RRK_START_TIME

       iint = rkint_start 


       !==============================================================================!
       !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
       !==============================================================================!
       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "============== STARTING MAIN LOOP AT:==============="
          if(DBG_SET(dbg_log)) &
               & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
          write(IPT,*) "===================================================="
       end if

       CALL REPORT('INITIAL CONDITIONS')

       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
       end if

       ICYC = 1
       RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)

       !////////////////////////// MAIN LOOP //////////////////////////////////////////
       DO WHILE(IntTime < ENDTIME)

          IntTime=IntTime + IMDTI
          iint = iint +1

          CALL INTERNAL_STEP    

          !==============================================================================!
          !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
          !==============================================================================!
          if(DBG_SET(dbg_log)) &
               & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)

          !      IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS: Assimilation Mode')

          !==============================================================================!
          !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
          !==============================================================================!
          !      CALL ARCHIVE

          !      CALL DUMP_PROBE_DATA
          !==============================================================================!
          !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
          !==============================================================================!
          !      CALL SHUTDOWN_CHECK(D1)

          !==============================================================================!
          !    LAGRANGIAN PARTICLE TRACKING                                              |
          !==============================================================================!
          !      CALL LAG_UPDATE

          !==============================================================================!
          !    NESTING OUTPUT                                                            |
          !==============================================================================!
          !      IF(NESTING_ON .and. NESTING_MODE == "large_domain") CALL ARCHIVE_NEST

          !==============================================================================!
          !    KALMAN FILTERING                                                          |
          !==============================================================================!
          IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          if (msr) print *, 'finish the filter_obs_data'
          IF(IntTime>=RRK_START_TIME .AND. IntTime<=RRK_END_TIME) THEN
             IF(IntTime == RRK_CYC) THEN
                CALL PRINT_TIME(IntTime,IPT,"DOING RRK_RRKF")
#if defined (MULTIPROCESSOR)
                IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif 
                CALL filter_obs_data
#if defined (MULTIPROCESSOR)
                IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif
                CALL set_rrkf_assim_data
                print *, 'kf observation location is:', nloc
                if (nloc>0) then
                   print *, 'rrk_rrkf'
                   CALL RRK_RRKF
                end if
                call deallocate_obs_data
                IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
                print *, 'finish rrkf'
                ICYC = ICYC+1
                RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)	    
             END IF
          END IF


       END DO
       !////////////////////////// END MAIN LOOP ////////////////////////////////////// 


       ! #####################################################################################

# endif
    CASE(FVCOM_RRKF_WITH_SSA)
# if defined (RRKF)
# if defined (DATA_ASSIM) 
       IF (SSTGRD_ASSIM .or. SSHGRD_ASSIM) THEN
       ELSE
          PRINT *, 'STOP, RRKF ONLY WITH SSTGRD_ASSIM OR SSHGRD_ASSIM'
          CALL PSTOP
       END IF
       IF(MSR) WRITE(IPT,*) 'Starting RRKF input data processing......'
       CALL RRK_REF
       print *, 'finish rrk_ref'
       IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       print *, 'call setup_rrk_restart'   

       NCF_RRKRE => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE))

       NCF_RRKFCT => SETUP_RESTART(TRIM(OUTPUT_DIR)//TRIM(GROUP_FILE1))
       print *, 'call rrk_rrk' 
       CALL RRK_RRK(1)    !  perturbation


       CALL NC_OPEN(NCF_RRKRE)
       print *, 'finish nc_open'
       DO RRK_EOFCONTR =1, RRK_NEOF
          print *, 'before readrestart' 
          IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          print *, NCF_RRKRE%FNAME
          CALL READRESTART(NCF_RRKRE,RRK_EOFCONTR)
          IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          print *, 'finish readrestart'
          IntTime = REF_START_TIME
          ExtTime = REF_START_TIME

          iint = rkint_start         

          !==============================================================================!
          !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
          !==============================================================================!
          if(DBG_SET(dbg_log)) THEN
             write(IPT,*) "===================================================="
             write(IPT,*) "============   ", 'RRK_NEOF =', RRK_EOFCONTR,"=============="
             write(IPT,*) "===================================================="
             if(DBG_SET(dbg_log)) &
                  & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
             write(IPT,*) "===================================================="
          end if

          CALL REPORT('INITIAL CONDITIONS')

          if(DBG_SET(dbg_log)) THEN
             write(IPT,*) "===================================================="
             write(IPT,*) "===================================================="
             write(IPT,*) "===================================================="
          end if
          IF (starttime%MuSOD ==0 )  THEN
             IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
             CALL ALLOC_BUFFER

             ! SET THE ASSIMILATION/SIMULATION RESET TIME
             ASSIM_RESTART_TIME = StartTime


             ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
             INT_START = ISTART
             INT_COUNT = ISTART
             INT_END   = ISTART + 2*(IEND-ISTART)

             !==============================================================================!
             !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
             !==============================================================================!
             if(DBG_SET(dbg_log)) THEN
                write(IPT,*) "===================================================="
                write(IPT,*) "===================================================="
                write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
                Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
                write(IPT,*) "===================================================="
             end if

             CALL REPORT('INITIAL CONDITIONS')

             if(DBG_SET(dbg_log)) THEN
                write(IPT,*) "===================================================="
                write(IPT,*) "===================================================="
                write(IPT,*) "===================================================="
             end if

          end if  !(starttime%MuSOD==0)

          DO WHILE(IntTime < EndTime)
             IF (starttime%MuSOD ==0 )  THEN
                IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME

                ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY

                IF(SSTGRD_ASSIM) THEN
                   SST_SAVE_INDEX = 0
                   SST_SAVE_TIME  = IntTime + sst_save_interval
                ENDIF


                IF(SSHGRD_ASSIM) THEN
                   SSH_SAVE_INDEX = 0
                   SSH_SAVE_TIME  = IntTime + ssh_save_interval
                ENDIF

                CALL SAVE_STATE

                if(DBG_SET(dbg_log)) THEN
                   Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
                   write(IPT,*) "======= Start 1 Day Simulation  ===================="
                end if

                !==============================================================================!
                !    RUN PURE SIMULATION MODE:
                !==============================================================================!

                ASSIM_FLAG = 0

                IF(SSHGRD_ASSIM) THEN
                   CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
                ENDIF
                IF(SSTGRD_ASSIM) THEN
                   CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
                ENDIF
                RECALC_RHO_MEAN%mjd = Inttime%mjd
                RECALC_RHO_MEAN%musod=0
                DO WHILE(IntTime < ASSIM_RESTART_TIME)

                   IINT = IINT + 1
                   IntTime=IntTime + IMDTI

                   CALL INTERNAL_STEP

                   !==============================================================================!
                   !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
                   !==============================================================================!
                   if(DBG_SET(dbg_log)) &
                        & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

                   IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

                   !==============================================================================!
                   !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
                   !==============================================================================!
                   CALL SHUTDOWN_CHECK(D1)

                END DO

                !==============================================================================!
                !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
                !==============================================================================!
                if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
                CALL RESTORE_STATE
                CALL RHO_PMEAN

                IF(SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE

                IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

                if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="
             ELSE  !(starttime%MuSOD  .NOT 0 )
                IF (MSR) PRINT *, 'DO NOT NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
             END IF
             if(DBG_SET(dbg_log)) THEN
                Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
                write(IPT,*) "======= Start 1 Day Assimilation  ===================="
             end if

             !==============================================================================!
             !    RUN DATA ASSIMILATION MODE:
             !==============================================================================!
             ASSIM_FLAG = 1
             IF(SSHGRD_ASSIM) THEN
                CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
             ENDIF
             IF(SSTGRD_ASSIM) THEN
                CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
             ENDIF
             RECALC_RHO_MEAN%mjd = Inttime%mjd
             RECALC_RHO_MEAN%musod=0

             DO WHILE(IntTime < ASSIM_RESTART_TIME)

                IINT = IINT + 1
                IntTime=IntTime + IMDTI

                CALL INTERNAL_STEP

                !==============================================================================!
                !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
                !==============================================================================!
                if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
                if(DBG_SET(dbg_log)) &
                     & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

                IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')


                !==============================================================================!
                !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
                !==============================================================================!
                CALL SHUTDOWN_CHECK(D1)

             END DO ! ASSIMILATION STAGE LOOP


          END DO ! ASSIMILATION TOTAL DAY CONTROL LOOP

          CALL WRITERESTART(NCF_RRKFCT,RRK_EOFCONTR)  

       END DO  ! RRK_NEOF


       print *, 'before rrk_rrk(2)'  
       CALL RRK_RRK(2)   ! calculate linearized model matrix
       print *, 'before rrk_rrk(4)'

       IF(MSR) CALL RRK_RRK(4) ! doubling algorithm


       !JQI4000 CONTINUE   

       IF(MSR) WRITE(IPT,*)
       IF(MSR) WRITE(IPT,*) 'Finish Kalman filter preparation! '
       IF(MSR) WRITE(IPT,*) 'Start data assimilation! '
       IF(MSR) WRITE(IPT,*)

       CALL NC_OPEN(NC_START)  


       CALL READRESTART(NC_START,REF_START_TIME)


       IntTime = RRK_START_TIME
       ExtTime = RRK_START_TIME

       iint = rkint_start 


       !==============================================================================!
       !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
       !==============================================================================!
       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "============== STARTING MAIN LOOP AT:==============="
          if(DBG_SET(dbg_log)) &
               & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
          write(IPT,*) "===================================================="
       end if

       CALL REPORT('INITIAL CONDITIONS')

       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
       end if

       ICYC = 1
       RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)
       CALL ALLOC_BUFFER

       ! SET THE ASSIMILATION/SIMULATION RESET TIME
       ASSIM_RESTART_TIME = StartTime
       IF(TSGRD_ASSIM) THEN
          CALL Now_2_month_days(StartTime,Pyear,Pmonth,Pmdays)
          TSC_SAVE_TIME2  = IntTime + days2time(real(Pmdays,kind=DP))
       ENDIF

       ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
       INT_START = ISTART
       INT_COUNT = ISTART
       INT_END   = ISTART + 2*(IEND-ISTART)

       !==============================================================================!
       !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
       !==============================================================================!
       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
          Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
          write(IPT,*) "===================================================="
       end if

       CALL REPORT('INITIAL CONDITIONS')

       if(DBG_SET(dbg_log)) THEN
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
          write(IPT,*) "===================================================="
       end if

       IF(TSGRD_ASSIM) THEN
          IF(DBG_SET(DBG_LOG)) WRITE(IPT,*) "SETUP ASSIMILATION RSTFILE"
          CALL SETUP_RSTFILE_ASSIM
       ENDIF

       DO WHILE(IntTime < EndTime)

          ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY

          IF(SSTGRD_ASSIM) THEN
             SST_SAVE_INDEX = 0
             SST_SAVE_TIME  = IntTime + sst_save_interval
          ENDIF

          IF(TSGRD_ASSIM) THEN
             TSC_SAVE_INDEX = 0
             TSC_SAVE_TIME = IntTime + TSC_SAVE_INTERVAL 

             CALL Now_2_month_days(IntTime+IMDTI,Pyear,Pmonth,Pmdays)

             IF(ASSIM_RESTART_TIME > TSC_SAVE_TIME2) THEN
                TSC_SAVE_INDEX2 = 1
                TSC_SAVE_TIME2 = TSC_SAVE_TIME2 + days2time(real(Pmdays,kind=DP))
             ELSE
                TSC_SAVE_INDEX2 = 0
             ENDIF
          ENDIF

          IF(SSHGRD_ASSIM) THEN
             SSH_SAVE_INDEX = 0
             SSH_SAVE_TIME  = IntTime + ssh_save_interval
          ENDIF

          CALL SAVE_STATE

          if(DBG_SET(dbg_log)) THEN
             Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
             write(IPT,*) "======= Start 1 Day Simulation  ===================="
          end if

          !==============================================================================!
          !    RUN PURE SIMULATION MODE:
          !==============================================================================!

          ASSIM_FLAG = 0

          IF(SSHGRD_ASSIM) THEN
             CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
          ENDIF
          IF(SSTGRD_ASSIM) THEN
             CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
          ENDIF
          RECALC_RHO_MEAN%mjd = Inttime%mjd
          RECALC_RHO_MEAN%musod=0
          DO WHILE(IntTime < ASSIM_RESTART_TIME)

             IINT = IINT + 1
             IntTime=IntTime + IMDTI

             CALL INTERNAL_STEP

             !==============================================================================!
             !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
             !==============================================================================!
             if(DBG_SET(dbg_log)) &
                  & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

             IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

             !==============================================================================!
             !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
             !==============================================================================!
             CALL SHUTDOWN_CHECK(D1)

          END DO

          !==============================================================================!
          !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
          !==============================================================================!
          if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
          CALL RESTORE_STATE
          CALL RHO_PMEAN

          IF(SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE
          IF(TSGRD_ASSIM) CALL TSGRD_OBSERVATION_UPDATE(IntTime+(days2time(1.0_DP))/2)
          IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

          if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="

          if(DBG_SET(dbg_log)) THEN
             Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
             write(IPT,*) "======= Start 1 Day Assimilation  ===================="
          end if

          !==============================================================================!
          !    RUN DATA ASSIMILATION MODE:
          !==============================================================================!
          ASSIM_FLAG = 1
          IF(SSHGRD_ASSIM) THEN
             CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
          ENDIF
          IF(SSTGRD_ASSIM) THEN
             CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
          ENDIF
          RECALC_RHO_MEAN%mjd = Inttime%mjd
          RECALC_RHO_MEAN%musod=0

          DO WHILE(IntTime < ASSIM_RESTART_TIME)

             IINT = IINT + 1
             IntTime=IntTime + IMDTI

             CALL INTERNAL_STEP

             !==============================================================================!
             !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
             !==============================================================================!
             if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
             if(DBG_SET(dbg_log)) &
                  & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

             IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

             !==============================================================================!
             !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
             !==============================================================================!
             !           CALL ARCHIVE

             !           CALL DUMP_PROBE_DATA
             !==============================================================================!
             !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
             !==============================================================================!
             CALL SHUTDOWN_CHECK(D1)

             !==============================================================================!
             !    LAGRANGIAN PARTICLE TRACKING                                              |
             !==============================================================================!
             !# if defined (LAG_PARTICLE)
             !           CALL LAG_UPDATE
             !# endif

             !==============================================================================!
             !    NESTING OUTPUT                                                            |
             !==============================================================================!
             !        IF(NCNEST_ON) CALL ARCHIVE_NEST

             !==============================================================================!
             !    KALMAN FILTERING                                                          |
             !==============================================================================!
             IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
             if (msr) print *, 'finish the filter_obs_data'
             IF(IntTime>=RRK_START_TIME .AND. IntTime<=RRK_END_TIME) THEN
                IF(IntTime == RRK_CYC) THEN
                   CALL PRINT_TIME(IntTime,IPT,"DOING RRK_RRKF")
#if defined (MULTIPROCESSOR)
                   IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif 
                   CALL filter_obs_data
#if defined (MULTIPROCESSOR)
                   IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif
                   CALL set_rrkf_assim_data
                   print *, 'kf observation location is:', nloc
                   if (nloc>0) then
                      print *, 'rrk_rrkf'
                      CALL RRK_RRKF
                   end if
                   call deallocate_obs_data
                   IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
                   print *, 'finish rrkf'
                   ICYC = ICYC+1
                   RRK_CYC = RRK_START_TIME + (RRK_INTERVAL*ICYC)	    
                END IF
             END IF


          END DO

          ! RUN THE NEXT DAY

       END DO
       !qxu{
# endif
# endif  

       !============================================END RRKF_WITH_SSA================================| 





    CASE(FVCOM_ENKF_WITHOUT_SSA)
# if defined (ENKF)
       IF(MSR) WRITE(IPT,*) 'Starting Kalman filter input data processing......'
       ICYC = 1
       enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
       DO WHILE(IntTime < ENDTIME) ! biggest control runtime
          if(msr)   print *,'start a forecast cycle'

          kf_times=kf_times+1 ! 


          DO ENKF_memberCONTR =1, ENKF_NENS 

             CALL ENKF_RESTART
             if(DBG_SET(dbg_log)) THEN
                write(IPT,*) "===================================================="
                write(IPT,*) "============   ", 'ENKF_NENS =', ENKF_memberCONTR,"=============="
                write(IPT,*) "===================================================="
                if(DBG_SET(dbg_log)) &
                     & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
                write(IPT,*) "===================================================="
             end if
             !     istart=iint+1
             !     NSTEPS = CALCULATE_NUMBER_OF_TIMESTEPS(StartTime,EndTime)
             !     IEND = ISTART + NSTEPS
             CALL REPORT('INITIAL CONDITIONS')

             if(DBG_SET(dbg_log)) THEN
                write(IPT,*) "===================================================="
                write(IPT,*) "===================================================="
                write(IPT,*) "===================================================="
             end if
             !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
             !////////////////////////// MAIN LOOP //////////////////////////////////////////
             DO while (IntTime < EndTime)

                IntTime=IntTime + IMDTI
                IINT = IINT + 1
                CALL INTERNAL_STEP

                !==============================================================================!
                !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
                !==============================================================================!
                if(DBG_SET(dbg_log)) &
                     & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)

                IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

                !==============================================================================!
                !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
                !==============================================================================!

                !for Eulerian velocity output
# if defined (WAVE_CURRENT_INTERACTION)
                !U = U-U_STOKES_3D
                !V = V-V_STOKES_3D
                !UA=UA-U_STOKES_2D
                !VA=VA-V_STOKES_2D
# endif


                IF (ENKF_memberCONTR ==1 ) THEN
                   if (msr) print *, "call archive for ensemble member: ", ENKF_memberCONTR
                   CALL ARCHIVE

                END IF

                IF(NC_ON)THEN
                   ! bounds checking
                   !IF(NC_DAT%FTIME%NEXT_IO == IntTime .or. FORCE_ARCHIVE) THEN
                   if (inttime > enkf_out_Time+0.1_SP*IMDTI) then
                      print *, 'inttime',inttime
                      print *, 'enkf_out_time',enkf_out_time
                      Call Fatal_Error("The time of the output &
                           &must be greater than or equal to the current time")
                   end if
                   IF(abs(enkf_out_Time -IntTime)<0.1_SP*IMDTI) THEN
                      call    ENKF_SETUP_NC(NC_DAT)
                      IF (ENKF_memberCONTR == ENKF_NENS) then
                         ARCHIVE_TIMES=ARCHIVE_TIMES+1
                         enkf_out_Time=enkf_out_Time+NC_DAT%FTIME%INTERVAL
                      end if

                   END IF
                END IF


#  if defined (WAVE_CURRENT_INTERACTION)
                !Convert Eulerian velocity back to Lagrangian velocity to maintain consistent with Mellor's equation
                !U = U+U_STOKES_3D
                !V = V+V_STOKES_3D
                !UA=UA+U_STOKES_2D
                !VA=VA+V_STOKES_2D
#  endif

                CALL DUMP_PROBE_DATA
                IF(OUT_STATION_TIMESERIES_ON)CALL OUT_STATION_TIMESERIES
# if defined (WAVE_CURRENT_INTERACTION)
                IF(OUT_WAVE_SPARSE_TIMESERIES_ON)CALL OUT_WAVE_SPARSE_TIMESERIES
# endif
                !==============================================================================!
                !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
                !==============================================================================!
                CALL SHUTDOWN_CHECK(D1)

                !==============================================================================!
                !  CALL BOUNDS CHECK TO SEE IF VARIABLES EXCEED USER-DEFINED THRESHOLDS 
                !==============================================================================!
                CALL BOUNDSCHK  !bounds checking

                !==============================================================================!
                !    LAGRANGIAN PARTICLE TRACKING                                              |
                !==============================================================================!
# if defined (LAG_PARTICLE)
                CALL LAG_UPDATE
# endif

                !==============================================================================!
                !    NESTING OUTPUT                                                            |
                !==============================================================================!
                IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
                IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif
                if (msr ) print *,  'inttime',inttime, 'enkf_cyc',enkf_cyc
                IF(IntTime == enkf_cyc) EXIT ! stop run when KF launch , exit from assimilation stage
             END DO
             !////////////////////////// END MAIN LOOP //////////////////////////////////////
             write(enkf_num,'(i2.2)')  ENKF_memberCONTR
             NCF_ENKFfct => SETUP_RESTART(TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_"//enkf_num//".nc")


             CALL ENKF_WRITERESTART(NCF_ENKFfct)
          end do ! enkf_nens
          IF(IntTime == enkf_cyc) THEN
             CALL PRINT_TIME(IntTime,IPT,"DOING ENKF")
             IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
             if (enkf_test) then
                call getnloc
             else 
                CALL filter_obs_data
                IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
                if (msr) print *, 'finish the filter_obs_data'
                CALL set_enkf_assim_data
             end if
             !           IF (MSR) THEN
             print *, 'kf observation location is:', nloc
             if (nloc>0) then
                CALL ENKF_ASS_EVE
             end if
             !           END IF
             IF (MSR) THEN
                if (ICYC==1 .AND. NLOC==0) THEN
                   open(435,file='idobs.dat')
                   write(435,*) -31,-711
                   close(435)

                END IF
                PRINT *, 'FINISH ENKF, ',inttime
             END IF
             if (.not. enkf_test) then
                call deallocate_obs_data
             end if
             IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
             ICYC = ICYC+1
             enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
          END IF

       END DO  ! end biggest run control


       IF(MSR) WRITE(IPT,*)
       IF(MSR) WRITE(IPT,*) 'Finish Ensemble Kalman filter! '
       IF(MSR) WRITE(IPT,*)
       IF(MSR) WRITE(IPT,*)

# endif

       !============================================END ENKF_WITHOUT_SSA================================|     

    CASE(FVCOM_ENKF_WITH_SSA)
# if defined (ENKF)
# if defined (DATA_ASSIM) 
       IF (SSTGRD_ASSIM .or. SSHGRD_ASSIM) THEN
       ELSE
          PRINT *, 'STOP, ENKF ONLY WITH SSTGRD_ASSIM OR SSHGRD_ASSIM'
          CALL PSTOP
       END IF
       IF(MSR) WRITE(IPT,*) 'Starting Kalman filter input data processing......'
       ICYC = 1
       enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
       if(msr)  print *, "inttime",inttime
       if(msr)  print *,  "endtime",endtime
       DO WHILE(IntTime < ENDTIME) ! biggest control runtime
          if(msr)   print *,'start a forecast cycle'

          kf_times=kf_times+1 ! 


          DO ENKF_memberCONTR =1, ENKF_NENS 

             CALL ENKF_RESTART
             if(DBG_SET(dbg_log)) THEN
                write(IPT,*) "===================================================="
                write(IPT,*) "============   ", 'ENKF_NENS =', ENKF_memberCONTR,"=============="
                write(IPT,*) "===================================================="
                if(DBG_SET(dbg_log)) &
                     & Call REPORT_TIME(IINT,ISTART,IEND,IntTime)
                write(IPT,*) "===================================================="
             end if

             CALL REPORT('INITIAL CONDITIONS')

             if(DBG_SET(dbg_log)) THEN
                write(IPT,*) "===================================================="
                write(IPT,*) "===================================================="
                write(IPT,*) "===================================================="
             end if
             !============================SSH AND SST ASSIMILATION CODE BELOW=============== 
             IF (starttime%MuSOD ==0 )  THEN
                IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
                CALL ALLOC_BUFFER

                ! SET THE ASSIMILATION/SIMULATION RESET TIME
                ASSIM_RESTART_TIME = StartTime


                ! INITIALIZE THE COUNT VARIABLES FOR THE SST LOOP
                INT_START = ISTART
                INT_COUNT = ISTART
                INT_END   = ISTART + 2*(IEND-ISTART)

                !==============================================================================!
                !  PREPARE TO START FVCOM'S MAIN LOOP                                          !
                !==============================================================================!
                if(DBG_SET(dbg_log)) THEN
                   write(IPT,*) "===================================================="
                   write(IPT,*) "===================================================="
                   write(IPT,*) "======STARTING MAIN LOOP ASSIMILATION MODE:========="
                   Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
                   write(IPT,*) "===================================================="
                end if

                CALL REPORT('INITIAL CONDITIONS')

                if(DBG_SET(dbg_log)) THEN
                   write(IPT,*) "===================================================="
                   write(IPT,*) "===================================================="
                   write(IPT,*) "===================================================="
                end if

             end if  !(starttime%MuSOD==0)

             DO WHILE(IntTime < EndTime)
                IF (starttime%MuSOD ==0 )  THEN
                   IF (MSR) PRINT *, 'NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME

                   ASSIM_RESTART_TIME = ASSIM_RESTART_TIME + days2time(1.0_DP)  ! ADD ONE MORE DAY

                   IF(SSTGRD_ASSIM) THEN
                      SST_SAVE_INDEX = 0
                      SST_SAVE_TIME  = IntTime + sst_save_interval
                   ENDIF


                   IF(SSHGRD_ASSIM) THEN
                      SSH_SAVE_INDEX = 0
                      SSH_SAVE_TIME  = IntTime + ssh_save_interval
                   ENDIF

                   CALL SAVE_STATE

                   if(DBG_SET(dbg_log)) THEN
                      Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
                      write(IPT,*) "======= Start 1 Day Simulation  ===================="
                   end if

                   !==============================================================================!
                   !    RUN PURE SIMULATION MODE:
                   !==============================================================================!

                   ASSIM_FLAG = 0

                   IF(SSHGRD_ASSIM) THEN
                      CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
                   ENDIF
                   IF(SSTGRD_ASSIM) THEN
                      CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
                   ENDIF
                   RECALC_RHO_MEAN%mjd = Inttime%mjd
                   RECALC_RHO_MEAN%musod=0
                   DO WHILE(IntTime < ASSIM_RESTART_TIME)

                      IINT = IINT + 1
                      IntTime=IntTime + IMDTI

                      CALL INTERNAL_STEP

                      !==============================================================================!
                      !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
                      !==============================================================================!
                      if(DBG_SET(dbg_log)) &
                           & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

                      IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

                      !==============================================================================!
                      !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
                      !==============================================================================!
                      CALL SHUTDOWN_CHECK(D1)

                   END DO

                   !==============================================================================!
                   !    CALL RESTORE TO RUN ASSIMILATION CODE                                     |
                   !==============================================================================!
                   if(DBG_SET(dbg_log)) write(IPT,*) "=====Start restore the state for assimilation===="
                   CALL RESTORE_STATE
                   CALL RHO_PMEAN

                   IF(SSTGRD_ASSIM) CALL SSTGRD_OBSERVATION_UPDATE

                   IF(SSHGRD_ASSIM) CALL SSHGRD_OBSERVATION_UPDATE

                   if(DBG_SET(dbg_log)) write(IPT,*) "======finish update the obs state  ============="
                ELSE  !(starttime%MuSOD  .NOT 0 )
                   IF (MSR) PRINT *, 'DO NOT NEED SIMULATION STAGE BECAUSE STARTTIME IS',STARTTIME
                END IF
                if(DBG_SET(dbg_log)) THEN
                   Call REPORT_TIME(INT_COUNT,INT_START,INT_END,IntTime)
                   write(IPT,*) "======= Start 1 Day Assimilation  ===================="
                end if

                !==============================================================================!
                !    RUN DATA ASSIMILATION MODE:
                !==============================================================================!
                ASSIM_FLAG = 1
                IF(SSHGRD_ASSIM) THEN
                   CALL SET_SSH_SAVE_TIME_INDEX(ASSIM_FLAG)
                ENDIF
                IF(SSTGRD_ASSIM) THEN
                   CALL SET_SST_SAVE_TIME_INDEX(ASSIM_FLAG)
                ENDIF
                RECALC_RHO_MEAN%mjd = Inttime%mjd
                RECALC_RHO_MEAN%musod=0

                DO WHILE(IntTime < ASSIM_RESTART_TIME)

                   IINT = IINT + 1
                   IntTime=IntTime + IMDTI

                   CALL INTERNAL_STEP

                   !==============================================================================!
                   !    OUTPUT SCREEN REPORT/TIME SERIES DATA/OUTPUT FILES                        |
                   !==============================================================================!
                   if(DBG_SET(dbg_log)) write(IPT,*) "=======  Assimilation  ================"
                   if(DBG_SET(dbg_log)) &
                        & Call REPORT_TIME(IINT,INT_START,INT_END,IntTime)

                   IF(REPORT_NOW(IINT,IREPORT)) CALL REPORT('FLOW FIELD STATS')

                   !==============================================================================!
                   !  CALL ARCHIVE TO WRITE THE OUTPUT (SELECTED BASED ON INTTIME)                !
                   !==============================================================================!
                   IF (ENKF_memberCONTR ==1 ) THEN
                      if (msr) print *, "call archive for ensemble member: ", ENKF_memberCONTR
                      CALL ARCHIVE

                   END IF
                   !==============================================================================!
                   !  CALL SHUTDOWN CHECK TO LOOK FOR BAD VALUES                                  !
                   !==============================================================================!
                   CALL SHUTDOWN_CHECK(D1)

                   !==============================================================================!
                   !    LAGRANGIAN PARTICLE TRACKING                                              |
                   !==============================================================================!
# if defined (LAG_PARTICLE)
                   CALL LAG_UPDATE
# endif

                   !==============================================================================!
                   !    NESTING OUTPUT                                                            |
                   !==============================================================================!
                   IF(NCNEST_ON)      CALL ARCHIVE_NEST
#       if defined (WAVE_CURRENT_INTERACTION)
                   IF(NCNEST_ON_WAVE) CALL ARCHIVE_NEST_WAVE
#       endif	
                   if (msr ) print *,  'inttime',inttime, 'enkf_cyc',enkf_cyc
                   IF(IntTime == enkf_cyc) EXIT ! stop run when KF launch , exit from assimilation stage
                END DO ! ASSIMILATION STAGE LOOP
                IF(IntTime == enkf_cyc) EXIT ! not run one day, exit from assimilation loop when KF should laungh


             END DO ! ASSIMILATION TOTAL DAY CONTROL LOOP
             write(enkf_num,'(i2.2)')  ENKF_memberCONTR
             NCF_ENKFfct => SETUP_RESTART(TRIM(OUTPUT_DIR)//"enkftemp/enkf_forecast_"//enkf_num//".nc")


             CALL ENKF_WRITERESTART(NCF_ENKFfct)
          end do ! enkf_nens
          IF(IntTime == enkf_cyc) THEN

             CALL PRINT_TIME(IntTime,IPT,"DOING ENKF")
             IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
             !               call getobsloc
             CALL filter_obs_data
             IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
             if (msr) print *, 'finish the filter_obs_data'
             CALL set_enkf_assim_data
             IF (MSR) THEN
                print *, 'kf observation location is:', nloc
                if (nloc>0) then
                   CALL ENKF_ASS_EVE
                end if
                PRINT *, 'FINISH ENKF, ',inttime
             END IF
             call deallocate_obs_data
             IF(PAR)CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
             ICYC = ICYC+1
             enkf_cyc = ENKF_START_TIME + (ENKF_INTERVAL*ICYC)
          END IF

       END DO  ! end biggest run control


       IF(MSR) WRITE(IPT,*)
       IF(MSR) WRITE(IPT,*) 'Finish Ensemble Kalman filter! '
       IF(MSR) WRITE(IPT,*)
       IF(MSR) WRITE(IPT,*)

# endif
# endif

       !============================================END ENKF_WITH_SSA================================|     
    CASE DEFAULT
       CALL FATAL_ERROR("UNKNOWN FVCOM_RUN_MODE :'"//TRIM(FVCOM_RUN_MODE),&
            & "Options are the following: '"//TRIM(FVCOM_PURE_SIM)//"&
            &' OR '"//TRIM(FVCOM_NUDGE_OI_ASSIM) )
    END SELECT


    RETURN
  END SUBROUTINE FVCOM_run

  !=========================================================================!
  !                        FINALIZE                                          !
  !=========================================================================!

  SUBROUTINE    FVCOM_finalize
#  if defined (SEMI_IMPLICIT) || defined (NH) || (defined (WAVE_CURRENT_INTERACTION) && !defined (EXPLICIT))
     CALL PETSc_CLEANUP
#  endif


#  if defined (DATA_ASSIM)   
    IF(TSGRD_ASSIM) THEN
       TSC_SAVE_INDEX2 = 1
       CALL SAVE_STATE  !! save the information netcdf file
    ENDIF
#  endif

#  if defined (WAVE_CURRENT_INTERACTION) && !defined (WAVE_OFFLINE)
     CALL SWMAIN_CLOSE
#  endif


    if(DBG_SET(dbg_log)) write(IPT,*)"TADA!"
    CALL PSHUTDOWN

    RETURN
  END SUBROUTINE FVCOM_finalize

END MODULE mod_ocean_control
