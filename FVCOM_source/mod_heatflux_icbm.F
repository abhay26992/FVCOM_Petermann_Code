!/===========================================================================/
! Copyright (c) 2007, The University of Massachusetts Dartmouth 
! Produced at the School of Marine Science & Technology 
! Marine Ecosystem Dynamics Modeling group
! All rights reserved.
!
! FVCOM has been developed by the joint UMASSD-WHOI research team. For 
! details of authorship and attribution of credit please see the FVCOM
! technical manual or contact the MEDM group.
!
! 
! This file is part of FVCOM. For details, see http://fvcom.smast.umassd.edu 
! The full copyright notice is contained in the file COPYRIGHT located in the 
! root directory of the FVCOM code. This original header must be maintained
! in all distributed versions.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED.  
!
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

!==============================================================================|
!   CONTROL VARIABLES                                                          |
!==============================================================================|

# if defined (HEATING_CALCULATED)
MODULE MOD_HEATFLUX

   USE ALL_VARS
   USE MOD_UTILS
   USE MOD_PREC

   IMPLICIT NONE
   SAVE

   LOGICAL HEATING_CALCULATE_ON
   CHARACTER(LEN=80) HEATING_CALCULATE_TYPE
   CHARACTER(LEN=80) HEATING_CALCULATE_FILE
   CHARACTER(LEN=80) HEATING_CALCULATE_KIND
   REAL(SP) :: ZUU                            !!HEIGHT OF WIND SPEED (M) 
   REAL(SP) :: ZTT                            !!HEIGHT OF AIR TEMPERATURE (M)
   REAL(SP) :: ZQQ                            !!HEIGHT OF RELATIVE HUMIDITY (M)
   REAL(SP) :: AIR_TEMPERATURE
   REAL(SP) :: RELATIVE_HUMIDITY
   REAL(SP) :: SURFACE_PRESSURE
   REAL(SP) :: LONGWAVE_RADIATION
   REAL(SP) :: SHORTWAVE_RADIATION
   REAL(SP) :: HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX
   REAL(SP) :: HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX
   REAL(SP) :: HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX
   
   NAMELIST /NML_HEATING_CALCULATED/           &
        & HEATING_CALCULATE_ON,                &
	& HEATING_CALCULATE_TYPE,              &
	& HEATING_CALCULATE_FILE,              &
	& HEATING_CALCULATE_KIND,              &
	& ZUU,                                 &
	& ZTT,                                 &
	& ZQQ,                                 &
	& AIR_TEMPERATURE,                     &
	& RELATIVE_HUMIDITY,                   &
	& SURFACE_PRESSURE,                    &
	& LONGWAVE_RADIATION,                  &
	& SHORTWAVE_RADIATION,                 &
        & HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX,    &
        & HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX, &
        & HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX
   
   REAL(SP), ALLOCATABLE :: CORRG(:),CORR(:)           !!LATITUDE OF NODES
   REAL(SP), PARAMETER :: emmiss = 0.97_SP           ! ejw 8/16/2006 based on ROMS
   REAL(SP), PARAMETER :: StefBo = 5.67e-8_SP        ! ejw 8/16/2006 based on ROMS

!==============================================================================|


   CONTAINS !---------------------------------------------------------------------|
            ! HEATING_CALCULATE_NAMELIST_INITIALIZE : Initialize the values in    |
	    !                                     namelist NML_HEATING_CALCULATED |
            ! HEATING_CALCULATE_NAMELIST_PRINT      : Print the values of namelist|
	    !                                         NML_HEATING_CALCULATED      |
            ! HEATING_CALCULATE_NAMELIST_READ       : Read the values of namelist |
	    !                                         NML_HEATING_CALCULATED      |
	    ! COARE26Z     : Calculates Sensible Heat Flux and Latent Heat Flux   |
            ! PSIT_26      : Computes Temperature Structure Function              |
	    ! PSIU_26      : Computes Velocity Structure Function                 |
	    ! QSAT26       : Computes Saturation Specific Humidity                |
	    ! GRV          : Computes g Given Lat in Deg                          |
	    
!==============================================================================|
!   Input Parameters Which Control the Calculation of Heat Flux                |
!==============================================================================|

   SUBROUTINE HEATING_CALCULATE_NAMELIST_INITIALIZE
   USE control, only : casename
   IMPLICIT NONE
   
   HEATING_CALCULATE_ON   = .FALSE.                
   HEATING_CALCULATE_TYPE = "'flux' or 'body'"              
   HEATING_CALCULATE_FILE = trim(casename)//"_hfx.nc"              
   HEATING_CALCULATE_KIND = "Options:"//TRIM(CNSTNT)//","//TRIM(STTC)//","//TRIM(TMDPNDNT)//","//TRIM(PRDC)//","//TRIM(VRBL)              
   ZUU                    = 2  ! Unit = m                                 
   ZTT                    = 2  ! Unit = m                                 
   ZQQ                    = 2  ! Unit = m                                 
   AIR_TEMPERATURE        = 0.0_SP                     
   RELATIVE_HUMIDITY      = 0.0_SP
   SURFACE_PRESSURE       = 0.0_SP                    
   LONGWAVE_RADIATION     = 0.0_SP                  
   SHORTWAVE_RADIATION    = 0.0_SP
   HEATING_LONGWAVE_PERCTAGE_IN_HEATFLUX = 0.78_SP
   HEATING_LONGWAVE_LENGTHSCALE_IN_HEATFLUX = 1.4_SP
   HEATING_SHORTWAVE_LENGTHSCALE_IN_HEATFLUX= 6.3_SP
   
   RETURN
   END SUBROUTINE HEATING_CALCULATE_NAMELIST_INITIALIZE
   
!------------------------------------------------------------------------------|
   SUBROUTINE HEATING_CALCULATE_NAMELIST_PRINT
   USE CONTROL, ONLY : IPT
   
   IMPLICIT NONE
   
   WRITE(UNIT=IPT,NML=NML_HEATING_CALCULATED)
   
   RETURN
   END SUBROUTINE HEATING_CALCULATE_NAMELIST_PRINT  
   
!------------------------------------------------------------------------------|
   SUBROUTINE HEATING_CALCULATE_NAMELIST_READ    
   USE CONTROL, ONLY : casename,NMLUNIT
   
   IMPLICIT NONE
   
   INTEGER :: IOS, I
   CHARACTER(LEN=120) :: FNAME
   CHARACTER(LEN=160) :: PATHNFILE
   
   IF(DBG_SET(DBG_SBR)) &
         & WRITE(IPT,*) "Subroutine Begins: Read_Heating_Calculate_Namelist;"

   IOS = 0

   FNAME = "./"//trim(casename)//"_run.nml"

   CALL FOPEN(NMLUNIT,trim(FNAME),'cfr')

   !READ NAME LIST FILE
    REWIND(NMLUNIT)

   ! Read IO Information
   READ(UNIT=NMLUNIT, NML=NML_HEATING_CALCULATED,IOSTAT=ios)
   if(ios .NE. 0 ) Then
     if(DBG_SET(dbg_log)) write(UNIT=IPT,NML=NML_HEATING_CALCULATED)
     Call Fatal_error("Can Not Read NameList NML_HEATING_CALCULATED from file: "//trim(FNAME))
   end if
   CLOSE(NMLUNIT)
   
   IF(HEATING_CALCULATE_ON .AND. .NOT. HEATING_ON)THEN
     HEATING_ON = .TRUE.
   ELSE IF(HEATING_CALCULATE_ON .AND. HEATING_ON)THEN
     CALL FATAL_ERROR("IN NAMELIST CONTROL FILE, IF HEATING_CALCULATE_ON = TRUE, ", &
                      "THEN SET HEATING_ON = FALSE") 
   END IF 
   IF(HEATING_CALCULATE_TYPE == 'flux') HEATING_TYPE = 'flux'
   IF(HEATING_CALCULATE_TYPE == 'body') HEATING_TYPE = 'body'

   RETURN
   END SUBROUTINE HEATING_CALCULATE_NAMELIST_READ

!==============================================================================|

!==============================================================================|
   SUBROUTINE COARE26Z(UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT)
!  A=coare26(u,zu,Ta,zt,rh,zq,Pa,Ts,dlw,dsw)
! Simplified non-vectorized version of coare2.6 code
! with cool skin option retained but warm layer and 
! surface wave options removed, and rain set to zero. 
! Assumes input are single scalars and that u is the
! magitude of the difference between wind and surface 
! current vectors, if latter available. Output:
! A = [tau qsen qlat Cd Ch Ce Cdn_10 Chn_10 Cen_10].
 
   REAL(SP) :: UR,ZU,TA,ZT,RH,ZQ,PA,TS,DLW,DSW,TAU,HSB,HLB,LAT,LE,L,JCOOL,L10
   REAL(SP) :: U,US,T,P,RL,RS,QSAT26S,QS,Q,ZI,BETA,VON,FDG,TDK,GRVS,GRAV,RGAS
   REAL(SP) :: CPA,CPV,RHOA,VISA,AL,BE,CPW,RHOW,VISW,TCW,BIGC,WETC,RNS,RNL,DU
   REAL(SP) :: DT,DQ,UG,DTER,DQER,UT,U10,USR,ZO10,CD10,CH10,CT10,ZOT10,CD,CT,CC
   REAL(SP) :: RIBCU,RIBU,ZETU,PSIU_26S,PSIT_26S,TSR,QSR,TKT,CHARN,ZET
   REAL(SP) :: ZO,ZOQ,ZOT,BF,QOUT,DELS,QCOL,ALQ,XLAMX,CH,CE,CDN_10,CHN_10,CEN_10,RR
   
   INTEGER  :: I
   INTEGER  :: NITS  ! new, Karsten Lettmann
   
   ! ---- new : Karsten Lettmann , Jan. 2018 --------
   REAL(SP) :: EPSI = 1.0e-16_SP   ! small value to stabilize the code
   ! -------- end new -------------------------------

   
! Set jcool=0 if Ts is surface, =1 if Ts is bulk.
! rcb checked 6/9/04
! set jcool=1 if Ts is bulk, 0 if Ts is true skin jcool=1;

! Changes to code:
! - changed line for zetu to original coare3 code
!   ZETU=CC*RIBU*(1._SP + 27._SP/9._SP*RIBU/CC)
! - Changed number of iteration in Bulk loop to NITS = 3
! - Set wind to have at least 0.1 m/s
! - introduce a small epsi to stabilize similar to ROMS implementation in file bulk_fluxes.F




! --------- New. karsten Lettmann -----
#  if !defined (NO_COOL_SKIN_HEAT_FLUX)
   JCOOL=1. ! original line
#  else
   !  Switch off the cool skin effect, as it might cause problems
   ! Cool skin is not used
   JCOOL=0.0_SP
#  endif   
! ---------end new ---------------------------- 

! rename variables from fairall et al coare3 code
 

!  wind speed (m/s) at height zr (m)
   U=UR
   
!  --------- New. Karsten Lettmann ---------
!  Allways have some small wind to make the allgorithm more stable
!  At the moment U > 0.1 m/s is set
#  if !defined(DOUBLE_PRECISION)
   U=AMAX1(.1_SP,U)
#  else
   U=DMAX1(.1_SP,U)
#  endif
!  -----------------------------------------
   
!  surface current speed in the wind direction(m/s)
   US=0*UR
!  water temperature (deg C)
   TS=TS
!  BULK AIR TEMPERATURE (C) AT HEIGHT ZT(m)
   T=TA
!  RELATIVE HUMIDITY (%) AT HEIGHT zq(M)
   RH=RH
! ----- new: Karsten Lettmann, Feb. 2016 --------------------------
!  SURFACE PRESSURE (mb)
   P=PA    ! switched on again, in Feb. 2016
!  The following was the orginal line to convert pressure from Pa -> mbar
!  This was commented out, as it is assumed that surface pressure for
!  heat flux is provided in mbar as desired by the COARE allgortihm
!   P=PA/100.  !ejw 8/15/2011 from pa to mb in this subroutine only, original line

   ! ------ new K. Lettmann , 2016.02.16 -------
   !print*,'------- Pressure in COAR26 P: ',  P
   !print*,'------- Pressure in COAR26 Pa: ', PA
! -------  end new             ------------------------------------

!  DOWNWARD LONGWAVE RADIATION (W/m2)
   RL=DLW
!  DOWNWARD SHORTWAVE RADIATION (W/m2)
   RS=DSW
!  CONVERT RH TO SPECIFIC HUMIDITY at sea surface(G/KG)
   CALL QSAT26(TS,P,QSAT26S)
   QS=0.98_SP*QSAT26S/1000._SP   ! and convert to KG/KG
   
   ! ------ new J Kirchner , Oct 2017 --------------------
   ! Switch off the calculation of specific humidity and use values
   ! for specific humidity that are stored in heatflux input file
   !    
#  if !defined (USE_SPEC_HUM_FOR_HEAT_FLUX)   
   !  SPECIFIC HUMIDITY OF AIR (G/KG)  
   CALL QSAT26(T,P,QSAT26S)                      ! original line
   Q=(0.01_SP*RH)*QSAT26S/1000._SP               ! and convert to KG/KG / original line
#  else
   ! here it is assumed RH already contains specific humidity 
   ! in units KG/KG
   Q = RH   
#  endif
! ----------end new -------------------------------------


! -------- new : some output ---------------
!   print*,'===================== some output ================'
!   print*,'------- U10     : ', U
!   print*,'------- P       : ', P
!   print*,'------- QS      : ', QS
!   print*,'------- Q       : ', Q
!   print*,'------- Tair    : ', T
!   print*,'------- Twater  : ', TS
!   print*,'------- SW      : ', RS
!   print*,'------- LW      : ', DLW
   ! ------------------------------------------

!            print *,'qq',Q
!   SET RAIN TO ZERO RAIN=0*U
!   SET RAIN RATE (MM/HR) - KEEP AS OPTION



!  ***********SET LOCAL CONSTANTS *********
   !ZI:  PBL HEIGHT (M)
   !Zi:  Height (m) of atmospheric boundary layer (in ROMS: blk_Zabl = 600.m)
   ZI = 600.0_SP

   ! Number of iterations of Coare Allgorithm default
   ! Orignal 6 iterations
   ! ROMS uses 3 iterations
   NITS = 3
   
   
!  LATITUDE (DEG,N=+)- GEORGES BANK
!  LAT=42.
!  *************************************


!  ************SET CONSTANTS **************
   BETA=1.2_SP
   VON=0.4_SP
   FDG=1.00_SP
   TDK=273.16_SP
   CALL GRV(LAT,GRVS)
   GRAV=GRVS
! ************AIR CONSTANTS **************
   RGAS=287.1_SP
   LE=(2.501_SP-0.00237_SP*TS)*1000000._SP
   CPA=1004.67_SP
   CPV=CPA*(1._SP+0.84_SP*Q)
   RHOA=P*100._SP/(RGAS*(T+TDK)*(1+0.61_SP*Q))
   VISA=1.326_SP*0.00001_SP*(1+6.542_SP*0.001_SP*T+  &
        8.301_SP*0.000001_SP*T*T-4.84_SP*0.000000001_SP*T*T*T)
! ***********COOL SKIN CONSTANTS******************
   AL=2.1_SP*0.00001_SP*((TS+3.2_SP)**0.79_SP)
   BE=0.026_SP
   CPW=4000._SP
   RHOW=1022._SP
   VISW=0.000001_SP
   TCW=0.6_SP
   BIGC=16._SP*GRAV*CPW*((RHOW*VISW)**3)/(TCW*TCW*RHOA*RHOA)
   WETC=0.622_SP*LE*QS/(RGAS*((TS+TDK)**2))
! ***************COMPUTE AUX STUFF *********
   RNS=RS*0.945_SP
   RNL=0.97_SP*(5.67_SP*0.00000001_SP*((TS-0.3_SP*JCOOL+TDK)**4)-RL)


! **************BEGIN BULK LOOP ***********

! **************FIRST GUESS ***************
   DU=U-US
   DT=TS-T-0.0098_SP*ZT
   DQ=QS-Q
   TA=T+TDK
   UG=0.5_SP
   DTER=0.3_SP
   DQER=WETC*DTER
   UT=SQRT(DU*DU+UG*UG) 
   U10=UT*LOG(10._SP/1e-4)/LOG(ZU/1e-4)
   USR=0.035_SP*U10

   ZO10=0.011_SP*USR*USR/GRAV+0.11_SP*VISA/USR
   CD10=(VON/LOG(10._SP/ZO10))**2
   CH10=0.00115_SP
   CT10=Ch10/SQRT(CD10)
   ZOT10=10._SP/EXP(VON/CT10)
   CD=(VON/LOG(ZU/ZO10))**2
   CT=VON/LOG(ZT/ZOT10)
   CC=VON*CT/CD
   RIBCU=-ZU/(ZI*0.004_SP*(BETA**3))
   RIBU=-GRAV*ZU*((DT-DTER*JCOOL)+.61_SP*TA*DQ)/(TA*(UT**2))
!       same as edson
!   NITS=6._SP    ! original line 
   IF(RIBU < 0)THEN
     ZETU=CC*RIBU/(1._SP+RIBU/RIBCU)
   ELSE
     ! ------ new : Karsten Lettmann, 2017 Dec --------
     ! The following line was different from COARE 3.0 code
     ! This inconsistency was already found in mod_heatflux_gl.F
     ! MDR, corrected from F90 version
     !
     ! ZETU=CC*RIBU/(1._SP+27._SP/(9*RIBU*CC))      ! original FVCOM line     
     ! ZETU=CC*RIBU*(1._SP + 27._SP/9._SP*RIBU/CC)  ! original coare3
     ZETU  =CC*RIBU*(1._SP + 3._SP       *RIBU/CC)  ! slightly changed
     ! ---- end new -----------------------------------
   END IF

   L10=ZU/ZETU
   IF(ZETU > 50)THEN
     NITS=1
   END IF
   CALL PSIU_26(ZU/L10,PSIU_26S)
   USR=UT*VON/(LOG(ZU/ZO10)-PSIU_26S)
   CALL PSIT_26(ZT/L10,PSIT_26S)
   TSR=-(DT-DTER*JCOOL)*VON*FDG/(LOG(ZT/ZOT10)-PSIT_26S)
   CALL PSIT_26(ZQ/L10,PSIT_26S)
   
   QSR=-(DQ-WETC*DTER*JCOOL)*VON*FDG/(LOG(ZQ/ZOT10)-PSIT_26S)
   TKT=.001_SP
   CHARN=0.011_SP
   
   IF(UT > 10)THEN 
     CHARN=0.011_SP+(UT-10)/(18-10)*(0.018_SP-0.011_SP)
   END IF
   
   IF(UT > 18)THEN
     CHARN=0.018_SP
   END IF

!*************** bulk loop ************
    do I=1,nits
!   DO I=1,1
     ZET=VON*GRAV*ZU*(TSR*(1._SP+0.61_SP*Q)+0.61_SP*TA*QSR)/   &
         (TA*(USR*USR)*(1._SP+0.61_SP*Q) + EPSI)
     ZO=CHARN*USR*USR/GRAV+0.11_SP*VISA/(USR+EPSI)
     RR=ZO*USR/VISA
     L=ZU/(ZET+EPSI)
#    if !defined(DOUBLE_PRECISION)
     ZOQ=AMIN1(1.15e-4_SP,5.5e-5_SP/(RR**0.6_SP)) 
#    else
     ZOQ=DMIN1(1.15e-4_SP,5.5e-5_SP/(RR**0.6_SP)) 
#    endif
     ZOT=ZOQ
     CALL PSIU_26(ZU/L,PSIU_26S)
     USR=UT*VON/(LOG(ZU/ZO)-PSIU_26S)
     ! ----- new, KL -------
#    if !defined(DOUBLE_PRECISION)     
     USR=AMAX1(EPSI,USR) 
#    else     
     USR=DMAX1(EPSI,USR) 
#    endif
     ! ---- end new --------
     CALL PSIT_26(ZT/L,PSIT_26S)
     TSR=-(DT-DTER*JCOOL)*VON*FDG/(LOG(ZT/ZOT)-PSIT_26S)
     CALL PSIT_26(ZQ/L,PSIT_26S)
     QSR=-(DQ-WETC*DTER*JCOOL)*VON*FDG/(LOG(ZQ/ZOQ)-PSIT_26S) 
     BF=-GRAV/TA*USR*(TSR+.61_SP*TA*QSR)               
!      print *,usr,tsr,qsr,Bf
     IF(BF > 0)THEN
       UG=BETA*((BF*ZI)**0.333_SP)
     ELSE
       UG=0.2_SP
     END IF
     
     UT=SQRT(DU*DU+UG*UG) 
     RNL=0.97_SP*(5.67_SP*0.00000001_SP*((TS-DTER*JCOOL+TDK)**4)-RL)
     HSB=-RHOA*CPA*USR*TSR
     HLB=-RHOA*LE*USR*QSR
     QOUT=RNL+HSB+HLB
 
!JQI DELS=RNS*(.065_SP+11*TKT-6.6_SP*0.00001_SP/(TKT*(1-EXP(-TKT/8.0_SP*0.0001_SP))))
     DELS=RNS*(.065_SP+11.0_SP*TKT-6.6_SP*0.00001_SP/TKT*(1.0_SP-EXP(-TKT/8.0e-4_SP)))

     QCOL=QOUT-DELS
     ALQ=AL*QCOL+BE*HLB*CPW/LE
     IF(ALQ > 0)THEN
       XLAMX=6._SP/(1._SP+(BIGC*ALQ/(USR+EPSI)**4)**0.75_SP)**0.333_SP
       TKT=XLAMX*VISW/(SQRT(RHOA/RHOW)*(USR+EPSI))
     ELSE
       XLAMX=6.0_SP
#    if !defined(DOUBLE_PRECISION)
       TKT=AMIN1(.01_SP,XLAMX*VISW/(SQRT(RHOA/RHOW)*(USR+EPSI)))
#    else
       TKT=DMIN1(.01_SP,XLAMX*VISW/(SQRT(RHOA/RHOW)*(USR+EPSI)))
#    endif
     END IF
!           cool skin
     DTER=QCOL*TKT/TCW
     DQER=WETC*DTER
   END DO
!  of  end bulk iter loop

!****** compute fluxes ******************************************
!             wind stress
   TAU=RHOA*USR*USR*DU/(UT+EPSI)
!              sensible heat flux
   HSB=RHOA*CPA*USR*TSR
!              latent heat flux
   HLB=RHOA*LE*USR*QSR

! -------- new : some output ---------------
!   print*,' ' 
!   print*,'the fluxes:'
!   print*,'------- Tau     : ', TAU
!   print*,'------- HSB     : ', HSB
!   print*,'------- HLB     : ', HLB
! ------------------------------------------
   
   

!****** compute transfer coeffs relative to ut @ meas. ht ********
#  if !defined(DOUBLE_PRECISION)
   CD=TAU/RHOA/(UT+EPSI)/AMAX1(.1_SP,DU)
#  else
   CD=TAU/RHOA/(UT+EPSI)/DMAX1(.1_SP,DU)
#  endif
   CH=-USR*TSR/(UT+EPSI)/(DT-DTER*JCOOL)
   CE=-USR*QSR/(DQ-DQER*JCOOL)/(UT+EPSI)

!****** compute 10-m neutral coeff relative to ut ****************
   CDN_10=VON*VON/LOG(10._SP/ZO)/LOG(10._SP/ZO)
   CHN_10=VON*VON*FDG/LOG(10._SP/ZO)/LOG(10._SP/ZOT)
   CEN_10=VON*VON*FDG/LOG(10._SP/ZO)/LOG(10._SP/ZOQ)

!******** rain heat flux (save to use if desired) *************
! dwat=2.11e-5*((t+tdk)/tdk)^1.94; %! water vapour diffusivity
! dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); %!heat diffusivity
! alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp)); %! wet bulb factor
! RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600;
!**************************************************************
!----------------------------------------------------------
   RETURN
   END SUBROUTINE COARE26Z
   

!-----------------------------------------------------------------------|
   SUBROUTINE PSIT_26(ZET,PSI)
! computes temperature structure function
   IMPLICIT NONE
   REAL(SP) :: ZET,PSI,X,PSIK,PSIC,F,C
   
   IF(ZET < 0.0_SP)THEN
     X=(1._SP-15._SP*ZET)**.5_SP
     PSIK=2._SP*LOG((1._SP+X)/2._SP)
     X=(1._SP-34.15_SP*ZET)**.3333_SP
     PSIC=1.5_SP*LOG((1._SP+X+X*X)/3._SP)-SQRT(3._SP)    &
          *ATAN((1._SP+2._SP*X)/SQRT(3._SP))+4._SP*ATAN(1._SP)/SQRT(3._SP)
     F=ZET*ZET/(1._SP+ZET*ZET)
     PSI=(1._SP-F)*PSIK+F*PSIC
   ELSE
#    if !defined(DOUBLE_PRECISION)
     C=AMIN1(50._SP,.35_SP*ZET)
#    else
     C=DMIN1(50._SP,.35_SP*ZET)
#    endif
     PSI=-((1._SP+2._SP/3._SP*ZET)**1.5_SP+.6667_SP     &
         *(ZET-14.28_SP)/EXP(C)+8.525_SP)
   END IF
   
   RETURN
   END SUBROUTINE PSIT_26
  

!----------------------------------------------------------
   SUBROUTINE PSIU_26(ZET,PSI)
! computes velocity structure function
   IMPLICIT NONE
   REAL(SP) :: ZET,PSI,X,PSIK,PSIC,F,C

   IF(ZET < 0.0_SP)THEN
     X=(1._SP-15._SP*ZET)**.25_SP
     PSIK=2._SP*LOG((1._SP+X)/2._SP)+LOG((1._SP+X*X)/2._SP)   &
          -2._SP*ATAN(X)+2._SP*ATAN(1._SP)
     X=(1._SP-10.15_SP*ZET)**.3333_SP
     PSIC=1.5_SP*LOG((1._SP+X+X*X)/3._SP)-SQRT(3._SP)   &
          *ATAN((1._SP+2._SP*X)/SQRT(3._SP))+4._SP*ATAN(1._SP)/SQRT(3._SP)
     F=ZET*ZET/(1._SP+ZET*ZET)
     PSI=(1._SP-F)*PSIK+F*PSIC
   ELSE
#    if !defined(DOUBLE_PRECISION)
     C=AMIN1(50._SP,.35_SP*ZET)
#    else
     C=DMIN1(50._SP,.35_SP*ZET)
#    endif
     PSI=-((1._SP+1.0_SP*ZET)**1.0_SP+.667_SP*(ZET-14.28_SP)/EXP(C)+8.525_SP)
   END IF
   
   RETURN
   END SUBROUTINE PSIU_26


!-----------------------------------------------------------
   SUBROUTINE QSAT26(T,P,QS)
   ! computes saturation specific humidity in g/kg
   !
   ! For Literature see e.g.:
   !     - https://en.wikipedia.org/wiki/Clausius%E2%80%93Clapeyron_relation
   !     - Jacobsen, Fundamentals of atmospheric modelling, page 28-33, 45
   !       QS follows from eq (2.32) and p_d = P - ES
   
   ! Input variables are:
   ! T  : Temperature in °C
   ! P  : total air pressure in hPa = mbar
   
   IMPLICIT NONE
   REAL(SP) :: T,P,QS,ES
   
   ! first, calulate the saturation water pressure ES in hPa
   ! Using a modified form of the August-Roche-Magnus formula
   ES=6.112_SP*EXP(17.502_SP*T/(T+241.0_SP))     &
      *(1.0007_SP+3.46_SP*0.000001_SP*P)
      
   !second, specific humidity in (g water vapour) / (kg moist air)
   !
   ! 622 = 1000 * eps with eps = 0.622 = R_d / R_v (ratio of gas constants)
   ! 0.378 = 1 - eps
   !
   QS=ES*622/(P-.378_SP*ES)
   
   RETURN
   END SUBROUTINE QSAT26


!------function g=grv(lat)----------------------------------------------|
   SUBROUTINE GRV(LAT,G)
   IMPLICIT NONE
   REAL(SP) :: PI,GAMMA,C1,C2,C3,C4,PHI,X,G,LAT

   PI=3.1415926_SP
! computes g given lat in deg
   GAMMA=9.7803267715_SP
   C1=0.0052790414_SP
   C2=0.0000232718_SP
   C3=0.0000001262_SP
   C4=0.0000000007_SP
   PHI=LAT*PI/180._SP
   X=SIN(PHI)
   G=GAMMA*(1._SP+C1*X**2+C2*X**4+C3*X**6+C4*X**8)
   
   RETURN
   END SUBROUTINE GRV

   
   
! ============ New Bulk Heat Flux by ICBM ==============================
! Subroutine for simple bulk Heat Flux implementation
!
! It uses the bulk fluxes after ?????
!
! Implemented by Julia Kirchner, ICBM, Oldenburg, Germany, 2017 Oct.
!

!   SUBROUTINE HEAT_BULK_ICBM(UR,ZU,TA,ZT,RH,ZQ,PA,TS,TAU,HSB,HLB)

   
   ! ------ Input variables -----------
!   REAL(SP), INTENT(IN) :: UR      ! Wind speed (m/s) at height ZU (m)
!   REAL(SP), INTENT(IN) :: ZU      ! height for wind speed (m), usually 10m
!   REAL(SP), INTENT(IN) :: TA      ! air temperature (°C) at height ZT (m)
!   REAL(SP), INTENT(IN) :: ZT      ! height for air temperature (m), usually 2m
!   REAL(SP), INTENT(IN) :: RH      ! RELATIVE HUMIDITY (%) AT HEIGHT zq(M)
!   REAL(SP), INTENT(IN) :: ZQ      ! height for relative humidity (m), usually 2m
!   REAL(SP), INTENT(IN) :: PA      ! SURFACE PRESSURE (mbar)
!   REAL(SP), INTENT(IN) :: TS      ! SURFACE water temperature (°C) (in top layer)
   
   
   ! ----- Output variables -----------
!   REAL(SP), INTENT(OUT) :: TAU      ! amplitude of Wind stress (N/m^2) at height ZU (m)
!   REAL(SP), INTENT(OUT) :: HSB      ! sensible heat flux (W/m^2) (positive if ocean is heated)
!   REAL(SP), INTENT(OUT) :: HLB      ! latent   heat flux (W/m^2) (positive if ocean is heated)
   
   
   
   ! ---- local variables --------------
   
   
!   INTEGER  :: I

 


!  CONVERT RH TO SPECIFIC HUMIDITY (G/KG)
!   CALL QSAT26(TS,P,QSAT26S)
!   QS=0.98_SP*QSAT26S/1000._SP
!  SPECIFIC HUMIDITY OF AIR (G/KG)  
!   CALL QSAT26(T,P,QSAT26S)
!   Q=(0.01_SP*RH)*QSAT26S/1000._SP

!   RETURN
!   END SUBROUTINE HEAT_BULK_ICBM
! ======== end new ======================================================
   
   
END MODULE MOD_HEATFLUX
# endif


