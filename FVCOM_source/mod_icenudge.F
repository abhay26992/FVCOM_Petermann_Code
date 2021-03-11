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
! THIS MODULE is developed by QIN ZHOU (APN) & TORE HATTERMANN (NPI)
!/---------------------------------------------------------------------------/
! CVS VERSION INFORMATION
! $Id$
! $Name$
! $Revision$
!/===========================================================================/

MODULE MOD_ICENUDGE
# if defined (ICENUDGE)   
   USE ALL_VARS
   USE MOD_PREC
   USE MOD_PAR
   IMPLICIT NONE
   SAVE

  CONTAINS

  SUBROUTINE ICENUDGE_OCN
!------------------------------------------------------------------------------|
!          DETERMINE NET HEAT, FRESHWATER AND SALT FLUX TO OCEAN
!------------------------------------------------------------------------------|
  
! read input data from ocean

   USE MOD_PREC
   USE ALL_VARS
   USE MOD_PAR


 IMPLICIT NONE
    
    INTEGER :: I,K,J
    REAL(SP):: vrel,fh,dhdt,Tf,SPRO ! relative ice-ocean velocity & tmp variable for albedo
    REAL (SP), DIMENSION(:),ALLOCATABLE, TARGET :: ustar,ustar_tmp,fbot,albv,transbot


    !-------SET CONSTANT-------------------------------------
    REAL(DP),PARAMETER :: cp_ocn  = 4218._DP  !! Specific heat capacity of ocean(J/kg/C)
    REAL(DP),PARAMETER :: c_h   = 0.006_DP  !! heat transfer coefficient 
    REAL(DP),PARAMETER :: Lfresh = 3.34E5_DP !! Latent heat of melting of fresh ice (J/kg)
    REAL(DP),PARAMETER :: ustar_min =0.005_DP  !! Minimum friction velocity for ocean heat flux (m/s)
    REAL(DP),PARAMETER :: rhoi = 917.0_DP !!  Density of ice (kg/m^3)
    REAL(DP),PARAMETER :: rhow = 1026.0_DP !!  Density of seawater (kg/m^3)
    REAL(DP),PARAMETER :: rhofresh = 1000.0_DP !! Density of fresh water (kg/m^3)
    REAL(DP),PARAMETER :: depressT = 0.054_DP !! Tf:brine salinity ratio (C/ppt)
    REAL(DP),PARAMETER :: SPCP = 4.2173E3_DP !! Specific heat capacity, keep the same value as in bcond 
    ! drag force
    REAL(DP),PARAMETER :: Cw = 0.0055_DP !! ocean-ice neutral drag coefficient
    REAL(SP), PARAMETER :: COSW = 1.0_SP  !cos(ocean turning angle) !turning angle = 0.
    REAL(SP), PARAMETER :: SINW = 0.0_SP  !sin(ocean turning angle) !turning angle = 0.
    ! shortwave
    REAL(SP), PARAMETER :: i0vis     = 0.70_DP   ! fraction of penetrating solar rad (visible)
    REAL(SP), PARAMETER :: albicev   = 0.78_DP   ! visible ice albedo for h > ahmax
    REAL(SP), PARAMETER :: ahmax     = 0.3_DP    ! thickness above which ice albedo constant (m)
    REAL(SP), PARAMETER :: c4i       = 4.0_DP
    REAL(SP), PARAMETER :: c1i       = 1.0_DP
    REAL(SP), PARAMETER :: p001      = 0.001_DP
    REAL(SP), PARAMETER :: albocn    = 0.06_DP   !ocean albedo
    REAL(SP), PARAMETER :: kappav    = 1.4_DP   ! vis extnctn coef in ice, wvlngth<700nm (1/m)


! in the ice-shelf interface

  ALLOCATE(ustar_tmp(0:NT));ustar_tmp = ZERO  ! friction velocity on cells(m/s)
  ALLOCATE(ustar(0:MT));   ustar = ZERO  ! friction velocity (m/s)
  ALLOCATE(fbot(0:MT));    fbot = ZERO  ! heat flux for basal melting 
  ALLOCATE(transbot(0:MT));    transbot = ZERO  !  sw transmitted to the ice bottom
  ALLOCATE(albv(0:MT));   albv = ZERO  ! albedo for visiable components absorbed in the ice interior or transmitted to the ocean. 
  

! compute ocean-ice stress

  CALL N2E2D(AICE,AICEU)

#   if defined (MULTIPROCESSOR)
        IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,AICEU)
#   endif


  DO I = 1,NT
     IF(AICEU(I)>=1.0E-5_SP)  THEN	 
       vrel = rhow*Cw*sqrt((U(I,1) - UICE(I))**2 + (V(I,1) - VICE(I))**2)  ! m/s
       ! turning angle is oppsite in the northen and southern hemisphere
       STROCNX(I) =  vrel*((U(I,1) - UICE(I))*COSW - (V(I,1) - VICE(I))*SINW*SIGN(1.0_SP,COR(I)))   
       STROCNY(I) =  vrel*((V(I,1) - VICE(I))*COSW + (U(I,1) - UICE(I))*SINW*SIGN(1.0_SP,COR(I)))
       ustar_tmp(I) = sqrt(sqrt(STROCNX(I)**2+STROCNY(I)**2)/rhow)
     ENDIF
  END DO

#   if defined (MULTIPROCESSOR)
        IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,STROCNX,STROCNY)
        IF(PAR) CALL AEXCHANGE(EC,MYID,NPROCS,ustar_tmp)
#   endif

    CALL E2N2D(ustar_tmp,ustar)

#   if defined (MULTIPROCESSOR)
        IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,ustar)
        IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,ustar)
#   endif    

! compute flux to ocean
    SPRO = SPCP*rhow  
    DO I = 1,MT
        IF(AICE(I)>=1.0E-5_SP)  THEN

          ! SW penetrating thru ice into ocean, only visible components (28% direct and 24% diffuse) absorbed 
          ! in the ice interior or transmitted to the ocean & assume frsnow = 0
          fh = min(atan(HICE(i)*c4i)/atan(ahmax*c4i),c1i)
          albv(I)  = albicev*fh  
          transbot(I) = exp(-kappav*HICE(I))
          FSWTHRU(I) = SWRAD(I)*0.52*(1-albv(I))*i0vis*transbot(I)
        
          ! basal melt/freeze heat flux
          Tf  = -depressT*S1(I,1)
          FHEAT(I) = -rhow*cp_ocn*c_h*MAX(ustar(I),ustar_min)*(T1(I,1)-Tf)  !  <0 melting, ocean loses heat
          dhdt    = FHEAT(I)/(rhoi*Lfresh)              ! dhdt < 0, melting
 
          ! salt flux and fresh water flux, ocean gets water & salt when dhdt<0 melting
          FSALT(I) = (-1)*rhoi*dhdt*ISALT(I)*p001
          FRESH(I) = (-1)*rhoi*dhdt

          !===============================================
          ! add the net heat flux and sw flux,! WTSURF & SWRAD is multiplied with (-1) 
          ! before it is called by vdif_ts,  that is why we times (-1)* here
          WTSURF(I) = (1-AICE(I))*WTSURF(I)+(-1)*AICE(I)*(FHEAT(I)+FSWTHRU(I))/SPRO 
          SWRAD(I)  = (1-AICE(I))*SWRAD(I) +(-1)*AICE(I)*FSWTHRU(I)/SPRO

        ENDIF

   END DO

# if defined (MULTIPROCESSOR)
       IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,WTSURF,SWRAD)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,WTSURF,SWRAD)
       IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,FSALT,FRESH)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,FSALT,FRESH)
       IF(PAR) CALL NODE_MATCH(1,NBN,BN_MLT,BN_LOC,BNC,MT,1,MYID,NPROCS,FHEAT,FSWTHRU)
       IF(PAR) CALL AEXCHANGE(NC,MYID,NPROCS,FHEAT,FSWTHRU)

# endif

 
!!======================================================
!!======================================================
 END SUBROUTINE ICENUDGE_OCN

!==============================================================================|
!==============================================================================|
# endif
END MODULE  MOD_ICENUDGE
