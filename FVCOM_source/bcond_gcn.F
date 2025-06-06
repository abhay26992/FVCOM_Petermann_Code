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
!     SET BOUNDARY CONDITIONS FOR ALMOST ALL VARIABLES                         |
!                                                                              |
!         idx: identifies which variables are considered                       |
!              1=tidal forcing                                                 |
!              2=solid bcs for external mode uaf and vaf                       |
!              3=solid bcs for internal mode uf and vf                         |
!              4=open bcs for s and t                                          |
!              5=solid bcs for internal mode u and v                           |
!              6=unused                                                        |
!              7=unused                                                        |
!              8=the surface forcings for internal mode                        |
!              9=the surface forcings for external mode                        |
!                                                                              |
!==============================================================================|

   SUBROUTINE BCOND_GCN(IDX,K_RK)

!==============================================================================|
   USE ALL_VARS
   USE BCS
   USE MOD_OBCS
   USE MOD_FORCE
   USE MOD_PAR
   USE MOD_WD
   USE MOD_BULK

# if defined(FABM)
   USE MOD_FABM_DATA
!   USE fabm
# endif
#  if defined (WATER_QUALITY)
   USE MOD_WQM
#  endif
#  if defined (SEDIMENT)
   USE MOD_SED
#  endif

#  if defined (THIN_DAM)
   USE MOD_DAM, ONLY : NBE_DAM,KDAM1,CLP_CELL,CLP_ALPHA,CLP_KDAM   !Jadon
#  endif

#  if defined (HEATING_CALCULATED) || (HEATING_CALCULATED_GL)
   USE MOD_HEATFLUX, ONLY : HEATING_CALCULATE_ON
#  endif
#  if defined (HEATING_SOLAR)
   USE MOD_SOLAR, ONLY : HEATING_SOLAR_ON
#  endif

   IMPLICIT NONE
   INTEGER, INTENT(IN) :: IDX
   REAL(SP) :: ZREF(KBM1),ZREFJ(KBM1),TSIGMA(KBM1),SSIGMA(KBM1)
   REAL(SP) :: TTMP(KBM1),STMP(KBM1),TREF(KBM1),SREF(KBM1)
   REAL(SP) :: PHY_Z(KBM1),PHY_Z1(KBM1)
   REAL(SP) :: TT1(KBM1),TT2(KBM1),SS1(KBM1),SS2(KBM1)
!   REAL(SP) :: TIME1,FACT,UFACT,FORCE,QPREC,QEVAP,UI,VI,UNTMP,VNTMP,TX,TY,HFLUX
   REAL(SP) :: TIME1,FACT,UFACT,FORCE,UI,VI,UNTMP,VNTMP,TX,TY,HFLUX
!   REAL(SP) :: DTXTMP,DTYTMP,QPREC2,QEVAP2,SPRO,WDS,CD,SPCP,ROSEA
   REAL(SP) :: DTXTMP,DTYTMP,SPRO,WDS,CD,SPCP,ROSEA,ROSEA1(MT),SPRO1(MT)
   REAL(SP) :: PHAI_IJ,ALPHA1,DHFLUXTMP,DHSHORTTMP,HSHORT,TIMERK1

   ! VARIABLES FOR LONG SHORE FLOW STUFF
   REAL(SP) :: ANG_WND,WNDALONG,RHOINTNXT,RHOINTCUR,CUMEL
   REAL(SP) :: TAU_X,TAU_Y, mag_wnd
   REAL(SP), POINTER,DIMENSION(:) :: lcl_dat, gbl_dat


   REAL(SP),POINTER :: eta_lcl(:), eta_gbl(:) ,elfgeo_gbl(:),elfgeo_lcl(:)

   INTEGER  I,J,K,I1,I2,J1,J2,II,L1,L2,IERR
   INTEGER  cdx,ndx,KDAM_TMP,K_RK
# if defined(FABM)
   INTEGER JJ,IFABM
# endif

!!$!=================|DEBUG|============================
!!$   LOGICAL, save :: INIT=.false.
!!$   INTEGER, save :: icount
!!$   character(len=5):: ccount
!!$   REAL(SP),POINTER :: temp_lcl(:),temp_gl(:) 
!!$!=================|DEBUG|============================


  if(dbg_set(dbg_sbr)) write(ipt,*) "Start: bcond_gcn: ",idx


   SELECT CASE(IDX)
!==============================================================================|
   CASE(1) !Surface Elevation Boundary Conditions (Tidal Forcing)              !
!==============================================================================|
   IF(OBC_ELEVATION_FORCING_ON) CALL BCOND_ASL
   CALL BCOND_ASL_CLP
   CALL BCOND_GWI(K_RK)
   CALL BCOND_BKI(K_RK)
   CALL BCOND_ORE

!
!--Allow setup/down on north boundary in response to longshore wind
!--Corrects for Wind-Driven Barotropic Response (See Schwing 1989)
!--Implemented by Jamie Pringle - Rewritten for parallel by D.Stuebe
!

   IF (OBC_LONGSHORE_FLOW_ON) THEN
      
      ! CALCULATE THE ADJUSTMENT DUE TO WIND DRIVEN FLOW ALLONG THE
      ! COAST OF NOVA SCOTIA
      DO I = 1,nobclsf
         TAU_X = 0.0_sp
         TAU_y = 0.0_sp

         I1= ibclsf(I)
         DO J = 1,NTVE(I1)
            K = NBVE(I1,J)
            TAU_X = TAU_X + WUSURF2(K)
            TAU_Y = TAU_Y + WVSURF2(K)
         END DO
         
         TAU_X = TAU_X/real(NTVE(I1),SP)
         TAU_Y = TAU_Y/real(NTVE(I1),SP)

         TAU_X = TAU_X * 1000.0_SP ! Approximate scale factor
         TAU_Y = TAU_Y * 1000.0_SP ! Approximate scale factor

         MAG_WND = SQRT(TAU_X**2+TAU_Y**2)
         
         
         IF (MAG_WND .GT. 0.0_sp) THEN
            ANG_WND=ATAN2(TAU_Y,TAU_X) - WDF_ANG(I)
         ELSE
            ANG_WND=0.0_sp
         END IF
         

         WNDALONG=SIN(ANG_WND)*MAG_WND

         ! SUBTRACT THE COMPONET ALONG SHORE AND INTO THE DOMAIN FROM
         ! THE SEA SURFACE HEIGHT
         
         ELF(I1)=ELF(I1) - WNDALONG*RBC_WDF(I)

!!$         !=================|DEBUG|============================
!!$ THIS IS NOT PARALLEL SAFE!
!!$            if(.not. INIT)then
!!$               init =.true.
!!$               icount=0
!!$            end if
!!$            
!!$            IF(I==1) THEN
!!$               write(ccount,'(I5.5)') icount
!!$               
!!$               call FOPEN(157,"WDF_"//ccount,'ofr')
!!$               write(157,*) "WDF DEBUG: ",icount
!!$               write(157,*) "Node#, ANG,  ELFGEO"
!!$            END IF
!!$
!!$            write(157,*) IBCLSF(I), WDF_ANG(I)/deg2rad, RBC_WDF(I),&
!!$                 &-1.0_SP*WNDALONG*RBC_WDF(I), MAG_WND*SIN(ANG_WND)
!!$
!!$            IF(I==nobclsf) THEN
!!$               icount = icount +1
!!$               
!!$               close(157)
!!$            END IF
!!$!=================|DEBUG|============================

      END DO
      
      

      ! NOW ADJUSTMENT ELF UNDER THERMAL WIND SUCH THAT THE BOTTOM
      ! VELOCITY IS ZERO 
      allocate(eta_lcl(NOBCLSF)); eta_lcl=0.0_sp
! DEBUG
!      allocate(temp_lcl(NOBCLSF)); temp_lcl=0.0_sp

      RHOINTCUR= 0.0_SP
      RHOINTNXT = 0.0_SP

      DO I=NOBCLSF,1,-1  ! COUNT BACKWARDS - The onshore direction
         ndx  = NBCLSF(I)
         cdx = IBCLSF(i)

         !INTEGRATE RHO IN DEPTH, CONVERT TO MKS UNITS DAMIT
         RHOINTCUR=RHOINTNXT
         RHOINTNXT=0.0_SP
         DO K=1,KBM1
            RHOINTNXT=RHOINTNXT+(1.0_SP+RHO1(CDX,K))*1.0E3_SP*DZ(CDX,K)
         END DO
! DEBUG         
!         temp_lcl(I)=RHOINTNXT

!ESTIMATE DENSITY GRADIENT, AND MODIFY BOUNDARY ELEVATION
!NOTE THE FACTOR OF 1000 AND 2 TO COMPENSATE FOR THE
!FACT THAT THE MODEL STORES RHO1 AS SIGMA, AND IN CGS.

! ADJUST THE SEA SURFACE HEIGHT FOR A MEAN FLOW DUE TO DENSITY GRADIENT
! CALCULATE THE THERMAL WIND AND ADJUST THE SEA SURFACE HEIGHT SO THAT THE BOTTOM
! BOUNDARY IS A LEVEL OF NO MOTION:
!
! ETATAN=[ -1/(Rho_bar)*(del(h*rho)/del(s) - Rho_bot*(del(h)/del(s) ] *del(s)
!
         IF (I /= nobclsf) THEN
            ETA_LCL(I)=-(1.0_SP/(0.5_SP*(RHOINTNXT+RHOINTCUR))) &
                 *((H(CDX)*RHOINTNXT-H(NDX)*RHOINTCUR) &
                 -0.5_SP*1.0e3_SP*(2.0_SP+RHO1(CDX,KBM1)+RHO1(NDX,KBM1)) &
                 *(H(CDX)-H(NDX)))
         END IF

      END DO

      CUMEL = 0.0_SP
      
      IF(SERIAL) THEN
         DO I=NOBCLSF,1,-1 ! COUNT BACKWARD FROM THE SHELF TOWARD SHORE
            NDX=ibclsf(I)
            CUMEL=CUMEL+ETA_LCL(I)
            ELF(NDX)=ELF(NDX)+ CUMEL *RBC_GEO(I)
!            write(ipt,*) "NDX=",NDX,"ELF(NDX)=",ELF(NDX)
         END DO
      END IF

# if defined (MULTIPROCESSOR)
      IF(PAR)THEN
         IF(MSR) THEN
            ALLOCATE(ETA_GBL(NOBCLSF_GL)); ETA_GBL = 0.0_SP
            ALLOCATE(ELFGEO_GBL(NOBCLSF_GL)); ELFGEO_GBL = 0.0_SP
! DEBUG
!            allocate(temp_gl(NOBCLSF_GL)); temp_gl=0.0_sp
         END IF

         CALL PCOLLECT(MYID,MSRID,NPROCS,LSFMAP,ETA_LCL,ETA_GBL)
! DEBUG
!         CALL PCOLLECT(MYID,MSRID,NPROCS,LSFMAP,temp_lcl,temp_gl)
         
         IF (MSR) THEN
                        
            DO I=NOBCLSF_GL,1,-1 ! COUNT BACKWARD FROM THE SHELF TOWARD SHORE
               CUMEL=CUMEL+ETA_GBL(I)
               ELFGEO_GBL(I)=CUMEL
            END DO

!!$!=================|DEBUG|============================
!!$
!!$            write(ipt,*) "writing file"
!!$
!!$            if(.not. INIT)then
!!$               init =.true.
!!$               icount=0
!!$            end if
!!$            
!!$            write(ccount,'(I5.5)') icount
!!$
!!$            call FOPEN(157,"LSF_"//ccount,'ofr')
!!$
!!$            write(157,*) "LSF DEBUG: ",icount
!!$            write(157,*) "one dens =",rho1(1,:)
!!$            write(157,*) "Node#, ETA,  ELFGEO"
!!$            do I=NOBCLSF_GL,1,-1
!!$               write(157,*) I, temp_gl(I), ETA_GBL(I), ELFGEO_GBL(I)
!!$            end do
!!$            
!!$            icount = icount +1
!!$
!!$            close(157)
!!$            
!!$!=================|DEBUG|============================

         END IF

         ALLOCATE(ELFGEO_LCL(NOBCLSF)); ELFGEO_LCL = 0.0_SP

         CALL PDEAL(MYID,MSRID,NPROCS,LSFMAP,ELFGEO_GBL,ELFGEO_LCL)

         DO I=NOBCLSF,1,-1 ! COUNT BACKWARD FROM THE SHELF TOWARD SHORE
            NDX=ibclsf(I)
            ELF(NDX)=ELF(NDX)+ELFGEO_LCL(I)*RBC_GEO(I)            
         END DO

         IF (MSR) THEN
            DEALLOCATE(ETA_GBL)
            DEALLOCATE(ELFGEO_GBL)
! DEBUG
!            deallocate(temp_gl)
         END IF

         DEALLOCATE(elfgeo_lcl)
      END IF

# endif
! DEBUG
!      deallocate(temp_lcl)
      DEALLOCATE(eta_lcl)

   END IF


!==============================================================================|
   CASE(2) !External Mode Velocity Boundary Conditions                         |
!==============================================================================|

   DO I=1,N

!
!--2 SOLID BOUNDARY EDGES------------------------------------------------------|
!
   IF(ISBCE(I) == 3) THEN
     UAF(I)=0.0_SP
     VAF(I)=0.0_SP
   END IF

!  GWC SPEED REPLACE ABOVE 4 LINES
!   UAF(LISBCE_3(1:NISBCE_3)) = 0.  
!   VAF(LISBCE_3(1:NISBCE_3)) = 0.  

!
!--1 SOLID BOUNDARY EDGE-------------------------------------------------------|
!
   IF(ISBCE(I) == 1) THEN
     ALPHA1=ALPHA(I)
     IF(NUMQBC > 0) THEN
       IF(RIVER_INFLOW_LOCATION == 'node') THEN
         DO J=1,NUMQBC
           I1=INODEQ(J)
           J1=NBVE(I1,1)
           J2=NBVE(I1,NTVE(I1))
           IF((I == J1).OR.(I == J2)) THEN
             UNTMP=UAF(I)*COS(ANGLEQ(J))+VAF(I)*SIN(ANGLEQ(J))
             VNTMP=-UAF(I)*SIN(ANGLEQ(J))+VAF(I)*COS(ANGLEQ(J))
             UNTMP=MAX(UNTMP,0.0_SP)
             UAF(I)=UNTMP*COS(ANGLEQ(J))-VNTMP*SIN(ANGLEQ(J))
             VAF(I)=UNTMP*SIN(ANGLEQ(J))+VNTMP*COS(ANGLEQ(J))
               GOTO 21
               END IF 
             END DO
            ELSE IF(RIVER_INFLOW_LOCATION == 'edge') THEN
             DO J=1,NUMQBC
               J1=ICELLQ(J)
               IF(I == J1) THEN
                UNTMP=UAF(I)*COS(ANGLEQ(J))+VAF(I)*SIN(ANGLEQ(J))
                VNTMP=-UAF(I)*SIN(ANGLEQ(J))+VAF(I)*COS(ANGLEQ(J))
                UNTMP=MAX(UNTMP,0.0_SP)
                UAF(I)=UNTMP*COS(ANGLEQ(J))-VNTMP*SIN(ANGLEQ(J))
                VAF(I)=UNTMP*SIN(ANGLEQ(J))+VNTMP*COS(ANGLEQ(J))
               GOTO 21
               END IF
             END DO
           END IF
            
           END IF 

#          if defined (THIN_DAM)
           IF(NBE_DAM(I) == 0)THEN
#          endif

           UI= UAF(I)*(SIN(ALPHA1))**2-VAF(I)*SIN(ALPHA1)*COS(ALPHA1)
           VI=-UAF(I)*SIN(ALPHA1)*COS(ALPHA1)+VAF(I)*(COS(ALPHA1))**2
           UAF(I)=UI
           VAF(I)=VI

!!#          if defined (WET_DRY)
!!           UI= UAS(I)*(SIN(ALPHA1))**2-VAS(I)*SIN(ALPHA1)*COS(ALPHA1)
!!           VI=-UAS(I)*SIN(ALPHA1)*COS(ALPHA1)+VAS(I)*(COS(ALPHA1))**2
!!           UAS(I)=UI
!!           VAS(I)=VI
!!#          endif

#          if defined (THIN_DAM)
           ELSE
               UNTMP=-UAF(I)*COS(ALPHA1)-VAF(I)*SIN(ALPHA1)
               VNTMP=UAF(I)*SIN(ALPHA1)-VAF(I)*COS(ALPHA1)
               IF(KDAM1(I)<=KBM1/2)UNTMP=UNTMP*KDAM1(I)/KBM1
               UAF(I)=-UNTMP*COS(ALPHA1)+VNTMP*SIN(ALPHA1)
               VAF(I)=-UNTMP*SIN(ALPHA1)-VNTMP*COS(ALPHA1)
           END IF
#          endif

21        CONTINUE
          END IF
       END DO

!!# if defined (THIN_DAM)
!!    DO I=1,NT
!!      IF(CLP_CELL(I)==1)THEN
!!        ALPHA1=CLP_ALPHA(I)
!!        UNTMP=-UAF(I)*COS(ALPHA1)-VAF(I)*SIN(ALPHA1)
!!        VNTMP=UAF(I)*SIN(ALPHA1)-VAF(I)*COS(ALPHA1)
!!        IF(CLP_KDAM(I,1)==1)KDAM_TMP=KDAM1(CLP_KDAM(I,2))
!!        IF(CLP_KDAM(I,1)==2)KDAM_TMP=(KDAM1(CLP_KDAM(I,2))+KDAM1(CLP_KDAM(I,2)))/2
!!        IF(KDAM_TMP<=KBM1/2)UNTMP=UNTMP*KDAM_TMP/KBM1
!!        UAF(I)=-UNTMP*COS(ALPHA1)+VNTMP*SIN(ALPHA1)
!!        VAF(I)=-UNTMP*SIN(ALPHA1)-VNTMP*COS(ALPHA1)
!!      END IF
!!    END DO
!!# endif

!==============================================================================|
   CASE(3) !3-D Velocity Boundary Conditions                                   !
!==============================================================================|

!  GWC SPEED REPLACE NEXT 4 LINES
!   UF(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.  
!   VF(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.  

      ! Suggest move the kbm1 loop inside 'if((i.eq.j1).or.(i.eq.j2)) then'

      ! THIS SHOULD BE SAFE FOR MORE THAN ONE RIVER AT THE SAME NODE
      ! OR EDGE

       do i= 1, n
         do k =1, kbm1
          if(isbce(i).eq.3) then
           uf(i,k)=0.0_SP
           vf(i,k)=0.0_SP
          end if

          if(isbce(i).eq.1) then
          alpha1=alpha(i)
          if(numqbc.ge.1) then
          if(river_inflow_location.eq.'node') then
            do j=1,numqbc
              i1=inodeq(j)
              j1=nbve(i1,1)
              j2=nbve(i1,ntve(i1))
              if((i.eq.j1).or.(i.eq.j2)) then
               untmp=uf(i,k)*cos(angleq(j))+vf(i,k)*sin(angleq(j))
               vntmp=-uf(i,k)*sin(angleq(j))+vf(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               uf(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               vf(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 31
              end if
            end do
           else if(river_inflow_location.eq.'edge') then
            do j=1,numqbc
              j1=icellq(j)
              if(i.eq.j1) then
               untmp=uf(i,k)*cos(angleq(j))+vf(i,k)*sin(angleq(j))
               vntmp=-uf(i,k)*sin(angleq(j))+vf(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               uf(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               vf(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 31
              end if
            end do
           else
             print*, 'river_inflow_location not correct'
             call pstop
          end if
          end if

#         if defined (THIN_DAM)
          IF(NBE_DAM(I) == 0.OR.K > KDAM1(I))THEN
!          IF(NBE_DAM(I) == 0)THEN
#         endif

          ui= uf(i,k)*(sin(alpha1))**2-vf(i,k)*sin(alpha1)*cos(alpha1)
          vi=-uf(i,k)*sin(alpha1)*cos(alpha1)+vf(i,k)*(cos(alpha1))**2
          uf(i,k)=ui
          vf(i,k)=vi

!!#         if defined (WET_DRY)
!!          ui= us(i,k)*(sin(alpha1))**2-vs(i,k)*sin(alpha1)*cos(alpha1)
!!          vi=-us(i,k)*sin(alpha1)*cos(alpha1)+vs(i,k)*(cos(alpha1))**2
!!          us(i,k)=ui
!!          vs(i,k)=vi
!!#         endif

#         if defined (THIN_DAM)
          end if
#         endif

31        continue
          end if
         end do
       end do



!==============================================================================|
   CASE(4) !Blank                                                              !
!==============================================================================|

!==============================================================================|
   CASE(5) !!SOLID BOUNDARY CONDITIONS ON U AND V                              !
!==============================================================================|


!  GWC SPEED REPLACE NEXT 4 LINES
!   U(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.  
!   V(LISBCE_3(1:NISBCE_3),1:KBM1) = 0.  

       do i= 1, n
         do k =1, kbm1
          if(isbce(i).eq.3) then
          u(i,k)=0.0_SP
          v(i,k)=0.0_SP
          end if

          if(isbce(i).eq.1) then
          alpha1=alpha(i)
          if(numqbc.ge.1) then
          if(river_inflow_location.eq.'node') then
            do j=1,numqbc
              i1=inodeq(j)
              j1=nbve(i1,1)
              j2=nbve(i1,ntve(i1))
              if((i.eq.j1).or.(i.eq.j2)) then
               untmp=u(i,k)*cos(angleq(j))+v(i,k)*sin(angleq(j))
               vntmp=-u(i,k)*sin(angleq(j))+v(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               u(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               v(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 51
              end if
            end do
           else if(river_inflow_location.eq.'edge') then
            do j=1,numqbc
              j1=icellq(j)
              if(i.eq.j1) then
               untmp=u(i,k)*cos(angleq(j))+v(i,k)*sin(angleq(j))
               vntmp=-u(i,k)*sin(angleq(j))+v(i,k)*cos(angleq(j))
               untmp=max(untmp,0.0_SP)
               u(i,k)=untmp*cos(angleq(j))-vntmp*sin(angleq(j))
               v(i,k)=untmp*sin(angleq(j))+vntmp*cos(angleq(j))
              goto 51
              end if
            end do
           else
             print*, 'river_inflow_location not correct'
             call pstop
          end if
          end if

#         if defined (THIN_DAM)
          IF(NBE_DAM(I) == 0.OR.K > KDAM1(I))THEN          
!          IF(NBE_DAM(I) == 0)THEN
#         endif

          ui= u(i,k)*(sin(alpha1))**2-v(i,k)*sin(alpha1)*cos(alpha1)
          vi=-u(i,k)*sin(alpha1)*cos(alpha1)+v(i,k)*(cos(alpha1))**2
          u(i,k)=ui
          v(i,k)=vi

!!#         if defined (WET_DRY)
!!          ui= us(i,k)*(sin(alpha1))**2-vs(i,k)*sin(alpha1)*cos(alpha1)
!!          vi=-us(i,k)*sin(alpha1)*cos(alpha1)+vs(i,k)*(cos(alpha1))**2
!!          us(i,k)=ui
!!          vs(i,k)=vi
!!#         endif

#         if defined (THIN_DAM)
          END IF
#         endif

51        continue
          end if
         end do
       end do



!==============================================================================|
   CASE(6) !Blank                                                              !
!==============================================================================|

!==============================================================================|
   CASE(7) !Blank                                                              !
!==============================================================================|

!==============================================================================|
   CASE(8) !!SURFACE FORCING FOR INTERNAL MODE                                 !
!==============================================================================|

!
!--Fresh Water Discharge-------------------------------------------------------|
!

      IF(NUMQBC_GL .GT. 0) THEN
# if defined (WATER_QUALITY) && !defined(SEDIMENT)
         CALL UPDATE_RIVERS(IntTime,QDIS,TEMP=TDIS,SALT=SDIS,WQM=WQMDIS)
# endif
# if defined (WATER_QUALITY) && defined(SEDIMENT)
         CALL UPDATE_RIVERS(IntTime,QDIS,TEMP=TDIS,SALT=SDIS,WQM=WQMDIS,SED=SEDDIS)
# endif
# if defined(SEDIMENT) && !defined(WATER_QUALITY)
         CALL UPDATE_RIVERS(IntTime,QDIS,TEMP=TDIS,SALT=SDIS,SED=SEDDIS)
# endif
# if defined(BioGen)
         IF(BIOLOGICAL_MODEL)CALL UPDATE_RIVERS(IntTime,QDIS,BIO=BIODIS)
# endif
# if defined(FABM)
         IF (FABM_MODEL) CALL UPDATE_RIVERS(IntTime,QDIS,FABMR=FABMDIS)
# endif
# if !defined(SEDIMENT) && !defined(WATER_QUALITY) 
         CALL UPDATE_RIVERS(IntTime,QDIS,TEMP=TDIS,SALT=SDIS)
# endif	 
                  
         QDIS    = QDIS*RAMP
         
!#  if defined (WATER_QUALITY)
!         DO N1 = 1, NB
!            ! NOT YET OPERATIONAL!!!
!            WDIS(:,N1) = UFACT*DWDIS(:,N1,L1) + FACT*DWDIS(:,N1,L2)
!         END DO
!#  endif
         
      END IF
      
      
      
      IF (WIND_ON) THEN
         IF (WIND_TYPE == SPEED)THEN
            CALL UPDATE_WIND(IntTime,UUWIND,VVWIND)
            
            CALL ASIMPLE_DRAG(UUWIND,VVWIND,WUSURF,WVSURF)

         ELSEIF(WIND_TYPE == STRESS)THEN
            CALL UPDATE_WIND(IntTime,WUSURF,WVSURF)
         END IF
         

         ! For output only - don't divide by density
         WUSURF_save = WUSURF * RAMP
         WVSURF_save = WVSURF * RAMP

         ! Divide by density 
         WUSURF = -WUSURF * RAMP *0.001_SP
         WVSURF = -WVSURF * RAMP *0.001_SP
         
         ! MAJOR MISTAKE: NEED THE NEGATIVE SIGN FOR THE INTERNAL
         ! STEP WIND STRESS

      END IF
      
# if !defined (HEATING_CALCULATED) && !defined (HEATING_CALCULATED_GL) && !defined (HEATING_SOLAR)
      IF (HEATING_ON) THEN
         CALL UPDATE_HEAT(IntTime,SWRAD_WATTS,WTSURF_WATTS)
         
         ! LOAD INTO THES VARIABLE TO SAVE THE FORCING IN THE CORRECT UNITS
         SWRAD_WATTS  = SWRAD_WATTS * RAMP
         WTSURF_WATTS = WTSURF_WATTS * RAMP

         SPCP=4.2174E3_SP
         ROSEA = 1.023E3_SP
         SPRO = SPCP*ROSEA
         
         WTSURF    = -WTSURF_WATTS/SPRO
         SWRAD     = -SWRAD_WATTS/SPRO

# if defined (MULTIPROCESSOR)
         IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SWRAD_WATTS,WTSURF_WATTS)         
# endif

      END IF
# endif      
      
# if defined (HEATING_CALCULATED) || (HEATING_CALCULATED_GL)
      IF (HEATING_CALCULATE_ON) THEN
         CALL UPDATE_HEAT_CALCULATED(IntTime,SWRAD_WATTS,WTSURF_WATTS)
         
         ! LOAD INTO THES VARIABLE TO SAVE THE FORCING IN THE CORRECT UNITS
         SWRAD_WATTS  = SWRAD_WATTS * RAMP
         WTSURF_WATTS = WTSURF_WATTS * RAMP

# if defined (HEATING_CALCULATED)
         SPCP=4.2174E3_SP
         ROSEA = 1.023E3_SP
         SPRO = SPCP*ROSEA
# else
         SPCP= 4.186E3_SP
         ROSEA = 1.000E3_SP
         SPRO = SPCP*ROSEA
# endif
         
         WTSURF    = -WTSURF_WATTS/SPRO
         SWRAD     = -SWRAD_WATTS/SPRO
	 
# if defined (MULTIPROCESSOR)
         IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SWRAD_WATTS,WTSURF_WATTS)         
# endif

      END IF
# endif      
# if defined (HEATING_SOLAR)
      IF (HEATING_SOLAR_ON) THEN
         CALL UPDATE_HEAT_SOLAR(IntTime,SWRAD_WATTS,WTSURF_WATTS)
         
         ! LOAD INTO THES VARIABLE TO SAVE THE FORCING IN THE CORRECT UNITS
         SWRAD_WATTS  = SWRAD_WATTS * RAMP
         WTSURF_WATTS = WTSURF_WATTS * RAMP

         ROSEA1 = (-.0057)*T1(:,1)*T1(:,1) + 0.0249*T1(:,1) + 1000
         !IF(TWAT .LT. 3.9) WRITE(IPT,*)'WATER TEMP BELOW 4 C, BEWARE DENSITY CALCULATION!, TWAT = ',TWAT
         !IF(TWAT .GT. 30.1) WRITE(IPT,*)'WATER TEMP ABOVE 30 C, BEWARE DENSITY CALCULATION!, TWAT = ',TWAT
         SPCP  = 4186
         !SPRO  = SPCP*ROSEA
         !SPRO  = SPCP*ROSEA
         SPRO1  = SPCP*ROSEA1
                    
!         SPCP=4.2174E3_SP
!         ROSEA = 1.023E3_SP
!         SPRO = SPCP*ROSEA
         
         WTSURF    = -WTSURF_WATTS/SPRO1
         SWRAD     = -SWRAD_WATTS/SPRO1
	 
# if defined (MULTIPROCESSOR)
         IF(PAR)CALL AEXCHANGE(NC,MYID,NPROCS,SWRAD_WATTS,WTSURF_WATTS)         
# endif

      END IF
# endif      
      
      IF (PRECIPITATION_ON) THEN
         CALL UPDATE_PRECIPITATION(IntTime,Qprec,Qevap)
         ! NO RAMP FOR PRECIP/EVAP
      END IF


# if defined (WAVE_CURRENT_INTERACTION) && (WAVE_OFFLINE)
      IF (WAVE_ON) THEN
         CALL UPDATE_WAVE(IntTime,WHS,WDIR,WPER,WLENGTH,WPER_BOT,WUB_BOT)
      END IF
# endif

# if defined (ICENUDGE)
     IF (ICENUDGE_ON)  THEN
        CALL UPDATE_ICENUDGE(IntTime,AICE,HICE,ISALT,UICE,VICE)
     END IF
# endif

# if defined (OFFLINE_SEDIMENT)
      CALL UPDATE_OFFLINE_FORCING(IntTime,offline_u,  offline_v, offline_w&
                                        &,offline_s1, offline_t1,offline_el,offline_kh&
# if defined (GOTM)
                                        &,offline_teps &
# else
                                        &,offline_q2, offline_q2l &
# endif
                                        &,offline_wn, offline_wc)
# endif

# if defined (FABM) && defined (OFFLINE_FABM)
      IF (FABM_MODEL) CALL UPDATE_OFFLINE_FABM_FORCING(IntTime,offlineb_u,  offlineb_v, offlineb_w&
                                        &,offlineb_s1,offlineb_t1, offlineb_el,offlineb_kh&
                                        &,offlineb_wn, offlineb_wc,offlineb_uard_obcn)
# endif

# if defined (OFFLINE_BIOLOGY)
      CALL UPDATE_OFFLINE_BIOLOGY_FORCING(IntTime,offline_u,  offline_v, offline_w&
                                        &,offline_s1, offline_t1,offline_el,offline_kh&
# if defined (GOTM)
                                        &,offline_teps &
# else
                                        &,offline_q2, offline_q2l &
# endif
                                        &,offline_wn, offline_wc)
# endif

!
!-- Set Groundwater flux ------------------------------------------------------|
!
      
      IF (GROUNDWATER_ON) THEN
         CALL UPDATE_GROUNDWATER(IntTime,BFWDIS,GW_TEMP=BFWTMP,GW_SALT=BFWSLT)
         BFWDIS = RAMP * BFWDIS
         ! DO NOT RAMP TEMP AND SALINITY
      END IF
      
!==============================================================================|
   CASE(9) !External Mode Surface BCs (River Flux/Wind Stress/Heat/Moist)      !
!==============================================================================|
# if !defined (SEMI_IMPLICIT)
!
!-Freshwater Flux: Set  Based on Linear Interpolation Between Two Data Times---|
!

      IF(NUMQBC_GL .GT. 0) THEN
         CALL UPDATE_RIVERS(ExtTime,QDIS2)
         
         QDIS2 = QDIS2*RAMP
         
      END IF
      
!
!-- Set Precipitation/Evaporation/Surface Wind ---------------------------------|
!

      IF (PRECIPITATION_ON) THEN
         CALL UPDATE_PRECIPITATION(ExtTime,Qprec2,Qevap2)
         ! NO RAMP FOR PRECIP/EVAP
      END IF
      
      IF (WIND_ON) THEN

         IF (WIND_TYPE == SPEED)THEN
            CALL UPDATE_WIND(ExtTime,UUWIND,VVWIND)
            
            CALL ASIMPLE_DRAG(UUWIND,VVWIND,WUSURF2,WVSURF2)

         ELSEIF(WIND_TYPE == STRESS)THEN
            CALL UPDATE_WIND(ExtTime,WUSURF2,WVSURF2)
         END IF

# if defined (TWO_D_MODEL)
         ! For output only - don't divide by density
         WUSURF_save = WUSURF2 * RAMP
         WVSURF_save = WVSURF2 * RAMP
# endif	 

         WUSURF2 = WUSURF2 * RAMP *0.001_SP
         WVSURF2 = WVSURF2 * RAMP *0.001_SP
         
      END IF
   
!
!-- Set Groundwater flux ------------------------------------------------------|
!

      IF (GROUNDWATER_ON) THEN
         CALL UPDATE_GROUNDWATER(ExtTime,BFWDIS2)
         BFWDIS2 = RAMP * BFWDIS2
      END IF
# endif      
   END SELECT
   
   if(dbg_set(dbg_sbr)) write(ipt,*)&
        & "End: bcond_gcn"
   
 END SUBROUTINE BCOND_GCN
!==============================================================================|


