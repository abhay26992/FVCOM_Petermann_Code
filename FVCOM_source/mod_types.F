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

MODULE MOD_TYPES
  USE MOD_PREC
  IMPLICIT NONE
  TYPE MAP
     INTEGER :: NSIZE =0 ! SIZE OF LOCAL DATA (Size of Local data)
     INTEGER :: GSIZE =0 ! GLOBAL SIZE (SIZE OF GLOBAL DATA)
     INTEGER :: LSIZE =0 ! SIZE OF LOCAL ARRAY (Usually ==nsize or MT,NT)
     ! INDEX FROM LOCAL ARRAY TO GLOBAL ARRAY
     INTEGER, POINTER,  DIMENSION(:) :: LOC_2_GL => NULL()
     ! INDEX FROM LOCAL ARRAY TO LOCAL GRID
     INTEGER, POINTER,  DIMENSION(:) :: LOC_2_Grid => NULL()
  END TYPE MAP

  TYPE COMM
     !----------------------------------------------------------
     ! SND: TRUE IF YOU ARE TO SEND TO PROCESSOR               |
     ! RCV: TRUE IF YOU ARE TO RECEIVE FROM PROCESSOR          |
     ! NSND: NUMBER OF DATA TO SEND TO PROCESSOR               |
     ! NRCV: NUMBER OF DATA TO RECEIVE FROM PROCESSOR          |
     ! SNDP: ARRAY POINTING TO LOCATIONS TO SEND TO PROCESSOR  | 
     ! RCVP: ARRAY POINTING TO LOCATIONS RECEIVED FROM PROCESS |
     ! RCPT: POINTER TO LOCATION IN RECEIVE BUFFER             |
     !----------------------------------------------------------

     !  LOGICAL :: SND,RCV
     INTEGER  :: NSND =0,NRCV=0,RCPT=0
     INTEGER, POINTER,  DIMENSION(:) :: SNDP => NULL(),RCVP => NULL()
     REAL(SP), POINTER,   DIMENSION(:) :: MLTP => NULL()
  END TYPE COMM

  TYPE BC
     INTEGER :: NTIMES =0
     REAL(SP), POINTER,  DIMENSION(:) :: TIMES => NULL()
     CHARACTER(LEN=80) :: LABEL
  END TYPE BC


  ! THIS IS NOT REALLY A GRID - IT ONLY HOLDS INTEGER INDEX DATA THAT
  ! IS GRID DEPENDENT
  TYPE GRID
     CHARACTER(len=80):: UNITS
     CHARACTER(len=160):: NAME
     INTEGER :: M   =0
     INTEGER :: N   =0
     INTEGER :: MT  =0
     INTEGER :: NT  =0
     INTEGER :: MGL =0
     INTEGER :: NGL =0
     INTEGER :: KB  =0
     INTEGER :: KBM1=0
     INTEGER :: KBM2=0

     
     ! PARALLEL STUFF
!     INTEGER, POINTER :: EL_PID(:) !!PROCESSOR OWNER OF GLOBAL ELEMENT

     INTEGER, POINTER :: egid(:) => NULL() !!GLOBAL ID OF LOCAL ELEMENT
     INTEGER, POINTER :: ngid(:) => NULL() !!GLOBAL ID OF LOCAL NODE 
     INTEGER, POINTER :: elid(:) => NULL() !!LOCAL ID OF GLOBAL ELEMENT
     INTEGER, POINTER :: nlid(:) => NULL() !!LOCAL ID OF GLOBAL NODE 

     ! TRANSFER INDEX INCLUDING HALO
     INTEGER, POINTER :: egid_X(:) => NULL()
     INTEGER, POINTER :: ngid_X(:) => NULL()
     INTEGER, POINTER :: elid_X(:) => NULL()
     INTEGER, POINTER :: nlid_X(:) => NULL()

     ! GRID DATA POINTERS

     INTEGER, POINTER :: NV(:,:) => NULL()


     REAL(SP),POINTER :: XM(:) => NULL()
     REAL(SP),POINTER :: YM(:) => NULL()
     REAL(SP),POINTER :: LON(:) => NULL()
     REAL(SP),POINTER :: LAT(:) => NULL()

     REAL(SP),POINTER :: XMC(:) => NULL()
     REAL(SP),POINTER :: YMC(:) => NULL()
     REAL(SP),POINTER :: LONC(:) => NULL()
     REAL(SP),POINTER :: LATC(:) => NULL()

     ! NOTES: SHOULD MAKE AN ARRAY TO STORE 1/ART, 1/ART2 and 1/ART2
     ! IT is faster and safer

!     REAL(SP), POINTER :: ART(:)              !!AREA OF ELEMENT
!     REAL(SP), POINTER :: ART1(:)             !!AREA OF NODE-BASE CONTROl VOLUME
!     REAL(SP), POINTER :: ART2(:)             !!AREA OF ELEMENTS AROUND NODE

     !----------------Node, Boundary Condition, and Control Volume-----------------------!


     INTEGER, POINTER :: NBE(:,:) => NULL()            !!INDICES OF ELMNT NEIGHBORS
 !    INTEGER, POINTER :: NTVE(:)             !! NUMBER OF ELEMENTS SURROUNDING EACH NODE
 !    INTEGER, POINTER :: NTSN(:)             !! NUMBER OF NODES SURROUNDING EACH NODE
 !    INTEGER, POINTER :: ISONB(:)            !!NODE MARKER = 0,1,2
 !    INTEGER, POINTER :: ISBC(:)     
 !    INTEGER, POINTER :: ISBCE(:)     
 !    INTEGER, POINTER :: IEC(:,:)
 !    INTEGER, POINTER :: IENODE(:,:)
     INTEGER, POINTER :: NBSN(:,:) => NULL()            !! INDICES OF NODES SURROUNDING EACH NODE
!     INTEGER, POINTER :: NIEC(:,:)
!     INTEGER, POINTER :: NTRG(:)
     INTEGER, POINTER :: NBVE(:,:) => NULL()            !! INDICIES OF ELEMENTS SURROUNDING EACH NODE
!     INTEGER, POINTER :: NBVT(:,:)

     !----------------2-d arrays for the general vertical coordinate -------------------------------!

     REAL(SP), POINTER :: H(:) => NULL()                      !!BATHYMETRY 
     REAL(SP), POINTER :: Z(:,:) => NULL()                    !!SIGMA COORDINATE VALUE 
     REAL(SP), POINTER :: ZZ(:,:) => NULL()                   !!INTRA LEVEL SIGMA VALUE
     REAL(SP), POINTER :: H1(:) => NULL()                      !!BATHYMETRY 
     REAL(SP), POINTER :: Z1(:,:) => NULL()                    !!SIGMA COORDINATE VALUE 
     REAL(SP), POINTER :: ZZ1(:,:) => NULL()                   !!INTRA LEVEL SIGMA VALUE

     !------------shape coefficient arrays and control volume metrics--------------------!
     
!     REAL(SP), POINTER :: A1U(:,:)      
!     REAL(SP), POINTER :: A2U(:,:)     
!     REAL(SP), POINTER :: AWX(:,:)   
!     REAL(SP), POINTER :: AWY(:,:)  
!     REAL(SP), POINTER :: AW0(:,:) 
!     REAL(SP), POINTER :: ALPHA(:)
     

     ! OBC DATA

!     INTEGER :: IOBCN =0
!     INTEGER, POINTER :: I_OBC_N(:)
!     INTEGER, POINTER :: TYPE_OBC(:)

     ! END DATA POINTERS FOR GRID TYPE

# if defined (WAVE_CURRENT_INTERACTION)
     INTEGER :: MSC =0
     INTEGER :: MDC =0
# endif
  END TYPE GRID

  CONTAINS

    SUBROUTINE KILL_GRID(G)
      IMPLICIT NONE
      TYPE(GRID), POINTER :: G
      

      IF(.NOT. ASSOCIATED(G)) RETURN

!      IF(ASSOCIATED(G%EL_PID)) DEALLOCATE(G%EL_PID)

      IF(ASSOCIATED(G%egid)) DEALLOCATE(G%egid)
      IF(ASSOCIATED(G%ngid)) DEALLOCATE(G%ngid)
      IF(ASSOCIATED(G%elid)) DEALLOCATE(G%elid)
      IF(ASSOCIATED(G%nlid)) DEALLOCATE(G%nlid)

      IF(ASSOCIATED(G%egid_X)) DEALLOCATE(G%egid_X)
      IF(ASSOCIATED(G%elid_X)) DEALLOCATE(G%elid_X)
      IF(ASSOCIATED(G%ngid_X)) DEALLOCATE(G%ngid_X)
      IF(ASSOCIATED(G%nlid_X)) DEALLOCATE(G%nlid_X)

      IF(ASSOCIATED(G%NV)) DEALLOCATE(G%NV)

      IF(ASSOCIATED(G%XM)) DEALLOCATE(G%XM)
      IF(ASSOCIATED(G%YM)) DEALLOCATE(G%YM)
      IF(ASSOCIATED(G%LON)) DEALLOCATE(G%LON)
      IF(ASSOCIATED(G%LAT)) DEALLOCATE(G%LAT)

      IF(ASSOCIATED(G%XMC)) DEALLOCATE(G%XMC)
      IF(ASSOCIATED(G%YMC)) DEALLOCATE(G%YMC)
      IF(ASSOCIATED(G%LONC)) DEALLOCATE(G%LONC)
      IF(ASSOCIATED(G%LATC)) DEALLOCATE(G%LATC)

!      IF(ASSOCIATED(G%ART)) DEALLOCATE(G%ART)
!      IF(ASSOCIATED(G%ART1)) DEALLOCATE(G%ART1)
!      IF(ASSOCIATED(G%ART2)) DEALLOCATE(G%ART2)

      IF(ASSOCIATED(G%NBE)) DEALLOCATE(G%NBE)
!      IF(ASSOCIATED(G%NTVE)) DEALLOCATE(G%NTVE)
!      IF(ASSOCIATED(G%NTSN)) DEALLOCATE(G%NTSN)
!      IF(ASSOCIATED(G%ISONB)) DEALLOCATE(G%ISONB)
!      IF(ASSOCIATED(G%ISBC)) DEALLOCATE(G%ISBC)
!      IF(ASSOCIATED(G%IEC)) DEALLOCATE(G%IEC)
!      IF(ASSOCIATED(G%IENODE)) DEALLOCATE(G%IENODE)
      IF(ASSOCIATED(G%NBSN)) DEALLOCATE(G%NBSN)
!      IF(ASSOCIATED(G%NIEC)) DEALLOCATE(G%NIEC)
!      IF(ASSOCIATED(G%NTRG)) DEALLOCATE(G%NTRG)
      IF(ASSOCIATED(G%NBVE)) DEALLOCATE(G%NBVE)
!      IF(ASSOCIATED(G%NBVT)) DEALLOCATE(G%NBVT)

!      IF(ASSOCIATED(G%Z)) DEALLOCATE(G%Z)
!      IF(ASSOCIATED(G%ZZ)) DEALLOCATE(G%ZZ)
 
!      IF(ASSOCIATED(G%A1U)) DEALLOCATE(G%A1U)
!      IF(ASSOCIATED(G%A2U)) DEALLOCATE(G%A2U)
!      IF(ASSOCIATED(G%AWX)) DEALLOCATE(G%AWX)
!      IF(ASSOCIATED(G%AWY)) DEALLOCATE(G%AWY)
!      IF(ASSOCIATED(G%AW0)) DEALLOCATE(G%AW0)
!      IF(ASSOCIATED(G%ALPHA)) DEALLOCATE(G%ALPHA)

!      IF(ASSOCIATED(G%I_OBC_N)) DEALLOCATE(G%I_OBC_N)
!      IF(ASSOCIATED(G%TYPE_OBC)) DEALLOCATE(G%TYPE_OBC)

      
      DEALLOCATE(G)


    END SUBROUTINE KILL_GRID


END MODULE MOD_TYPES

