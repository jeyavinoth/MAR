!     ######spl
    SUBROUTINE CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,             &
                           PDTCONV, ODEEP, OSHAL, OREFRESH_ALL, ODOWN, KICE,   &
                           OSETTADJ, PTADJD, PTADJS,                           &
                           KENSM,                                              &
                           PPABS, PZZ, PDXDY,                                  &
                           PT, PRV, PRC, PRI, PU, PV, PW,                      &
                           KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,              &
                           PPRTEN, PPRSTEN,                                    &
                           PUMF, PDMF, PPRLFLX, PPRSFLX, PCAPE, KCLTOP, KCLBAS,& 
                           OCHTRANS, KCH1, PCH1, PCH1TEN                       )
!   ############################################################################
!
!!**** Interface routine to the fast MNH convection code developed for ECMWF/ARPEGE IFS
!!     having a structure typical for operational routines
!!
!!
!!    PURPOSE
!!    -------
!!      The routine interfaces the MNH convection code as developed for operational
!!      forecast models like ECMWF, ARPEGE or HIRLAM with the typical MNH array structure
!!      Calls the deep and/or shallow convection routine
!!
!!
!!**  METHOD
!!    ------
!!     Returns one tendency for shallow+deep convection but each part can
!!     be activated/desactivated separately.
!!     For deep convection one can enable up to 3 additional ensemble members
!!     - this substantially improves the smoothness of the scheme and 
!!       allows for runs with different cloud radii (entrainment rates) and 
!!       reduces the arbitrariness inherent to convective trigger condition
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_DEEP
!!    CONVECT_SHALLOW
!!    INI_CONVPAR, INI_CONVPAR1
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/12/98 
!!      Modif       11/04/O2 allow for ensemble of deep updrafts/downdrafts
!!
!!    REFERENCE
!!    ---------
!!    Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc., Vol 127, pp 869-886: 
!!           A mass flux convection scheme for regional and global models.
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                  INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                  INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                  INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                  INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                  INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                ! KBDIA that is at least 1
INTEGER,                  INTENT(IN) :: KTDIA    ! vertical computations can be
                                                 ! limited to KLEV + 1 - KTDIA
                                                 ! default=1
REAL,                     INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                 ! calls of the deep convection
                                                 ! scheme
LOGICAL,                  INTENT(IN) :: ODEEP    ! switch for deep convection
LOGICAL,                  INTENT(IN) :: OSHAL    ! switch for shallow convection
LOGICAL,                  INTENT(IN) :: OREFRESH_ALL ! refresh or not all 
                                                 ! tendencies  at every call
LOGICAL,                  INTENT(IN) :: ODOWN    ! take or not convective
                                                 ! downdrafts into account
INTEGER,                  INTENT(IN) :: KICE     ! flag for ice ( 1 = yes, 
                                                 !                0 = no ice )
INTEGER,                  INTENT(IN) :: KENSM    ! number of additional deep convection calls
                                                 ! for ensemble (presently limited to 3)
                                                 ! KENSM=0 corresponds to base run with
                                                 ! 1 deep and 1 shallow call
LOGICAL,                  INTENT(IN) :: OSETTADJ ! logical to set convective
                                                 ! adjustment time by user
REAL,                     INTENT(IN) :: PTADJD   ! user defined deep adjustment time (s)
REAL,                     INTENT(IN) :: PTADJS   ! user defined shal. adjustment time (s)
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PT     ! grid scale T at time t  (K)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRV    ! grid scale water vapor  (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRC    ! grid scale r_c  (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRI    ! grid scale r_i  (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PU     ! grid scale horiz. wind u (m/s) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PV     ! grid scale horiz. wind v (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PW     ! grid scale vertical velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABS  ! grid scale pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ    ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY  ! grid area (m2)
!   
INTEGER, DIMENSION(KLON), INTENT(INOUT):: KCOUNT   ! convective counter(recompute
                                                   ! tendency or keep it
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperat. tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRTEN ! total surf precipitation tendency (m/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf precipitation tendency (m/s)
!
! Chemical Tracers:
LOGICAL,                    INTENT(IN)        :: OCHTRANS ! flag to compute convective
                                                          ! transport for chemical tracer
INTEGER,                    INTENT(IN)        :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN)   :: PCH1     ! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN  ! chemical convective tendency
                                                          ! (1/s)
!					 
! Diagnostic variables:
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux   (kg/s m2)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PPRLFLX! liquid precip flux  (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PPRSFLX! solid precip flux   (m/s)
REAL, DIMENSION(KLON),   INTENT(INOUT)    :: PCAPE  ! CAPE (J/kg)
INTEGER, DIMENSION(KLON),INTENT(INOUT)    :: KCLTOP ! cloud top level (number of model level)
INTEGER, DIMENSION(KLON),INTENT(INOUT)    :: KCLBAS ! cloud base level(number of model level)
                                                    ! they are given a value of
                                                    ! 0 if no convection
!						 
!*       0.2   Declarations of local variables :
!
INTEGER  :: JI, JK, JN  ! loop index
!
REAL, DIMENSION(KLON)               :: ZTIMEC, ZPRLTEN
!
! special for shallow convection
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,  &
                                       ZUMFS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCH1TENS
INTEGER, DIMENSION(:), ALLOCATABLE  :: ICLBASS, ICLTOPS
!
!*       0.3   Declarations of additional Ensemble fields:
!
INTEGER                             :: KENS     ! number of allowed additional deep convection calls
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTTENE   ! convective temperat. tendency (K/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRVTENE  ! convective r_v tendency (1/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRCTENE  ! convective r_c tendency (1/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRITENE  ! convective r_i tendency (1/s)
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZPRLTENE ! liquid surf precipitation tendency (m/s)
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZPRSTENE ! solid surf precipitation tendency (m/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZUMFE    ! updraft mass flux   (kg/s m2)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDMFE    ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZPRLFLXE ! liquid precip flux  (m/s)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZPRSFLXE ! solid precip flux   (m/s)
REAL, DIMENSION(:,:,:,:),ALLOCATABLE:: ZCH1TENE ! chemical convective tendency
INTEGER, DIMENSION(:,:),ALLOCATABLE :: ICLTOPE  ! cloud top level (number of model level)
INTEGER, DIMENSION(:,:),ALLOCATABLE :: ICLBASE  ! cloud base level(number of model level)
REAL, DIMENSION(:),     ALLOCATABLE :: ZEDUMMY  ! field not to be recomputed by ensemble
INTEGER, DIMENSION(:),  ALLOCATABLE :: IEDUMMY  ! field not to be recomputed by ensemble
REAL, DIMENSION(:),     ALLOCATABLE :: ZWEIGHT  ! weighting factor for ensemble members
REAL                                :: ZSUM     ! sum of weighting factors
!
!-------------------------------------------------------------------------------
!
!
!*       0.5  Allocate 2D (horizontal, vertical) arrays and additional ensemble arrays
!             ------------------------------------------------------------------------
!
    ALLOCATE( ZTTENS(KLON,KLEV) ); ALLOCATE( ZRVTENS(KLON,KLEV) ) 
    ALLOCATE( ZRCTENS(KLON,KLEV) ); ALLOCATE( ZRITENS(KLON,KLEV) ) 
    ALLOCATE( ZCH1TENS(KLON,KLEV,KCH1) ) 
    ALLOCATE( ZUMFS(KLON,KLEV) )
    ALLOCATE( ICLBASS(KLON) ); ALLOCATE( ICLTOPS(KLON) )
!
    KCLTOP(:)  = 1 ! set default value when no convection
    KCLBAS(:)  = 1 ! can be changed  depending on user
    ICLTOPS(:) = 1
    ICLBASS(:) = 1
!
KENS = MIN( KENSM, 3 )
IF ( KENS > 0 ) THEN
    ALLOCATE( ZTTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZRVTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZRCTENE(KLON,KLEV,KENS) )
    ALLOCATE( ZRITENE(KLON,KLEV,KENS) )
    ALLOCATE( ZUMFE(KLON,KLEV,KENS) )
    ALLOCATE( ZDMFE(KLON,KLEV,KENS) )
    ALLOCATE( ZCH1TENE(KLON,KLEV,KCH1,KENS) )
    ALLOCATE( ZPRLFLXE(KLON,KLEV,KENS) )
    ALLOCATE( ZPRSFLXE(KLON,KLEV,KENS) )
    ALLOCATE( ZPRLTENE(KLON,KENS) )
    ALLOCATE( ZPRSTENE(KLON,KENS) )
    ALLOCATE( ICLTOPE(KLON,KENS) )
    ALLOCATE( ICLBASE(KLON,KENS) )
    ALLOCATE( ZEDUMMY(KLON) )
    ALLOCATE( IEDUMMY(KLON) )
    ALLOCATE( ZWEIGHT(KENS) )
END IF
!
!*       4.a  Call deep convection routine
!             ----------------------------
!
IF ( ODEEP ) THEN
!
! 1. Base version
!
    CALL INI_CONVPAR
!
    IF ( OSETTADJ ) ZTIMEC(:) = PTADJD

!
    CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,           &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,  & 
                          PPABS, PZZ, PDXDY, ZTIMEC,                     &
                          PT, PRV, PRC, PRI, PU, PV, PW,                 &
                          KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                          ZPRLTEN, PPRSTEN,                              &
                          KCLTOP, KCLBAS, PPRLFLX, PPRSFLX,              & 
                          PUMF, PDMF, PCAPE,                             &
                          OCHTRANS, KCH1, PCH1, PCH1TEN                  )
!
!  2. Additional Ensemble members
!
    IF ( KENS > 0 ) THEN
!
    CALL INI_CONVPAR1
!
!* first member - changes in MODD_CONVPAR (cloud radius of 500 m or so)
!                                          specified in INI_CONVPAR1
!
    CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,                                   &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,                          & 
                          PPABS, PZZ, PDXDY, ZTIMEC,                                             &
                          PT, PRV, PRC, PRI, PU, PV, PW,                                         &
                          IEDUMMY, ZTTENE(:,:,1), ZRVTENE(:,:,1), ZRCTENE(:,:,1), ZRITENE(:,:,1),&
                          ZPRLTENE(:,1), ZPRSTENE(:,1),                                          &
                          ICLTOPE(:,1), ICLBASE(:,1), ZPRLFLXE(:,:,1), ZPRSFLXE(:,:,1),          & 
                          ZUMFE(:,:,1), ZDMFE(:,:,1), ZEDUMMY,                                   &
                          OCHTRANS, KCH1, PCH1, ZCH1TENE(:,:,:,1)                                )
    END IF
!
    IF ( KENS > 1 ) THEN
!
    CALL INI_CONVPAR
!
!* second member (positive vertical velocity perturb for Trigger)
!
    CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,                                   &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,                          & 
                          PPABS, PZZ, PDXDY, ZTIMEC,                                             &
                          PT, PRV, PRC, PRI, PU, PV, PW*1.5+1.e-4,                               &
                          IEDUMMY, ZTTENE(:,:,2), ZRVTENE(:,:,2), ZRCTENE(:,:,2), ZRITENE(:,:,2),&
                          ZPRLTENE(:,2), ZPRSTENE(:,2),                                          &
                          ICLTOPE(:,2), ICLBASE(:,2), ZPRLFLXE(:,:,2), ZPRSFLXE(:,:,2),          & 
                          ZUMFE(:,:,2), ZDMFE(:,:,2), ZEDUMMY,                                   &
                          OCHTRANS, KCH1, PCH1, ZCH1TENE(:,:,:,2)                                )
    END IF
!
    IF ( KENS > 2 ) THEN
!
!* third member (negative vertical velocity perturb for Trigger)
!
    CALL CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,                                   &
                          PDTCONV, KICE, OREFRESH_ALL, ODOWN, OSETTADJ,                          & 
                          PPABS, PZZ, PDXDY, ZTIMEC,                                             &
                          PT, PRV, PRC, PRI, PU, PV, PW*.5-1.e-4,                                &
                          IEDUMMY, ZTTENE(:,:,3), ZRVTENE(:,:,3), ZRCTENE(:,:,3), ZRITENE(:,:,3),&
                          ZPRLTENE(:,3), ZPRSTENE(:,3),                                          &
                          ICLTOPE(:,3), ICLBASE(:,3), ZPRLFLXE(:,:,3), ZPRSFLXE(:,:,3),          & 
                          ZUMFE(:,:,3), ZDMFE(:,:,3), ZEDUMMY,                                   &
                          OCHTRANS, KCH1, PCH1, ZCH1TENE(:,:,:,3)                                )
    END IF
!
ENDIF
IF ( .NOT. ODEEP ) THEN
     KCOUNT(:)  =0
     PTTEN(:,:) =0.
     PRVTEN(:,:)=0.
     PRCTEN(:,:)=0.
     PRITEN(:,:)=0.
     PUMF(:,:)  =0.
     PDMF(:,:)  =0.
   ! KCLTOP(:)  =1
   ! KCLBAS(:)  =1
     PCH1TEN(:,:,:)=0.
     ZPRLTEN(:)  =0.
     PPRSTEN(:)  =0.
     PPRLFLX(:,:)=0.
     PPRSFLX(:,:)=0.
     PCAPE(:)    =0.
END IF
!
!*       4.b  Call shallow convection routine
!             -------------------------------
!
IF ( OSHAL ) THEN
!
    IF ( .NOT. ODEEP ) CALL INI_CONVPAR 
    CALL INI_CONVPAR_SHAL
!
    CALL CONVECT_SHALLOW( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                             PDTCONV, KICE, OSETTADJ, PTADJS,            &
                             PPABS, PZZ,                                 &
                             PT, PRV, PRC, PRI, PW,                      &
                             ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,          &
                             ICLTOPS, ICLBASS, ZUMFS,                    &
                             OCHTRANS, KCH1, PCH1, ZCH1TENS              )
ENDIF
IF ( .NOT. OSHAL ) THEN
     ZTTENS(:,:) =0.
     ZRVTENS(:,:)=0.
     ZRCTENS(:,:)=0.
     ZRITENS(:,:)=0.
     ZUMFS(:,:)  =0
   ! ICLTOPS(:)  =1
   ! ICLBASS(:)  =1
     ZCH1TENS(:,:,:)=0.
END IF
!
!*       5.  Add  - if activated - ensemble average values for deep 
!            and then shallow convective tendencies
!            ---------------------------------------------------------
!
ZSUM = 1.
IF ( KENS > 0 ) THEN
    IF ( KENS == 1 ) ZWEIGHT(:) = .5
    IF ( KENS >  1 ) ZWEIGHT(:) = 1.
    DO JN = 1, KENS
       PTTEN(:,:)  = PTTEN(:,:)  + ZWEIGHT(JN) * ZTTENE(:,:,JN)
       PRVTEN(:,:) = PRVTEN(:,:) + ZWEIGHT(JN) * ZRVTENE(:,:,JN)
       PRCTEN(:,:) = PRCTEN(:,:) + ZWEIGHT(JN) * ZRCTENE(:,:,JN)
       PRITEN(:,:) = PRITEN(:,:) + ZWEIGHT(JN) * ZRITENE(:,:,JN)
       PPRLFLX(:,:)= PPRLFLX(:,:)+ ZWEIGHT(JN) * ZPRLFLXE(:,:,JN)
       PPRSFLX(:,:)= PPRSFLX(:,:)+ ZWEIGHT(JN) * ZPRSFLXE(:,:,JN)
       PUMF(:,:)   = PUMF(:,:)   + ZWEIGHT(JN) * ZUMFE(:,:,JN)
       PDMF(:,:)   = PDMF(:,:)   + ZWEIGHT(JN) * ZDMFE(:,:,JN)
       ZPRLTEN(:)  = ZPRLTEN(:)  + ZWEIGHT(JN) * ZPRLTENE(:,JN)
       PPRSTEN(:)  = PPRSTEN(:)  + ZWEIGHT(JN) * ZPRSTENE(:,JN)
       KCLTOP(:)   = MAX(KCLTOP(:), ICLTOPE(:,JN))
       KCLBAS(:)   = MAX(KCLBAS(:), ICLBASE(:,JN))
       IF ( OCHTRANS )  &
         & PCH1TEN(:,:,:) = PCH1TEN(:,:,:) + ZWEIGHT(JN) * ZCH1TENE(:,:,:,JN)
    END DO
!
    ZSUM = 1. / ( 1. + SUM( ZWEIGHT(:) ) )
END IF
!
       PTTEN(:,:)  = PTTEN(:,:)  * ZSUM + ZTTENS(:,:) 
       PRVTEN(:,:) = PRVTEN(:,:) * ZSUM + ZRVTENS(:,:)
       PRCTEN(:,:) = PRCTEN(:,:) * ZSUM + ZRCTENS(:,:)
       PRITEN(:,:) = PRITEN(:,:) * ZSUM + ZRITENS(:,:)
       PPRLFLX(:,:)= PPRLFLX(:,:)* ZSUM
       PPRSFLX(:,:)= PPRSFLX(:,:)* ZSUM
       PUMF(:,:)   = PUMF(:,:)   * ZSUM + ZUMFS(:,:)
       PDMF(:,:)   = PDMF(:,:)   * ZSUM
       PPRTEN(:)   = ( ZPRLTEN(:) + PPRSTEN(:) ) * ZSUM
       PPRSTEN(:)  = PPRSTEN(:)  * ZSUM
       KCLTOP(:)   = MAX(KCLTOP(:), ICLTOPS(:))
       KCLBAS(:)   = MAX(KCLBAS(:), ICLBASS(:))
    IF ( OCHTRANS ) THEN
          PCH1TEN(:,:,:) = PCH1TEN(:,:,:) * ZSUM + ZCH1TENS(:,:,:)
    END IF
!
!*       6.  Deallocate local arrays
!
       DEALLOCATE( ICLBASS ); DEALLOCATE( ICLTOPS )
       DEALLOCATE( ZUMFS )
       DEALLOCATE( ZCH1TENS ) 
       DEALLOCATE( ZRCTENS ); DEALLOCATE( ZRITENS ) 
       DEALLOCATE( ZTTENS ); DEALLOCATE( ZRVTENS ) 

IF ( KENS > 0 ) THEN
       DEALLOCATE( ZTTENE );   DEALLOCATE( ZRVTENE )
       DEALLOCATE( ZRCTENE );  DEALLOCATE( ZRITENE )
       DEALLOCATE( ZUMFE );    DEALLOCATE( ZDMFE )
       DEALLOCATE( ZCH1TENE )
       DEALLOCATE( ZPRLFLXE ); DEALLOCATE( ZPRSFLXE )
       DEALLOCATE( ZPRLTENE ); DEALLOCATE( ZPRSTENE )
       DEALLOCATE( ZEDUMMY );  DEALLOCATE( IEDUMMY )
       DEALLOCATE( ZWEIGHT )
END IF
!
!
END SUBROUTINE CONVECTION
!     ######spl
      MODULE MODD_CONVPAREXT
!     ######################
!
IMPLICIT NONE
!
INTEGER, SAVE :: JCVEXB ! start vertical computations at
                        ! 1 + JCVEXB = 1 + ( KBDIA - 1 )
INTEGER, SAVE :: JCVEXT ! limit vertical computations to
                        ! KLEV - JCVEXT = KLEV - ( KTDIA - 1 )

!$OMP threadprivate(JCVEXB,JCVEXT)
!
END MODULE MODD_CONVPAREXT
!     ######spl
      MODULE MODD_CST
!     ###############
!
IMPLICIT NONE
!
REAL, SAVE :: XP00   ! reference pressure
REAL, SAVE :: XPI    ! Pi
REAL, SAVE ::  XG    ! gravity constant
REAL, SAVE :: XMD    ! molecular weight of dry air
REAL, SAVE :: XMV    ! molecular weight of water vapor
REAL, SAVE :: XRD    ! gaz constant for dry air
REAL, SAVE :: XRV    ! gaz constant for water vapor
REAL, SAVE :: XCPD   ! specific heat of dry air
REAL, SAVE :: XCPV   ! specific heat of water vapor
REAL, SAVE :: XRHOLW ! density of liquid water
REAL, SAVE :: XCL    ! specific heat of liquid water
REAL, SAVE :: XCI    ! specific heat of ice
REAL, SAVE :: XTT    ! triple point temperature
REAL, SAVE :: XLVTT  ! latent heat of vaporisation at XTT
REAL, SAVE :: XLSTT  ! latent heat of sublimation at XTT 
REAL, SAVE :: XLMTT  ! latent heat of melting at XTT
REAL, SAVE :: XESTT  ! saturation pressure at XTT
REAL, SAVE :: XALPW  ! constants in saturation pressure over liquid water
REAL, SAVE :: XBETAW 
REAL, SAVE :: XGAMW 
REAL, SAVE :: XALPI  ! constants in saturation pressure over ice
REAL, SAVE :: XBETAI 
REAL, SAVE :: XGAMI 
!$OMP threadprivate(XP00,XPI,XG,XMD,XMV,XRD,XRV,XCPD,XCPV,XRHOLW,XCL, &
!$OMP XCI,XTT,XLVTT,XLSTT,XLMTT,XESTT,XALPW,XBETAW,XGAMW,XALPI, &
!$OMP XBETAI,XGAMI)

!
END MODULE MODD_CST
!     ######spl
      MODULE MODD_CONVPAR
!     ###################
!
!!****  *MODD_CONVPAR* - Declaration of convection constants 
!!
!!    PURPOSE
!!    -------
!      The purpose of this declarative module is to declare  the 
!      constants in the deep convection parameterization.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CONVPAR)
!!          
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96                      
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
!
REAL, SAVE :: XA25        ! 25 km x 25 km reference grid area
!
REAL, SAVE :: XCRAD       ! cloud radius 
REAL, SAVE :: XCDEPTH     ! minimum necessary cloud depth
REAL, SAVE :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)  
!
REAL, SAVE :: XZLCL       ! maximum allowed allowed height 
                          ! difference between departure level and surface
REAL, SAVE :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL, SAVE :: XWTRIG      ! constant in vertical velocity trigger
REAL, SAVE :: XDTHPBL     ! temperature perturbation in PBL for trigger
REAL, SAVE :: XDRVPBL     ! moisture perturbation in PBL for trigger
!
!
REAL, SAVE :: XNHGAM      ! accounts for non-hydrost. pressure 
                          ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
REAL, SAVE :: XTFRZ1      ! begin of freezing interval
REAL, SAVE :: XTFRZ2      ! end of freezing interval
!
REAL, SAVE :: XRHDBC      ! relative humidity below cloud in downdraft
!
REAL, SAVE :: XRCONV      ! constant in precipitation conversion 
REAL, SAVE :: XSTABT      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
REAL, SAVE :: XSTABC      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE
REAL, SAVE :: XUSRDPTH    ! pressure thickness used to compute updraft
                          ! moisture supply rate for downdraft
REAL, SAVE :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                          ! allowed below  melting level
REAL, SAVE :: XUVDP       ! constant for pressure perturb in momentum transport

!$OMP threadprivate(XA25,XCRAD,XCDEPTH,XENTR,XZLCL,XZPBL,XWTRIG, &
!$OMP XDTHPBL,XDRVPBL,XNHGAM,XTFRZ1,XTFRZ2,XRHDBC,XRCONV,XSTABT, &
!$OMP XSTABC,XUSRDPTH,XMELDPTH,XUVDP)
!
END MODULE MODD_CONVPAR 
!     ######spl
      SUBROUTINE INI_CONVPAR
!     ######################
!
!!****  *INI_CONVPAR * - routine to initialize the constants modules 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the constants
!     stored in  modules MODD_CONVPAR, MODD_CST, MODD_CONVPAREXT.
!      
!
!!**  METHOD
!!    ------
!!      The deep convection constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR   : contains deep convection constants
!!      Module MODD_CST       : contains physical constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CONVPAR, routine INI_CONVPAR)
!!      
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAR
USE MODD_CST
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------
!
!
XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
XCRAD    = 1500.     ! cloud radius 
XCDEPTH  = 3.E3      ! minimum necessary cloud depth
XENTR    = 0.03      ! entrainment constant (m/Pa) = 0.2 (m)  
!
XZLCL    = 3.5E3     ! maximum allowed allowed height 
                     ! difference between the surface and the LCL
XZPBL    = 60.E2     ! minimum mixed layer depth to sustain convection
XWTRIG   = 6.00      ! constant in vertical velocity trigger
XDTHPBL  = .3        ! Temp. perturbation in PBL for trigger
XDRVPBL  = 1.e-4     ! moisture  perturbation in PBL for trigger
!
!
XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XTFRZ1   = 273.16    ! begin of freezing interval
XTFRZ2   = 250.16    ! end of freezing interval
!
XRHDBC   = 0.9       ! relative humidity below cloud in downdraft

XRCONV   = 0.015     ! constant in precipitation conversion 
XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
XUSRDPTH = 165.E2    ! pressure thickness used to compute updraft
                     ! moisture supply rate for downdraft
XMELDPTH = 200.E2    ! layer (Pa) through which precipitation melt is
                     ! allowed below downdraft
XUVDP    = 0.7       ! constant for pressure perturb in momentum transport
!
!
!*       2.    Set the fundamental thermodynamical constants
!              these have the same values (not names) as in ARPEGE IFS 
!              -------------------------------------------------------
!
!
XP00   = 1.E5        ! reference pressure
XPI    = 3.141592654 ! Pi
 XG    = 9.80665     ! gravity constant
XMD    = 28.9644E-3  ! molecular weight of dry air
XMV    = 18.0153E-3  ! molecular weight of water vapor
XRD    = 287.05967   ! gaz constant for dry air
XRV    = 461.524993  ! gaz constant for water vapor
XCPD   = 1004.708845 ! specific heat of dry air
XCPV   = 1846.1      ! specific heat of water vapor
XRHOLW = 1000.       ! density of liquid water
XCL    = 4218.       ! specific heat of liquid water
XCI    = 2106.       ! specific heat of ice
XTT    = 273.16      ! triple point temperature
XLVTT  = 2.5008E6    ! latent heat of vaporisation at XTT
XLSTT  = 2.8345E6    ! latent heat of sublimation at XTT 
XLMTT  = 0.3337E6    ! latent heat of melting at XTT
XESTT  = 611.14      ! saturation pressure at XTT
XALPW  = 60.22416    ! constants in saturation pressure over liquid water
XBETAW = 6822.459384
XGAMW  = 5.13948
XALPI  = 32.62116    ! constants in saturation pressure over ice
XBETAI = 6295.421
XGAMI  = 0.56313
!
!
END SUBROUTINE INI_CONVPAR 
!     ######spl
      SUBROUTINE INI_CONVPAR1
!     #######################
!
!!****  *INI_CONVPAR * - routine to initialize the convective constants modules 
!!                       with modifications for ensemble run.
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the constants
!     stored in  modules MODD_CONVPAR, MODD_CST, MODD_CONVPAREXT.
!      
!
!!**  METHOD
!!    ------
!!      The deep convection constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CONVPAR, routine INI_CONVPAR)
!!      
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAR
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------
!
!
XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
XCRAD    =  500.     ! cloud radius 
XCDEPTH  = 3.E3      ! minimum necessary cloud depth
XENTR    = 0.03      ! entrainment constant (m/Pa) = 0.2 (m)  
!
XZLCL    = 3.5E3     ! maximum allowed allowed height 
                     ! difference between the surface and the LCL
XZPBL    = 60.E2     ! minimum mixed layer depth to sustain convection
XWTRIG   = 6.00      ! constant in vertical velocity trigger
XDTHPBL  = .3        ! Temp. perturbation in PBL for trigger
XDRVPBL  = 1.e-4     ! moisture  perturbation in PBL for trigger
!
!
XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XTFRZ1   = 273.16    ! begin of freezing interval
XTFRZ2   = 250.16    ! end of freezing interval
!
XRHDBC   = 0.9       ! relative humidity below cloud in downdraft

XRCONV   = 0.015     ! constant in precipitation conversion 
XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
XUSRDPTH = 165.E2    ! pressure thickness used to compute updraft
                     ! moisture supply rate for downdraft
XMELDPTH = 200.E2    ! layer (Pa) through which precipitation melt is
                     ! allowed below downdraft
XUVDP    = 0.7       ! constant for pressure perturb in momentum transport
!
!
END SUBROUTINE INI_CONVPAR1 
!     ######spl
    SUBROUTINE CONVECT_DEEP( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,           &
                                PDTCONV, KICE, OREFRESH, ODOWN, OSETTADJ,      &
                                PPABST, PZZ, PDXDY, PTIMEC,                    &
                                PTT, PRVT, PRCT, PRIT, PUT, PVT, PWT,          &
                                KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                                PPRLTEN, PPRSTEN,                              &
                                KCLTOP, KCLBAS, PPRLFLX, PPRSFLX,              &
                                PUMF, PDMF, PCAPE,                             &
                                OCH1CONV, KCH1, PCH1, PCH1TEN                  )
!   ############################################################################
!
!!**** Monitor routine to compute all convective tendencies by calls
!!     of several subroutines.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      tendencies. The routine first prepares all necessary grid-scale
!!      variables. The final convective tendencies are then computed by
!!      calls of different subroutines.
!!
!!
!!**  METHOD
!!    ------
!!      We start by selecting convective columns in the model domain through
!!      the call of routine TRIGGER_FUNCT. Then, we allocate memory for the
!!      convection updraft and downdraft variables and gather the grid scale
!!      variables in convective arrays. 
!!      The updraft and downdraft computations are done level by level starting
!!      at the  bottom and top of the domain, respectively.
!!      All computations are done on MNH thermodynamic levels. The depth
!!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_TRIGGER_FUNCT 
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_TSTEP_PREF
!!    CONVECT_DOWNDRAFT
!!    CONVECT_PRECIP_ADJUST
!!    CONVECT_CLOSURE
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                   ! gravity constant
!!          XPI                  ! number Pi
!!          XP00                 ! reference pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XRHOLW               ! density of liquid water
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!          XCL, XCI             ! specific heat for liquid water and ice
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
!!
!!      Module MODD_CONVPAR
!!          XA25                 ! reference grid area
!!          XCRAD                ! cloud radius
!!
!!         
!!    REFERENCE
!!    ---------
!!
!!      Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc. : 
!!           A mass flux convection scheme for regional and global models.
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Peter Bechtold 04/10/97 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAREXT
USE MODD_CONVPAR
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                    INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                    INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                  ! KBDIA that is at least 1
INTEGER,                    INTENT(IN) :: KTDIA    ! vertical computations can be
                                                   ! limited to KLEV + 1 - KTDIA
                                                   ! default=1
REAL,                       INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                   ! calls of the deep convection
                                                   ! scheme
INTEGER,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes, 
                                                   !                0 = no ice )
LOGICAL,                    INTENT(IN) :: OREFRESH ! refresh or not tendencies
						   ! at every call
LOGICAL,                    INTENT(IN) :: ODOWN    ! take or not convective
						   ! downdrafts into account
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUT      ! grid scale horiz. wind u "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PVT      ! grid scale horiz. wind v "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical 
						   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY    ! horizontal grid area (m-a2)
REAL, DIMENSION(KLON),      INTENT(IN) :: PTIMEC   ! value of convective adjustment
                                                   ! time if OSETTADJ=.TRUE.
!   
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCOUNT ! convective counter (recompute
                                                   ! tendency or keep it)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature
						   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRLTEN! liquid surf. precipitation
                                                   ! tendency (m/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf. precipitation
                                                   ! tendency (m/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRLFLX! liquid precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRSFLX! solid  precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF   ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PCAPE  ! maximum CAPE (J/kg)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: ITEST, ICONV, ICONV1    ! number of convective columns
INTEGER  :: IIB, IIE                ! horizontal loop bounds
INTEGER  :: IKB, IKE                ! vertical loop bounds
INTEGER  :: IKS                     ! vertical dimension
INTEGER  :: JI, JL                  ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKP, JKM            ! vertical loop index
INTEGER  :: IFTSTEPS                ! only used for chemical tracers
REAL     :: ZEPS, ZEPSA, ZEPSB      ! R_d / R_v, R_v / R_d, XCPV / XCPD - ZEPSA
REAL     :: ZCPORD, ZRDOCP          ! C_p/R_d,  R_d/C_p
!
LOGICAL, DIMENSION(KLON, KLEV)     :: GTRIG3 ! 3D logical mask for convection 
LOGICAL, DIMENSION(KLON)           :: GTRIG  ! 2D logical mask for trigger test
REAL, DIMENSION(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, 
							   ! theta_v, theta_es
REAL, DIMENSION(KLON)              :: ZTIME  ! convective time period
REAL, DIMENSION(KLON)              :: ZWORK2, ZWORK2B ! work array
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(:),ALLOCATABLE  :: IDPL    ! index for parcel departure level
INTEGER, DIMENSION(:),ALLOCATABLE  :: IPBL    ! index for source layer top
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILCL    ! index for lifting condensation level 
INTEGER, DIMENSION(:),ALLOCATABLE  :: IETL    ! index for zero buoyancy level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILFS    ! index for level of free sink
INTEGER, DIMENSION(:),ALLOCATABLE  :: IDBL    ! index for downdraft base level  
INTEGER, DIMENSION(:),ALLOCATABLE  :: IML     ! melting level  
!
INTEGER, DIMENSION(:), ALLOCATABLE :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(:), ALLOCATABLE :: ISLCL   ! index for lifting condensation level 
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTHLCL ! updraft theta at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSRVLCL ! updraft rv at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSWLCL  ! updraft w at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSZLCL  ! LCL height
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTHVELCL! envir. theta_v at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSDXDY  ! grid area (m^2)
!
! grid scale variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZZ      ! height of model layer (m) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPRES   ! grid scale pressure
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDPRES  ! pressure difference between 
                                              ! bottom and top of layer (Pa)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZU      ! grid scale horiz. u component on theta grid
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZV      ! grid scale horiz. v component on theta grid
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZW      ! grid scale vertical velocity on theta grid
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTT     ! temperature
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTH     ! grid scale theta     
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHV    ! grid scale theta_v     
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHL    ! grid scale enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHES, ZTHEST ! grid scale saturated theta_e
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRW     ! grid scale total water (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRV     ! grid scale water vapor (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRC     ! grid scale cloud water (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRI     ! grid scale cloud ice (kg/kg) 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZDXDY   ! grid area (m^2)
!
! updraft variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUMF    ! updraft mass flux (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUER    ! updraft entrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUDR    ! updraft detrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUPR    ! updraft precipitation in
                                              ! flux units (kg water / s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHL   ! updraft enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHV   ! updraft theta_v (K)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURW    ! updraft total water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURC    ! updraft cloud water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURI    ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURR    ! liquid precipit. (kg/kg)
                                              ! produced in  model layer
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURS    ! solid precipit. (kg/kg)
                                              ! produced in  model layer
REAL, DIMENSION(:),   ALLOCATABLE  :: ZUTPR   ! total updraft precipitation (kg/s)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZMFLCL  ! cloud base unit mass flux(kg/s) 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZCAPE   ! available potent. energy     
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTHLCL  ! updraft theta at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTLCL   ! updraft temp. at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZRVLCL  ! updraft rv at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZWLCL   ! updraft w at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZZLCL   ! LCL height
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTHVELCL! envir. theta_v at LCL
!
! downdraft variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDMF    ! downdraft mass flux (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDER    ! downdraft entrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDDR    ! downdraft detrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDTHL   ! downdraft enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDRW    ! downdraft total water (kg/kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZMIXF   ! mixed fraction at LFS        
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTPR    ! total surf precipitation (kg/s)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZSPR    ! solid surf precipitation (kg/s)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZDTEVR  ! donwndraft evapor. (kg/s)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZPREF   ! precipitation efficiency
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDTEVRF ! donwndraft evapor. (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPRLFLX ! liquid precip flux
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPRSFLX ! solid precip flux
!
! closure variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZLMASS  ! mass of model layer (kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTIMEA  ! advective time period
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTIMEC, ZTIMED! time during which convection is
                                              ! active at grid point (as ZTIME)
!
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRVC    ! conv. adj. grid scale r_w 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRCC    ! conv. adj. grid scale r_c 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRIC    ! conv. adj. grid scale r_i 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZWSUB   ! envir. compensating subsidence (Pa/s)
!
LOGICAL, DIMENSION(:),ALLOCATABLE  :: GTRIG1  ! logical mask for convection    
LOGICAL, DIMENSION(:),ALLOCATABLE  :: GWORK   ! logical work array
INTEGER, DIMENSION(:),ALLOCATABLE  :: IINDEX, IJINDEX, IJSINDEX, IJPINDEX!hor.index
REAL, DIMENSION(:),   ALLOCATABLE  :: ZCPH    ! specific heat C_ph 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZLV, ZLS! latent heat of vaporis., sublim.
REAL                               :: ZES     ! saturation vapor mixng ratio
!
! Chemical Tracers:
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1C   ! conv. adjust. chemical specy 1
REAL, DIMENSION(:,:),   ALLOCATABLE:: ZWORK3  ! work arrays 
LOGICAL, DIMENSION(:,:,:),ALLOCATABLE::GTRIG4 ! logical mask
!
!-------------------------------------------------------------------------------
!
!
!*       0.3    Compute loop bounds
!               -------------------
!
IIB = KIDIA
IIE = KFDIA
JCVEXB = MAX( 0, KBDIA - 1 )
IKB = 1 + JCVEXB 
IKS = KLEV
JCVEXT = MAX( 0, KTDIA - 1 )
IKE = IKS - JCVEXT 
!
!
!*       0.5    Update convective counter ( where KCOUNT > 0 
!               convection is still active ).
!               ---------------------------------------------
!
KCOUNT(IIB:IIE) = KCOUNT(IIB:IIE) - 1 
!
IF ( OREFRESH ) THEN
KCOUNT(:) = 1
KCOUNT(IIB:IIE) = 0 ! refresh or not at every call
END IF
!
GTRIG(:)  = KCOUNT(:) <= 0
ITEST = COUNT( GTRIG(:) )
IF ( ITEST == 0 ) RETURN  ! if convection is already active at every grid point
                          ! exit CONVECT_DEEP
!
!
!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------
!
GTRIG3(:,:) = SPREAD( GTRIG(:), DIM=2, NCOPIES=IKS )
WHERE ( GTRIG3(:,:) ) 
    PTTEN(:,:)  = 0.
    PRVTEN(:,:) = 0.
    PRCTEN(:,:) = 0.
    PRITEN(:,:) = 0.
    PPRLFLX(:,:)= 0.
    PPRSFLX(:,:)= 0.
!   PUTEN(:,:)  = 0.
!   PVTEN(:,:)  = 0.
    PUMF(:,:)   = 0.
    PDMF(:,:)   = 0.
END WHERE
WHERE ( GTRIG(:) ) 
   PPRLTEN(:) = 0.
   PPRSTEN(:) = 0.
!  KCLTOP(:)  = 0 ! already initialized in CONVECTION
!  KCLBAS(:)  = 0
   PCAPE(:)   = 0.
END WHERE
IF ( OCH1CONV ) THEN
   ALLOCATE( GTRIG4(KLON,KLEV,KCH1) )
   GTRIG4(:,:,:) = SPREAD( GTRIG3(:,:), DIM=3, NCOPIES=KCH1 )
   WHERE( GTRIG4(:,:,:) ) PCH1TEN(:,:,:) = 0.
   DEALLOCATE( GTRIG4 )
END IF
!
!
!*       1.     Initialize  local variables
!               ----------------------------
!
ZEPS   = XRD / XRV
ZEPSA  = XRV / XRD 
ZEPSB  = XCPV / XCPD - ZEPSA
ZCPORD = XCPD / XRD
ZRDOCP = XRD / XCPD
!
!
!*       1.1    Set up grid scale theta, theta_v, theta_es 
!               ------------------------------------------
!
ZTHT(:,:) = 300.
ZSTHV(:,:)= 300.
ZSTHES(:,:) = 400.
DO JK = IKB, IKE
DO JI = IIB, IIE
   IF ( PPABST(JI,JK) > 40.E2 ) THEN
      ZTHT(JI,JK)  = PTT(JI,JK) * ( XP00 / PPABST(JI,JK) ) ** ZRDOCP
      ZSTHV(JI,JK) = ZTHT(JI,JK) * ( 1. + ZEPSA * PRVT(JI,JK) ) /              &
                     ( 1. + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )
!
          ! use conservative Bolton (1980) formula for theta_e
          ! it is used to compute CAPE for undilute parcel ascent
          ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here
!
      ZES = EXP( XALPW - XBETAW / PTT(JI,JK) - XGAMW * LOG( PTT(JI,JK) ) )
      ZES = MIN( 1., ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
      ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **             &
                ( 1. - 0.28 * ZES ) * EXP( MIN(500.,                           &
                                           ( 3374.6525 / PTT(JI,JK) - 2.5403 ) &
                                          * ZES * ( 1. + 0.81 * ZES ) ) )
   END IF
END DO
END DO
!
!
!
!*       2.     Test for convective columns and determine properties at the LCL 
!               --------------------------------------------------------------
!
!*       2.1    Allocate arrays depending on number of model columns that need
!               to be tested for convection (i.e. where no convection is present
!               at the moment.
!               --------------------------------------------------------------
!
     ALLOCATE( ZPRES(ITEST,IKS) )
     ALLOCATE( ZZ(ITEST,IKS) )
     ALLOCATE( ZW(ITEST,IKS) )
     ALLOCATE( ZTH(ITEST,IKS) )
     ALLOCATE( ZTHV(ITEST,IKS) )
     ALLOCATE( ZTHEST(ITEST,IKS) )
     ALLOCATE( ZRV(ITEST,IKS) )
     ALLOCATE( ZSTHLCL(ITEST) )
     ALLOCATE( ZSTLCL(ITEST) )
     ALLOCATE( ZSRVLCL(ITEST) )
     ALLOCATE( ZSWLCL(ITEST) )
     ALLOCATE( ZSZLCL(ITEST) )
     ALLOCATE( ZSTHVELCL(ITEST) )
     ALLOCATE( ISDPL(ITEST) )
     ALLOCATE( ISPBL(ITEST) )
     ALLOCATE( ISLCL(ITEST) )
     ALLOCATE( ZSDXDY(ITEST) )
     ALLOCATE( GTRIG1(ITEST) )
     ALLOCATE( ZCAPE(ITEST) )
     ALLOCATE( IINDEX(KLON) )
     ALLOCATE( IJSINDEX(ITEST) )
     DO JI = 1, KLON
        IINDEX(JI) = JI
     END DO
     IJSINDEX(:) = PACK( IINDEX(:), MASK=GTRIG(:) )
!
  DO JK = IKB, IKE
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZPRES(JI,JK)  = PPABST(JL,JK)
     ZZ(JI,JK)     = PZZ(JL,JK)
     ZTH(JI,JK)    = ZTHT(JL,JK)
     ZTHV(JI,JK)   = ZSTHV(JL,JK)
     ZTHEST(JI,JK) = ZSTHES(JL,JK)
     ZRV(JI,JK)    = MAX( 0., PRVT(JL,JK) )
     ZW(JI,JK)     = PWT(JL,JK)
  END DO
  END DO
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZSDXDY(JI)    = PDXDY(JL)
  END DO
!
!*       2.2    Compute environm. enthalpy and total water = r_v + r_i + r_c 
!               and envir. saturation theta_e
!               ------------------------------------------------------------
!
!
!*       2.3    Test for convective columns and determine properties at the LCL 
!               --------------------------------------------------------------
!
     ISLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL 
     ISDPL(:) = IKB
     ISPBL(:) = IKB
!
!
     CALL CONVECT_TRIGGER_FUNCT( ITEST, KLEV,                              &
                                 ZPRES, ZTH, ZTHV, ZTHEST,                 &
                                 ZRV, ZW, ZZ, ZSDXDY,                      &
                                 ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                                 ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1,   &
                                 ZCAPE )
!
     DO JI = 1, ITEST
        JL = IJSINDEX(JI)
        PCAPE(JL) = ZCAPE(JI)
     END DO
!
     DEALLOCATE( ZPRES )
     DEALLOCATE( ZZ )
     DEALLOCATE( ZTH )
     DEALLOCATE( ZTHV )
     DEALLOCATE( ZTHEST )
     DEALLOCATE( ZRV )
     DEALLOCATE( ZW )
     DEALLOCATE( ZCAPE )
!
!
!*       3.     After the call of TRIGGER_FUNCT we allocate all the dynamic
!               arrays used in the convection scheme using the mask GTRIG, i.e.
!               we do calculus only in convective columns. This corresponds to
!               a GATHER operation.
!               --------------------------------------------------------------
!
     ICONV = COUNT( GTRIG1(:) )
     IF ( ICONV == 0 )  THEN 
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( GTRIG1 )
         DEALLOCATE( IINDEX )
         DEALLOCATE( IJSINDEX )
         RETURN   ! no convective column has been found, exit CONVECT_DEEP
     ENDIF
!
     ! vertical index variables
!
         ALLOCATE( IDPL(ICONV) )
         ALLOCATE( IPBL(ICONV) )
         ALLOCATE( ILCL(ICONV) )
         ALLOCATE( ICTL(ICONV) )
         ALLOCATE( IETL(ICONV) )
!
	 ! grid scale variables
!
         ALLOCATE( ZZ(ICONV,IKS) )
         ALLOCATE( ZPRES(ICONV,IKS) )
         ALLOCATE( ZDPRES(ICONV,IKS) )
         ALLOCATE( ZU(ICONV,IKS) )
         ALLOCATE( ZV(ICONV,IKS) )
         ALLOCATE( ZTT(ICONV, IKS) )
         ALLOCATE( ZTH(ICONV,IKS) )
         ALLOCATE( ZTHV(ICONV,IKS) )
         ALLOCATE( ZTHL(ICONV,IKS) )
         ALLOCATE( ZTHES(ICONV,IKS) )
         ALLOCATE( ZRV(ICONV,IKS) )
         ALLOCATE( ZRC(ICONV,IKS) )
         ALLOCATE( ZRI(ICONV,IKS) )
         ALLOCATE( ZRW(ICONV,IKS) )
         ALLOCATE( ZDXDY(ICONV) )
!
         ! updraft variables
!
         ALLOCATE( ZUMF(ICONV,IKS) )
         ALLOCATE( ZUER(ICONV,IKS) )
         ALLOCATE( ZUDR(ICONV,IKS) )
         ALLOCATE( ZUPR(ICONV,IKS) )
         ALLOCATE( ZUTHL(ICONV,IKS) )
         ALLOCATE( ZUTHV(ICONV,IKS) )
         ALLOCATE( ZURW(ICONV,IKS) )
         ALLOCATE( ZURC(ICONV,IKS) )
         ALLOCATE( ZURI(ICONV,IKS) )
         ALLOCATE( ZURR(ICONV,IKS) )
         ALLOCATE( ZURS(ICONV,IKS) )
         ALLOCATE( ZUTPR(ICONV) )
         ALLOCATE( ZTHLCL(ICONV) )
         ALLOCATE( ZTLCL(ICONV) )
         ALLOCATE( ZRVLCL(ICONV) )
         ALLOCATE( ZWLCL(ICONV) )
         ALLOCATE( ZMFLCL(ICONV) )
         ALLOCATE( ZZLCL(ICONV) )
         ALLOCATE( ZTHVELCL(ICONV) )
         ALLOCATE( ZCAPE(ICONV) )
!
         ! work variables
!
         ALLOCATE( IJINDEX(ICONV) )
         ALLOCATE( IJPINDEX(ICONV) )
         ALLOCATE( ZCPH(ICONV) )
         ALLOCATE( ZLV(ICONV) )
         ALLOCATE( ZLS(ICONV) )
!
!
!*           3.1    Gather grid scale and updraft base variables in
!                   arrays using mask GTRIG
!                   ---------------------------------------------------
!
         GTRIG(:)      = UNPACK( GTRIG1(:), MASK=GTRIG(:), FIELD=.FALSE. )  
         IJINDEX(:)    = PACK( IINDEX(:), MASK=GTRIG(:) )
!
    DO JK = IKB, IKE
    DO JI = 1, ICONV
         JL = IJINDEX(JI)
         ZZ(JI,JK)     = PZZ(JL,JK)
         ZPRES(JI,JK)  = PPABST(JL,JK)
         ZTT(JI,JK)    = PTT(JL,JK)
         ZTH(JI,JK)    = ZTHT(JL,JK)
         ZTHES(JI,JK)  = ZSTHES(JL,JK)
         ZRV(JI,JK)    = MAX( 0., PRVT(JL,JK) )
         ZRC(JI,JK)    = MAX( 0., PRCT(JL,JK) )
         ZRI(JI,JK)    = MAX( 0., PRIT(JL,JK) )
         ZTHV(JI,JK)   = ZSTHV(JL,JK)
         ZU(JI,JK)     = PUT(JL,JK)
         ZV(JI,JK)     = PVT(JL,JK)
    END DO
    END DO
    IF ( OSETTADJ ) THEN
         ALLOCATE( ZTIMED(ICONV) )
         DO JI = 1, ICONV
            JL = IJINDEX(JI)
            ZTIMED(JI) = PTIMEC(JL)
         END DO
    END IF
!
    DO JI = 1, ITEST
       IJSINDEX(JI) = JI	
    END DO
    IJPINDEX(:) = PACK( IJSINDEX(:), MASK=GTRIG1(:) )
    DO JI = 1, ICONV
         JL = IJPINDEX(JI)
         IDPL(JI)      = ISDPL(JL)
         IPBL(JI)      = ISPBL(JL)
         ILCL(JI)      = ISLCL(JL)
         ZTHLCL(JI)    = ZSTHLCL(JL)
         ZTLCL(JI)     = ZSTLCL(JL)
         ZRVLCL(JI)    = ZSRVLCL(JL)
         ZWLCL(JI)     = ZSWLCL(JL)
         ZZLCL(JI)     = ZSZLCL(JL)
         ZTHVELCL(JI)  = ZSTHVELCL(JL)
         ZDXDY(JI)     = ZSDXDY(JL)
    END DO
         ALLOCATE( GWORK(ICONV) )
         GWORK(:)      = PACK( GTRIG1(:),  MASK=GTRIG1(:) ) 
         DEALLOCATE( GTRIG1 )
         ALLOCATE( GTRIG1(ICONV) )
         GTRIG1(:)     = GWORK(:)
!                 
         DEALLOCATE( GWORK )
         DEALLOCATE( IJPINDEX )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )
!
!
!*           3.2    Compute pressure difference 
!                   ---------------------------------------------------
!
        ZDPRES(:,IKB) = 0.
        DO JK = IKB + 1, IKE
            ZDPRES(:,JK)  = ZPRES(:,JK-1) - ZPRES(:,JK)
        END DO
!
!*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c 
!                  ----------------------------------------------------------
!
        DO JK = IKB, IKE, 1
            ZRW(:,JK)  = ZRV(:,JK) + ZRC(:,JK) + ZRI(:,JK)
            ZCPH(:)    = XCPD + XCPV * ZRW(:,JK)
            ZLV(:)     = XLVTT + ( XCPV - XCL ) * ( ZTT(:,JK) - XTT ) ! compute L_v
            ZLS(:)     = XLSTT + ( XCPV - XCI ) * ( ZTT(:,JK) - XTT ) ! compute L_i
            ZTHL(:,JK) = ZCPH(:) * ZTT(:,JK) + ( 1. + ZRW(:,JK) ) * XG * ZZ(:,JK) &
                         - ZLV(:) * ZRC(:,JK) - ZLS(:) * ZRI(:,JK)
        END DO
!
!
!*           4.     Compute updraft properties 
!                   ----------------------------
!
!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s ) 
!                   -------------------------------------------------------------
!
         DO JI = 1, ICONV
               JK = ILCL(JI) - 1
               ZMFLCL(JI) = ZPRES(JI,JK) / ( XRD * ZTT(JI,JK) *                &
                         ( 1. + ZEPS * ZRVLCL(JI) ) ) * XPI * XCRAD * XCRAD 
         END DO
!
         DEALLOCATE( ZCPH )
         DEALLOCATE( ZLV )
         DEALLOCATE( ZLS )
!
!
     CALL CONVECT_UPDRAFT( ICONV, KLEV,                                     &
                           KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
                           ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   & 
                           ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
                           ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                           ZURC, ZURI, ZURR, ZURS, ZUPR,                    &
                           ZUTPR, ZCAPE, ICTL, IETL                         )
!
!
!
!*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud 
!                   thickness is smaller than 3 km
!                   -----------------------------------------------------------
!
!
     ICONV1 = COUNT(GTRIG1) 
!
     IF ( ICONV1 > 0 )  THEN
!
!*       4.3    Allocate memory for downdraft variables
!               ---------------------------------------
!
! downdraft variables
!
        ALLOCATE( ILFS(ICONV) )
        ALLOCATE( IDBL(ICONV) )
        ALLOCATE( IML(ICONV) )
        ALLOCATE( ZDMF(ICONV,IKS) )
        ALLOCATE( ZDER(ICONV,IKS) )
        ALLOCATE( ZDDR(ICONV,IKS) )
        ALLOCATE( ZDTHL(ICONV,IKS) )
        ALLOCATE( ZDRW(ICONV,IKS) )
        ALLOCATE( ZLMASS(ICONV,IKS) )
        DO JK = IKB, IKE
           ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / XG  ! mass of model layer
        END DO
	ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
        ALLOCATE( ZMIXF(ICONV) )
        ALLOCATE( ZTPR(ICONV) )
        ALLOCATE( ZSPR(ICONV) )
        ALLOCATE( ZDTEVR(ICONV) )
        ALLOCATE( ZPREF(ICONV) )
        ALLOCATE( ZDTEVRF(ICONV,IKS) )
        ALLOCATE( ZPRLFLX(ICONV,IKS) )
        ALLOCATE( ZPRSFLX(ICONV,IKS) )
!
! closure variables
!
        ALLOCATE( ZTIMEA(ICONV) )
        ALLOCATE( ZTIMEC(ICONV) )
        ALLOCATE( ZTHC(ICONV,IKS) )
        ALLOCATE( ZRVC(ICONV,IKS) )
        ALLOCATE( ZRCC(ICONV,IKS) )
        ALLOCATE( ZRIC(ICONV,IKS) )
        ALLOCATE( ZWSUB(ICONV,IKS) )
!
!
!*           5.     Compute downdraft properties 
!                   ----------------------------
!
!*           5.1    Compute advective time period and precipitation 
!                   efficiency as a function of mean ambient wind (shear) 
!                   --------------------------------------------------------
!
        CALL CONVECT_TSTEP_PREF( ICONV, KLEV,                          &
                                 ZU, ZV, ZPRES, ZZ, ZDXDY, ILCL, ICTL, &
                                 ZTIMEA, ZPREF )
!
          ! exclude convective downdrafts if desired
        IF ( .NOT. ODOWN ) ZPREF(:) = 1.
!
          ! Compute the period during which convection is active
        ZTIMEC(:) = MAX( 1800., MIN( 3600., ZTIMEA(:) ) )
        ZTIMEC(:) = REAL( INT( ZTIMEC(:) / PDTCONV ) ) * PDTCONV
        ZTIMEC(:) = MAX( PDTCONV, ZTIMEC(:) ) ! necessary if PDTCONV > 1800
        IF ( OSETTADJ ) THEN
             ZTIMEC(:) = MAX( PDTCONV, ZTIMED(:) )
        END IF
!
!
!*           5.2    Compute melting level
!                   ----------------------
!
        IML(:) = IKB
        DO JK = IKE, IKB, -1
          WHERE( ZTT(:,JK) <= XTT )  IML(:) = JK
        END DO
!
        CALL CONVECT_DOWNDRAFT( ICONV, KLEV,                               &
                                KICE, ZPRES, ZDPRES, ZZ, ZTH, ZTHES,       & 
                                ZRW, ZRC, ZRI,                             &
                                ZPREF, ILCL, ICTL, IETL,                   &
                                ZUTHL, ZURW, ZURC, ZURI,                   &
                                ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,             &
                                ZMIXF, ZDTEVR, ILFS, IDBL, IML,            &
                                ZDTEVRF                                    )
!
!
!*           6.     Adjust up and downdraft mass flux to be consistent
!                   with precipitation efficiency relation.
!                   --------------------------------------------------- 
!
       CALL CONVECT_PRECIP_ADJUST( ICONV, KLEV,                              &
                                   ZPRES,ZUMF, ZUER, ZUDR, ZUPR, ZUTPR, ZURW,&
                                   ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,            &
                                   ZPREF, ZTPR, ZMIXF, ZDTEVR,               &
                                   ILFS, IDBL, ILCL, ICTL, IETL,             &
                                   ZDTEVRF                                   )
!
!
!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------
!
       CALL CONVECT_CLOSURE( ICONV, KLEV,                                &
                             ZPRES, ZDPRES, ZZ, ZDXDY, ZLMASS,           &
                             ZTHL, ZTH, ZRW, ZRC, ZRI, GTRIG1,           &
                             ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,              &
                             ILCL, IDPL, IPBL, ILFS, ICTL, IML,          &
                             ZUMF, ZUER, ZUDR, ZUTHL, ZURW,              &
                             ZURC, ZURI, ZUPR,                           &
                             ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,              &
                             ZTPR, ZSPR, ZDTEVR,                         &
                             ZCAPE, ZTIMEC,                              &
                             IFTSTEPS,                                   &
                             ZDTEVRF, ZPRLFLX, ZPRSFLX )
!

 
!
!*           8.     Determine the final grid-scale (environmental) convective 
!                   tendencies and set convective counter
!                   --------------------------------------------------------
!
!
!*           8.1    Grid scale tendencies
!                   ---------------------
!
          ! in order to save memory, the tendencies are temporarily stored
          ! in the tables for the adjusted grid-scale values
!
      DO JK = IKB, IKE
         ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
           * ( ZPRES(:,JK) / XP00 ) ** ZRDOCP ! change theta in temperature
         ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) ) &
					         / ZTIMEC(:) 

         ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
         ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:) 
!
         ZPRLFLX(:,JK) = ZPRLFLX(:,JK) / ( XRHOLW * ZDXDY(:) )
         ZPRSFLX(:,JK) = ZPRSFLX(:,JK) / ( XRHOLW * ZDXDY(:) )
!
      END DO
!
      ZPRLFLX(:,IKB) = ZPRLFLX(:,IKB+1)
      ZPRSFLX(:,IKB) = ZPRSFLX(:,IKB+1)
!
!
!*           8.2    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = 0.
       ZWORK2B(:) = 0.
       DO JK = JKM, IKB+1, -1
	 JKP = JK + 1
         DO JI = 1, ICONV
           ZWORK2(JI) = ZWORK2(JI) + ( ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK) ) *   & ! moisture
                                             (ZPRES(JI,JK) - ZPRES(JI,JKP)) / XG
      !    ZWORK2B(JI) = ZWORK2B(JI) + (                                             & ! energy
      !                                ( XCPD + XCPV * ZRW(JI,JK) )* ZTHC(JI,JK)   - &
      !          ( XLVTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRCC(JI,JK)   - & 
      !          ( XLSTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRIC(JI,JK) ) * & 
      !                                      (ZPRES(JI,JK) - ZPRES(JI,JKP)) / XG
         END DO
       END DO
!
          ! Budget error (compare integral to surface precip.)
!
       DO JI = 1, ICONV
         IF ( ZTPR(JI) > 0.) THEN
           JKP = ICTL(JI) + 1
           ZWORK2(JI) = ( ZTPR(JI) / ZDXDY(JI) + ZWORK2(JI) ) * XG /                 &
                                        ( ZPRES(JI,IKB+1) - ZPRES(JI,JKP) )
      !    ZWORK2B(JI) = ( ZTPR(JI) / ZDXDY(JI) *                                    &
      !       ( XLVTT + ( XCPV - XCL ) * ( ZTT(JI,IKB) - XTT ) ) - ZWORK2B(JI) )     &
      !                                * XG / ( ZPRES(JI,IKB+1) - ZPRES(JI,JKP) )
         END IF
       END DO
!
          ! Apply uniform correction
!
       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > 0. .AND. JK <= ICTL(JI) ) THEN
                ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)                                ! moisture
       !        ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2B(JI) / ( XCPD + XCPV * ZRW(JI,JK) )! energy
           END IF
         END DO
       END DO
!
!

              ! execute a "scatter"= pack command to store the tendencies in
              ! the final 2D tables
!
      DO JK = IKB, IKE
      DO JI = 1, ICONV
         JL = IJINDEX(JI)
         PTTEN(JL,JK)   = ZTHC(JI,JK)
         PRVTEN(JL,JK)  = ZRVC(JI,JK)
         PRCTEN(JL,JK)  = ZRCC(JI,JK)
         PRITEN(JL,JK)  = ZRIC(JI,JK)
!
         PPRLFLX(JL,JK) = ZPRLFLX(JI,JK)
         PPRSFLX(JL,JK) = ZPRSFLX(JI,JK)
      END DO
      END DO
!
!
!
!
!*           8.3    Convective rainfall tendency
!                   ----------------------------
!
	         ! liquid and solid surface rainfall tendency in m/s
       ZTPR(:)   = ZTPR(:) / ( XRHOLW * ZDXDY(:) ) ! total surf precip
       ZSPR(:)   = ZSPR(:) / ( XRHOLW * ZDXDY(:) ) ! solid surf precip
       ZTPR(:)   = ZTPR(:) - ZSPR(:) ! compute liquid part
!
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        PPRLTEN(JL) = ZTPR(JI)
        PPRSTEN(JL) = ZSPR(JI)
     END DO
!
!
!                   Cloud base and top levels
!                   -------------------------
!
     ILCL(:) = MIN( ILCL(:), ICTL(:) )
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        KCLTOP(JL) = ICTL(JI)
        KCLBAS(JL) = ILCL(JI)
     END DO
!
!
!*           8.4    Set convective counter
!                   ----------------------
!
	 ! compute convective counter for just activated convective
         ! grid points
         ! If the advective time period is less than specified
         ! minimum for convective period, allow feedback to occur only
         ! during advective time
!
     ZTIME(:) = 1.
     ZWORK2(:) = 0.
     DO JI = 1, ICONV
       JL = IJINDEX(JI)
       ZTIME(JL)  =  ZTIMEC(JI)
       ZWORK2(JL) =  ZTIMEA(JI)
       ZWORK2(JL) =  MIN( ZWORK2(JL), ZTIME(JL) )
       ZWORK2(JL) =  MAX( ZWORK2(JL), PDTCONV )
       IF ( GTRIG(JL) )  KCOUNT(JL) = INT( ZWORK2(JL) / PDTCONV )
       IF ( GTRIG(JL) .AND. PPRLTEN(JL)<1.E-14 ) KCOUNT(JL) = 0
     END DO
!
!
!
!*           8.7    Compute convective tendencies for Tracers
!                   ------------------------------------------
!
  IF ( OCH1CONV ) THEN
!
       ALLOCATE( ZCH1(ICONV,IKS,KCH1) )
       ALLOCATE( ZCH1C(ICONV,IKS,KCH1) )
       ALLOCATE( ZWORK3(ICONV,KCH1) )
!
       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          ZCH1(JI,JK,:) = PCH1(JL,JK,:)
       END DO
       END DO
!
      CALL CONVECT_CHEM_TRANSPORT( ICONV, KLEV, KCH1, ZCH1, ZCH1C,      &
                                   IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,  &
                                   ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,  &
                                   ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB, &
				   IFTSTEPS )
!
       DO JK = IKB, IKE
       DO JN = 1, KCH1
          ZCH1C(:,JK,JN) = ( ZCH1C(:,JK,JN)- ZCH1(:,JK,JN) ) / ZTIMEC(:)
       END DO
       END DO
!
!*           8.8    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
       JKM = MAXVAL( ICTL(:) )
       ZWORK3(:,:) = 0.
       DO JK = JKM, IKB+1, -1
	 JKP = JK + 1
         DO JI = 1, ICONV
           ZWORK3(JI,:) = ZWORK3(JI,:) + ZCH1C(JI,JK,:) *                    &
                              (ZPRES(JI,JK) - ZPRES(JI,JKP)) / XG
         END DO
       END DO
!
          ! Mass error (integral must be zero)
!
       DO JI = 1, ICONV
         IF ( ZTPR(JI) > 0.) THEN
           JKP = ICTL(JI) + 1
           ZWORK3(JI,:) = ZWORK3(JI,:) *                                     &
                                    XG / ( ZPRES(JI,IKB+1) - ZPRES(JI,JKP) )
         END IF
       END DO
!
          ! Apply uniform correction but assure positive mass at each level
!
       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ZTPR(JI) > 0. .AND. JK <= ICTL(JI) ) THEN
                ZCH1C(JI,JK,:) = ZCH1C(JI,JK,:) - ZWORK3(JI,:)
                ZCH1C(JI,JK,:) = MAX( ZCH1C(JI,JK,:), -ZCH1(JI,JK,:)/ZTIMEC(JI) )
           END IF
         END DO
       END DO
!
       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PCH1TEN(JL,JK,:) = ZCH1C(JI,JK,:)
       END DO
       END DO
  END IF
!
!
!*           9.     Write up- and downdraft mass fluxes
!                   ------------------------------------
!
    DO JK = IKB, IKE
       ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
       ZDMF(:,JK)  = ZDMF(:,JK) / ZDXDY(:)
    END DO
    ZWORK2(:) = 1.
    WHERE ( PPRLTEN(:)<1.E-14 ) ZWORK2(:) = 0.
    DO JK = IKB, IKE
    DO JI = 1, ICONV
       JL = IJINDEX(JI)
       PUMF(JL,JK) = ZUMF(JI,JK) * ZWORK2(JL)
       PDMF(JL,JK) = ZDMF(JI,JK) * ZWORK2(JL)
    END DO
    END DO
!
!
!*           10.    Deallocate all local arrays
!                   ---------------------------
!
! downdraft variables
!
      DEALLOCATE( ZDMF )
      DEALLOCATE( ZDER )
      DEALLOCATE( ZDDR )
      DEALLOCATE( ZDTHL )
      DEALLOCATE( ZDRW )
      DEALLOCATE( ZLMASS )
      DEALLOCATE( ZMIXF )
      DEALLOCATE( ZTPR )
      DEALLOCATE( ZSPR )
      DEALLOCATE( ZDTEVR )
      DEALLOCATE( ZPREF )
      DEALLOCATE( IML )
      DEALLOCATE( ILFS )
      DEALLOCATE( IDBL )
      DEALLOCATE( ZDTEVRF )
      DEALLOCATE( ZPRLFLX )
      DEALLOCATE( ZPRSFLX )
!
!   closure variables
!
      DEALLOCATE( ZTIMEA )
      DEALLOCATE( ZTIMEC )
      DEALLOCATE( ZTHC )
      DEALLOCATE( ZRVC )
      DEALLOCATE( ZRCC )
      DEALLOCATE( ZRIC )
      DEALLOCATE( ZWSUB )
!
       IF ( OCH1CONV ) THEN
           DEALLOCATE( ZCH1 )
           DEALLOCATE( ZCH1C )
           DEALLOCATE( ZWORK3 )
       END IF
!
    ENDIF
!
!    vertical index
!
    DEALLOCATE( IDPL )
    DEALLOCATE( IPBL )
    DEALLOCATE( ILCL )
    DEALLOCATE( ICTL )
    DEALLOCATE( IETL )
!
! grid scale variables
!
    DEALLOCATE( ZZ )
    DEALLOCATE( ZPRES )
    DEALLOCATE( ZDPRES )
    DEALLOCATE( ZU )
    DEALLOCATE( ZV )
    DEALLOCATE( ZTT )
    DEALLOCATE( ZTH )
    DEALLOCATE( ZTHV )
    DEALLOCATE( ZTHL )
    DEALLOCATE( ZTHES )
    DEALLOCATE( ZRW )
    DEALLOCATE( ZRV )
    DEALLOCATE( ZRC )
    DEALLOCATE( ZRI )
    DEALLOCATE( ZDXDY )
!
! updraft variables
!
    DEALLOCATE( ZUMF )
    DEALLOCATE( ZUER )
    DEALLOCATE( ZUDR )
    DEALLOCATE( ZUTHL )
    DEALLOCATE( ZUTHV )
    DEALLOCATE( ZURW )
    DEALLOCATE( ZURC )
    DEALLOCATE( ZURI )
    DEALLOCATE( ZURR )
    DEALLOCATE( ZURS )
    DEALLOCATE( ZUPR )
    DEALLOCATE( ZUTPR )
    DEALLOCATE( ZTHLCL )
    DEALLOCATE( ZTLCL )
    DEALLOCATE( ZRVLCL )
    DEALLOCATE( ZWLCL )
    DEALLOCATE( ZZLCL )
    DEALLOCATE( ZTHVELCL )
    DEALLOCATE( ZMFLCL )
    DEALLOCATE( ZCAPE )
    IF ( OSETTADJ ) DEALLOCATE( ZTIMED )
!
! work arrays
!
    DEALLOCATE( IINDEX )
    DEALLOCATE( IJINDEX )
    DEALLOCATE( IJSINDEX )
    DEALLOCATE( GTRIG1 )
!
!
END SUBROUTINE CONVECT_DEEP
!     ######spl
      SUBROUTINE CONVECT_TRIGGER_FUNCT( KLON, KLEV,                           &
                                        PPRES, PTH, PTHV, PTHES,              &
                                        PRV, PW, PZ, PDXDY,                   &
                                        PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL,  &
                                        PTHVELCL, KLCL, KDPL, KPBL, OTRIG,    &
                                        PCAPE )
!     ######################################################################
!
!!**** Determine convective columns as well as the cloudy values of theta,
!!     and qv at the lifting condensation level (LCL) 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine convective columns
!!   
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      What we look for is the undermost unstable level at each grid point.
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! Reference pressure
!!          XRD, XRV           ! Gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XTT                ! triple point temperature
!!          XBETAW, XGAMW      ! constants for vapor saturation pressure
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! maximum height difference between
!!                             ! the surface and the DPL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  20/03/97  Select first departure level
!!                            that produces a cloud thicker than XCDEPTH
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                   :: KLON      ! horizontal loop index
INTEGER, INTENT(IN)                   :: KLEV      ! vertical loop index
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY     ! grid area
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH, PTHV ! theta, theta_v
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHES     ! envir. satur. theta_e
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRV       ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PPRES     ! pressure
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PZ        ! height of grid point (m)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PW        ! vertical velocity
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL    ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL     ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PWLCL     ! parcel velocity at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL     ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(OUT):: OTRIG     ! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KLCL    ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KDPL    ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KPBL    ! contains index of source layer top
REAL, DIMENSION(KLON),     INTENT(OUT):: PCAPE     ! CAPE (J/kg) for diagnostics
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JKK, JK, JKP, JKM, JKDL, JL, JKT, JT! vertical loop index
INTEGER :: JI                                  ! horizontal loop index 
INTEGER :: IIE, IKB, IKE                       ! horizontal + vertical loop bounds
REAL    :: ZEPS, ZEPSA                         ! R_d / R_v, R_v / R_d 
REAL    :: ZCPORD, ZRDOCP                      ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON) :: ZTHLCL, ZTLCL, ZRVLCL, & ! locals for PTHLCL,PTLCL
                               ZWLCL,  ZZLCL, ZTHVELCL  ! PRVLCL, ....
INTEGER, DIMENSION(KLON) :: IDPL, IPBL, ILCL      ! locals for KDPL, ...
REAL, DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL, DIMENSION(KLON) :: ZZDPL    ! height of DPL 
REAL, DIMENSION(KLON) :: ZTHVLCL  ! theta_v at LCL = mixed layer value
REAL, DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL, DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure 
REAL, DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL, DIMENSION(KLON) :: ZCAPE    ! convective available energy (m^2/s^2/g)
REAL, DIMENSION(KLON) :: ZTHEUL   ! updraft equiv. pot. temperature (K)
REAL, DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL, DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(KLON) :: ZTOP     ! estimated cloud top (m)
REAL, DIMENSION(KLON,KLEV):: ZCAP ! CAPE at every level for diagnostics
!INTEGER, DIMENSION(KLON) :: ITOP  ! work array to store highest test layer
REAL, DIMENSION(KLON) :: ZWORK1, ZWORK2, ZWORK3    ! work arrays
LOGICAL, DIMENSION(KLON) :: GTRIG, GTRIG2          ! local arrays for OTRIG
LOGICAL, DIMENSION(KLON) :: GWORK1                 ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Initialize local variables
!               --------------------------
!
ZEPS       = XRD / XRV
ZEPSA      = XRV / XRD 
ZCPORD     = XCPD / XRD
ZRDOCP     = XRD / XCPD
OTRIG(:)   = .FALSE.
IDPL(:)    = KDPL(:)
IPBL(:)    = KPBL(:)
ILCL(:)    = KLCL(:)
!ITOP(:)    = IKB
PWLCL(:)   = 0.
ZWLCL(:)   = 0.
PTHLCL(:)  = 1.
PTHVELCL(:)= 1.
PTLCL(:)   = 1.
PRVLCL(:)  = 0.
PWLCL(:)   = 0.
PZLCL(:)   = PZ(:,IKB)
ZZDPL(:)   = PZ(:,IKB)
GTRIG2(:)  = .TRUE.
ZCAP(:,:)  = 0.
!
!
!
!       1.     Determine highest necessary loop test layer
!              -------------------------------------------
!
JT = IKE - 2
DO JK = IKB + 1, IKE - 2
 ! DO JI = 1, IIE
 !    IF ( PZ(JI,JK) - PZ(JI,IKB) <= XZLCL ) ITOP(JI) = JK
 ! END DO
   IF ( PZ(1,JK) - PZ(1,IKB) < 12.E3 ) JT = JK 
END DO
!
!
!*       2.     Enter loop for convection test
!               ------------------------------
!
JKP = MINVAL( IDPL(:) ) + 1
!JKT = MAXVAL( ITOP(:) )
JKT = JT
DO JKK = JKP, JKT
!
     GWORK1(:) = ZZDPL(:) - PZ(:,IKB) < XZLCL
          ! we exit the trigger test when the center of the mixed layer is more
          ! than 3500 m  above soil level.
     WHERE ( GWORK1(:) )
        ZDPTHMIX(:) = 0.
        ZPRESMIX(:) = 0.
        ZTHLCL(:)   = 0.
        ZRVLCL(:)   = 0.
        ZZDPL(:)    = PZ(:,JKK)
        IDPL(:)     = JKK
     END WHERE
!
!
!*       3.     Construct a mixed layer of at least 60 hPa (XZPBL)
!               ------------------------------------------
!
     DO JK = JKK, IKE - 1
       JKM = JK + 1
       DO JI = 1, IIE     
         IF ( GWORK1(JI) .AND. ZDPTHMIX(JI) < XZPBL ) THEN
            IPBL(JI)     = JK
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            ZTHLCL(JI)   = ZTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            ZRVLCL(JI)   = ZRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
         END IF
       END DO
        IF ( MINVAL ( ZDPTHMIX(:) ) >= XZPBL ) EXIT
     END DO
!
!
     WHERE ( GWORK1(:) )
!
        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
      ! ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) 
      ! ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:) 
        ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) + XDTHPBL
        ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:) + XDRVPBL
        ZTHVLCL(:)  = ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                 &
				/ ( 1. + ZRVLCL(:) )
!
!*       4.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL. 
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               ----------------------------------------------------
!
! 
        ZTMIX(:)  = ZTHLCL(:) * ( ZPRESMIX(:) / XP00 ) ** ZRDOCP 
        ZEVMIX(:) = ZRVLCL(:) * ZPRESMIX(:) / ( ZRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8, ZEVMIX(:) )
        ZWORK1(:) = LOG( ZEVMIX(:) / 613.3 )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8 - 32.19 * ZWORK1(:) ) / ( 17.502 - ZWORK1(:) ) 
              ! adiabatic saturation temperature
        ZTLCL(:)  = ZWORK1(:) - ( .212 + 1.571E-3 * ( ZWORK1(:) - XTT )      &
                   - 4.36E-4 * ( ZTMIX(:) - XTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        ZTLCL(:)  = MIN( ZTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = XP00 * ( ZTLCL(:) / ZTHLCL(:) ) ** ZCPORD
!
     END WHERE
!
!
!*       4.2    Correct ZTLCL in order to be completely consistent
!               with MNH saturation formula
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, ZTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTLCL(:) * ( XBETAW / ZTLCL(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                        ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        ZTLCL(:)  = ZTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
!
     END WHERE
!
!
!*       4.3    If ZRVLCL = PRVMIX is oversaturated set humidity 
!               and temperature to saturation values. 
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) .AND. ZRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( XBETAW / ZTMIX(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                       ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) 
        ZTLCL(:)  = ZTMIX(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
        ZRVLCL(:) = ZRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        ZTHLCL(:) = ZTLCL(:) * ( XP00 / ZPLCL(:) ) ** ZRDOCP
        ZTHVLCL(:)= ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                   &
                              / ( 1. + ZRVLCL(:) )
     END WHERE
!
!
!*        5.1   Determine  vertical loop index at the LCL and DPL
!               --------------------------------------------------
!
    DO JK = JKK, IKE - 1
       DO JI = 1, IIE
         IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. GWORK1(JI) ) ILCL(JI) = JK + 1
       END DO
    END DO
!
!
!*        5.2   Estimate height and environm. theta_v at LCL
!               --------------------------------------------------
!
    DO JI = 1, IIE
        JK   = ILCL(JI)
        JKM  = JK - 1
        ZDP(JI)    = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                     LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI) 
           ! we compute the precise value of the LCL
           ! The precise height is between the levels ILCL and ILCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    END DO
    WHERE( GWORK1(:) )
        ZTHVELCL(:) = ZWORK1(:) 
        ZZLCL(:)    = ZWORK2(:)
    END WHERE
!        
!
!*       6.     Check to see if cloud is bouyant 
!               --------------------------------
!
!*      6.1    Compute grid scale vertical velocity perturbation term ZWORK1
!               -------------------------------------------------------------
! 
             !  normalize w grid scale to a 25 km refer. grid
     DO JI = 1, IIE
        JK  = ILCL(JI)
        JKM = JK - 1 
        JKDL = IDPL(JI)
       !ZWORK1(JI) =  ( PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) )  &
        ZWORK1(JI) =  ( PW(JI,JK)  +  PW(JI,JKDL)*ZZLCL(JI)/PZ(JI,JKDL) ) * .5 &  
                           * SQRT( PDXDY(JI) / XA25 )
!                         - 0.02 * ZZLCL(JI) / XZLCL ! avoid spurious convection
     END DO
             ! compute sign of normalized grid scale w
        ZWORK2(:) = SIGN( 1., ZWORK1(:) ) 
        ZWORK1(:) = XWTRIG * ZWORK2(:) * ABS( ZWORK1(:) ) ** 0.333       &
                           * ( XP00 / ZPLCL(:) ) ** ZRDOCP
!
!*       6.2    Compute parcel vertical velocity at LCL
!               ---------------------------------------
!                   
     DO JI = 1, IIE
        JKDL = IDPL(JI)
        ZWORK3(JI) = XG * ZWORK1(JI) * ( ZZLCL(JI) - PZ(JI,JKDL) )       &
                       / ( PTHV(JI,JKDL) + ZTHVELCL(JI) )
     END DO
     WHERE( GWORK1(:) )
       ZWLCL(:)  = 1. + .5 * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) ) 
       GTRIG(:)  = ZTHVLCL(:) - ZTHVELCL(:) + ZWORK1(:) > 0. .AND.       &
                   ZWLCL(:) > 0. 
     END WHERE
!
!
!*       6.3    Look for parcel that produces sufficient cloud depth.
!               The cloud top is estimated as the level where the CAPE 
!               is smaller  than a given value (based on vertical velocity eq.)
!               --------------------------------------------------------------
!
     ZTHEUL(:) = ZTLCL(:) * ( ZTHLCL(:) / ZTLCL(:) )                       &
                                             ** ( 1. - 0.28 * ZRVLCL(:) )  &
                          * EXP( ( 3374.6525 / ZTLCL(:) - 2.5403 ) *       &
                               ZRVLCL(:) * ( 1. + 0.81 * ZRVLCL(:) ) )
!
     ZCAPE(:) = 0.
     ZTOP(:)  = 0.
     ZWORK3(:)= 0.
     JKM = MINVAL( ILCL(:) )
     DO JL = JKM, JT
        JK = JL + 1
        DO JI = 1, IIE
           ZWORK1(JI) = ( 2. * ZTHEUL(JI) /                                &
            ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1. ) * ( PZ(JI,JK) - PZ(JI,JL) )
           IF ( JL < ILCL(JI) ) ZWORK1(JI) = 0.
           ZCAPE(JI)  = ZCAPE(JI) + ZWORK1(JI)
           ZCAP(JI,JKK) = ZCAP(JI,JKK) + XG * MAX( 0., ZWORK1(JI) ) ! actual CAPE
           ZWORK2(JI) = XNHGAM * XG * ZCAPE(JI) + 1.05 * ZWLCL(JI) * ZWLCL(JI)
               ! the factor 1.05 takes entrainment into account
           ZWORK2(JI) = SIGN( 1., ZWORK2(JI) )
           ZWORK3(JI) = ZWORK3(JI) + MIN(0., ZWORK2(JI) )
           ZWORK3(JI) = MAX( -1., ZWORK3(JI) )
               ! Nota, the factors ZWORK2 and ZWORK3 are only used to avoid
               ! if and goto statements, the difficulty is to extract only
               ! the level where the criterium is first fullfilled
           ZTOP(JI)   = PZ(JI,JL) * .5 * ( 1. + ZWORK2(JI) ) * ( 1. + ZWORK3(JI) ) + &
                        ZTOP(JI) * .5 * ( 1. - ZWORK2(JI) )
         END DO
     END DO
!
!
     WHERE( ZTOP(:) - ZZLCL(:) .GE. XCDEPTH  .AND. GTRIG(:) .AND. GTRIG2(:) )
        GTRIG2(:)   = .FALSE.
        OTRIG(:)    = GTRIG(:)     ! we  select the first departure level
        PTHLCL(:)   = ZTHLCL(:)    ! that gives sufficient cloud depth
        PRVLCL(:)   = ZRVLCL(:)
        PTLCL(:)    = ZTLCL(:)
        PWLCL(:)    = ZWLCL(:)
        PZLCL(:)    = ZZLCL(:)
        PTHVELCL(:) = ZTHVELCL(:)
        KDPL(:)     = IDPL(:)
        KPBL(:)     = IPBL(:)
        KLCL(:)     = ILCL(:)
     END WHERE
!
END DO
!
     DO JI = 1, IIE
       PCAPE(JI) = MAXVAL( ZCAP(JI,:) ) ! maximum CAPE for diagnostics
     END DO
!
!
END SUBROUTINE CONVECT_TRIGGER_FUNCT
!     ######spl
      SUBROUTINE CONVECT_UPDRAFT( KLON, KLEV,                                      &
				  KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW, &
                                  PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,   &
                                  PMFLCL, OTRIG, KLCL, KDPL, KPBL,                 &
                                  PUMF, PUER, PUDR, PUTHL, PUTHV, PURW,            &
                                  PURC, PURI, PURR, PURS, PUPR,                    &
                                  PUTPR, PCAPE, KCTL, KETL )
!     #############################################################################
!
!!**** Compute updraft properties from DPL to CTL. 
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine updraft properties
!!      ( mass flux, thermodynamics, precipitation ) 
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_MIXING_FUNCT
!!     Routine CONVECT_CONDENS
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV, XCL    ! Cp of dry air, water vapor and liquid water
!!          XTT                ! triple point temperature
!!          XLVTT              ! vaporisation heat at XTT
!!        
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XCRAD              ! cloud radius
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XENTR              ! entrainment constant
!!          XRCONV             ! constant in precipitation conversion 
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!          XTFRZ1             ! begin of freezing interval
!!          XTFRZ2             ! begin of freezing interval
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_UPDRAFT)
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  10/12/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                    :: KLON  ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV  ! vertical dimension
INTEGER, INTENT(IN)                    :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHL  ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHV  ! grid scale theta_v     
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between 
                                                ! bottom and top of layer (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PTLCL  ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMFLCL ! cloud  base unit mass flux
                                                ! (kg/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(INOUT):: OTRIG! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! contains vert. index of DPL 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   !  " vert. index of source layertop
!
!
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KCTL   ! contains vert. index of CTL 
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KETL   ! contains vert. index of        &
                                                !equilibrium (zero buoyancy) level 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHL ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHV ! updraft theta_v (K)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURW  ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURC  ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURI  ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURR  ! liquid precipit. (kg/kg)
                                                ! produced in  model layer
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::PURS ! solid precipit. (kg/kg)
                                                ! produced in  model layer
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::PUPR ! updraft precipitation in
                                                ! flux units (kg water / s)
REAL, DIMENSION(KLON),     INTENT(OUT):: PUTPR  ! total updraft precipitation
                                                ! in flux units (kg water / s)
REAL, DIMENSION(KLON),     INTENT(OUT):: PCAPE  ! available potent. energy
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP, JKM, JK1, JK2, JKMIN  ! vertical loop index
REAL    :: ZEPSA, ZCVOCD  ! R_v / R_d, C_pv / C_pd 
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON)    :: ZUT             ! updraft temperature (K)
REAL, DIMENSION(KLON)    :: ZUW1, ZUW2      ! square of updraft vert.
                                            ! velocity at levels k and k+1
REAL, DIMENSION(KLON)    :: ZE1,ZE2,ZD1,ZD2 ! fractional entrainm./detrain
                                            ! rates at levels k and k+1
REAL, DIMENSION(KLON)    :: ZMIXF           ! critical mixed fraction  
REAL, DIMENSION(KLON)    :: ZCPH            ! specific heat C_ph 
REAL, DIMENSION(KLON)    :: ZLV, ZLS        ! latent heat of vaporis., sublim.       
REAL, DIMENSION(KLON)    :: ZURV            ! updraft water vapor at level k+1
REAL, DIMENSION(KLON)    :: ZPI             ! Pi=(P0/P)**(Rd/Cpd)  
REAL, DIMENSION(KLON)    :: ZTHEUL          ! theta_e for undilute ascent
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5,   &
                            ZWORK6          ! work arrays
INTEGER, DIMENSION(KLON) :: IWORK           ! wok array
LOGICAL, DIMENSION(KLON) :: GWORK1, GWORK2, GWORK4, GWORK5 
                                            ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK6     ! work array
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
IIE = KLON
!
!
!*       1.     Initialize updraft properties and local variables
!               -------------------------------------------------
!
ZEPSA      = XRV / XRD 
ZCVOCD     = XCPV / XCPD
ZCPORD     = XCPD / XRD
ZRDOCP     = XRD / XCPD
!
PUMF(:,:)  = 0.
PUER(:,:)  = 0.
PUDR(:,:)  = 0.
PUTHL(:,:) = 0.
PUTHV(:,:) = 0.
PURW(:,:)  = 0.
PURC(:,:)  = 0.
PURI(:,:)  = 0.
PUPR(:,:)  = 0.
PURR(:,:)  = 0.
PURS(:,:)  = 0.
PUTPR(:)   = 0.
ZUW1(:)    = PWLCL(:) * PWLCL(:)
ZUW2(:)    = 0.
ZE1(:)     = 1.
ZD1(:)     = 0.
PCAPE(:)   = 0.
KCTL(:)    = IKB
KETL(:)    = KLCL(:)
GWORK2(:)  = .TRUE.
GWORK5(:)  = .TRUE.
ZPI(:)     = 1.
ZWORK3(:)  = 0.
ZWORK4(:)  = 0.
ZWORK5(:)  = 0.
ZWORK6(:)  = 0.
GWORK1(:)  = .FALSE.
GWORK4(:)  = .FALSE.
!
!
!*       1.1    Compute undilute updraft theta_e for CAPE computations
!               Bolton (1980) formula.
!               Define accurate enthalpy for updraft
!               -----------------------------------------------------
!
ZTHEUL(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** ( 1. - 0.28 * PRVLCL(:) )  &
            * EXP( ( 3374.6525 / PTLCL(:) - 2.5403 ) *                        &
                                   PRVLCL(:) * ( 1. + 0.81 * PRVLCL(:) ) )
!
!
ZWORK1(:) = ( XCPD + PRVLCL(:) * XCPV ) * PTLCL(:)                            &
            + ( 1. + PRVLCL(:) ) * XG * PZLCL(:)
!
!
!*       2.     Set updraft properties between DPL and LCL
!               ------------------------------------------
!
JKP = MAXVAL( KLCL(:) )
JKM = MINVAL( KDPL(:) )
DO JK = JKM, JKP
   DO JI = 1, IIE
    IF ( JK >= KDPL(JI) .AND. JK < KLCL(JI) ) THEN
        PUMF(JI,JK)  = PMFLCL(JI)
        PUTHL(JI,JK) = ZWORK1(JI) 
        PUTHV(JI,JK) = PTHLCL(JI) * ( 1. + ZEPSA * PRVLCL(JI) ) /             &
                                  ( 1. + PRVLCL(JI) )
        PURW(JI,JK)  = PRVLCL(JI) 
   END IF
   END DO
END DO                        
!
!
!*       3.     Enter loop for updraft computations
!               ------------------------------------
!
JKMIN = MINVAL( KLCL(:) - 1 )
DO JK = MAX( IKB + 1, JKMIN ), IKE - 1
  ZWORK6(:) = 1.
  JKP = JK + 1  
!
  GWORK4(:) = JK >= KLCL(:) - 1 
  GWORK1(:) = GWORK4(:) .AND. GWORK2(:) ! this mask is used to confine
                           ! updraft computations between the LCL and the CTL
!                                                         
  WHERE( JK == KLCL(:) - 1 ) ZWORK6(:) = 0. ! factor that is used in buoyancy
                                        ! computation at first level above LCL
!
!
!*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1   
!               ----------------------------------------------------------
!
    ZWORK1(:) = PURC(:,JK) + PURR(:,JK)
    ZWORK2(:) = PURI(:,JK) + PURS(:,JK)
    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), PUTHL(:,JK), PURW(:,JK),&
                          ZWORK1, ZWORK2, PZ(:,JKP), GWORK1, ZUT, ZURV,     &
                          PURC(:,JKP), PURI(:,JKP), ZLV, ZLS, ZCPH )
!
!
  ZPI(:) = ( XP00 / PPRES(:,JKP) ) ** ZRDOCP   
  WHERE ( GWORK1(:) )
!
    PUTHV(:,JKP) = ZPI(:) * ZUT(:) * ( 1. + ZEPSA * ZURV(:) )           &  
                         / ( 1. + PURW(:,JK) )     
!
!
!*       5.     Compute square of vertical velocity using entrainment   
!               at level k
!               -----------------------------------------------------
!    
    ZWORK3(:) = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -         &
                     ( 1. - ZWORK6(:) ) * PZLCL(:)          ! level thickness  
    ZWORK4(:) = PTHV(:,JK) * ZWORK6(:) +                   &
                 ( 1. - ZWORK6(:) ) * PTHVELCL(:)
    ZWORK5(:) = 2. * ZUW1(:) * PUER(:,JK) / MAX( .1, PUMF(:,JK) )
    ZUW2(:)   = ZUW1(:) + ZWORK3(:) * XNHGAM * XG *        & 
                  ( ( PUTHV(:,JK) + PUTHV(:,JKP) ) /       &
                  ( ZWORK4(:) + PTHV(:,JKP) ) - 1. )       & ! buoyancy term
                - ZWORK5(:)                                  ! entrainment term
!
!
!*       6.     Update total precipitation: dr_r=(r_c+r_i)*exp(-rate*dz)  
!               --------------------------------------------------------
!
!                    compute level mean vertical velocity
    ZWORK2(:)   = 0.5 *                                                    &
                       ( SQRT( MAX( 1.E-2, ZUW2(:) ) ) +                   &
                         SQRT( MAX( 1.E-2, ZUW1(:) ) ) )          
    PURR(:,JKP) = 0.5 * ( PURC(:,JK) + PURC(:,JKP) + PURI(:,JK) + PURI(:,JKP) )&
                      * ( 1. - EXP( - XRCONV  * ZWORK3(:) / ZWORK2(:) ) )
    PUPR(:,JKP) = PURR(:,JKP) * PUMF(:,JK) ! precipitation rate ( kg water / s)
    PUTPR(:)    = PUTPR(:) + PUPR(:,JKP)   ! total precipitation rate
    ZWORK2(:)   = PURR(:,JKP) / MAX( 1.E-8, PURC(:,JKP) + PURI(:,JKP) )
    PURR(:,JKP) = ZWORK2(:) * PURC(:,JKP)          ! liquid precipitation
    PURS(:,JKP) = ZWORK2(:) * PURI(:,JKP)          ! solid precipitation
!
!
!*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation 
!               -------------------------------------------------------
!
    PURW(:,JKP)  = PURW(:,JK) - PURR(:,JKP) - PURS(:,JKP) 
    PURC(:,JKP)  = PURC(:,JKP) - PURR(:,JKP)
    PURI(:,JKP)  = PURI(:,JKP) - PURS(:,JKP)       
    PUTHL(:,JKP) = ( XCPD + PURW(:,JKP) * XCPV ) * ZUT(:)                     &
                   + ( 1. + PURW(:,JKP) ) * XG * PZ(:,JKP)                    &
                   - ZLV(:) * PURC(:,JKP) - ZLS(:) * PURI(:,JKP)             
!    
    ZUW1(:)      = ZUW2(:)       
!
  END WHERE
!
!
!*       8.     Compute entrainment and detrainment using conservative
!               variables adjusted for precipitation ( not for entrainment)
!               -----------------------------------------------------------
!
!*       8.1    Compute critical mixed fraction by estimating unknown  
!               T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
!               We determine the zero crossing of the linear curve
!               evaluating the derivative using ZMIXF=0.1.
!               -----------------------------------------------------
!    
    ZMIXF(:)  = 0.1   ! starting value for critical mixed fraction
    ZWORK1(:) = ZMIXF(:) * PTHL(:,JKP)                                     &
                     + ( 1. - ZMIXF(:) ) * PUTHL(:,JKP) ! mixed enthalpy
    ZWORK2(:) = ZMIXF(:) * PRW(:,JKP)                                      &
                     + ( 1. - ZMIXF(:) ) * PURW(:,JKP)  ! mixed r_w
!
    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), ZWORK1, ZWORK2,        &
                          PURC(:,JKP), PURI(:,JKP), PZ(:,JKP), GWORK1, ZUT,&
                          ZWORK3, ZWORK4, ZWORK5, ZLV, ZLS, ZCPH )
!        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)
!        
     ! compute theta_v of mixture
    ZWORK3(:) = ZUT(:) * ZPI(:) * ( 1. + ZEPSA * (                         &
                ZWORK2(:) - ZWORK4(:) - ZWORK5(:) ) ) / ( 1. + ZWORK2(:) )
     ! compute final value of critical mixed fraction using theta_v
     ! of mixture, grid-scale and updraft
    ZMIXF(:) = MAX( 0., PUTHV(:,JKP) - PTHV(:,JKP) ) * ZMIXF(:) /          &
                              ( PUTHV(:,JKP) - ZWORK3(:) + 1.E-10 )
    ZMIXF(:) = MAX( 0., MIN( 1., ZMIXF(:) ) )
!    
!
!*       8.2     Compute final midlevel values for entr. and detrainment    
!                after call of distribution function
!                -------------------------------------------------------
!    
!
    CALL CONVECT_MIXING_FUNCT ( KLON, ZMIXF, 1, ZE2, ZD2 )
!       Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates
!
! ZWORK1(:) = XENTR * PMFLCL(:) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  zwork1(:) = xentr * xg / xcrad * pumf(:,jk) * ( pz(:,jkp) - pz(:,jk) )
! ZWORK1(:) = XENTR * pumf(:,jk) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  ZWORK2(:) = 0.
  WHERE ( GWORK1(:) ) ZWORK2(:) = 1.
  WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) )
    ze2=.5; zd2=.5  ! modif entrainment=detrainment, this avoids
		    ! too large mass flux values at upper levels
    PUER(:,JKP) = 0.5 * ZWORK1(:) * ( ZE1(:) + ZE2(:) ) * ZWORK2(:)
    PUDR(:,JKP) = 0.5 * ZWORK1(:) * ( ZD1(:) + ZD2(:) ) * ZWORK2(:)
  ELSEWHERE
    PUER(:,JKP) = 0.
    PUDR(:,JKP) = ZWORK1(:) * ZWORK2(:)
  END WHERE
!
!*       8.3     Determine equilibrium temperature level
!                --------------------------------------
!
   WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) .AND. JK > KLCL(:) + 1 &   
           .AND. GWORK1(:) )
         KETL(:) = JKP            ! equilibrium temperature level 
   END WHERE
!
!*       8.4     If the calculated detrained mass flux is greater than    
!                the total updraft mass flux, or vertical velocity is
!                negative, all cloud mass detrains at previous model level,
!                exit updraft calculations - CTL is attained
!                -------------------------------------------------------
!
  WHERE( GWORK1(:) )                                                   &
        GWORK2(:) = PUMF(:,JK) - PUDR(:,JKP) > 10. .AND. ZUW2(:) > 0.        
  WHERE ( GWORK2(:) ) KCTL(:) = JKP   ! cloud top level
  GWORK1(:) = GWORK2(:) .AND. GWORK4(:)
!
  IF ( COUNT( GWORK2(:) ) == 0 ) EXIT           
!
!
!*       9.   Compute CAPE for undilute ascent using theta_e and 
!             theta_es instead of theta_v. This estimation produces 
!             a significantly larger value for CAPE than the actual one.
!             ----------------------------------------------------------
!
  WHERE ( GWORK1(:) )
!
    ZWORK3(:)   = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -                      &
                  ( 1. - ZWORK6(:) ) *  PZLCL(:)              ! level thickness
    ZWORK2(:)   = PTHES(:,JK) + ( 1. - ZWORK6(:) ) *                      &
     ( PTHES(:,JKP) - PTHES(:,JK) ) / ( PZ(:,JKP) - PZ(:,JK) ) *          &
     ( PZLCL(:) - PZ(:,JK) ) ! linear interpolation for theta_es at LCL
                            ! ( this is only done for model level just above LCL
!
    ZWORK1(:) = ( 2. * ZTHEUL(:) ) / ( ZWORK2(:) + PTHES(:,JKP) ) - 1.   
    PCAPE(:)  = PCAPE(:) + XG * ZWORK3(:) * MAX( 0., ZWORK1(:) )
!
!
!*       10.   Compute final values of updraft mass flux, enthalpy, r_w 
!              at level k+1    
!              --------------------------------------------------------
!    
    PUMF(:,JKP)  = PUMF(:,JK) - PUDR(:,JKP) + PUER(:,JKP) 
    PUMF(:,JKP)  = MAX( PUMF(:,JKP), 0.1 )
    PUTHL(:,JKP) = ( PUMF(:,JK) * PUTHL(:,JK) +                              &
                     PUER(:,JKP) * PTHL(:,JK) - PUDR(:,JKP) * PUTHL(:,JK) )  &
                    / PUMF(:,JKP) + PUTHL(:,JKP) - PUTHL(:,JK)
    PURW(:,JKP)  = ( PUMF(:,JK) * PURW(:,JK) +                               &
                     PUER(:,JKP) * PRW(:,JK) - PUDR(:,JKP) * PURW(:,JK) )    &
                    / PUMF(:,JKP) - PURR(:,JKP) - PURS(:,JKP)
!    
    ZE1(:) = ZE2(:) ! update fractional entrainment/detrainment
    ZD1(:) = ZD2(:)
!
  END WHERE
!
END DO
!
!*       12.1    Set OTRIG to False if cloud thickness < XCDEPTH
!                or CAPE < 1
!                ------------------------------------------------
!
    DO JI = 1, IIE
          JK  = KCTL(JI)
          OTRIG(JI) = PZ(JI,JK) - PZLCL(JI) >= XCDEPTH               &
                     .AND. PCAPE(JI) > 1. 
    END DO
    WHERE( .NOT. OTRIG(:) )
          KCTL(:) = IKB 
    END WHERE
KETL(:) = MAX( KETL(:), KLCL(:) + 2 )
KETL(:) = MIN( KETL(:), KCTL(:) )
!
!
!*       12.2    If the ETL and CTL are the same detrain updraft mass   
!                flux at this level
!                ------------------------------------------------------- 
!
ZWORK1(:) = 0.
WHERE ( KETL(:) == KCTL(:) ) ZWORK1(:) = 1.
!
DO JI = 1, IIE
    JK = KETL(JI) 
    PUDR(JI,JK)   = PUDR(JI,JK) +                                    &
                          ( PUMF(JI,JK) - PUER(JI,JK) )  * ZWORK1(JI)  
    PUER(JI,JK)   = PUER(JI,JK) * ( 1. - ZWORK1(JI) )
    PUMF(JI,JK)   = PUMF(JI,JK) * ( 1. - ZWORK1(JI) )
    JKP = KCTL(JI) + 1
    PUER(JI,JKP)  = 0. ! entrainm/detr rates have been already computed
    PUDR(JI,JKP)  = 0. ! at level KCTL+1, set them to zero
END DO
!    
!*       12.3    Adjust mass flux profiles, detrainment rates, and   
!                precipitation fallout rates to reflect linear decrease
!                in mass flux between the ETL and CTL
!                -------------------------------------------------------        
! 
ZWORK1(:) = 0.
JK1 = MINVAL( KETL(:) )
JK2 = MAXVAL( KCTL(:) )
DO JK = JK1, JK2
    DO JI = 1, IIE
    IF( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
        ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
    END IF
    END DO
END DO
!
DO JI = 1, IIE
    JK = KETL(JI) 
    ZWORK1(JI) = PUMF(JI,JK) / MAX( 1., ZWORK1(JI) )
END DO
!
DO JK = JK1 + 1, JK2
    JKP = JK - 1
    DO JI = 1, IIE
    IF ( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
      ! PUTPR(JI)    = PUTPR(JI) - ( PURR(JI,JK) + PURS(JI,JK) ) * PUMF(JI,JKP)      
        PUTPR(JI)    = PUTPR(JI) - PUPR(JI,JK)
        PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
        PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
        PUPR(JI,JK)  = PUMF(JI,JKP) * ( PURR(JI,JK) + PURS(JI,JK) )
        PUTPR(JI)    = PUTPR(JI) + PUPR(JI,JK)
    END IF
    END DO
END DO
!
!         12.4   Set mass flux and entrainment in the source layer.
!                Linear increase throughout the source layer.
!                -------------------------------------------------------
!
!IWORK(:) = MIN( KPBL(:), KLCL(:) - 1 )
IWORK(:) = KPBL(:)
DO JI = 1, IIE
     JK  = KDPL(JI)
     JKP = IWORK(JI)
!          mixed layer depth
     ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP) + PDPRES(JI,JK)
END DO
!
JKP = MAXVAL( IWORK(:) )
DO JK = JKM, JKP
   DO JI = 1, IIE
   IF ( JK >= KDPL(JI)  .AND. JK <= IWORK(JI) ) THEN
       PUER(JI,JK) = PUER(JI,JK) + PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
       PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
   END IF
   END DO
END DO
!
!
!*       13.   If cloud thickness is smaller than  3 km, no
!              convection is allowed
!              Nota: For technical reasons, we stop the convection
!                    computations in this case and do not go back to
!                    TRIGGER_FUNCT to look for the next unstable LCL
!                    which could produce a thicker cloud.
!              ---------------------------------------------------
!
GWORK6(:,:) = SPREAD( OTRIG(:), DIM=2, NCOPIES=KLEV )
WHERE ( .NOT. OTRIG(:) ) PUTPR(:) = 0.
WHERE ( .NOT. GWORK6(:,:) )
    PUMF(:,:)  = 0.
    PUDR(:,:)  = 0.
    PUER(:,:)  = 0.
    PUTHL(:,:) = PTHL(:,:)
    PURW(:,:)  = PRW(:,:)
    PUPR(:,:)  = 0.
    PURC(:,:)  = 0.
    PURI(:,:)  = 0.
    PURR(:,:)  = 0.
    PURS(:,:)  = 0.
END WHERE
!
END SUBROUTINE CONVECT_UPDRAFT
!     ######spl
      SUBROUTINE CONVECT_CONDENS( KLON,                                           &
                                  KICE, PPRES, PTHL, PRW, PRCO, PRIO, PZ, OWORK1, &
                                  PT, PEW, PRC, PRI, PLV, PLS, PCPH   )
!     ###########################################################################
!
!!**** Compute temperature cloud and ice water content from enthalpy and r_w 
!!
!!
!!    PURPOSE
!!    -------
!!     The purpose of this routine is to determine cloud condensate
!!     and to return values for L_v, L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!     Condensate is extracted iteratively 
!!     
!!
!!    EXTERNAL
!!    --------
!!     None
!!     
!!
!!    IMPLICIT ARGUMENTS     
!!    ------------------
!!
!!      Module MODD_CST
!!          XG                   ! gravity constant
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XALPI, XBETAI, XGAMI ! constants for ice saturation pressure
!!          XP00                 ! reference pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XCL, XCI             ! specific heat for liquid water and ice
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR
!!          XTFRZ1               ! begin of freezing interval
!!          XTFRZ2               ! end of freezing interval
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CONDENS)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                :: KLON    ! horizontal loop index
INTEGER, INTENT(IN)                :: KICE    ! flag for ice ( 1 = yes,
                                              !                0 = no ice )
REAL, DIMENSION(KLON),   INTENT(IN) :: PPRES  ! pressure
REAL, DIMENSION(KLON),   INTENT(IN) :: PTHL   ! enthalpy (J/kg)
REAL, DIMENSION(KLON),   INTENT(IN) :: PRW    ! total water mixing ratio  
REAL, DIMENSION(KLON),   INTENT(IN) :: PRCO   ! cloud water estimate (kg/kg)
REAL, DIMENSION(KLON),   INTENT(IN) :: PRIO   ! cloud ice   estimate (kg/kg)
REAL, DIMENSION(KLON),   INTENT(IN) :: PZ     ! level height (m)
LOGICAL, DIMENSION(KLON),INTENT(IN) :: OWORK1 ! logical mask         
!
!
REAL, DIMENSION(KLON),   INTENT(OUT):: PT     ! temperature   
REAL, DIMENSION(KLON),   INTENT(OUT):: PRC    ! cloud water mixing ratio(kg/kg)
REAL, DIMENSION(KLON),   INTENT(OUT):: PRI    ! cloud ice mixing ratio  (kg/kg)
REAL, DIMENSION(KLON),   INTENT(OUT):: PLV    ! latent heat L_v    
REAL, DIMENSION(KLON),   INTENT(OUT):: PLS    ! latent heat L_s  
REAL, DIMENSION(KLON),   INTENT(OUT):: PCPH   ! specific heat C_ph   
REAL, DIMENSION(KLON),   INTENT(OUT):: PEW    ! water saturation mixing ratio  
!
!*       0.2   Declarations of local variables KLON
!
INTEGER :: JITER          ! iteration index
REAL    :: ZEPS, ZEPSA    ! R_d / R_v, 1 / ZEPS
REAL    :: ZCVOCD         ! XCPV / XCPD
REAL    :: ZRDOCP         ! R_d / C_pd
!
REAL, DIMENSION(KLON)    :: ZEI           ! ice saturation mixing ratio
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZT ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!*       1.     Initialize temperature and Exner function
!               -----------------------------------------
!
ZRDOCP      = XRD / XCPD  
ZEPS        = XRD / XRV
ZEPSA       = 1. / ZEPS
ZCVOCD      = XCPV / XCPD
!
!
    ! Make a first temperature estimate, based e.g. on values of
    !  r_c and r_i at lower level
!
      !! Note that the definition of ZCPH is not the same as used in
      !! routine CONVECT_SATMIXRATIO
     PCPH(:)   = XCPD + XCPV * PRW(:)
     ZWORK1(:) = ( 1. + PRW(:) ) * XG * PZ(:)
     PT(:)     = ( PTHL(:) + PRCO(:) * XLVTT + PRIO(:) * XLSTT - ZWORK1(:) )   &
                 / PCPH(:)
     PT(:)     = MAX(180., MIN( 330., PT(:) ) ) ! set overflow bounds in
                                                    ! case that PTHL=0     
!
!
!*       2.     Enter the iteration loop
!               ------------------------
!    
DO JITER = 1,6
     PEW(:) = EXP( XALPW - XBETAW / PT(:) - XGAMW * LOG( PT(:) ) )
     ZEI(:) = EXP( XALPI - XBETAI / PT(:) - XGAMI * LOG( PT(:) ) )
     PEW(:) = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )
     ZEI(:) = ZEPS * ZEI(:) / ( PPRES(:) - ZEI(:) )    
!
     PLV(:)    = XLVTT + ( XCPV - XCL ) * ( PT(:) - XTT ) ! compute L_v
     PLS(:)    = XLSTT + ( XCPV - XCI ) * ( PT(:) - XTT ) ! compute L_i
!    
     ZWORK2(:) = ( PT(:) - XTFRZ2 ) / ( XTFRZ1 - XTFRZ2 ) ! freezing interval
     ZWORK2(:) = MAX( 0., MIN(1., ZWORK2(:) ) ) 
     ZWORK2(:) = ZWORK2(:) * ZWORK2(:)
     IF ( KICE == 0 ) ZWORK2(:) = 1.
     ZWORK3(:) = ( 1. - ZWORK2(:) ) * ZEI(:) + ZWORK2(:) * PEW(:)
     PRC(:)    = MAX( 0., ZWORK2(:) * ( PRW(:) - ZWORK3(:) ) )
     PRI(:)    = MAX( 0., ( 1. - ZWORK2(:) ) * ( PRW(:) - ZWORK3(:) ) )
     ZT(:)     = ( PTHL(:) + PRC(:) * PLV(:) + PRI(:) * PLS(:) - ZWORK1(:) )   &
                 / PCPH(:)
     PT(:) = PT(:) + ( ZT(:) - PT(:) ) * 0.4  ! force convergence
     PT(:) = MAX( 175., MIN( 330., PT(:) ) )
END DO
!
!
END SUBROUTINE CONVECT_CONDENS
!     ######spl
      SUBROUTINE CONVECT_SATMIXRATIO( KLON,                          &
                                      PPRES, PT, PEW, PLV, PLS, PCPH )      
!     ################################################################
!
!!**** Compute vapor saturation mixing ratio over liquid water
!!
!!
!!    PDRPOSE
!!    -------
!!     The purpose of this routine is to determine saturation mixing ratio
!!     and to return values for L_v L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!     None
!!
!!
!!    IMPLICIT ARGUMENTS    
!!    ------------------
!!      Module MODD_CST
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XCL, XCI             ! specific heat for liquid water and ice
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_SATMIXRATIO)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!------------------------- ------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                INTENT(IN) :: KLON    ! horizontal loop index
REAL, DIMENSION(KLON),  INTENT(IN) :: PPRES   ! pressure
REAL, DIMENSION(KLON),  INTENT(IN) :: PT      ! temperature   
!
REAL, DIMENSION(KLON),  INTENT(OUT):: PEW     ! vapor saturation mixing ratio
REAL, DIMENSION(KLON),  INTENT(OUT):: PLV     ! latent heat L_v    
REAL, DIMENSION(KLON),  INTENT(OUT):: PLS     ! latent heat L_s  
REAL, DIMENSION(KLON),  INTENT(OUT):: PCPH    ! specific heat C_ph   
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(KLON)              :: ZT      ! temperature   
REAL    :: ZEPS           ! R_d / R_v
!
!
!-------------------------------------------------------------------------------
!
    ZEPS      = XRD / XRV
!
    ZT(:)     = MIN( 400., MAX( PT(:), 10. ) ) ! overflow bound
    PEW(:)    = EXP( XALPW - XBETAW / ZT(:) - XGAMW * LOG( ZT(:) ) )
    PEW(:)    = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )
!
    PLV(:)    = XLVTT + ( XCPV - XCL ) * ( ZT(:) - XTT ) ! compute L_v
    PLS(:)    = XLSTT + ( XCPV - XCI ) * ( ZT(:) - XTT ) ! compute L_i
!    
    PCPH(:)   = XCPD + XCPV * PEW(:)                     ! compute C_ph 
!
END SUBROUTINE CONVECT_SATMIXRATIO
!     ######spl
      SUBROUTINE CONVECT_MIXING_FUNCT( KLON,                &
                                       PMIXC, KMF, PER, PDR ) 
!     #######################################################
!
!!**** Determine the area under the distribution function
!!     KMF = 1 : gaussian  KMF = 2 : triangular distribution function
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the entrainment and
!!      detrainment rate by evaluating the are under the distribution 
!!      function. The integration interval is limited by the critical
!!      mixed fraction PMIXC
!!   
!!
!!
!!**  METHOD
!!    ------
!!      Use handbook of mathemat. functions by Abramowitz and Stegun, 1968
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine MIXING_FUNCT)
!!      Abramovitz and Stegun (1968), handbook of math. functions 
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,               INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,               INTENT(IN) :: KMF    ! switch for dist. function
REAL, DIMENSION(KLON), INTENT(IN) :: PMIXC  ! critical mixed fraction
!
REAL, DIMENSION(KLON), INTENT(OUT):: PER    ! normalized entrainment rate
REAL, DIMENSION(KLON), INTENT(OUT):: PDR    ! normalized detrainment rate
!
!*       0.2   Declarations of local variables :
!
REAL    :: ZSIGMA = 0.166666667                   ! standard deviation 
REAL    :: ZFE    = 4.931813949                   ! integral normalization 
REAL    :: ZSQRTP = 2.506628,  ZP  = 0.33267      ! constants
REAL    :: ZA1    = 0.4361836, ZA2 =-0.1201676    ! constants
REAL    :: ZA3    = 0.9372980, ZT1 = 0.500498     ! constants
REAL    :: ZE45   = 0.01111                       ! constant
!
REAL, DIMENSION(KLON) :: ZX, ZY, ZW1, ZW2         ! work variables
REAL    :: ZW11
!
!
!-------------------------------------------------------------------------------
!
!       1.     Use gaussian function for KMF=1
!              -------------------------------
!
IF( KMF == 1 ) THEN 
    ! ZX(:)  = ( PMIXC(:) - 0.5 ) / ZSIGMA
      ZX(:)  = 6. * PMIXC(:) - 3.
      ZW1(:) = 1. / ( 1.+ ZP * ABS ( ZX(:) ) )
      ZY(:)  = EXP( -0.5 * ZX(:) * ZX(:) )
      ZW2(:) = ZA1 * ZW1(:) + ZA2 * ZW1(:) * ZW1(:) +                   &
               ZA3 * ZW1(:) * ZW1(:) * ZW1(:)
      ZW11   = ZA1 * ZT1 + ZA2 * ZT1 * ZT1 + ZA3 * ZT1 * ZT1 * ZT1
ENDIF 
!
WHERE ( KMF == 1 .AND. ZX(:) >= 0. )
	PER(:) = ZSIGMA * ( 0.5 * ( ZSQRTP - ZE45 * ZW11                 &
		 - ZY(:) * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )        &
		 - 0.5 * ZE45 * PMIXC(:) * PMIXC(:)
	PDR(:) = ZSIGMA*( 0.5 * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )       &
		 + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
		 - ZE45 * ( 0.5 + 0.5 * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE
WHERE ( KMF == 1 .AND. ZX(:) < 0. ) 
	PER(:) = ZSIGMA*( 0.5 * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )       &
		 + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
		 - 0.5 * ZE45 * PMIXC(:) * PMIXC(:)
	PDR(:) = ZSIGMA * ( 0.5 * ( ZSQRTP - ZE45 * ZW11 - ZY(:)         &
		 * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )                &
		 - ZE45 * ( 0.5 + 0.5 * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE
!
      PER(:) = PER(:) * ZFE
      PDR(:) = PDR(:) * ZFE
!
!
!       2.     Use triangular function KMF=2
!              -------------------------------
!
!     not yet released
!
!
END SUBROUTINE CONVECT_MIXING_FUNCT
!     ######spl
      SUBROUTINE CONVECT_TSTEP_PREF( KLON, KLEV,                           &
				     PU, PV, PPRES, PZ, PDXDY, KLCL, KCTL, &
                                     PTIMEA, PPREF )
!     ######################################################################
!
!!**** Routine to compute convective advection time step and precipitation 
!!     efficiency 
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      advection time step PTIMEC as a function of the mean ambient 
!!      wind as well as the precipitation efficiency as a function
!!      of wind shear and cloud base height.
!!
!!
!!**  METHOD
!!    ------
!!     
!!
!!    EXTERNAL
!!    --------
!!     None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation 
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                    :: KLON   ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV   ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PU     ! grid scale horiz. wind u 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PV     ! grid scale horiz. wind v
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY  ! grid area (m^2)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL   ! lifting condensation level index
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL   ! cloud top level index
!
REAL, DIMENSION(KLON),      INTENT(OUT):: PTIMEA ! advective time period
REAL, DIMENSION(KLON),      INTENT(OUT):: PPREF  ! precipitation efficiency 
!
!
!*       0.2   Declarations of local variables KLON
!
INTEGER :: IIE, IKB, IKE                      ! horizontal + vertical loop bounds
INTEGER :: JI                                 ! horizontal loop index
INTEGER :: JK, JKLC, JKP5, JKCT               ! vertical loop index
!
INTEGER, DIMENSION(KLON)  :: IP500       ! index of 500 hPa levels
REAL, DIMENSION(KLON)     :: ZCBH        ! cloud base height 
REAL, DIMENSION(KLON)     :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Determine vertical index for 500 hPa levels 
!               ------------------------------------------
!
!
IP500(:) = IKB
DO JK = IKB, IKE
    WHERE ( PPRES(:,JK) >= 500.E2 ) IP500(:) = JK
END DO
!
!
!*       2.     Compute convective time step 
!               ----------------------------
!
	    ! compute wind speed at LCL, 500 hPa, CTL

DO JI = 1, IIE
   JKLC = KLCL(JI)
   JKP5 = IP500(JI)
   JKCT = KCTL(JI)
   ZWORK1(JI) = SQRT( PU(JI,JKLC) * PU(JI,JKLC) +           &
		      PV(JI,JKLC) * PV(JI,JKLC)  ) 
   ZWORK2(JI) = SQRT( PU(JI,JKP5) * PU(JI,JKP5) +           &
		      PV(JI,JKP5) * PV(JI,JKP5)  ) 
   ZWORK3(JI) = SQRT( PU(JI,JKCT) * PU(JI,JKCT) +           &
		      PV(JI,JKCT) * PV(JI,JKCT)  ) 
END DO
!
ZWORK2(:) = MAX( 0.1, 0.5 * ( ZWORK1(:) + ZWORK2(:) ) )
!
PTIMEA(:) = SQRT( PDXDY(:) ) / ZWORK2(:) 
!
!
!*       3.     Compute precipitation efficiency 
!               -----------------------------------
!
!*       3.1    Precipitation efficiency as a function of wind shear
!               ----------------------------------------------------
!
ZWORK2(:) = SIGN( 1., ZWORK3(:) - ZWORK1(:) )
DO JI = 1, IIE
    JKLC = KLCL(JI)
    JKCT = KCTL(JI)
    ZWORK1(JI) = ( PU(JI,JKCT) - PU(JI,JKLC) )  *          &
                 ( PU(JI,JKCT) - PU(JI,JKLC) )  +          &
                 ( PV(JI,JKCT) - PV(JI,JKLC) )  *          &
                 ( PV(JI,JKCT) - PV(JI,JKLC) )  
    ZWORK1(JI) = 1.E3 * ZWORK2(JI) * SQRT( ZWORK1(JI) ) /  &
	         MAX( 1.E-2, PZ(JI,JKCT) - PZ(JI,JKLC) )
END DO
!
PPREF(:)  = 1.591 + ZWORK1(:) * ( -.639 + ZWORK1(:) * (        &
				9.53E-2 - ZWORK1(:) * 4.96E-3 ) ) 
PPREF(:)  = MAX( .4, MIN( PPREF(:), .92 ) )
!
!*       3.2    Precipitation efficiency as a function of cloud base height 
!               ----------------------------------------------------------
!
DO JI = 1, IIE
   JKLC = KLCL(JI)
   ZCBH(JI)   = MAX( 3., ( PZ(JI,JKLC) - PZ(JI,IKB) ) * 3.281E-3 ) 
END DO
ZWORK1(:) = .9673 + ZCBH(:) * ( -.7003 + ZCBH(:) * ( .1622 + &
	      ZCBH(:) *  ( -1.2570E-2 + ZCBH(:) * ( 4.2772E-4 -  &
              ZCBH(:) * 5.44E-6 ) ) ) )
ZWORK1(:) = MAX( .4, MIN( .92, 1./ ( 1. + ZWORK1(:) ) ) )
!
!*       3.3    Mean precipitation efficiency is used to compute rainfall 
!               ----------------------------------------------------------
!
PPREF(:) = 0.5 * ( PPREF(:) + ZWORK1(:) )
!
!
END SUBROUTINE CONVECT_TSTEP_PREF
!     ######spl
     SUBROUTINE CONVECT_DOWNDRAFT( KLON, KLEV,                                &
                                   KICE, PPRES, PDPRES, PZ, PTH, PTHES,       &
                                   PRW, PRC, PRI,                             &
                                   PPREF, KLCL, KCTL, KETL,                   &
                                   PUTHL, PURW, PURC, PURI,                   &
                                   PDMF, PDER, PDDR, PDTHL, PDRW,             &
                                   PMIXF, PDTEVR, KLFS, KDBL, KML,            &
                                   PDTEVRF )
!    ########################################################################
!
!!**** Compute downdraft properties from LFS to DBL. 
!!
!!
!!    PDRPOSE                                                       
!!    -------
!!      The purpose of this routine is to determine downdraft properties
!!      ( mass flux, thermodynamics ) 
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from top.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO        
!!                
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XPI                ! Pi
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XCPV, XCL, XCI     ! Cp of water vapor, liquid water and ice
!!          XTT                ! triple point temperature
!!          XLVTT, XLSTT       ! vaporisation/sublimation heat at XTT
!!
!!      Module MODD_CONVPAR
!!          XCRAD              ! cloud radius
!!          XZPBL              ! thickness of downdraft detrainment layer
!!          XENTR              ! entrainment constant in pressure coordinates
!!          XRHDBC             ! relative humidity in downdraft below cloud
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_DOWNDRAFT)
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
INTEGER,                    INTENT(IN) :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! grid scale theta        
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRC   ! grid scale r_c (cloud water) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRI   ! grid scale r_i (cloud ice) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between 
						! bottom and top of layer (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! level height (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL 
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of 
						! equilibrium (zero buoyancy) level 
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KML   ! " vert. index of melting level
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUTHL ! updraft enthalpy (J/kg)      
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURC  ! updraft r_c (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURI  ! updraft r_i (kg/kg)
REAL, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDMF   ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDER   ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDDR   ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTHL  ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDRW   ! downdraft total water (kg/kg)
REAL, DIMENSION(KLON),      INTENT(OUT):: PMIXF  ! mixed fraction at LFS
REAL, DIMENSION(KLON),      INTENT(OUT):: PDTEVR ! total downdraft evaporation
                                                 ! rate at LFS (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PDTEVRF! downdraft evaporation rate
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLFS    ! contains vert. index of LFS 
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KDBL    ! contains vert. index of DBL   
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE     ! horizontal + vertical loop bounds
INTEGER :: JK, JKP, JKM, JKT ! vertical loop index
INTEGER :: JI, JL            ! horizontal loop index
INTEGER :: JITER          ! iteration loop index
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
REAL    :: ZEPS           ! R_d / R_v
REAL    :: ZEPSA, ZCVOCD  ! R_v / R_d, C_pv / C_pd
!
INTEGER, DIMENSION(KLON) :: IDDT      ! top level of detrainm. layer
REAL, DIMENSION(KLON)    :: ZTHE      ! environm. theta_e (K)
REAL, DIMENSION(KLON)    :: ZDT, ZDTP ! downdraft temperature (K)
REAL, DIMENSION(KLON)    :: ZCPH      ! specific heat C_ph 
REAL, DIMENSION(KLON)    :: ZLV, ZLS  ! latent heat of vaporis., sublim.       
REAL, DIMENSION(KLON)    :: ZDDT      ! thickness (hPa) of detrainm. layer
REAL, DIMENSION(KLON)    :: ZPI       ! Pi=(P0/P)**(Rd/Cpd)  
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZWORK4,  &
                                   ZWORK5                  ! work arrays 
LOGICAL, DIMENSION(KLON) :: GWORK1                         ! work array
!
!
!-------------------------------------------------------------------------------
!
!        0.3    Set loop bounds
!               ---------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Initialize downdraft properties
!               -------------------------------
!
ZCPORD     = XCPD / XRD
ZRDOCP     = XRD / XCPD
ZEPS       = XRD / XRV
ZEPSA      = XRV / XRD
ZCVOCD     = XCPV / XCPD
PDMF(:,:)  = 0.
PDER(:,:)  = 0.
PDDR(:,:)  = 0.
PDRW(:,:)  = 0.
PDTHL(:,:) = 0.
PDTEVR(:)  = 0.
PMIXF(:)   = 0.
ZTHE(:)    = 0.
ZDDT(:)    = PDPRES(:,IKB+2)
KDBL(:)    = IKB + 1
KLFS(:)    = IKB + 1
IDDT(:)    = KDBL(:) + 1
!  
!
!*       2.     Determine the LFS by looking for minimum of environmental 
!               saturated theta_e 
!               ----------------------------------------------------------
!
ZWORK1(:) = 900.   ! starting value for search of minimum envir. theta_e
DO JK = MINVAL( KLCL(:) ) + 2, MAXVAL( KETL(:) )
   DO JI = 1, IIE
      GWORK1(JI) = JK >= KLCL(JI) + 2 .AND. JK < KETL(JI)  
      IF ( GWORK1(JI) .AND. ZWORK1(JI) > PTHES(JI,JK) ) THEN
         KLFS(JI)   = JK
         ZWORK1(JI) = MIN( ZWORK1(JI), PTHES(JI,JK) )
      END IF
   END DO
END DO      
!
!
!*       3.     Determine the mixed fraction using environmental and updraft
!               values of theta_e at LFS
!               ---------------------------------------------------------   
!
DO JI = 1, IIE
    JK = KLFS(JI)
    ZPI(JI)    = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
      ! compute updraft theta_e
    ZWORK3(JI) = PURW(JI,JK) - PURC(JI,JK) - PURI(JI,JK)
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI) 
    ZLV(JI)    = XLVTT + ( XCPV - XCL ) * ( ZDT(JI) - XTT )                   
    ZLS(JI)    = XLSTT + ( XCPV - XCI ) * ( ZDT(JI) - XTT )                   
    ZCPH(JI)   = XCPD + XCPV * PURW(JI,JK)
    ZDT(JI)    = ( PUTHL(JI,JK) - ( 1. + PURW(JI,JK) ) * XG * PZ(JI,JK)       &
                 + ZLV(JI) * PURC(JI,JK) + ZLS(JI) * PURI(JI,JK) ) / ZCPH(JI)           
    ZWORK1(JI) = ZDT(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )              &
                                  * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )     &
                                  * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
      ! compute environmental theta_e
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI)
    ZLV(JI)    = XLVTT + ( XCPV - XCL ) * ( ZDT(JI) - XTT )                   
    ZLS(JI)    = XLSTT + ( XCPV - XCI ) * ( ZDT(JI) - XTT )                   
    ZWORK3(JI) = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
    ZCPH(JI)   = XCPD + XCPV * PRW(JI,JK)
    ZWORK2(JI) = ZDT(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )              &
                                  * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )     &
                                  * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
      ! compute mixed fraction
    PMIXF(JI)  = MAX( 0., ( ZWORK1(JI) - PTHES(JI,JK) ) )                   &
                  / ( ZWORK1(JI) - ZWORK2(JI) + 1.E-10 )
    PMIXF(JI)  = MAX(0., MIN( 1., PMIXF(JI) ) )
    ZWORK4(JI) = PPRES(JI,JK)
END DO
!
!
!*       4.     Estimate the effect of melting on the downdraft  
!               ---------------------------------------------
!
ZWORK1(:) = 0.
      ! use total solid precipitation
!DO JK = IKB + 1, IKE
!    ZWORK1(:) = ZWORK1(:) + PURS(:,JK) ! total snow/hail content
!END DO
!
DO JI = 1, IIE
     JK  = KLCL(JI)
     JKP = KCTL(JI)
     ZWORK1(JI) = 0.5 * ( PURW(JI,JK) - PURW(JI,JKP) )
END DO
!
      ! temperature perturbation due to melting at LFS
ZWORK3(:) = 0.
WHERE( KML(:) > IKB + 2 )
	  ZWORK3(:) = ZWORK1(:) * ( ZLS(:) - ZLV(:) ) / ZCPH(:)
	  ZDT(:)    = ZDT(:) - ZWORK3(:) * REAL(KICE)
END WHERE
!
!
!*       5.     Initialize humidity at LFS as a saturated mixture of
!               updraft and environmental air
!               -----------------------------------------------------    
!
DO JI = 1, IIE
     JK = KLFS(JI)
     PDRW(JI,JK)  = PMIXF(JI) * PRW(JI,JK) + ( 1. - PMIXF(JI) ) * PURW(JI,JK)
     ZWORK2(JI)   = PDRW(JI,JK) - ( 1. - PMIXF(JI) )                          &
                                     * ( PURC(JI,JK) + PURI(JI,JK) )
END DO
!
!
!*       6.1    Determine the DBL by looking for level where the envir.
!               theta_es at the LFS corrected by melting effects  becomes
!               larger than envir. value
!               ---------------------------------------------------------
!
      ! compute satur. mixing ratio for melting corrected temperature
CALL CONVECT_SATMIXRATIO( KLON, ZWORK4, ZDT, ZWORK3, ZLV, ZLS, ZCPH )  
!
      ! compute envir. saturated theta_e for melting corrected temperature
    ZWORK1(:) = MIN( ZWORK2(:), ZWORK3(:) )
    ZWORK3(:) = ZWORK3(:) * ZWORK4(:) / ( ZWORK3(:) + ZEPS ) ! sat. pressure
    ZWORK3(:) = LOG( ZWORK3(:) / 613.3 )
              ! dewp point temperature
    ZWORK3(:) = ( 4780.8 - 32.19 * ZWORK3(:) ) / ( 17.502 - ZWORK3(:) )
              ! adiabatic saturation temperature
    ZWORK3(:) = ZWORK3(:) - ( .212 + 1.571E-3 * ( ZWORK3(:) - XTT )          &
                  - 4.36E-4 * ( ZDT(:) - XTT ) ) * ( ZDT(:) - ZWORK3(:) )
    ZWORK4(:) = SIGN(0.5, ZWORK2(:) - ZWORK3(:) )
    ZDT(:)    = ZDT(:) * ( .5 + ZWORK4(:) ) + ( .5 - ZWORK4(:) ) * ZWORK3(:) 
    ZWORK2(:) = ZDT(:) * ZPI(:) ** ( 1. - 0.28 * ZWORK2(:) )                 &
                                  * EXP( ( 3374.6525 / ZDT(:) - 2.5403 )     &
                                  * ZWORK1(:) * ( 1. + 0.81 * ZWORK1(:) ) )
!
GWORK1(:) = .TRUE.
JKM = MAXVAL( KLFS(:) )
DO JK = JKM - 1, IKB + 1, -1
  DO JI = 1, IIE
     IF ( JK < KLFS(JI) .AND. ZWORK2(JI) > PTHES(JI,JK) .AND. GWORK1(JI) ) THEN
	  KDBL(JI) = JK
          GWORK1(JI) = .FALSE.
     END IF
  END DO
END DO
!
!
!*       7.     Define mass flux and entr/detr. rates at LFS
!               -------------------------------------------
!
DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI)  = PPRES(JI,JK) /                                            &
                   ( XRD * ZDT(JI) * ( 1. + ZEPS * ZWORK1(JI) ) ) ! density
     PDMF(JI,JK) = - ( 1. - PPREF(JI) ) * ZWORK1(JI) * XPI * XCRAD * XCRAD
     PDTHL(JI,JK)= ZWORK2(JI)   ! theta_l is here actually theta_e
     ZWORK2(JI)  = PDMF(JI,JK)
     PDDR(JI,JK) = 0.
     PDER(JI,JK) = - PMIXF(JI) * PDMF(JI,JK)
END DO
!
!
!         7.1   Downdraft detrainment is assumed to occur in a layer
!               of 60 hPa, determine top level IDDT of this layer
!               ---------------------------------------------------------
!
ZWORK1(:) = 0.
DO JK = IKB + 2, JKM
      ZWORK1(:) = ZWORK1(:) + PDPRES(:,JK)
      WHERE ( JK > KDBL(:) .AND. ZWORK1(:) <= XZPBL )
           ZDDT(:) = ZWORK1(:) 
           IDDT(:) = JK
      END WHERE
END DO
!
!
!*       8.     Enter loop for downdraft computations. Make a first guess
!               of initial downdraft mass flux. 
!               In the downdraft computations we use theta_es instead of 
!               enthalpy as it allows to better take into account evaporation
!               effects. As the downdraft detrainment rate is zero apart 
!               from the detrainment layer, we just compute enthalpy 
!               downdraft from theta_es in this layer.
!               ----------------------------------------------------------
!
!
ZWORK5(:) = 0.
!
DO JK =  JKM - 1, IKB + 1, -1
  JKP = JK + 1
  DO JI = 1, IIE
    IF ( JK < KLFS(JI) .AND. JK >= IDDT(JI) )  THEN
      PDER(JI,JK)  = - ZWORK2(JI) * XENTR * PDPRES(JI,JKP) / XCRAD 
                                               ! DER and DPRES are positive
      PDMF(JI,JK)  = PDMF(JI,JKP) - PDER(JI,JK) 
      ZPI(JI)      = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
      ZDT(JI)      = PTH(JI,JK) / ZPI(JI)
      ZWORK1(JI)   = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
      ZTHE(JI)     = ZDT(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK1(JI) )           &
                               * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )         &
                               * ZWORK1(JI) * ( 1. + 0.81 * ZWORK1(JI) ) )
         ! PDTHL is here theta_es, later on in this routine this table is
         ! reskipped to enthalpy 
      PDTHL(JI,JK) = ( PDTHL(JI,JKP) * PDMF(JI,JKP) - ZTHE(JI) * PDER(JI,JK)    &
                    ) / ( PDMF(JI,JK) - 1.E-7 )      
      PDRW(JI,JK)  = ( PDRW(JI,JKP) * PDMF(JI,JKP) - PRW(JI,JK) * PDER(JI,JK)   &
                    ) / ( PDMF(JI,JK) - 1.E-7 )       
    END IF
    IF ( JK < IDDT(JI) .AND. JK >= KDBL(JI) )   THEN
      JL = IDDT(JI)
      PDDR(JI,JK)  = - PDMF(JI,JL) * PDPRES(JI,JKP) / ZDDT(JI) 
      PDMF(JI,JK)  = PDMF(JI,JKP) + PDDR(JI,JK) 
      PDTHL(JI,JK) = PDTHL(JI,JKP)
      PDRW(JI,JK)  = PDRW(JI,JKP)
    END IF
  END DO
END DO
!
!
!*       9.     Calculate total downdraft evaporation 
!               rate for given mass flux (between DBL and IDDT)
!               -----------------------------------------------
!
PDTEVRF(:,:) = 0.
!
JKT = MAXVAL( IDDT(:) )
DO JK = IKB + 1, JKT
!
       ZPI(:) = ( XP00 / PPRES(:,JK) ) ** ZRDOCP
       ZDT(:) = PTH(:,JK) / ZPI(:)
!
!*       9.1    Determine wet bulb temperature at DBL from theta_e.
!               The iteration algoritm is similar to that used in
!               routine CONVECT_CONDENS
!               --------------------------------------------------
!
   DO JITER = 1, 4
       CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK1, ZLV, ZLS, ZCPH )  
       ZDTP(:) = PDTHL(:,JK) / ( ZPI(:) ** ( 1. - 0.28 * ZWORK1(:) )         &
                      * EXP( ( 3374.6525 / ZDT(:) - 2.5403 )                 &
                             * ZWORK1(:) * ( 1. + 0.81 * ZWORK1(:) ) ) )
       ZDT(:)  = 0.4 * ZDTP(:) + 0.6 * ZDT(:) ! force convergence
   END DO
!
!
!*       9.2    Sum total downdraft evaporation rate. No evaporation
!               if actual humidity is larger than specified one.
!               -----------------------------------------------------
!
   ZWORK2(:) = ZWORK1(:) / ZDT(:) * ( XBETAW / ZDT(:) - XGAMW ) ! dr_sat/dT
   ZWORK2(:) = ZLV(:) / ZCPH(:) * ZWORK1(:) * ( 1. - XRHDBC ) /              &
                    ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) ! temperature perturb                                                           ! due to evaporation
   ZDT(:)    = ZDT(:) + ZWORK2(:)
!
   CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK3, ZLV, ZLS, ZCPH )
!
   ZWORK3(:)    = ZWORK3(:) * XRHDBC
   ZWORK1(:)    = MAX( 0., ZWORK3(:) - PDRW(:,JK) ) 
   PDTEVR(:)    = PDTEVR(:) + ZWORK1(:) * PDDR(:,JK) 
   PDTEVRF(:,JK)= PDTEVRF(:,JK) + ZWORK1(:) * PDDR(:,JK) 
      ! compute enthalpie and humidity in the detrainment layer
   PDRW(:,JK)   = MAX( PDRW(:,JK), ZWORK3(:) ) 
   PDTHL(:,JK)  = ( ( XCPD + PDRW(:,JK) * XCPV ) * ZDT(:)                    &
                    + ( 1. + PDRW(:,JK) ) * XG * PZ(:,JK) ) 
!
END DO
!
!
!*      12.     If downdraft does not evaporate any water for specified 
!               relative humidity, no downdraft is allowed
!               ---------------------------------------------------------
!
ZWORK2(:) = 1.
WHERE ( PDTEVR(:) < 1. .OR. KLFS(:) == IKB + 1 ) ZWORK2(:) = 0.
DO JK = IKB, JKM
      KDBL(:)     = KDBL(:) * INT( ZWORK2(:) ) + ( 1 - INT( ZWORK2(:) ) ) * IKB
      KLFS(:)     = KLFS(:) * INT( ZWORK2(:) ) + ( 1 - INT( ZWORK2(:) ) ) * IKB
      PDMF(:,JK)  = PDMF(:,JK)  * ZWORK2(:)
      PDER(:,JK)  = PDER(:,JK)  * ZWORK2(:) 
      PDDR(:,JK)  = PDDR(:,JK)  * ZWORK2(:) 
      ZWORK1(:)   = REAL( KLFS(:) - JK )         ! use this to reset thl_d
      ZWORK1(:)   = MAX( 0.,MIN(1.,ZWORK1(:) ) ) ! and rv_d to zero above LFS
      PDTHL(:,JK) = PDTHL(:,JK) * ZWORK2(:) * ZWORK1(:)
      PDRW(:,JK)  = PDRW(:,JK)  * ZWORK2(:) * ZWORK1(:)
      PDTEVR(:)   = PDTEVR(:)   * ZWORK2(:)
      PDTEVRF(:,JK)= PDTEVRF(:,JK) * ZWORK2(:)
END DO
!
END SUBROUTINE CONVECT_DOWNDRAFT
!     ######spl
      SUBROUTINE CONVECT_PRECIP_ADJUST( KLON, KLEV,                        &
                                        PPRES, PUMF, PUER, PUDR,           &
                                        PUPR, PUTPR, PURW,                 &
                                        PDMF, PDER, PDDR, PDTHL, PDRW,     &
                                        PPREF, PTPR, PMIXF, PDTEVR,        &
                                        KLFS, KDBL, KLCL, KCTL, KETL,      &
                                        PDTEVRF )
!     ######################################################################
!
!!**** Adjust up- and downdraft mass fluxes to be consistent with the
!!     mass transport at the LFS given by the precipitation efficiency
!!     relation. 
!!
!!
!!    PURPOSE                                                       
!!    -------
!!      The purpose of this routine is to adjust up- and downdraft mass
!!      fluxes below the LFS to be consistent with the precipitation
!!      efficiency relation
!!
!!
!!
!!**  METHOD
!!    ------
!!      
!!
!!    EXTERNAL
!!    --------
!!     None
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!     Module MODD_CONVPAR
!!        XUSRDPTH             ! pressure depth to compute updraft humidity
!!                             ! supply rate for downdraft
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_PRECIP_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT
USE MODD_CONVPAR
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PUTPR ! updraft  total precipit. (kg/s
REAL, DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
REAL, DIMENSION(KLON),      INTENT(IN) :: PMIXF ! critical mixed fraction at LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of equilibrium 
						! (zero buoyancy) level 
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KLFS ! contains vert. index of LFS
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KDBL ! contains vert. index of DBL
!
REAL, DIMENSION(KLON),      INTENT(INOUT) :: PDTEVR ! total downdraft evaporation
                                                    ! rate at LFS   
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTEVRF! downdraft evaporation rate
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER   ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR   ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUPR   ! updraft  precipit. (kg/s)     
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER   ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR   ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTHL  ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDRW   ! downdraft total water (kg/kg)
!
REAL, DIMENSION(KLON),     INTENT(OUT)   :: PTPR    ! total precipitation (kg/s) 
                                                 ! = downdraft precipitation
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE        ! horizontal + vertical loop bounds
INTEGER :: JK, JKT1, JKT2, JKT3 ! vertical loop index
INTEGER :: JI                   ! horizontal loop index
!
INTEGER, DIMENSION(KLON) :: IPRL
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3,     &
				    ZWORK4, ZWORK5, ZWORK6 ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
IKB  = 1 + JCVEXB 
IKE  = KLEV - JCVEXT 
IIE  = KLON
JKT1 = MAXVAL( KLFS(:) )
JKT2 = MAXVAL( KCTL(:) )
JKT3 = MINVAL( KLCL(:) )
!
!
!        1.    Set some output variables for columns where no downdraft 
!              exists. Exit if there is no downdraft at all.
!              --------------------------------------------------------
!
IPRL(:) = IKB
PTPR(:) = 0.
!
WHERE ( PDTEVR(:) == 0. )
     PTPR(:)    = PUTPR(:)  ! no downdraft evaporation => no downdraft, all
			    ! precipitation occurs in updraft
END WHERE
IF ( COUNT( PDTEVR(:) > 0. ) == 0 ) RETURN ! exit routine if no downdraft exists
!
!*       2.     The total mass transported from the updraft to the down-  
!               draft at the LFS must be consistent with the three water
!               budget terms :
!               ---------------------------------------------------------
!
!*       2.1    Downdraft evaporation rate at the DBL. The evaporation
!               rate in downdraft must be consistent with precipitation
!               efficiency relation.
!               --------------------------------------------------------
!
!
DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI) = PDTEVR(JI) / MIN( -1.E-1, PDMF(JI,JK) )
     ZWORK6(JI) = PDMF(JI,JK)
END DO
!
!*       2.2    Some preliminar computations for downdraft = total 
!               precipitation rate. The precipitation is evaluated in 
!               a layer thickness DP=XUSRDPTH=165 hPa above the LCL.
!               The difference between updraft precipitation and downdraft
!               precipitation (updraft supply rate) is used to drive the
!               downdraft through evaporational cooling.
!               --------------------------------------------------------
!
DO JI = 1, IIE
     JK = KLCL(JI)
     ZWORK5(JI) = PPRES(JI,JK)
END DO
!
PTPR(:) = 0.
DO JK = JKT3, JKT2
    WHERE ( JK >= KLCL(:) .AND. PPRES(:,JK) >= ZWORK5(:) - XUSRDPTH )
	PTPR(:) = PTPR(:) + PUPR(:,JK)
	IPRL(:) = JK
    END WHERE
END DO
IPRL(:) = MIN( KETL(:), IPRL(:) )
!
DO JI = 1, IIE
     JK = IPRL(JI)
     PTPR(JI) = PUMF(JI,JK+1) * PURW(JI,JK+1) + PTPR(JI) 
END DO
!
PTPR(:) = PPREF(:) * MIN( PUTPR(:), PTPR(:) )
ZWORK4(:) = PUTPR(:) - PTPR(:) 
!
!
!*       2.3    Total amount of precipitation that falls out of the up-
!               draft between the LCL and the LFS.
!               Condensate transfer from up to downdraft at LFS
!               ---------------------------------------------------------
!
ZWORK5(:) = 0.
DO JK = JKT3, JKT1
     WHERE ( JK >= KLCL(:) .AND. JK <= KLFS(:) )
	   ZWORK5(:) = ZWORK5(:) +  PUPR(:,JK)
     END WHERE
END DO
!
DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK2(JI) = ( 1. - PPREF(JI) ) * ZWORK5(JI) *                     &
                  ( 1. - PMIXF(JI) ) / MAX( 1.E-1, PUMF(JI,JK) )
END DO
!
!
!*       2.4    Increase the first guess downdraft mass flux to satisfy
!               precipitation efficiency relation.
!               If downdraft does not evaporate any water at the DBL for  
!               the specified relative humidity, or if the corrected mass 
!               flux at the LFS is positive no downdraft is allowed
!               ---------------------------------------------------------
!    
!
ZWORK1(:) = ZWORK4(:) / ( ZWORK1(:) + ZWORK2(:) + 1.E-8 ) 
ZWORK2(:) = ZWORK1(:) / MIN( -1.E-1, ZWORK6(:) ) ! ratio of budget consistent to actual DMF
!
ZWORK3(:) = 1.
ZWORK6(:) = 1.
WHERE ( ZWORK1(:) > 0. .OR. PDTEVR(:) < 1. ) 
   KDBL(:)   = IKB
   KLFS(:)   = IKB
   PDTEVR(:) = 0. 
   ZWORK2(:) = 0.
   ZWORK3(:) = 0.
   ZWORK6(:) = 0.
END WHERE
!
DO JK = IKB, JKT1   
     PDMF(:,JK)  = PDMF(:,JK)  * ZWORK2(:)
     PDER(:,JK)  = PDER(:,JK)  * ZWORK2(:)  
     PDDR(:,JK)  = PDDR(:,JK)  * ZWORK2(:)  
   PDTEVRF(:,JK) = PDTEVRF(:,JK)* ZWORK2(:)  
     PDRW(:,JK)  = PDRW(:,JK)  * ZWORK3(:)  
     PDTHL(:,JK) = PDTHL(:,JK) * ZWORK3(:)  
END DO     
ZWORK4(:) = ZWORK2(:)
!
!
!*       3.     Increase updraft mass flux, mass detrainment rate, and water  
!               substance detrainment rates to be consistent with the transfer
!               of the estimated mass from the up- to the downdraft at the LFS
!               --------------------------------------------------------------
!
DO JI = 1, IIE
    JK = KLFS(JI)
    ZWORK2(JI) = ( 1. - ZWORK6(JI) ) + ZWORK6(JI) *                   &
		  ( PUMF(JI,JK) - ( 1. - PMIXF(JI) ) * ZWORK1(JI) ) / &
		  MAX( 1.E-1, PUMF(JI,JK) )
END DO
!
!
JKT1  = MAXVAL( KLFS(:) )  ! value of KLFS might have been reset to IKB above
DO JK = IKB, JKT1
    DO JI = 1, IIE
      IF ( JK <= KLFS(JI) ) THEN
	PUMF(JI,JK)  = PUMF(JI,JK)  * ZWORK2(JI) 
	PUER(JI,JK)  = PUER(JI,JK)  * ZWORK2(JI)
	PUDR(JI,JK)  = PUDR(JI,JK)  * ZWORK2(JI)
	PUPR(JI,JK)  = PUPR(JI,JK)  * ZWORK2(JI)
      END IF
    END DO
END DO
!
!
!*       4.     Increase total = downdraft precipitation and evaporation rate
!               -------------------------------------------------------------
!
WHERE ( PDTEVR(:) > 0. )
    PDTEVR(:)  = PDTEVR(:) * ZWORK4(:)
    PTPR(:)    = PTPR(:) + PPREF(:) * ZWORK5(:) * ( ZWORK2(:) - 1. )
ELSEWHERE
    PTPR(:)    = PUTPR(:)
END WHERE
!
!
END SUBROUTINE CONVECT_PRECIP_ADJUST
!     ######spl
     SUBROUTINE CONVECT_CLOSURE( KLON, KLEV,                                 &
                                 PPRES, PDPRES, PZ, PDXDY, PLMASS,           &
                                 PTHL, PTH, PRW, PRC, PRI, OTRIG1,           &
                                 PTHC, PRWC, PRCC, PRIC, PWSUB,              &
                                 KLCL, KDPL, KPBL, KLFS, KCTL, KML,          &
                                 PUMF, PUER, PUDR, PUTHL, PURW,              &
                                 PURC, PURI, PUPR,                           &
                                 PDMF, PDER, PDDR, PDTHL, PDRW,              &
                                 PTPR, PSPR, PDTEVR,                         &
                                 PCAPE, PTIMEC,                              &
                                 KFTSTEPS,                                   &
                                 PDTEVRF, PPRLFLX, PPRSFLX                   )
!    #######################################################################
!
!!**** Uses modified Fritsch-Chappell closure
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!     (over a time step PTIMEC) environmental values of THETA_l, R_w, R_c, R_i
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PTHC-PTH)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     
!!    CONVECT_CLOSURE_THRVLCL
!!    CONVECT_CLOSURE_ADJUST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV         ! specific heat for dry air and water vapor
!!          XCL, XCI           ! specific heat for liquid water and ice
!!          XTT                ! triple point temperature
!!          XLVTT, XLSTT       ! vaporization, sublimation heat constant
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration 
!!          XSTABC             ! stability factor in CAPE adjustment
!!          XMELDPTH           ! allow melting over specific pressure depth
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Peter Bechtold 04/10/97 change for enthalpie, r_c + r_i tendencies
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                   INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,                   INTENT(IN) :: KLEV   ! vertical dimension
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KML    ! index for melting level
REAL, DIMENSION(KLON),  INTENT(INOUT) :: PTIMEC ! convection time step 
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY  ! grid area (m^2)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH    ! grid scale theta        
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRW    ! grid scale total water  
			                        ! mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRC    ! grid scale r_c 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRI    ! grid scale r_i 
LOGICAL, DIMENSION(KLON),  INTENT(IN) :: OTRIG1 ! logical to keep trace of 
                                                ! convective arrays modified in UPDRAFT
!   
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES ! pressure difference between 
                                                 ! bottom and top of layer (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PLMASS ! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN)  :: PCAPE  ! available potent. energy
INTEGER,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                 ! only used for chemical tracers
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUPR  ! updraft precipitation in
                                                  ! flux units (kg water / s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF  ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDER  ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDDR  ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDTHL ! downdraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)   :: PDRW  ! downdraft total water (kg/kg)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PTPR  ! total surf precipitation (kg/s)
REAL, DIMENSION(KLON),      INTENT(OUT)  :: PSPR  ! solid surf precipitation (kg/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PDTEVR! donwndraft evapor. (kg/s)
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDTEVRF! downdraft evaporation rate
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRSFLX! solid  precip flux
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JK, JKP, JKMAX ! vertical loop index
INTEGER :: JI             ! horizontal loop index
INTEGER :: JITER          ! iteration loop index
INTEGER :: JSTEP          ! fractional time loop index
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
REAL    :: ZCVOCD, ZEPSA  ! C_pv / C_pd, R_v / R_d
!
REAL, DIMENSION(KLON,KLEV) :: ZTHLC       ! convectively adjusted 
                                          ! grid scale enthalpy
REAL, DIMENSION(KLON,KLEV) :: ZOMG        ! conv. environm. subsidence (Pa/s)
REAL, DIMENSION(KLON,KLEV) :: ZUMF        ! non-adjusted updraft mass flux
REAL, DIMENSION(KLON,KLEV) :: ZUER        !   "     updraft  entrainm. rate
REAL, DIMENSION(KLON,KLEV) :: ZUDR        !   "     updraft  detrainm. rate
REAL, DIMENSION(KLON,KLEV) :: ZDMF        !   "   downdraft mass flux
REAL, DIMENSION(KLON,KLEV) :: ZDER        !   "   downdraft  entrainm. rate
REAL, DIMENSION(KLON,KLEV) :: ZDDR        !   "   downdraft  detrainm. rate
REAL, DIMENSION(KLON)     :: ZTPR         !   "   total precipitation
REAL, DIMENSION(KLON)     :: ZDTEVR       !   "   total downdraft evapor. 
REAL, DIMENSION(KLON,KLEV):: ZPRLFLX      !   "   liquid precip flux
REAL, DIMENSION(KLON,KLEV):: ZPRSFLX      !   "   solid  precip flux
REAL, DIMENSION(KLON)     :: ZPRMELT      ! melting of precipitation
REAL, DIMENSION(KLON)     :: ZPRMELTO     ! non-adjusted  "
REAL, DIMENSION(KLON)     :: ZADJ         ! mass adjustment factor
REAL, DIMENSION(KLON)     :: ZADJMAX      ! limit value for ZADJ
REAL, DIMENSION(KLON)     :: ZCAPE        ! new CAPE after adjustment
REAL, DIMENSION(KLON)     :: ZTIMEC       ! fractional convective time step
REAL, DIMENSION(KLON,KLEV):: ZTIMC        ! 2D work array for ZTIMEC
!
REAL, DIMENSION(KLON)     :: ZTHLCL       ! new  theta at LCL
REAL, DIMENSION(KLON)     :: ZRVLCL       ! new  r_v at LCL
REAL, DIMENSION(KLON)     :: ZZLCL        ! height of LCL
REAL, DIMENSION(KLON)     :: ZTLCL        ! temperature at LCL
REAL, DIMENSION(KLON)     :: ZTELCL       ! envir. temper. at LCL
REAL, DIMENSION(KLON)     :: ZTHEUL       ! theta_e for undilute ascent
REAL, DIMENSION(KLON)     :: ZTHES1, ZTHES2! saturation environm. theta_e
REAL, DIMENSION(KLON,KLEV) :: ZTHMFIN, ZTHMFOUT, ZRWMFIN, ZRWMFOUT
REAL, DIMENSION(KLON,KLEV) :: ZRCMFIN, ZRCMFOUT, ZRIMFIN, ZRIMFOUT
                                    ! work arrays for environm. compensat. mass flux
REAL, DIMENSION(KLON)     :: ZPI          ! (P/P00)**R_d/C_pd 
REAL, DIMENSION(KLON)     :: ZLV          ! latent heat of vaporisation
REAL, DIMENSION(KLON)     :: ZLS          ! latent heat of sublimation 
REAL, DIMENSION(KLON)     :: ZLM          ! latent heat of melting
REAL, DIMENSION(KLON)     :: ZCPH         ! specific heat C_ph
REAL, DIMENSION(KLON)     :: ZMELDPTH     ! actual depth of melting layer 
INTEGER, DIMENSION(KLON)  :: ITSTEP       ! fractional convective time step
INTEGER, DIMENSION(KLON)  :: ICOUNT       ! timestep counter 
INTEGER, DIMENSION(KLON)  :: ILCL         ! index lifting condens. level
INTEGER, DIMENSION(KLON)  :: IWORK1       ! work array
REAL, DIMENSION(KLON)     :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5
REAL, DIMENSION(KLON,KLEV):: ZWORK6
LOGICAL, DIMENSION(KLON)  :: GWORK1, GWORK3! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK4    ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.2    Initialize  local variables
!               ----------------------------
!
!
PSPR(:)   = 0.
ZTIMC(:,:) = 0.
ZTHES2(:) = 0.
ZWORK1(:) = 0. 
ZWORK2(:) = 0. 
ZWORK3(:) = 0. 
ZWORK4(:) = 0. 
ZWORK5(:) = 0. 
GWORK1(:) = .FALSE.
GWORK3(:) = .FALSE.  
GWORK4(:,:) = .FALSE.  
ILCL(:)   = KLCL(:)
!
ZCPORD    = XCPD / XRD
ZRDOCP    = XRD / XCPD
ZCVOCD    = XCPV / XCPD 
ZEPSA     = XRV / XRD
!
ZADJ(:)   = 1.
ZWORK5(:) = 1.
WHERE( .NOT. OTRIG1(:) ) ZWORK5(:) = 0. 
!
!
!*       0.3   Compute loop bounds
!              ------------------- 
!
IIE    = KLON
IKB    = 1 + JCVEXB 
IKS    = KLEV
IKE    = KLEV - JCVEXT 
JKMAX  = MAXVAL( KCTL(:) )
!
!
!*       2.     Save initial mass flux values to be used in adjustment procedure
!               ---------------------------------------------------------------
!
ZUMF(:,:)  = PUMF(:,:)
ZUER(:,:)  = PUER(:,:)
ZUDR(:,:)  = PUDR(:,:)
ZDMF(:,:)  = PDMF(:,:)
ZDER(:,:)  = PDER(:,:)
ZDDR(:,:)  = PDDR(:,:)
ZTPR(:)    = PTPR(:)
ZDTEVR(:)  = PDTEVR(:)
ZOMG(:,:)  = 0.
PWSUB(:,:) = 0. 
ZPRMELT(:) = 0.
PPRLFLX(:,:) = 0.
ZPRLFLX(:,:) = 0.
PPRSFLX(:,:) = 0.
ZPRSFLX(:,:) = 0.
!
!
!*       2.1    Some preliminar computations for melting of precipitation
!               used later in section 9 and computation of precip fluxes
!               Precipitation fluxes are updated for melting and evaporation
!               ---------------------------------------------------------
!
!
ZWORK1(:) = 0.
ZMELDPTH(:) = 0.
ZWORK6(:,:) = 0.
DO JK = JKMAX + 1, IKB + 1, -1
   ! Nota: PUPR is total precipitation flux, but the solid, liquid
   !       precipitation is stored in units kg/kg; therefore we compute here
   !       the solid fraction of the total precipitation flux.
  DO JI = 1, IIE
     ZWORK2(JI)    = PUPR(JI,JK) / ( PURC(JI,JK) + PURI(JI,JK) + 1.E-8 )
     ZPRMELT(JI)   = ZPRMELT(JI) + PURI(JI,JK) * ZWORK2(JI)
     ZWORK1(JI)    = ZWORK1(JI) + PURC(JI,JK) * ZWORK2(JI) - PDTEVRF(JI,JK)
     ZPRLFLX(JI,JK)= MAX( 0., ZWORK1(JI) )
     ZPRMELT(JI)   = ZPRMELT(JI) + MIN( 0., ZWORK1(JI) )
     ZPRSFLX(JI,JK)= ZPRMELT(JI) 
     IF ( KML(JI) >= JK .AND. ZMELDPTH(JI) <= XMELDPTH ) THEN                 
          ZPI(JI)    = ( PPRES(JI,JK) / XP00 ) ** ZRDOCP 
          ZWORK3(JI) = PTH(JI,JK) * ZPI(JI)            ! temperature estimate
          ZLM(JI)    = XLSTT + ( XCPV - XCI ) * ( ZWORK3(JI) - XTT ) -       &
               ( XLVTT + ( XCPV - XCL ) * ( ZWORK3(JI) - XTT ) ) ! L_s - L_v
          ZCPH(JI)   = XCPD + XCPV * PRW(JI,JK)
          ZMELDPTH(JI) = ZMELDPTH(JI) + PDPRES(JI,JK)
          ZWORK6(JI,JK)= ZLM(JI) * PTIMEC(JI) / PLMASS(JI,JK) * PDPRES(JI,JK)
          ZOMG(JI,JK)= 1. ! at this place only used as work variable
     END IF
  END DO
!
END DO
!
ZWORK2(:) = 0.
DO JK = JKMAX, IKB + 1, -1
    ZWORK1(:) = ZPRMELT(:) * PDPRES(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
    ZWORK2(:) = ZWORK2(:) + ZWORK1(:) * ZOMG(:,JK)
    ZPRLFLX(:,JK) = ZPRLFLX(:,JK) + ZWORK2(:) 
    ZPRSFLX(:,JK) = ZPRSFLX(:,JK) - ZWORK2(:)
END DO 
WHERE( ZPRSFLX(:,:) < 1. ) ZPRSFLX(:,:)=0.
ZPRMELTO(:) = ZPRMELT(:)
!
!
!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger 
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------
!
ZADJMAX(:) = 1000.
IWORK1(:) = MAX( ILCL(:), KLFS(:) )
JKP = MINVAL( KDPL(:) )
DO JK = JKP, IKE
  DO JI = 1, IIE
    IF( JK > KDPL(JI) .AND. JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) /                                      &
                  ( ( PUER(JI,JK) + PDER(JI,JK) + 1.E-5 ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    END IF
  END DO
END DO
!
!
GWORK1(:) = OTRIG1(:)  ! logical array to limit adjustment to not definitively
                       ! adjusted columns
!
DO JK = IKB, IKE
  ZTHLC(:,:) = PTHL(:,:) ! initialize adjusted envir. values 
  PRWC(:,:)  = PRW(:,:)
  PRCC(:,:)  = PRC(:,:)
  PRIC(:,:)  = PRI(:,:)
  PTHC(:,:)  = PTH(:,:)
END DO
!
!
!
DO JITER = 1, 7  ! Enter adjustment loop to assure that all CAPE is
                 ! removed within the advective time interval TIMEC
!
     ZTIMEC(:) = PTIMEC(:)
     GWORK4(:,:)   = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKS )
     WHERE( GWORK4(:,:) ) PWSUB(:,:) = 0.
     ZOMG(:,:)=0.
!
     DO JK = IKB + 1, JKMAX
           JKP = MAX( IKB + 1, JK - 1 )
           WHERE ( GWORK1(:) .AND. JK <= KCTL(:) )
!
!
!*       4.     Determine vertical velocity at top and bottom of each layer
!               to satisfy mass continuity.
!               ---------------------------------------------------------------
              ! we compute here Domega/Dp = - g rho Dw/Dz = 1/Dt
!
             ZWORK1(:)   = - ( PUER(:,JKP) + PDER(:,JKP) -                   &
                           PUDR(:,JKP) - PDDR(:,JKP) ) / PLMASS(:,JKP)
!    
             PWSUB(:,JK) = PWSUB(:,JKP) - PDPRES(:,JK-1) * ZWORK1(:)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer
!
!   
!*       5.     Compute fractional time step. For stability or 
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------
!
             ZWORK1(:) = XSTABT * PDPRES(:,JKP) / ( ABS( PWSUB(:,JK) ) + 1.E-10 )
              ! the factor XSTABT is used for stability reasons
             ZTIMEC(:) = MIN( ZTIMEC(:), ZWORK1(:) ) 
!
              ! transform vertical velocity in mass flux units
             ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / XG 
         END WHERE
     END DO
!
!
     WHERE( GWORK4(:,:) )
           ZTHLC(:,:) = PTHL(:,:) ! reinitialize adjusted envir. values 
           PRWC(:,:)  = PRW(:,:)  ! when iteration criterium not attained
           PRCC(:,:)  = PRC(:,:)
           PRIC(:,:)  = PRI(:,:)
           PTHC(:,:)  = PTH(:,:)
     END WHERE
!
! 
!        6. Check for mass conservation, i.e. ZWORK1 > 1.E-2
!           If mass is not conserved, the convective tendencies
!           automatically become zero.
!           ----------------------------------------------------
!
    DO JI = 1, IIE
       JK=KCTL(JI)
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1 )    &
                                                            - PWSUB(JI,JK)
    END DO
    WHERE( GWORK1(:) .AND. ABS( ZWORK1(:) ) - .01 > 0. )
        GWORK1(:) = .FALSE.
        PTIMEC(:) = 1.E-1
        ZTPR(:)   = 0.
        ZWORK5(:) = 0.
    END WHERE
    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
        ZPRLFLX(:,JK) = ZPRLFLX(:,JK) * ZWORK5(:)
        ZPRSFLX(:,JK) = ZPRSFLX(:,JK) * ZWORK5(:)
    END DO
    GWORK4(:,1:IKB) = .FALSE.
    GWORK4(:,IKS)   = .FALSE.
!
    ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1 
    ZTIMEC(:) = PTIMEC(:) / REAL( ITSTEP(:) ) ! adjust  fractional time step
                                           ! to be an integer multiple of PTIMEC
    ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
    ICOUNT(:) = 0
!
!
!
    KFTSTEPS = MAXVAL( ITSTEP(:) )
    DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop here
!
     	 ICOUNT(:) = ICOUNT(:) + 1
!
	     GWORK3(:) =  ITSTEP(:) >= ICOUNT(:) .AND. GWORK1(:) 
!
!
!*       7.     Assign enthalpy and r_w values at the top and bottom of each
!               layer based on the sign of w
!               ------------------------------------------------------------
!
             ZTHMFIN(:,:)   = 0.
             ZRWMFIN(:,:)   = 0.
             ZRCMFIN(:,:)   = 0.
             ZRIMFIN(:,:)   = 0.
             ZTHMFOUT(:,:)  = 0.
             ZRWMFOUT(:,:)  = 0.
             ZRCMFOUT(:,:)  = 0.
             ZRIMFOUT(:,:)  = 0.
!
         DO JK = IKB + 1, JKMAX
         GWORK4(:,JK) = GWORK3(:) .AND. JK <= KCTL(:)
         JKP = MAX( IKB + 1, JK - 1 )
           DO JI = 1, IIE
           IF ( GWORK3(JI) ) THEN
!
               ZWORK1(JI)       = SIGN( 1., ZOMG(JI,JK) )
               ZWORK2(JI)       = 0.5 * ( 1. + ZWORK1(JI) )
               ZWORK1(JI)       = 0.5 * ( 1. - ZWORK1(JI) )
               ZTHMFIN(JI,JK)   = - ZOMG(JI,JK) * ZTHLC(JI,JKP) * ZWORK1(JI)
               ZTHMFOUT(JI,JK)  =   ZOMG(JI,JK) * ZTHLC(JI,JK)  * ZWORK2(JI)
               ZTHMFIN(JI,JKP)  = ZTHMFIN(JI,JKP)  + ZTHMFOUT(JI,JK) * ZWORK2(JI)
               ZTHMFOUT(JI,JKP) = ZTHMFOUT(JI,JKP) + ZTHMFIN(JI,JK)  * ZWORK1(JI)
               ZRWMFIN(JI,JK)   = - ZOMG(JI,JK) * PRWC(JI,JKP) * ZWORK1(JI)
               ZRWMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRWC(JI,JK)  * ZWORK2(JI)
               ZRWMFIN(JI,JKP)  = ZRWMFIN(JI,JKP)  + ZRWMFOUT(JI,JK) * ZWORK2(JI)
               ZRWMFOUT(JI,JKP) = ZRWMFOUT(JI,JKP) + ZRWMFIN(JI,JK)  * ZWORK1(JI)
               ZRCMFIN(JI,JK)   = - ZOMG(JI,JK) * PRCC(JI,JKP) * ZWORK1(JI)
               ZRCMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRCC(JI,JK)  * ZWORK2(JI)
               ZRCMFIN(JI,JKP)  = ZRCMFIN(JI,JKP)  + ZRCMFOUT(JI,JK) * ZWORK2(JI)
               ZRCMFOUT(JI,JKP) = ZRCMFOUT(JI,JKP) + ZRCMFIN(JI,JK)  * ZWORK1(JI)
               ZRIMFIN(JI,JK)   = - ZOMG(JI,JK) * PRIC(JI,JKP) * ZWORK1(JI)
               ZRIMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRIC(JI,JK)  * ZWORK2(JI)
               ZRIMFIN(JI,JKP)  = ZRIMFIN(JI,JKP)  + ZRIMFOUT(JI,JK) * ZWORK2(JI)
               ZRIMFOUT(JI,JKP) = ZRIMFOUT(JI,JKP) + ZRIMFIN(JI,JK)  * ZWORK1(JI)
!
           END IF
           END DO
         END DO
!
         WHERE ( GWORK4(:,:) )
!
!******************************************************************************
!
!*       8.     Update the environmental values of enthalpy and r_w at each level
!               NOTA: These are the MAIN EQUATIONS of the scheme
!               -----------------------------------------------------------------
!
!
           ZTHLC(:,:) = ZTHLC(:,:) + ZTIMC(:,:) / PLMASS(:,:) * (      &
                          ZTHMFIN(:,:) + PUDR(:,:) * PUTHL(:,:)  +     &
                          PDDR(:,:) * PDTHL(:,:) - ZTHMFOUT(:,:) -     &
                        ( PUER(:,:) + PDER(:,:) ) * PTHL(:,:)   )
           PRWC(:,:)  = PRWC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
                         ZRWMFIN(:,:) + PUDR(:,:) * PURW(:,:)  +       &
                         PDDR(:,:) * PDRW(:,:) - ZRWMFOUT(:,:) -       &
                        ( PUER(:,:) + PDER(:,:) ) * PRW(:,:)    )    
           PRCC(:,:)  = PRCC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
               ZRCMFIN(:,:) + PUDR(:,:) * PURC(:,:) - ZRCMFOUT(:,:) -  &
                        ( PUER(:,:) + PDER(:,:) ) * PRC(:,:)    )    
           PRIC(:,:)  = PRIC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
               ZRIMFIN(:,:) + PUDR(:,:) * PURI(:,:) - ZRIMFOUT(:,:) -  & 
                        ( PUER(:,:) + PDER(:,:) ) * PRI(:,:)    )    
!
!
!******************************************************************************
!
         END WHERE
!
    END DO ! Exit the fractional time step loop
!
! 
!*           9.    Allow frozen precipitation to melt over a 200 mb deep layer
!                  -----------------------------------------------------------
!
      DO JK = JKMAX, IKB + 1, -1
            ZTHLC(:,JK) = ZTHLC(:,JK) -                                &
               ZPRMELT(:) * ZWORK6(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
      END DO
!
!
!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------
!
      DO JK = IKB + 1, JKMAX
         DO JI = 1, IIE
         IF( GWORK1(JI) .AND. JK <= KCTL(JI) ) THEN
           ZPI(JI)    = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
           ZCPH(JI)   = XCPD + PRWC(JI,JK) * XCPV
           ZWORK2(JI) = PTH(JI,JK) / ZPI(JI)  ! first temperature estimate
           ZLV(JI)    = XLVTT + ( XCPV - XCL ) * ( ZWORK2(JI) - XTT )
           ZLS(JI)    = XLVTT + ( XCPV - XCI ) * ( ZWORK2(JI) - XTT )
             ! final linearized temperature
           ZWORK2(JI) = ( ZTHLC(JI,JK) + ZLV(JI) * PRCC(JI,JK) + ZLS(JI) * PRIC(JI,JK) &
                       - (1. + PRWC(JI,JK) ) * XG * PZ(JI,JK) ) / ZCPH(JI)
           ZWORK2(JI) = MAX( 180., MIN( 340., ZWORK2(JI) ) )
           PTHC(JI,JK)= ZWORK2(JI) * ZPI(JI) ! final adjusted envir. theta
         END IF
         END DO
      END DO
!
!
!*         11.     Compute new cloud ( properties at new LCL )
!                     NOTA: The computations are very close to
!                           that in routine TRIGGER_FUNCT
!                  ---------------------------------------------
!
      CALL CONVECT_CLOSURE_THRVLCL(  KLON, KLEV,                           &
                                     PPRES, PTHC, PRWC, PZ, GWORK1,        &
                                     ZTHLCL, ZRVLCL, ZZLCL, ZTLCL, ZTELCL, &
                                     ILCL, KDPL, KPBL )
!
!
       ZTLCL(:)  = MAX( 230., MIN( 335., ZTLCL(:) ) )  ! set some overflow bounds
       ZTELCL(:) = MAX( 230., MIN( 335., ZTELCL(:) ) )
       ZTHLCL(:) = MAX( 230., MIN( 345., ZTHLCL(:) ) )
       ZRVLCL(:) = MAX(   0., MIN(   1., ZRVLCL(:) ) )
!
!
!*         12.    Compute adjusted CAPE
!                 ---------------------
!
       ZCAPE(:)  = 0.
       ZPI(:)    = ZTHLCL(:) / ZTLCL(:)
       ZPI(:)    = MAX( 0.95, MIN( 1.5, ZPI(:) ) )
       ZWORK1(:) = XP00 / ZPI(:) ** ZCPORD ! pressure at LCL
!
       CALL CONVECT_SATMIXRATIO( KLON, ZWORK1, ZTELCL, ZWORK3, ZLV, ZLS, ZCPH )
       ZWORK3(:) = MIN(   .1, MAX(   0., ZWORK3(:) ) )
!
	        ! compute theta_e updraft undilute
       ZTHEUL(:) = ZTLCL(:) * ZPI(:) ** ( 1. - 0.28 * ZRVLCL(:) )            &
                                  * EXP( ( 3374.6525 / ZTLCL(:) - 2.5403 )   &
                                  * ZRVLCL(:) * ( 1. + 0.81 * ZRVLCL(:) ) )
!
	        ! compute theta_e saturated environment at LCL
       ZTHES1(:) = ZTELCL(:) * ZPI(:) ** ( 1. - 0.28 * ZWORK3(:) )           &
                                  * EXP( ( 3374.6525 / ZTELCL(:) - 2.5403 )  &
                                  * ZWORK3(:) * ( 1. + 0.81 * ZWORK3(:) ) )
!
      DO JK = MINVAL( ILCL(:) ), JKMAX
        JKP = JK - 1
        DO JI = 1, IIE
          ZWORK4(JI) = 1.
          IF ( JK == ILCL(JI) ) ZWORK4(JI) = 0.
!
           ! compute theta_e saturated environment and adjusted values
           ! of theta
!
          GWORK3(JI)  = JK >= ILCL(JI) .AND. JK <= KCTL(JI) .AND. GWORK1(JI) 
!
          ZPI(JI)     = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
          ZWORK2(JI)  = PTHC(JI,JK) / ZPI(JI)
        END DO
!
        CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZWORK2, ZWORK3, ZLV, ZLS, ZCPH )
!
!
        DO JI = 1, IIE
          IF ( GWORK3(JI) ) THEN
              ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )   &
                                   * EXP( ( 3374.6525 / ZWORK2(JI) - 2.5403 ) &
                                   * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
!
              ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                &
                           ( 1. - ZWORK4(JI) ) * ZZLCL(JI)    ! level thickness
              ZWORK1(JI)  = ( 2. * ZTHEUL(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - 1.
              ZCAPE(JI)   = ZCAPE(JI) + XG * ZWORK3(JI) * MAX( 0., ZWORK1(JI) )
              ZTHES1(JI)  = ZTHES2(JI)
          END IF
        END DO
      END DO
!
!                                                          
!*         13.     Determine mass adjustment factor knowing how much
!                  CAPE has been removed.
!                  -------------------------------------------------
!
       WHERE ( GWORK1(:) )
           ZWORK1(:) = MAX( PCAPE(:) - ZCAPE(:), 0.1 * PCAPE(:) )
           ZWORK2(:) = ZCAPE(:) / ( PCAPE(:) + 1.E-8 )
!       
           GWORK1(:) = ZWORK2(:) > 0.1 .OR. ZCAPE(:) == 0. ! mask for adjustment
       END WHERE
!
       WHERE ( ZCAPE(:) == 0. .AND. GWORK1(:) )  ZADJ(:) = ZADJ(:) * 0.5
       WHERE ( ZCAPE(:) /= 0. .AND. GWORK1(:) )                              &
               ZADJ(:) = ZADJ(:) * XSTABC * PCAPE(:) / ( ZWORK1(:) + 1.E-8 )
       ZADJ(:) = MIN( ZADJ(:), ZADJMAX(:) )  
!
!
!*         13.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization
!                 ----------------------------------------------------
!
       CALL CONVECT_CLOSURE_ADJUST( KLON, KLEV, ZADJ,                     &
                                    PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR,   &
                                    PDMF, ZDMF, PDER, ZDER, PDDR, ZDDR,   &
                                    ZPRMELT, ZPRMELTO, PDTEVR, ZDTEVR,    &
                                    PTPR, ZTPR,                           &
                                    PPRLFLX, ZPRLFLX, PPRSFLX, ZPRSFLX    )
!
!
      IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached 
                                          ! desired degree of stabilization.
!
END DO  ! end of big adjustment iteration loop
!
!
        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
  PRWC(:,JK) = MAX( 0., PRWC(:,JK) - PRCC(:,JK) - PRIC(:,JK) )
END DO
!
        ! compute surface solid (ice) precipitation 
PSPR(:) = ZPRMELT(:) * ( 1. - ZMELDPTH(:) / XMELDPTH )
PSPR(:) = MAX( 0., PSPR(:) )
!
!
END SUBROUTINE CONVECT_CLOSURE
!     ######spl
     SUBROUTINE CONVECT_CLOSURE_ADJUST( KLON, KLEV, PADJ,                      &
                                        PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR, &
                                        PDMF, PZDMF, PDER, PZDER, PDDR, PZDDR, &
                                        PPRMELT, PZPRMELT, PDTEVR, PZDTEVR,    &
                                        PTPR, PZTPR,                           &
                                        PPRLFLX, PZPRLFL, PPRSFLX, PZPRSFL     )

!    #########################################################################
!
!!**** Uses closure adjustment factor to adjust mass flux and to modify
!!     precipitation efficiency  when necessary. The computations are
!!     similar to routine CONVECT_PRECIP_ADJUST.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust the mass flux using the
!!      factor PADJ computed in CONVECT_CLOSURE
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!
!!    EXTERNAL
!!    --------
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!     
!!    None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    None
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
REAL, DIMENSION(KLON),      INTENT(IN) :: PADJ     ! mass adjustment factor
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF  ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER  ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR  ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDDR ! initial value of  "
REAL, DIMENSION(KLON),   INTENT(INOUT):: PTPR     ! total precipitation (kg/s)
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZTPR    ! initial value of "
REAL, DIMENSION(KLON),   INTENT(INOUT):: PDTEVR   ! donwndraft evapor. (kg/s)
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZDTEVR  ! initial value of " 
REAL, DIMENSION(KLON),   INTENT(INOUT):: PPRMELT  ! melting of precipitation
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZPRMELT ! initial value of " 
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRLFL! initial value "
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRSFLX! solid  precip flux
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRSFL! initial value "
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE                 ! horiz. + vert. loop bounds
INTEGER :: JK                            ! vertical loop index
!
!
!-------------------------------------------------------------------------------
!
!*       0.3   Compute loop bounds
!              -------------------
!
IIE  = KLON
IKB  = 1 + JCVEXB 
IKE  = KLEV - JCVEXT 
!
!
!*       1.     Adjust mass flux by the factor PADJ to converge to
!               specified degree of stabilization
!               ----------------------------------------------------
!
          PPRMELT(:)  = PZPRMELT(:)   * PADJ(:)
          PDTEVR(:)   = PZDTEVR(:)    * PADJ(:)
          PTPR(:)     = PZTPR(:)      * PADJ(:)
!
     DO JK = IKB + 1, IKE
	  PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
          PDMF(:,JK)  = PZDMF(:,JK)   * PADJ(:)
          PDER(:,JK)  = PZDER(:,JK)   * PADJ(:)
          PDDR(:,JK)  = PZDDR(:,JK)   * PADJ(:)
          PPRLFLX(:,JK) = PZPRLFL(:,JK) * PADJ(:)
          PPRSFLX(:,JK) = PZPRSFL(:,JK) * PADJ(:)
     END DO
!
END SUBROUTINE CONVECT_CLOSURE_ADJUST
!     ######spl
      SUBROUTINE CONVECT_CLOSURE_THRVLCL( KLON, KLEV,                         &
                                          PPRES, PTH, PRV, PZ, OWORK1,        &
                                         PTHLCL, PRVLCL, PZLCL, PTLCL, PTELCL,&
                                          KLCL, KDPL, KPBL )
!     ######################################################################
!
!!**** Determine thermodynamic properties at new LCL
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the thermodynamic
!!      properties at the new lifting condensation level LCL
!!   
!!
!!
!!**  METHOD
!!    ------
!!    see CONVECT_TRIGGER_FUNCT
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! Reference pressure
!!          XRD, XRV           ! Gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XTT                ! triple point temperature
!!          XBETAW, XGAMW      ! constants for vapor saturation pressure
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! lowest allowed pressure difference between
!!                             ! surface and LCL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! theta
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRV   ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of grid point (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KDPL  ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KPBL  ! " vert. index of source layer top
LOGICAL, DIMENSION(KLON),   INTENT(IN) :: OWORK1! logical mask 
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL  ! temperature at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTELCL ! environm. temp. at LCL (K)
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLCL   ! contains vert. index of LCL
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JK, JKM, JKMIN, JKMAX      ! vertical loop index
INTEGER :: JI                         ! horizontal loop index 
INTEGER :: IIE, IKB, IKE              ! horizontal + vertical loop bounds
REAL    :: ZEPS, ZEPSA    ! R_d / R_v, R_v / R_d 
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL, DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL, DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure 
REAL, DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL, DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL, DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(KLON) :: ZWORK1, ZWORK2     ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Initialize local variables
!               --------------------------
!
ZEPS      = XRD / XRV
ZEPSA     = XRV / XRD 
ZCPORD    = XCPD / XRD
ZRDOCP    = XRD / XCPD
!
ZDPTHMIX(:) = 0.
ZPRESMIX(:) = 0.
PTHLCL(:)   = 300.
PTLCL(:)    = 300.
PTELCL(:)   = 300.
PRVLCL(:)   = 0.
PZLCL(:)    = PZ(:,IKB)
ZTMIX(:)    = 230.
ZPLCL(:)    = 1.E4 
KLCL(:)     = IKB + 1
!
!
!*       2.     Construct a mixed layer as in TRIGGER_FUNCT
!               -------------------------------------------
!
     JKMAX = MAXVAL( KPBL(:) )
     JKMIN = MINVAL( KDPL(:) )
     DO JK = IKB + 1, JKMAX
        JKM = JK + 1
        DO JI = 1, IIE
        IF ( JK >= KDPL(JI) .AND. JK <= KPBL(JI) ) THEN
!           
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            PTHLCL(JI)   = PTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            PRVLCL(JI)   = PRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
!
        END IF
        END DO
     END DO
!
!
WHERE ( OWORK1(:) )
!
        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
        PTHLCL(:)   = PTHLCL(:)   / ZDPTHMIX(:)
        PRVLCL(:)   = PRVLCL(:)   / ZDPTHMIX(:)
!
!*       3.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL.
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               --------------------------------------------------
!
!
        ZTMIX(:)  = PTHLCL(:) * ( ZPRESMIX(:) / XP00 ) ** ZRDOCP
        ZEVMIX(:) = PRVLCL(:) * ZPRESMIX(:) / ( PRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8, ZEVMIX(:) )
        ZWORK1(:) = LOG( ZEVMIX(:) / 613.3 )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8 - 32.19 * ZWORK1(:) ) / ( 17.502 - ZWORK1(:) ) 
              ! adiabatic saturation temperature
        PTLCL(:)  = ZWORK1(:) - ( .212 + 1.571E-3 * ( ZWORK1(:) - XTT )      &
                  - 4.36E-4 * ( ZTMIX(:) - XTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        PTLCL(:)  = MIN( PTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = XP00 * ( PTLCL(:) / PTHLCL(:) ) ** ZCPORD
!
END WHERE
!
     ZPLCL(:) = MIN( 2.E5, MAX( 10., ZPLCL(:) ) ) ! bound to avoid overflow
!
!
!*       3.2    Correct PTLCL in order to be completely consistent
!               with MNH saturation formula
!               --------------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, PTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / PTLCL(:) * ( XBETAW / PTLCL(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - PRVLCL(:) ) /                              &
                        ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) 
        PTLCL(:)  = PTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
!
     END WHERE
!
!
!*       3.3    If PRVLCL is oversaturated set humidity and temperature
!               to saturation values.
!               -------------------------------------------------------
!
    CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) .AND. PRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( XBETAW / ZTMIX(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - PRVLCL(:) ) /                              &
                        ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        PTLCL(:)  = ZTMIX(:) + ZLV(:) / ZCPH(:) * ZWORK2(:)
        PRVLCL(:) = PRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        PTHLCL(:) = PTLCL(:) * ( XP00 / ZPLCL(:) ) ** ZRDOCP
     END WHERE
!
!
!*        4.1   Determine  vertical loop index at the LCL 
!               -----------------------------------------
!
     DO JK = JKMIN, IKE - 1
        DO JI = 1, IIE
        IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. OWORK1(JI) ) THEN
            KLCL(JI)  = JK + 1
            PZLCL(JI) = PZ(JI,JK+1)
        END IF
        END DO
     END DO
!
!
!*        4.2   Estimate height and environmental temperature at LCL
!               ----------------------------------------------------
!
    DO JI = 1, IIE
        JK   = KLCL(JI)
        JKM  = JK - 1
        ZDP(JI)     = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                      LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI)  = PTH(JI,JK)  * ( PPRES(JI,JK) / XP00 ) ** ZRDOCP
        ZWORK2(JI)  = PTH(JI,JKM) * ( PPRES(JI,JKM) / XP00 ) ** ZRDOCP
        ZWORK1(JI)  = ZWORK2(JI) + ( ZWORK1(JI) - ZWORK2(JI) ) * ZDP(JI) 
           ! we compute the precise value of the LCL
           ! The precise height is between the levels KLCL and KLCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    END DO
    WHERE( OWORK1(:) )
       PTELCL(:) = ZWORK1(:)
       PZLCL(:)  = ZWORK2(:)
    END WHERE
!        
!
!
END SUBROUTINE CONVECT_CLOSURE_THRVLCL
!     ######spl
      SUBROUTINE CONVECT_CHEM_TRANSPORT( KLON, KLEV, KCH, PCH1, PCH1C,       &
                                         KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                         PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                         PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                         KFTSTEPS )
!    #######################################################################
!
!!**** Compute  modified chemical tracer values due to convective event
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!      environmental values of the chemical tracers
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PCH1C-PCH1)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      main deep convection code
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    11/12/97
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAREXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                INTENT(IN) :: KCH      ! number of passive tracers
!
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(IN) :: PCH1 ! grid scale tracer concentr.
REAL,DIMENSION(KLON,KLEV,KCH),INTENT(OUT):: PCH1C! conv adjusted tracer concntr.
!
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)
!
REAL, DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: INCH1          ! number of chemical tracers
INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP        ! vertical loop index
INTEGER :: JN             ! chemical tracer loop index
INTEGER :: JSTEP          ! fractional time loop index
INTEGER :: JKLC, JKLD, JKLP, JKMAX ! loop index for levels
!
REAL, DIMENSION(KLON,KLEV)     :: ZOMG ! compensat. subsidence (Pa/s)
REAL, DIMENSION(KLON,KLEV,KCH) :: ZUCH1, ZDCH1 ! updraft/downdraft values
REAL, DIMENSION(KLON)          :: ZTIMEC  ! fractional convective time step
REAL, DIMENSION(KLON,KLEV)     :: ZTIMC! 2D work array for ZTIMEC
REAL, DIMENSION(KLON,KLEV,KCH) :: ZCH1MFIN, ZCH1MFOUT
                                   ! work arrays for environm. compensat. mass
REAL, DIMENSION(KLON,KCH)      :: ZWORK1, ZWORK2, ZWORK3
!
!-------------------------------------------------------------------------------
!
!*       0.3   Compute loop bounds
!              -------------------
!
INCH1  = KCH
IIE    = KLON
IKB    = 1 + JCVEXB 
IKS    = KLEV
IKE    = KLEV - JCVEXT 
JKMAX  = MAXVAL( KCTL(:) )
!
!
!*      2.      Updraft computations
!               --------------------
!
ZUCH1(:,:,:) = 0.
!
!*      2.1     Initialization  at LCL
!               ----------------------------------
!
DO JI = 1, IIE
    JKLC = KLCL(JI)
    JKLD = KDPL(JI)
    JKLP = KPBL(JI)
    ZWORK1(JI,:) = .5 * ( PCH1(JI,JKLD,:) + PCH1(JI,JKLP,:) )
END DO
!
!*      2.2     Final updraft loop
!               ------------------
!
DO JK = MINVAL( KDPL(:) ), JKMAX
JKP = JK + 1
!
    DO JN = 1, INCH1
     DO JI = 1, IIE
       IF ( KDPL(JI) <= JK .AND. KLCL(JI) > JK )                             &
            ZUCH1(JI,JK,JN) = ZWORK1(JI,JN)
!
       IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
                            !if you have reactive i.e. non-passive tracers
                            !add the corresponding sink term in the following equation
           ZUCH1(JI,JKP,JN) = ( PUMF(JI,JK) * ZUCH1(JI,JK,JN) +              &
                                PUER(JI,JKP) * PCH1(JI,JK,JN) )  /           & 
                              ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7 )
       END IF
     END DO
   END DO
!
END DO
!
!*      3.      Downdraft computations
!               ----------------------
!
ZDCH1(:,:,:) = 0.
!
!*      3.1     Initialization at the LFS
!               -------------------------
!
ZWORK1(:,:) = SPREAD( PMIXF(:), DIM=2, NCOPIES=INCH1 )
DO JI = 1, IIE
     JK = KLFS(JI)
     ZDCH1(JI,JK,:) = ZWORK1(JI,:) * PCH1(JI,JK,:) +                          &
                                       ( 1. - ZWORK1(JI,:) ) * ZUCH1(JI,JK,:)
END DO
!
!*      3.2     Final downdraft loop
!               --------------------
!
DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
JKP = JK - 1
    DO JN = 1, INCH1
    DO JI = 1, IIE
      IF ( JK <= KLFS(JI) .AND. JKP >= KDBL(JI) ) THEN
       ZDCH1(JI,JKP,JN) = ( ZDCH1(JI,JK,JN) * PDMF(JI,JK) -              &
                            PCH1(JI,JK,JN) *  PDER(JI,JKP) ) /           &
                          ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7 ) 
      END IF
    END DO
    END DO
END DO
!
!							   
!*      4.      Final closure (environmental) computations
!               ------------------------------------------
!
PCH1C(:,IKB:IKE,:) = PCH1(:,IKB:IKE,:) ! initialize adjusted envir. values
!
DO JK = IKB, IKE
   ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / XG ! environmental subsidence
END DO
!
ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
                                         ! to be an integer multiple of PTIMEC
WHERE ( PTIMEC(:) < 1. ) ZTIMEC(:) = 0.
ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
!
ZCH1MFIN(:,:,:)   = 0.
ZCH1MFOUT(:,:,:)  = 0.
!
DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop
!
      DO JK = IKB + 1, JKMAX
          JKP = MAX( IKB + 1, JK - 1 )
	  ZWORK3(:,:) = SPREAD( ZOMG(:,JK), DIM=2, NCOPIES=INCH1 )
          ZWORK1(:,:) = SIGN( 1., ZWORK3(:,:) )
          ZWORK2(:,:) = 0.5 * ( 1. + ZWORK1(:,:) )
          ZWORK1(:,:) = 0.5 * ( 1. - ZWORK1(:,:) )
          ZCH1MFIN(:,JK,:)  = - ZWORK3(:,:) * PCH1C(:,JKP,:) * ZWORK1(:,:)
          ZCH1MFOUT(:,JK,:) =   ZWORK3(:,:) * PCH1C(:,JK,:)  * ZWORK2(:,:)
         ZCH1MFIN(:,JKP,:) = ZCH1MFIN(:,JKP,:) + ZCH1MFOUT(:,JK,:) * ZWORK2(:,:)
         ZCH1MFOUT(:,JKP,:)= ZCH1MFOUT(:,JKP,:) + ZCH1MFIN(:,JK,:) * ZWORK1(:,:)
      END DO
!
      DO JN = 1, INCH1
       DO JK = IKB + 1, JKMAX
         PCH1C(:,JK,JN) = PCH1C(:,JK,JN) + ZTIMC(:,JK) / PLMASS(:,JK) *  (    &
                      ZCH1MFIN(:,JK,JN) + PUDR(:,JK) * ZUCH1(:,JK,JN) +       &
                      PDDR(:,JK) * ZDCH1(:,JK,JN) - ZCH1MFOUT(:,JK,JN) -      &
                      ( PUER(:,JK) + PDER(:,JK) ) * PCH1(:,JK,JN)    )
         PCH1C(:,JK,JN) = MAX( 0., PCH1C(:,JK,JN) )
       END DO
      END DO
!
END DO ! Exit the fractional time step loop
!
!
END SUBROUTINE CONVECT_CHEM_TRANSPORT
!     ######spl
      MODULE MODD_CONVPAR_SHAL
!     ########################
!
!!****  *MODD_CONVPAR_SHAL* - Declaration of convection constants 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare  the 
!!      constants in the deep convection parameterization.    
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CONVPAR_SHAL)
!!          
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96                      
!!   Last modified  04/10/98
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
!
REAL, SAVE :: XA25        ! 25 km x 25 km reference grid area
!
REAL, SAVE :: XCRAD       ! cloud radius 
REAL, SAVE :: XCTIME_SHAL ! convective adjustment time
REAL, SAVE :: XCDEPTH     ! minimum necessary cloud depth
REAL, SAVE :: XCDEPTH_D   ! maximum allowed cloud thickness
REAL, SAVE :: XDTPERT     ! add small Temp perturb. at LCL
REAL, SAVE :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)  
!
REAL, SAVE :: XZLCL       ! maximum allowed allowed height 
                          ! difference between departure level and surface
REAL, SAVE :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL, SAVE :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL, SAVE :: XNHGAM      ! accounts for non-hydrost. pressure 
			  ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
REAL, SAVE :: XTFRZ1      ! begin of freezing interval
REAL, SAVE :: XTFRZ2      ! end of freezing interval
!
!
REAL, SAVE :: XSTABT      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
REAL, SAVE :: XSTABC      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE
!$OMP threadprivate(XA25,XCRAD,XCTIME_SHAL,XCDEPTH,XCDEPTH_D,XDTPERT, &
!$OMP XENTR,XZLCL,XZPBL,XWTRIG,XNHGAM,XTFRZ1,XTFRZ2,XSTABT,XSTABC)
!
END MODULE MODD_CONVPAR_SHAL 
!     ######spl
      SUBROUTINE INI_CONVPAR_SHAL
!     ###########################
!
!!****  *INI_CONVPAR * - routine to initialize the constants modules 
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to initialize  the constants
!!     stored in  modules MODD_CONVPAR_SHAL
!!      
!!
!!**  METHOD
!!    ------
!!      The shallow convection constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR_SHAL   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CONVPAR_SHAL, routine INI_CONVPAR)
!!      
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAR_SHAL
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------
!
!
XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
XCRAD       =  50.   ! cloud radius 
XCTIME_SHAL = 10800. ! convective adjustment time
XCDEPTH     = 0.5E3  ! minimum necessary shallow cloud depth
XCDEPTH_D   = 3.0E3  ! maximum allowed shallow cloud depth
XDTPERT     = .2     ! add small Temp perturbation at LCL
!
XENTR    = 0.02      ! entrainment constant (m/Pa) = 0.2 (m)  
!
XZLCL    = 1.5E3     ! maximum allowed allowed height 
                     ! difference between the DPL and the LCL
XZPBL    = 50.E2     ! minimum mixed layer depth to sustain convection
!
!
XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XTFRZ1   = 268.16    ! begin of freezing interval
XTFRZ2   = 248.16    ! end of freezing interval
!

XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
!
!
END SUBROUTINE INI_CONVPAR_SHAL 
!     ######spl
    SUBROUTINE CONVECT_SHALLOW( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                                   PDTCONV, KICE, OSETTADJ, PTADJS,               &
                                   PPABST, PZZ,                                   &
                                   PTT, PRVT, PRCT, PRIT, PWT,                    &
                                   PTTEN, PRVTEN, PRCTEN, PRITEN,                 &
                                   KCLTOP, KCLBAS, PUMF,                          &
                                   OCH1CONV, KCH1, PCH1, PCH1TEN                  )
!   ############################################################################
!
!!**** Monitor routine to compute all convective tendencies by calls
!!     of several subroutines.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      tendencies. The routine first prepares all necessary grid-scale
!!      variables. The final convective tendencies are then computed by
!!      calls of different subroutines.
!!
!!
!!**  METHOD
!!    ------
!!      We start by selecting convective columns in the model domain through
!!      the call of routine TRIGGER_FUNCT. Then, we allocate memory for the
!!      convection updraft and downdraft variables and gather the grid scale
!!      variables in convective arrays. 
!!      The updraft and downdraft computations are done level by level starting
!!      at the  bottom and top of the domain, respectively.
!!      All computations are done on MNH thermodynamic levels. The depth
!!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_TRIGGER_SHAL
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT_SHAL
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_CLOSURE_SHAL
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST_SHAL
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                   ! gravity constant
!!          XPI                  ! number Pi
!!          XP00                 ! reference pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XRHOLW               ! density of liquid water
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!          XCL, XCI             ! specific heat for liquid water and ice
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
!!
!!      Module MODD_CONVPAR
!!          XA25                 ! reference grid area
!!          XCRAD                ! cloud radius
!!
!!         
!!    REFERENCE
!!    ---------
!!
!!      Bechtold, 1997 : Meso-NH scientific  documentation (31 pp)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci., Vol. 37, 1722-1761.
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Peter Bechtold 15/11/96 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAREXT
USE MODD_CONVPAR_SHAL
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                    INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                    INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                  ! KBDIA that is at least 1
INTEGER,                    INTENT(IN) :: KTDIA    ! vertical computations can be
						   ! limited to KLEV + 1 - KTDIA
                                                   ! default=1
REAL,                       INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                   ! calls of the deep convection
                                                   ! scheme
INTEGER,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes, 
                                                   !                0 = no ice )
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
						   ! adjustment time by user
REAL,                       INTENT(IN) :: PTADJS   ! user defined adjustment time
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical 
						   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m) 
!   
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature 
						   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: ITEST, ICONV, ICONV1    ! number of convective columns
INTEGER  :: IIB, IIE                ! horizontal loop bounds
INTEGER  :: IKB, IKE                ! vertical loop bounds
INTEGER  :: IKS                     ! vertical dimension
INTEGER  :: JI, JL                  ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKP, JKM            ! vertical loop index
INTEGER  :: IFTSTEPS                ! only used for chemical tracers
REAL     :: ZEPS, ZEPSA, ZEPSB      ! R_d / R_v, R_v / R_d, XCPV / XCPD - ZEPSA
REAL     :: ZCPORD, ZRDOCP          ! C_p/R_d,  R_d/C_p
!
LOGICAL, DIMENSION(KLON, KLEV)     :: GTRIG3 ! 3D logical mask for convection 
LOGICAL, DIMENSION(KLON)           :: GTRIG  ! 2D logical mask for trigger test
REAL, DIMENSION(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, theta_v
REAL, DIMENSION(KLON)              :: ZTIME  ! convective time period
REAL, DIMENSION(KLON)              :: ZWORK2, ZWORK2B ! work array 
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(:),ALLOCATABLE  :: IDPL    ! index for parcel departure level
INTEGER, DIMENSION(:),ALLOCATABLE  :: IPBL    ! index for source layer top
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILCL    ! index for lifting condensation level 
INTEGER, DIMENSION(:),ALLOCATABLE  :: IETL    ! index for zero buoyancy level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILFS    ! index for level of free sink
!
INTEGER, DIMENSION(:), ALLOCATABLE :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(:), ALLOCATABLE :: ISLCL   ! index for lifting condensation level 
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTHLCL ! updraft theta at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSRVLCL ! updraft rv at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSWLCL  ! updraft w at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSZLCL  ! LCL height
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTHVELCL! envir. theta_v at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSDXDY  ! grid area (m^2)
!
! grid scale variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZZ      ! height of model layer (m) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPRES   ! grid scale pressure
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDPRES  ! pressure difference between 
                                              ! bottom and top of layer (Pa)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZW      ! grid scale vertical velocity on theta grid
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTT     ! temperature
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTH     ! grid scale theta     
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHV    ! grid scale theta_v     
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHL    ! grid scale enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHES, ZTHEST ! grid scale saturated theta_e
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRW     ! grid scale total water (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRV     ! grid scale water vapor (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRC     ! grid scale cloud water (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRI     ! grid scale cloud ice (kg/kg) 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZDXDY   ! grid area (m^2)
!
! updraft variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUMF    ! updraft mass flux (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUER    ! updraft entrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUDR    ! updraft detrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHL   ! updraft enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHV   ! updraft theta_v (K)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURW    ! updraft total water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURC    ! updraft cloud water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURI    ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZMFLCL  ! cloud base unit mass flux(kg/s) 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZCAPE   ! available potent. energy     
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTHLCL  ! updraft theta at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTLCL   ! updraft temp. at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZRVLCL  ! updraft rv at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZWLCL   ! updraft w at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZZLCL   ! LCL height
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTHVELCL! envir. theta_v at LCL
!
! downdraft variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDMF    ! downdraft mass flux (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDER    ! downdraft entrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDDR    ! downdraft detrainment (kg/s)
!
! closure variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZLMASS  ! mass of model layer (kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTIMEC  ! advective time period
!
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRVC    ! conv. adj. grid scale r_w 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRCC    ! conv. adj. grid scale r_c 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRIC    ! conv. adj. grid scale r_i 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZWSUB   ! envir. compensating subsidence (Pa/s)
!
LOGICAL, DIMENSION(:),ALLOCATABLE  :: GTRIG1  ! logical mask for convection    
LOGICAL, DIMENSION(:),ALLOCATABLE  :: GWORK   ! logical work array
INTEGER, DIMENSION(:),ALLOCATABLE  :: IINDEX, IJINDEX, IJSINDEX, IJPINDEX!hor.index
REAL, DIMENSION(:),   ALLOCATABLE  :: ZCPH    ! specific heat C_ph 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZLV, ZLS! latent heat of vaporis., sublim.
REAL                               :: ZES     ! saturation vapor mixng ratio
REAl                               :: ZW1     ! work variable
!
! Chemical Tracers:
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1C   ! conv. adjust. chemical specy 1
REAL, DIMENSION(:,:),   ALLOCATABLE:: ZWORK3  ! conv. adjust. chemical specy 1
LOGICAL, DIMENSION(:,:,:),ALLOCATABLE::GTRIG4 ! logical mask
!
!-------------------------------------------------------------------------------
!
!
!*       0.3    Compute loop bounds
!               -------------------
!
IIB = KIDIA
IIE = KFDIA
JCVEXB = MAX( 0, KBDIA - 1 )
IKB = 1 + JCVEXB 
IKS = KLEV
JCVEXT = MAX( 0, KTDIA - 1)
IKE = IKS - JCVEXT 
!
!
!*       0.5    Update convective counter ( where KCOUNT > 0 
!               convection is still active ).
!               ---------------------------------------------
!
GTRIG(:) = .FALSE.
GTRIG(IIB:IIE) = .TRUE.
ITEST = COUNT( GTRIG(:) )
IF ( ITEST == 0 ) RETURN 
                        
!
!
!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------
!
GTRIG3(:,:) = SPREAD( GTRIG(:), DIM=2, NCOPIES=IKS )
WHERE ( GTRIG3(:,:) ) 
    PTTEN(:,:)  = 0.
    PRVTEN(:,:) = 0.
    PRCTEN(:,:) = 0.
    PRITEN(:,:) = 0.
!   PUTEN(:,:)  = 0.
!   PVTEN(:,:)  = 0.
    PUMF(:,:)   = 0.
END WHERE
WHERE ( GTRIG(:) ) 
 ! KCLTOP(:)  = 0 ! already initialized in CONVECTION
 ! KCLBAS(:)  = 0
END WHERE
IF ( OCH1CONV ) THEN
    ALLOCATE( GTRIG4(KLON,KLEV,KCH1) )
    GTRIG4(:,:,:) = SPREAD( GTRIG3(:,:), DIM=3, NCOPIES=KCH1 )
    WHERE( GTRIG4(:,:,:) ) PCH1TEN(:,:,:) = 0.
    DEALLOCATE( GTRIG4 )
END IF
!
!
!*       1.     Initialize  local variables
!               ----------------------------
!
ZEPS   = XRD / XRV
ZEPSA  = XRV / XRD 
ZEPSB  = XCPV / XCPD - ZEPSA
ZCPORD = XCPD / XRD
ZRDOCP = XRD / XCPD
!
!
!*       1.1    Set up grid scale theta, theta_v, theta_es 
!               ------------------------------------------
!
ZTHT(:,:) = 300.
ZSTHV(:,:)= 300.
ZSTHES(:,:) = 400.
DO JK = IKB, IKE
DO JI = IIB, IIE
   IF ( PPABST(JI,JK) > 40.E2 ) THEN
      ZTHT(JI,JK)  = PTT(JI,JK) * ( XP00 / PPABST(JI,JK) ) ** ZRDOCP
      ZSTHV(JI,JK) = ZTHT(JI,JK) * ( 1. + ZEPSA * PRVT(JI,JK) ) /              &
                     ( 1. + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )
!
          ! use conservative Bolton (1980) formula for theta_e
          ! it is used to compute CAPE for undilute parcel ascent
          ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here
!
      ZES = EXP( XALPW - XBETAW / PTT(JI,JK) - XGAMW * LOG( PTT(JI,JK) ) )
      ZES = MIN( 1., ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
      ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **             &
                ( 1. - 0.28 * ZES ) * EXP( MIN(500.,                           &
                                           ( 3374.6525 / PTT(JI,JK) - 2.5403 ) &
                                          * ZES * ( 1. + 0.81 * ZES ) ) )
   END IF
END DO
END DO
!
!
!
!*       2.     Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------
!
!*       2.1    Allocate arrays depending on number of model columns that need
!               to be tested for convection (i.e. where no convection is present
!               at the moment.
!               --------------------------------------------------------------
!
     ALLOCATE( ZPRES(ITEST,IKS) )
     ALLOCATE( ZZ(ITEST,IKS) )
     ALLOCATE( ZW(ITEST,IKS) )
     ALLOCATE( ZTH(ITEST,IKS) )
     ALLOCATE( ZTHV(ITEST,IKS) )
     ALLOCATE( ZTHEST(ITEST,IKS) )
     ALLOCATE( ZRV(ITEST,IKS) )
     ALLOCATE( ZSTHLCL(ITEST) )
     ALLOCATE( ZSTLCL(ITEST) )
     ALLOCATE( ZSRVLCL(ITEST) )
     ALLOCATE( ZSWLCL(ITEST) )
     ALLOCATE( ZSZLCL(ITEST) )
     ALLOCATE( ZSTHVELCL(ITEST) )
     ALLOCATE( ISDPL(ITEST) )
     ALLOCATE( ISPBL(ITEST) )
     ALLOCATE( ISLCL(ITEST) )
     ALLOCATE( ZSDXDY(ITEST) )
     ALLOCATE( GTRIG1(ITEST) )
     ALLOCATE( IINDEX(KLON) )
     ALLOCATE( IJSINDEX(ITEST) )
     DO JI = 1, KLON
        IINDEX(JI) = JI
     END DO
     IJSINDEX(:) = PACK( IINDEX(:), MASK=GTRIG(:) )
!
  DO JK = IKB, IKE
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZPRES(JI,JK)  = PPABST(JL,JK)
     ZZ(JI,JK)     = PZZ(JL,JK)
     ZTH(JI,JK)    = ZTHT(JL,JK)
     ZTHV(JI,JK)   = ZSTHV(JL,JK)
     ZTHEST(JI,JK) = ZSTHES(JL,JK)
     ZRV(JI,JK)    = MAX( 0., PRVT(JL,JK) )
     ZW(JI,JK)     = PWT(JL,JK)
  END DO
  END DO
  DO JI = 1, ITEST
     JL = IJSINDEX(JI)
     ZSDXDY(JI)    = XA25
  END DO
!
!*       2.2    Compute environm. enthalpy and total water = r_v + r_i + r_c 
!               and envir. saturation theta_e
!               ------------------------------------------------------------
!
!
!*       2.3    Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------
!
     ISLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL 
     ISDPL(:) = IKB
     ISPBL(:) = IKB
!
     CALL CONVECT_TRIGGER_SHAL(  ITEST, KLEV,                              &
                                 ZPRES, ZTH, ZTHV, ZTHEST,                 &
                                 ZRV, ZW, ZZ, ZSDXDY,                      &
                                 ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                                 ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1    )
!
     DEALLOCATE( ZPRES )
     DEALLOCATE( ZZ )
     DEALLOCATE( ZTH )
     DEALLOCATE( ZTHV )
     DEALLOCATE( ZTHEST )
     DEALLOCATE( ZRV )
     DEALLOCATE( ZW )
!
!
!*       3.     After the call of TRIGGER_FUNCT we allocate all the dynamic
!               arrays used in the convection scheme using the mask GTRIG, i.e.
!               we do calculus only in convective columns. This corresponds to
!               a GATHER operation.
!               --------------------------------------------------------------
!
     ICONV = COUNT( GTRIG1(:) )
     IF ( ICONV == 0 )  THEN 
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( GTRIG1 )
         DEALLOCATE( IINDEX )
         DEALLOCATE( IJSINDEX )
         RETURN   ! no convective column has been found, exit CONVECT_SHALLOW 
     ENDIF
!
     ! vertical index variables
!
         ALLOCATE( IDPL(ICONV) )
         ALLOCATE( IPBL(ICONV) )
         ALLOCATE( ILCL(ICONV) )
         ALLOCATE( ICTL(ICONV) )
         ALLOCATE( IETL(ICONV) )
!
	 ! grid scale variables
!
         ALLOCATE( ZZ(ICONV,IKS) )
         ALLOCATE( ZPRES(ICONV,IKS) )
         ALLOCATE( ZDPRES(ICONV,IKS) )
         ALLOCATE( ZTT(ICONV, IKS) )
         ALLOCATE( ZTH(ICONV,IKS) )
         ALLOCATE( ZTHV(ICONV,IKS) )
         ALLOCATE( ZTHL(ICONV,IKS) )
         ALLOCATE( ZTHES(ICONV,IKS) )
         ALLOCATE( ZRV(ICONV,IKS) )
         ALLOCATE( ZRC(ICONV,IKS) )
         ALLOCATE( ZRI(ICONV,IKS) )
         ALLOCATE( ZRW(ICONV,IKS) )
         ALLOCATE( ZDXDY(ICONV) )
!
         ! updraft variables
!
         ALLOCATE( ZUMF(ICONV,IKS) )
         ALLOCATE( ZUER(ICONV,IKS) )
         ALLOCATE( ZUDR(ICONV,IKS) )
         ALLOCATE( ZUTHL(ICONV,IKS) )
         ALLOCATE( ZUTHV(ICONV,IKS) )
         ALLOCATE( ZURW(ICONV,IKS) )
         ALLOCATE( ZURC(ICONV,IKS) )
         ALLOCATE( ZURI(ICONV,IKS) )
         ALLOCATE( ZTHLCL(ICONV) )
         ALLOCATE( ZTLCL(ICONV) )
         ALLOCATE( ZRVLCL(ICONV) )
         ALLOCATE( ZWLCL(ICONV) )
         ALLOCATE( ZMFLCL(ICONV) )
         ALLOCATE( ZZLCL(ICONV) )
         ALLOCATE( ZTHVELCL(ICONV) )
         ALLOCATE( ZCAPE(ICONV) )
!
         ! work variables
!
         ALLOCATE( IJINDEX(ICONV) )
         ALLOCATE( IJPINDEX(ICONV) )
         ALLOCATE( ZCPH(ICONV) )
         ALLOCATE( ZLV(ICONV) )
         ALLOCATE( ZLS(ICONV) )
!
!
!*           3.1    Gather grid scale and updraft base variables in
!                   arrays using mask GTRIG
!                   ---------------------------------------------------
!
         GTRIG(:)      = UNPACK( GTRIG1(:), MASK=GTRIG(:), FIELD=.FALSE. )  
         IJINDEX(:)    = PACK( IINDEX(:), MASK=GTRIG(:) )
!
    DO JK = IKB, IKE
    DO JI = 1, ICONV
         JL = IJINDEX(JI)
         ZZ(JI,JK)     = PZZ(JL,JK)
         ZPRES(JI,JK)  = PPABST(JL,JK)
         ZTT(JI,JK)    = PTT(JL,JK)
         ZTH(JI,JK)    = ZTHT(JL,JK)
         ZTHES(JI,JK)  = ZSTHES(JL,JK)
         ZRV(JI,JK)    = MAX( 0., PRVT(JL,JK) )
         ZRC(JI,JK)    = MAX( 0., PRCT(JL,JK) )
         ZRI(JI,JK)    = MAX( 0., PRIT(JL,JK) )
         ZTHV(JI,JK)   = ZSTHV(JL,JK)
    END DO
    END DO
!
    DO JI = 1, ITEST
       IJSINDEX(JI) = JI	
    END DO
    IJPINDEX(:) = PACK( IJSINDEX(:), MASK=GTRIG1(:) )
    DO JI = 1, ICONV
         JL = IJPINDEX(JI)
         IDPL(JI)      = ISDPL(JL)
         IPBL(JI)      = ISPBL(JL)
         ILCL(JI)      = ISLCL(JL)
         ZTHLCL(JI)    = ZSTHLCL(JL)
         ZTLCL(JI)     = ZSTLCL(JL)
         ZRVLCL(JI)    = ZSRVLCL(JL)
         ZWLCL(JI)     = ZSWLCL(JL)
         ZZLCL(JI)     = ZSZLCL(JL)
         ZTHVELCL(JI)  = ZSTHVELCL(JL)
         ZDXDY(JI)     = ZSDXDY(JL)
    END DO
         ALLOCATE( GWORK(ICONV) )
         GWORK(:)      = PACK( GTRIG1(:),  MASK=GTRIG1(:) ) 
         DEALLOCATE( GTRIG1 )
         ALLOCATE( GTRIG1(ICONV) )
         GTRIG1(:)     = GWORK(:)
!                 
         DEALLOCATE( GWORK )
         DEALLOCATE( IJPINDEX )
         DEALLOCATE( ISDPL )
         DEALLOCATE( ISPBL )
         DEALLOCATE( ISLCL )
         DEALLOCATE( ZSTHLCL )
         DEALLOCATE( ZSTLCL )
         DEALLOCATE( ZSRVLCL )
         DEALLOCATE( ZSWLCL )
         DEALLOCATE( ZSZLCL )
         DEALLOCATE( ZSTHVELCL )
         DEALLOCATE( ZSDXDY )
!
!
!*           3.2    Compute pressure difference 
!                   ---------------------------------------------------
!
        ZDPRES(:,IKB) = 0.
        DO JK = IKB + 1, IKE
            ZDPRES(:,JK)  = ZPRES(:,JK-1) - ZPRES(:,JK)
        END DO
!
!*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c
!                  ----------------------------------------------------------
!
        DO JK = IKB, IKE, 1
            ZRW(:,JK)  = ZRV(:,JK) + ZRC(:,JK) + ZRI(:,JK)
            ZCPH(:)    = XCPD + XCPV * ZRW(:,JK)
            ZLV(:)     = XLVTT + ( XCPV - XCL ) * ( ZTT(:,JK) - XTT ) ! compute L_v
            ZLS(:)     = XLSTT + ( XCPV - XCI ) * ( ZTT(:,JK) - XTT ) ! compute L_i
            ZTHL(:,JK) = ZCPH(:) * ZTT(:,JK) + ( 1. + ZRW(:,JK) ) * XG * ZZ(:,JK) &
                         - ZLV(:) * ZRC(:,JK) - ZLS(:) * ZRI(:,JK)
        END DO
!
        DEALLOCATE( ZCPH )
        DEALLOCATE( ZLV )
        DEALLOCATE( ZLS )
!
!
!*           4.     Compute updraft properties 
!                   ----------------------------
!
!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s ) 
!                   -------------------------------------------------------------
!
         ZDXDY(:)  = XA25
         ZMFLCL(:) = XA25 * 1.e-3
!
!
     CALL CONVECT_UPDRAFT_SHAL( ICONV, KLEV,                                     &
                                KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
                                ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   & 
                                ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
                                ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                                ZURC, ZURI, ZCAPE, ICTL, IETL                    )
!
!
!
!*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud 
!                   thickness is smaller than 3 km
!                   -----------------------------------------------------------
!
!
     ICONV1 = COUNT(GTRIG1) 
!
     IF ( ICONV1 > 0 )  THEN
!
!*       4.3    Allocate memory for downdraft variables
!               ---------------------------------------
!
! downdraft variables
!
        ALLOCATE( ZDMF(ICONV,IKS) )
        ALLOCATE( ZDER(ICONV,IKS) )
        ALLOCATE( ZDDR(ICONV,IKS) )
        ALLOCATE( ILFS(ICONV) )
        ALLOCATE( ZLMASS(ICONV,IKS) )
        ZDMF(:,:) = 0.
        ZDER(:,:) = 0.
        ZDDR(:,:) = 0.
        ILFS(:)   = IKB
        DO JK = IKB, IKE
           ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / XG  ! mass of model layer
        END DO
	ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
!
! closure variables
!
        ALLOCATE( ZTIMEC(ICONV) )
        ALLOCATE( ZTHC(ICONV,IKS) )
        ALLOCATE( ZRVC(ICONV,IKS) )
        ALLOCATE( ZRCC(ICONV,IKS) )
        ALLOCATE( ZRIC(ICONV,IKS) )
        ALLOCATE( ZWSUB(ICONV,IKS) )
!
!
!*           5.     Compute downdraft properties 
!                   ----------------------------
!
        ZTIMEC(:) = XCTIME_SHAL
        IF ( OSETTADJ ) ZTIMEC(:) = PTADJS
!
!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------
!
       CALL CONVECT_CLOSURE_SHAL( ICONV, KLEV,                         &
                                  ZPRES, ZDPRES, ZZ, ZDXDY, ZLMASS,    &
                                  ZTHL, ZTH, ZRW, ZRC, ZRI, GTRIG1,    &
                                  ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,       &
                                  ILCL, IDPL, IPBL, ICTL,              &
                                  ZUMF, ZUER, ZUDR, ZUTHL, ZURW,       &
                                  ZURC, ZURI, ZCAPE, ZTIMEC, IFTSTEPS  )

!
!
!
!*           8.     Determine the final grid-scale (environmental) convective 
!                   tendencies and set convective counter
!                   --------------------------------------------------------
!
!
!*           8.1    Grid scale tendencies
!                   ---------------------
!
          ! in order to save memory, the tendencies are temporarily stored
          ! in the tables for the adjusted grid-scale values
!
      DO JK = IKB, IKE
         ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
           * ( ZPRES(:,JK) / XP00 ) ** ZRDOCP ! change theta in temperature
         ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) ) &
                                                  / ZTIMEC(:)

         ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
         ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:)
!
      END DO
!
!
!*           8.2    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
       JKM = MAXVAL( ICTL(:) )
       ZWORK2(:) = 0.
       ZWORK2B(:) = 0.
       DO JK = JKM, IKB+1, -1
	 JKP = JK + 1
         DO JI = 1, ICONV
           ZW1 =  ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK)
           ZWORK2(JI) = ZWORK2(JI) +  ZW1 *          & ! moisture
                                             (ZPRES(JI,JK) - ZPRES(JI,JKP)) / XG
           ZW1 = ( XCPD + XCPV * ZRW(JI,JK) )* ZTHC(JI,JK)   - &
                 ( XLVTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRCC(JI,JK) - &
                 ( XLSTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRIC(JI,JK) 
           ZWORK2B(JI) = ZWORK2B(JI) + ZW1 *         & ! energy
                                             (ZPRES(JI,JK) - ZPRES(JI,JKP)) / XG
         END DO
       END DO
!
          ! Budget error (integral must be zero)
!
       DO JI = 1, ICONV
         IF ( ICTL(JI) > 2 ) THEN
           JKP = ICTL(JI) + 1
           ZWORK2(JI) =  ZWORK2(JI) * XG / ( ZPRES(JI,IKB+1) - ZPRES(JI,JKP) )
           ZWORK2B(JI) = ZWORK2B(JI) * XG / ( ZPRES(JI,IKB+1) - ZPRES(JI,JKP) )
         END IF
       END DO
!
          ! Apply uniform correction
!
       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( ICTL(JI) > 2 .AND. JK <= ICTL(JI) ) THEN
                 ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)                                ! moisture
                 ZTHC(JI,JK) = ZTHC(JI,JK) - ZWORK2B(JI) / ( XCPD + XCPV * ZRW(JI,JK) )! energy
           END IF
         END DO
       END DO
!
	      ! execute a "scatter"= pack command to store the tendencies in
	      ! the final 2D tables
!
      DO JK = IKB, IKE
      DO JI = 1, ICONV
         JL = IJINDEX(JI)
         PTTEN(JL,JK)   = ZTHC(JI,JK)
         PRVTEN(JL,JK)  = ZRVC(JI,JK)
         PRCTEN(JL,JK)  = ZRCC(JI,JK)
         PRITEN(JL,JK)  = ZRIC(JI,JK)
      END DO
      END DO
!
!
!                   Cloud base and top levels
!                   -------------------------
!
     ILCL(:) = MIN( ILCL(:), ICTL(:) )
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        KCLTOP(JL) = ICTL(JI)
        KCLBAS(JL) = ILCL(JI)
     END DO
!
!
!*           8.7    Compute convective tendencies for Tracers
!                   ------------------------------------------
!
  IF ( OCH1CONV ) THEN
!
       ALLOCATE( ZCH1(ICONV,IKS,KCH1) )
       ALLOCATE( ZCH1C(ICONV,IKS,KCH1) )
       ALLOCATE( ZWORK3(ICONV,KCH1) )
!
       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          ZCH1(JI,JK,:) = PCH1(JL,JK,:)
       END DO
       END DO
!
      CALL CONVECT_CHEM_TRANSPORT( ICONV, KLEV, KCH1, ZCH1, ZCH1C,          &
                                   IDPL, IPBL, ILCL, ICTL, ILFS, ILFS,      &
                                   ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,      &
                                   ZTIMEC, ZDXDY, ZDMF(:,1), ZLMASS, ZWSUB, &
                                   IFTSTEPS )
!
       DO JK = IKB, IKE
       DO JN = 1, KCH1
          ZCH1C(:,JK,JN) = ( ZCH1C(:,JK,JN)- ZCH1(:,JK,JN) ) / ZTIMEC(:)
       END DO
       END DO
!
!
!*           8.8    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
       JKM = MAXVAL( ICTL(:) )
       ZWORK3(:,:) = 0.
       DO JK = JKM, IKB+1, -1
	 JKP = JK + 1
         DO JI = 1, ICONV
           ZWORK3(JI,:) = ZWORK3(JI,:) + ZCH1C(JI,JK,:) *                    &
                              (ZPRES(JI,JK) - ZPRES(JI,JKP)) / XG
         END DO
       END DO
!
          ! Mass error (integral must be zero)
!
       DO JI = 1, ICONV
           JKP = ICTL(JI) + 1
           ZWORK3(JI,:) = ZWORK3(JI,:) *                                     &
                                    XG / ( ZPRES(JI,IKB+1) - ZPRES(JI,JKP) )
       END DO
!
          ! Apply uniform correction but assure positive mass at each level
!
       DO JK = JKM, IKB+1, -1
         DO JI = 1, ICONV
           IF ( JK <= ICTL(JI) ) THEN
                ZCH1C(JI,JK,:) = ZCH1C(JI,JK,:) - ZWORK3(JI,:)
       !        ZCH1C(JI,JK,:) = MAX( ZCH1C(JI,JK,:), -ZCH1(JI,JK,:)/ZTIMEC(JI) )
           END IF
         END DO
       END DO
!
       DO JK = IKB, IKE
       DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PCH1TEN(JL,JK,:) = ZCH1C(JI,JK,:)
       END DO
       END DO
  END IF
!
!
!*           9.     Write up- and downdraft mass fluxes
!                   ------------------------------------
!
    DO JK = IKB, IKE
       ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
    END DO
    ZWORK2(:) = 1.
    DO JK = IKB, IKE
    DO JI = 1, ICONV
       JL = IJINDEX(JI)
       PUMF(JL,JK) = ZUMF(JI,JK) * ZWORK2(JL)
    END DO
    END DO
!
!
!*           10.    Deallocate all local arrays
!                   ---------------------------
!
! downdraft variables
!
      DEALLOCATE( ZDMF )
      DEALLOCATE( ZDER )
      DEALLOCATE( ZDDR )
      DEALLOCATE( ILFS )
      DEALLOCATE( ZLMASS )
!
!   closure variables
!
      DEALLOCATE( ZTIMEC )
      DEALLOCATE( ZTHC )
      DEALLOCATE( ZRVC )
      DEALLOCATE( ZRCC )
      DEALLOCATE( ZRIC )
      DEALLOCATE( ZWSUB )
!
       IF ( OCH1CONV ) THEN
           DEALLOCATE( ZCH1 )
           DEALLOCATE( ZCH1C )
           DEALLOCATE( ZWORK3 )
       END IF
!
    ENDIF
!
!    vertical index
!
    DEALLOCATE( IDPL )
    DEALLOCATE( IPBL )
    DEALLOCATE( ILCL )
    DEALLOCATE( ICTL )
    DEALLOCATE( IETL )
!
! grid scale variables
!
    DEALLOCATE( ZZ )
    DEALLOCATE( ZPRES )
    DEALLOCATE( ZDPRES )
    DEALLOCATE( ZTT )
    DEALLOCATE( ZTH )
    DEALLOCATE( ZTHV )
    DEALLOCATE( ZTHL )
    DEALLOCATE( ZTHES )
    DEALLOCATE( ZRW )
    DEALLOCATE( ZRV )
    DEALLOCATE( ZRC )
    DEALLOCATE( ZRI )
    DEALLOCATE( ZDXDY )
!
! updraft variables
!
    DEALLOCATE( ZUMF )
    DEALLOCATE( ZUER )
    DEALLOCATE( ZUDR )
    DEALLOCATE( ZUTHL )
    DEALLOCATE( ZUTHV )
    DEALLOCATE( ZURW )
    DEALLOCATE( ZURC )
    DEALLOCATE( ZURI )
    DEALLOCATE( ZTHLCL )
    DEALLOCATE( ZTLCL )
    DEALLOCATE( ZRVLCL )
    DEALLOCATE( ZWLCL )
    DEALLOCATE( ZZLCL )
    DEALLOCATE( ZTHVELCL )
    DEALLOCATE( ZMFLCL )
    DEALLOCATE( ZCAPE )
!
! work arrays
!
    DEALLOCATE( IINDEX )
    DEALLOCATE( IJINDEX )
    DEALLOCATE( IJSINDEX )
    DEALLOCATE( GTRIG1 )
!
!
END SUBROUTINE CONVECT_SHALLOW
!     ######spl
      SUBROUTINE CONVECT_TRIGGER_SHAL(  KLON, KLEV,                           &
                                        PPRES, PTH, PTHV, PTHES,              &
                                        PRV, PW, PZ, PDXDY,                   &
                                        PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL,  &
                                        PTHVELCL, KLCL, KDPL, KPBL, OTRIG     )
!     ######################################################################
!
!!**** Determine convective columns as well as the cloudy values of theta,
!!     and qv at the lifting condensation level (LCL) 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine convective columns
!!   
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      What we look for is the undermost unstable level at each grid point.
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! Reference pressure
!!          XRD, XRV           ! Gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XTT                ! triple point temperature
!!          XBETAW, XGAMW      ! constants for vapor saturation pressure
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! maximum height difference between
!!                             ! the surface and the DPL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XCDEPTH_D          ! maximum allowed cloud depth
!!          XDTPERT            ! add small Temp peturbation
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  20/03/97  Select first departure level
!!                            that produces a cloud thicker than XCDEPTH
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR_SHAL
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                   :: KLON      ! horizontal loop index
INTEGER, INTENT(IN)                   :: KLEV      ! vertical loop index
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY     ! grid area
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH, PTHV ! theta, theta_v
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHES     ! envir. satur. theta_e
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRV       ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PPRES     ! pressure
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PZ        ! height of grid point (m)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PW        ! vertical velocity
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL    ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL     ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PWLCL     ! parcel velocity at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL     ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(OUT):: OTRIG     ! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KLCL    ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KDPL    ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KPBL    ! contains index of source layer top
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JKK, JK, JKP, JKM, JL, JKT, JT      ! vertical loop index
INTEGER :: JI                                  ! horizontal loop index 
INTEGER :: IIE, IKB, IKE                       ! horizontal + vertical loop bounds
REAL    :: ZEPS, ZEPSA                         ! R_d / R_v, R_v / R_d 
REAL    :: ZCPORD, ZRDOCP                      ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON) :: ZTHLCL, ZTLCL, ZRVLCL, & ! locals for PTHLCL,PTLCL
                               ZWLCL,  ZZLCL, ZTHVELCL  ! PRVLCL, ....
INTEGER, DIMENSION(KLON) :: IDPL, IPBL, ILCL      ! locals for KDPL, ...
REAL, DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL, DIMENSION(KLON) :: ZZDPL    ! height of DPL 
REAL, DIMENSION(KLON) :: ZTHVLCL  ! theta_v at LCL = mixed layer value
REAL, DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL, DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure 
REAL, DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL, DIMENSION(KLON) :: ZCAPE    ! convective available energy (m^2/s^2/g)
REAL, DIMENSION(KLON) :: ZTHEUL   ! updraft equiv. pot. temperature (K)
REAL, DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL, DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(KLON) :: ZTOP     ! estimated cloud top (m)
!INTEGER, DIMENSION(KLON) :: ITOP  ! work array to store highest test layer
REAL, DIMENSION(KLON) :: ZWORK1, ZWORK2, ZWORK3    ! work arrays
LOGICAL, DIMENSION(KLON) :: GTRIG, GTRIG2          ! local arrays for OTRIG
LOGICAL, DIMENSION(KLON) :: GWORK1                 ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Initialize local variables
!               --------------------------
!
ZEPS       = XRD / XRV
ZEPSA      = XRV / XRD 
ZCPORD     = XCPD / XRD
ZRDOCP     = XRD / XCPD
OTRIG(:)   = .FALSE.
IDPL(:)    = KDPL(:)
IPBL(:)    = KPBL(:)
ILCL(:)    = KLCL(:)
!ITOP(:)    = IKB
PWLCL(:)   = 0.
ZWLCL(:)   = 0.
PTHLCL(:)  = 1.
PTHVELCL(:)= 1.
PTLCL(:)   = 1.
PRVLCL(:)  = 0.
PWLCL(:)   = 0.
PZLCL(:)   = PZ(:,IKB)
ZZDPL(:)   = PZ(:,IKB)
GTRIG2(:)  = .TRUE.
!
!
!
!       1.     Determine highest necessary loop test layer
!              -------------------------------------------
!
JT = IKE - 2
DO JK = IKB + 1, IKE - 2
 ! DO JI = 1, IIE
 !    IF ( PZ(JI,JK) - PZ(JI,IKB) <= XZLCL ) ITOP(JI) = JK
 ! END DO
   IF ( PZ(1,JK) - PZ(1,IKB) < 5.E3 ) JT = JK 
END DO
!
!
!*       2.     Enter loop for convection test
!               ------------------------------
!
JKP = MINVAL( IDPL(:) ) + 1
!JKT = MAXVAL( ITOP(:) )
JKT = JT
DO JKK = JKP, JKT
!
     GWORK1(:) = ZZDPL(:) - PZ(:,IKB) < XZLCL
          ! we exit the trigger test when the center of the mixed layer is more
          ! than 1500 m  above soil level.
     WHERE ( GWORK1(:) )
        ZDPTHMIX(:) = 0.
        ZPRESMIX(:) = 0.
        ZTHLCL(:)   = 0.
        ZRVLCL(:)   = 0.
        ZZDPL(:)    = PZ(:,JKK)
        IDPL(:)     = JKK
     END WHERE
!
!
!*       3.     Construct a mixed layer of at least 50 hPa (XZPBL)
!               ------------------------------------------
!
     DO JK = JKK, IKE - 1
       JKM = JK + 1
       DO JI = 1, IIE     
         IF ( GWORK1(JI) .AND. ZDPTHMIX(JI) < XZPBL ) THEN
            IPBL(JI)     = JK
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            ZTHLCL(JI)   = ZTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            ZRVLCL(JI)   = ZRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
         END IF
       END DO
        IF ( MINVAL ( ZDPTHMIX(:) ) >= XZPBL ) EXIT
     END DO
!
!
     WHERE ( GWORK1(:) )
!
        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
        ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) + XDTPERT ! add small Temp Perturb.
        ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:)
        ZTHVLCL(:)  = ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                 &
				/ ( 1. + ZRVLCL(:) )
!
!*       4.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL. 
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               ----------------------------------------------------
!
! 
        ZTMIX(:)  = ZTHLCL(:) * ( ZPRESMIX(:) / XP00 ) ** ZRDOCP 
        ZEVMIX(:) = ZRVLCL(:) * ZPRESMIX(:) / ( ZRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8, ZEVMIX(:) )
        ZWORK1(:) = LOG( ZEVMIX(:) / 613.3 )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8 - 32.19 * ZWORK1(:) ) / ( 17.502 - ZWORK1(:) ) 
              ! adiabatic saturation temperature
        ZTLCL(:)  = ZWORK1(:) - ( .212 + 1.571E-3 * ( ZWORK1(:) - XTT )      &
                   - 4.36E-4 * ( ZTMIX(:) - XTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        ZTLCL(:)  = MIN( ZTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = XP00 * ( ZTLCL(:) / ZTHLCL(:) ) ** ZCPORD
!
     END WHERE
!
!
!*       4.2    Correct ZTLCL in order to be completely consistent
!               with MNH saturation formula
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, ZTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTLCL(:) * ( XBETAW / ZTLCL(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                        ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) 
        ZTLCL(:)  = ZTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
!
     END WHERE
!
!
!*       4.3    If ZRVLCL = PRVMIX is oversaturated set humidity 
!               and temperature to saturation values. 
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) .AND. ZRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( XBETAW / ZTMIX(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                       ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) 
        ZTLCL(:)  = ZTMIX(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
        ZRVLCL(:) = ZRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        ZTHLCL(:) = ZTLCL(:) * ( XP00 / ZPLCL(:) ) ** ZRDOCP
        ZTHVLCL(:)= ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                   &
                              / ( 1. + ZRVLCL(:) )
     END WHERE
!
!
!*        5.1   Determine  vertical loop index at the LCL and DPL
!               --------------------------------------------------
!
    DO JK = JKK, IKE - 1
       DO JI = 1, IIE
         IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. GWORK1(JI) ) ILCL(JI) = JK + 1
       END DO
    END DO
!
!
!*        5.2   Estimate height and environm. theta_v at LCL
!               --------------------------------------------------
!
    DO JI = 1, IIE
        JK   = ILCL(JI)
        JKM  = JK - 1
        ZDP(JI)    = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                     LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI) 
           ! we compute the precise value of the LCL
           ! The precise height is between the levels ILCL and ILCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    END DO
    WHERE( GWORK1(:) )
        ZTHVELCL(:) = ZWORK1(:) 
        ZZLCL(:)    = ZWORK2(:)
    END WHERE
!        
!
!*       6.     Check to see if cloud is bouyant 
!               --------------------------------
!
!*      6.1    Compute grid scale vertical velocity perturbation term ZWORK1
!               -------------------------------------------------------------
! 
!            !  normalize w grid scale to a 25 km refer. grid
!    DO JI = 1, IIE
!       JK  = ILCL(JI)
!       JKM = JK - 1 
!       ZWORK1(JI) =  ( PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) )  &
!                          * SQRT( PDXDY(JI) / XA25 )
!                         - 0.02 * ZZLCL(JI) / XZLCL ! avoid spurious convection
!    END DO
!            ! compute sign of normalized grid scale w
!       ZWORK2(:) = SIGN( 1., ZWORK1(:) ) 
!       ZWORK1(:) = XWTRIG * ZWORK2(:) * ABS( ZWORK1(:) ) ** 0.333       &
!                          * ( XP00 / ZPLCL(:) ) ** ZRDOCP
!
!*       6.2    Compute parcel vertical velocity at LCL
!               ---------------------------------------
!                   
!    DO JI = 1, IIE
!       JKDL = IDPL(JI)
!       ZWORK3(JI) = XG * ZWORK1(JI) * ( ZZLCL(JI) - PZ(JI,JKDL) )       &
!                      / ( PTHV(JI,JKDL) + ZTHVELCL(JI) )
!    END DO
!    WHERE( GWORK1(:) )
!      ZWLCL(:)  = 1. + .5 * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) ) 
!      GTRIG(:)  = ZTHVLCL(:) - ZTHVELCL(:) + ZWORK1(:) > 0. .AND.       &
!                  ZWLCL(:) > 0. 
!    END WHERE
     ZWLCL(:) = 1.
!
!
!*       6.3    Look for parcel that produces sufficient cloud depth.
!               The cloud top is estimated as the level where the CAPE 
!               is smaller  than a given value (based on vertical velocity eq.)
!               --------------------------------------------------------------
!
     ZTHEUL(:) = ZTLCL(:) * ( ZTHLCL(:) / ZTLCL(:) )                       &
                                             ** ( 1. - 0.28 * ZRVLCL(:) )  &
                          * EXP( ( 3374.6525 / ZTLCL(:) - 2.5403 ) *       &
                               ZRVLCL(:) * ( 1. + 0.81 * ZRVLCL(:) ) )
!
     ZCAPE(:) = 0.
     ZTOP(:)  = 0.
     ZWORK3(:)= 0.
     JKM = MINVAL( ILCL(:) )
     DO JL = JKM, JT
        JK = JL + 1
        DO JI = 1, IIE
           ZWORK1(JI) = ( 2. * ZTHEUL(JI) /                                &
            ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1. ) * ( PZ(JI,JK) - PZ(JI,JL) )
           IF ( JL < ILCL(JI) ) ZWORK1(JI) = 0.
           ZCAPE(JI)  = ZCAPE(JI) + ZWORK1(JI)
           ZWORK2(JI) = XNHGAM * XG * ZCAPE(JI) + 1.05 * ZWLCL(JI) * ZWLCL(JI)
               ! the factor 1.05 takes entrainment into account
           ZWORK2(JI) = SIGN( 1., ZWORK2(JI) )
           ZWORK3(JI) = ZWORK3(JI) + MIN(0., ZWORK2(JI) )
           ZWORK3(JI) = MAX( -1., ZWORK3(JI) )
               ! Nota, the factors ZWORK2 and ZWORK3 are only used to avoid
               ! if and goto statements, the difficulty is to extract only
               ! the level where the criterium is first fullfilled
           ZTOP(JI)   = PZ(JI,JL) * .5 * ( 1. + ZWORK2(JI) ) * ( 1. + ZWORK3(JI) ) + &
                        ZTOP(JI) * .5 * ( 1. - ZWORK2(JI) )
         END DO
     END DO
!
!
     ZWORK2(:) = ZTOP(:) - ZZLCL(:)
     WHERE( ZWORK2(:)  .GE. XCDEPTH  .AND. ZWORK2(:) < XCDEPTH_D .AND. GTRIG2(:) )
        GTRIG2(:)   = .FALSE.
        OTRIG(:)    = .TRUE.
      ! OTRIG(:)    = GTRIG(:)     ! we  select the first departure level
        PTHLCL(:)   = ZTHLCL(:)    ! that gives sufficient cloud depth
        PRVLCL(:)   = ZRVLCL(:)
        PTLCL(:)    = ZTLCL(:)
        PWLCL(:)    = ZWLCL(:)
        PZLCL(:)    = ZZLCL(:)
        PTHVELCL(:) = ZTHVELCL(:)
        KDPL(:)     = IDPL(:)
        KPBL(:)     = IPBL(:)
        KLCL(:)     = ILCL(:)
     END WHERE
!
END DO
!
!
END SUBROUTINE CONVECT_TRIGGER_SHAL
!     ######spl
    SUBROUTINE CONVECT_UPDRAFT_SHAL( KLON, KLEV,                                     &
     		                     KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW,&
                                     PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,  &
                                     PMFLCL, OTRIG, KLCL, KDPL, KPBL,                &
                                     PUMF, PUER, PUDR, PUTHL, PUTHV, PURW,           &
                                     PURC, PURI, PCAPE, KCTL, KETL )
!    ###############################################################################
!
!!**** Compute updraft properties from DPL to CTL. 
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine updraft properties
!!      ( mass flux, thermodynamics, precipitation ) 
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_MIXING_FUNCT
!!     Routine CONVECT_CONDENS
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV, XCL    ! Cp of dry air, water vapor and liquid water
!!          XTT                ! triple point temperature
!!          XLVTT              ! vaporisation heat at XTT
!!        
!!
!!      Module MODD_CONVPAR_SHAL
!!          XA25               ! reference grid area
!!          XCRAD              ! cloud radius
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XENTR              ! entrainment constant
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!          XTFRZ1             ! begin of freezing interval
!!          XTFRZ2             ! begin of freezing interval
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_UPDRAFT)
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  10/12/97
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR_SHAL
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                    :: KLON  ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV  ! vertical dimension
INTEGER, INTENT(IN)                    :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHL  ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHV  ! grid scale theta_v     
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water  
                                                ! mixing ratio 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between 
                                                ! bottom and top of layer (Pa) 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHLCL ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PTLCL  ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PRVLCL ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PMFLCL ! cloud  base unit mass flux
                                                ! (kg/s)
REAL, DIMENSION(KLON),     INTENT(IN) :: PZLCL  ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(INOUT):: OTRIG! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! contains vert. index of DPL 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   !  " vert. index of source layertop
!
!
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KCTL   ! contains vert. index of CTL 
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KETL   ! contains vert. index of        &
						!equilibrium (zero buoyancy) level 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHL ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHV ! updraft theta_v (K)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURW  ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURC  ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(OUT):: PURI  ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(KLON),     INTENT(OUT):: PCAPE  ! available potent. energy
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP, JKM, JK1, JK2, JKMIN   ! vertical loop index
REAL    :: ZEPSA, ZCVOCD  ! R_v / R_d, C_pv / C_pd 
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON)    :: ZUT             ! updraft temperature (K)
REAL, DIMENSION(KLON)    :: ZUW1, ZUW2      ! square of updraft vert.
                                            ! velocity at levels k and k+1
REAL, DIMENSION(KLON)    :: ZE1,ZE2,ZD1,ZD2 ! fractional entrainm./detrain
                                            ! rates at levels k and k+1
REAL, DIMENSION(KLON)    :: ZMIXF           ! critical mixed fraction  
REAL, DIMENSION(KLON)    :: ZCPH            ! specific heat C_ph 
REAL, DIMENSION(KLON)    :: ZLV, ZLS        ! latent heat of vaporis., sublim.       
REAL, DIMENSION(KLON)    :: ZURV            ! updraft water vapor at level k+1
REAL, DIMENSION(KLON)    :: ZPI             ! Pi=(P0/P)**(Rd/Cpd)  
REAL, DIMENSION(KLON)    :: ZTHEUL          ! theta_e for undilute ascent
REAL, DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5,   &
                            ZWORK6          ! work arrays
INTEGER, DIMENSION(KLON) :: IWORK           ! wok array
LOGICAL, DIMENSION(KLON) :: GWORK1, GWORK2, GWORK4, GWORK5 
                                            ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK6     ! work array
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
IIE = KLON
!
!
!*       1.     Initialize updraft properties and local variables
!               -------------------------------------------------
!
ZEPSA      = XRV / XRD 
ZCVOCD     = XCPV / XCPD
ZCPORD     = XCPD / XRD
ZRDOCP     = XRD / XCPD
!
PUMF(:,:)  = 0.
PUER(:,:)  = 0.
PUDR(:,:)  = 0.
PUTHL(:,:) = 0.
PUTHV(:,:) = 0.
PURW(:,:)  = 0.
PURC(:,:)  = 0.
PURI(:,:)  = 0.
ZUW1(:)    = PWLCL(:) * PWLCL(:)
ZUW2(:)    = 0.
ZE1(:)     = 0.
ZD1(:)     = 0.
PCAPE(:)   = 0.
KCTL(:)    = IKB
KETL(:)    = KLCL(:)
GWORK2(:)  = .TRUE.
GWORK5(:)  = .TRUE.
ZPI(:)     = 1.
ZWORK3(:)  = 0.
ZWORK4(:)  = 0.
ZWORK5(:)  = 0.
ZWORK6(:)  = 0.
GWORK1(:)  = .FALSE.
GWORK4(:)  = .FALSE.
!
!
!*       1.1    Compute undilute updraft theta_e for CAPE computations
!               Bolton (1980) formula.
!               Define accurate enthalpy for updraft
!               -----------------------------------------------------
!
ZTHEUL(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** ( 1. - 0.28 * PRVLCL(:) )  &
            * EXP( ( 3374.6525 / PTLCL(:) - 2.5403 ) *                        &
                                   PRVLCL(:) * ( 1. + 0.81 * PRVLCL(:) ) )
!
!
ZWORK1(:) = ( XCPD + PRVLCL(:) * XCPV ) * PTLCL(:)                            &
            + ( 1. + PRVLCL(:) ) * XG * PZLCL(:)
!
!
!*       2.     Set updraft properties between DPL and LCL
!               ------------------------------------------
!
JKP = MAXVAL( KLCL(:) )
JKM = MINVAL( KDPL(:) )
DO JK = JKM, JKP
   DO JI = 1, IIE
   IF ( JK >= KDPL(JI) .AND. JK < KLCL(JI) ) THEN
	PUMF(JI,JK)  = PMFLCL(JI)
	PUTHL(JI,JK) = ZWORK1(JI) 
	PUTHV(JI,JK) = PTHLCL(JI) * ( 1. + ZEPSA * PRVLCL(JI) ) /             &
                                  ( 1. + PRVLCL(JI) )
        PURW(JI,JK)  = PRVLCL(JI) 
   END IF
   END DO
END DO                        
!
!
!*       3.     Enter loop for updraft computations
!               ------------------------------------
!
JKMIN = MINVAL( KLCL(:) - 1 ) 
DO JK = MAX( IKB + 1, JKMIN ), IKE - 1
  ZWORK6(:) = 1.
  JKP = JK + 1  
!
  GWORK4(:) = JK >= KLCL(:) - 1 
  GWORK1(:) = GWORK4(:) .AND. GWORK2(:) ! this mask is used to confine
                           ! updraft computations between the LCL and the CTL
!                                                         
  WHERE( JK == KLCL(:) - 1 ) ZWORK6(:) = 0. ! factor that is used in buoyancy
                                        ! computation at first level above LCL
!
!
!*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1   
!               ----------------------------------------------------------
!
    ZWORK1(:) = PURC(:,JK) 
    ZWORK2(:) = PURI(:,JK) 
    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), PUTHL(:,JK), PURW(:,JK),&
                          ZWORK1, ZWORK2, PZ(:,JKP), GWORK1, ZUT, ZURV,     &
                          PURC(:,JKP), PURI(:,JKP), ZLV, ZLS, ZCPH )
!
!
  ZPI(:) = ( XP00 / PPRES(:,JKP) ) ** ZRDOCP   
  WHERE ( GWORK1(:) )
!
    PUTHV(:,JKP) = ZPI(:) * ZUT(:) * ( 1. + ZEPSA * ZURV(:) )           &  
                         / ( 1. + PURW(:,JK) )     
!
!
!*       5.     Compute square of vertical velocity using entrainment   
!               at level k
!               -----------------------------------------------------
!    
    ZWORK3(:) = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -         &
                     ( 1. - ZWORK6(:) ) * PZLCL(:)          ! level thickness  
    ZWORK4(:) = PTHV(:,JK) * ZWORK6(:) +                   &
                 ( 1. - ZWORK6(:) ) * PTHVELCL(:)
    ZWORK5(:) = 2. * ZUW1(:) * PUER(:,JK) / MAX( .1, PUMF(:,JK) )
    ZUW2(:)   = ZUW1(:) + ZWORK3(:) * XNHGAM * XG *        & 
                  ( ( PUTHV(:,JK) + PUTHV(:,JKP) ) /       &
                  ( ZWORK4(:) + PTHV(:,JKP) ) - 1. )       & ! buoyancy term
                - ZWORK5(:)                                  ! entrainment term
!
!
!*       6.     Update total precipitation: dr_r=(r_c+r_i)*exp(-rate*dz)  
!               --------------------------------------------------------
!
!                    compute level mean vertical velocity  
    ZWORK2(:)   = 0.5 *                                                    &
                       ( SQRT( MAX( 1.E-2, ZUW2(:) ) ) +                   &
                         SQRT( MAX( 1.E-2, ZUW1(:) ) ) )          
!
!
!*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation 
!               -------------------------------------------------------
!
    PURW(:,JKP)  = PURW(:,JK) 
    PURC(:,JKP)  = PURC(:,JKP)
    PURI(:,JKP)  = PURI(:,JKP)
    PUTHL(:,JKP) = PUTHL(:,JK)
!    
    ZUW1(:)      = ZUW2(:)       
!
  END WHERE
!
!
!*       8.     Compute entrainment and detrainment using conservative
!               variables adjusted for precipitation ( not for entrainment)
!               -----------------------------------------------------------
!
!*       8.1    Compute critical mixed fraction by estimating unknown  
!               T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
!               We determine the zero crossing of the linear curve
!               evaluating the derivative using ZMIXF=0.1.
!               -----------------------------------------------------
!    
    ZMIXF(:)  = 0.1   ! starting value for critical mixed fraction
    ZWORK1(:) = ZMIXF(:) * PTHL(:,JKP)                                     &
                     + ( 1. - ZMIXF(:) ) * PUTHL(:,JKP) ! mixed enthalpy
    ZWORK2(:) = ZMIXF(:) * PRW(:,JKP)                                      &
                     + ( 1. - ZMIXF(:) ) * PURW(:,JKP)  ! mixed r_w
!
    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), ZWORK1, ZWORK2,        &
                          PURC(:,JKP), PURI(:,JKP), PZ(:,JKP), GWORK1, ZUT,&
                          ZWORK3, ZWORK4, ZWORK5, ZLV, ZLS, ZCPH )
!        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)
!        
     ! compute theta_v of mixture
    ZWORK3(:) = ZUT(:) * ZPI(:) * ( 1. + ZEPSA * (                         &
                ZWORK2(:) - ZWORK4(:) - ZWORK5(:) ) ) / ( 1. + ZWORK2(:) )
     ! compute final value of critical mixed fraction using theta_v
     ! of mixture, grid-scale and updraft
    ZMIXF(:) = MAX( 0., PUTHV(:,JKP) - PTHV(:,JKP) ) * ZMIXF(:) /          &
                              ( PUTHV(:,JKP) - ZWORK3(:) + 1.E-10 )
    ZMIXF(:) = MAX( 0., MIN( 1., ZMIXF(:) ) )
!    
!
!*       8.2     Compute final midlevel values for entr. and detrainment    
!                after call of distribution function
!                -------------------------------------------------------
!    
!
    CALL CONVECT_MIXING_FUNCT ( KLON, ZMIXF, 1, ZE2, ZD2 )
!       Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates
!
! ZWORK1(:) = XENTR * PMFLCL(:) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  zwork1(:) = xentr * xg / xcrad * pumf(:,jk) * ( pz(:,jkp) - pz(:,jk) )
! ZWORK1(:) = XENTR * pumf(:,jk) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  ZWORK2(:) = 0.
  WHERE ( GWORK1(:) ) ZWORK2(:) = 1.
  WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) )
    PUER(:,JKP) = 0.5 * ZWORK1(:) * ( ZE1(:) + ZE2(:) ) * ZWORK2(:)
    PUDR(:,JKP) = 0.5 * ZWORK1(:) * ( ZD1(:) + ZD2(:) ) * ZWORK2(:)
  ELSEWHERE
    PUER(:,JKP) = 0.
    PUDR(:,JKP) = ZWORK1(:) * ZWORK2(:)
  END WHERE
!
!*       8.3     Determine equilibrium temperature level
!                --------------------------------------
!
   WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) .AND. JK > KLCL(:) + 1 &   
           .AND. GWORK1(:) )
         KETL(:) = JKP            ! equilibrium temperature level 
   END WHERE
!
!*       8.4     If the calculated detrained mass flux is greater than    
!                the total updraft mass flux, or vertical velocity is
!                negative, all cloud mass detrains at previous model level,
!                exit updraft calculations - CTL is attained
!                -------------------------------------------------------
!
  WHERE( GWORK1(:) )                                                   &
        GWORK2(:) = PUMF(:,JK) - PUDR(:,JKP) > 10. .AND. ZUW2(:) > 0.        
  WHERE ( GWORK2(:) ) KCTL(:) = JKP   ! cloud top level
  GWORK1(:) = GWORK2(:) .AND. GWORK4(:)
!
  IF ( COUNT( GWORK2(:) ) == 0 ) EXIT           
!
!
!*       9.   Compute CAPE for undilute ascent using theta_e and 
!             theta_es instead of theta_v. This estimation produces 
!             a significantly larger value for CAPE than the actual one.
!             ----------------------------------------------------------
!
  WHERE ( GWORK1(:) )
!
    ZWORK3(:)   = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -                      &
                  ( 1. - ZWORK6(:) ) *  PZLCL(:)              ! level thickness
    ZWORK2(:)   = PTHES(:,JK) + ( 1. - ZWORK6(:) ) *                      &
     ( PTHES(:,JKP) - PTHES(:,JK) ) / ( PZ(:,JKP) - PZ(:,JK) ) *          &
     ( PZLCL(:) - PZ(:,JK) ) ! linear interpolation for theta_es at LCL
                            ! ( this is only done for model level just above LCL
!
    ZWORK1(:) = ( 2. * ZTHEUL(:) ) / ( ZWORK2(:) + PTHES(:,JKP) ) - 1.   
    PCAPE(:)  = PCAPE(:) + XG * ZWORK3(:) * MAX( 0., ZWORK1(:) )
!
!
!*       10.   Compute final values of updraft mass flux, enthalpy, r_w 
!              at level k+1    
!              --------------------------------------------------------
!    
    PUMF(:,JKP)  = PUMF(:,JK) - PUDR(:,JKP) + PUER(:,JKP) 
    PUMF(:,JKP)  = MAX( PUMF(:,JKP), 0.1 )
    PUTHL(:,JKP) = ( PUMF(:,JK) * PUTHL(:,JK) +                              &
                     PUER(:,JKP) * PTHL(:,JK) - PUDR(:,JKP) * PUTHL(:,JK) )  &
                    / PUMF(:,JKP) 
    PURW(:,JKP)  = ( PUMF(:,JK) * PURW(:,JK) +                               &
                     PUER(:,JKP) * PRW(:,JK) - PUDR(:,JKP) * PURW(:,JK) )    &
                    / PUMF(:,JKP) 
!    
!
    ZE1(:) = ZE2(:) ! update fractional entrainment/detrainment
    ZD1(:) = ZD2(:)
!
  END WHERE
!
END DO
!
!*       12.1    Set OTRIG to False if cloud thickness < 0.5km
!                or > 3km (deep convection) or CAPE < 1
!                ------------------------------------------------
!
    DO JI = 1, IIE
          JK  = KCTL(JI)
          ZWORK1(JI) = PZ(JI,JK) - PZLCL(JI)
          OTRIG(JI) = ZWORK1(JI) >= XCDEPTH  .AND. ZWORK1(JI) < 3.E3        &
                     .AND. PCAPE(JI) > 1. 
    END DO
    WHERE( .NOT. OTRIG(:) )
          KCTL(:) = IKB 
    END WHERE
KETL(:) = MAX( KETL(:), KLCL(:) + 2 )
KETL(:) = MIN( KETL(:), KCTL(:) )
!
!
!*       12.2    If the ETL and CTL are the same detrain updraft mass   
!                flux at this level
!                ------------------------------------------------------- 
!
ZWORK1(:) = 0.
WHERE ( KETL(:) == KCTL(:) ) ZWORK1(:) = 1.
!
DO JI = 1, IIE
    JK = KETL(JI) 
    PUDR(JI,JK)   = PUDR(JI,JK) +                                    &
                          ( PUMF(JI,JK) - PUER(JI,JK) )  * ZWORK1(JI)  
    PUER(JI,JK)   = PUER(JI,JK) * ( 1. - ZWORK1(JI) )
    PUMF(JI,JK)   = PUMF(JI,JK) * ( 1. - ZWORK1(JI) )
    JKP = KCTL(JI) + 1
    PUER(JI,JKP)  = 0. ! entrainm/detr rates have been already computed
    PUDR(JI,JKP)  = 0. ! at level KCTL+1, set them to zero
END DO
!    
!*       12.3    Adjust mass flux profiles, detrainment rates, and   
!                precipitation fallout rates to reflect linear decrease
!                in mass flux between the ETL and CTL
!                -------------------------------------------------------        
! 
ZWORK1(:) = 0.
JK1 = MINVAL( KETL(:) )
JK2 = MAXVAL( KCTL(:) )

DO JK = JK1, JK2
    DO JI = 1, IIE
    IF( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
        ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
    END IF
    END DO
END DO
!
DO JI = 1, IIE
    JK = KETL(JI) 
    ZWORK1(JI) = PUMF(JI,JK) / MAX( 1., ZWORK1(JI) )
END DO
!
DO JK = JK1 + 1, JK2
    JKP = JK - 1
    DO JI = 1, IIE
    IF ( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
        PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
        PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
    END IF
    END DO
END DO
!
!         12.4   Set mass flux and entrainment in the source layer.
!                Linear increase throughout the source layer.
!                -------------------------------------------------------
!
!IWORK(:) = MIN( KPBL(:), KLCL(:) - 1 )
IWORK(:) = KPBL(:)
DO JI = 1, IIE
     JK  = KDPL(JI)
     JKP = IWORK(JI)
!          mixed layer depth
     ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP) + PDPRES(JI,JK)
END DO
!
JKP = MAXVAL( IWORK(:) )
DO JK = JKM, JKP
   DO JI = 1, IIE
   IF ( JK >= KDPL(JI)  .AND. JK <= IWORK(JI) ) THEN
       PUER(JI,JK) = PUER(JI,JK) + PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
       PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
   END IF
   END DO
END DO
!
!
!*       13.   If cloud thickness is smaller than  .5 km or > 3 km
!              no shallow convection is allowed
!              Nota: For technical reasons, we stop the convection
!                    computations in this case and do not go back to
!                    TRIGGER_FUNCT to look for the next unstable LCL
!                    which could produce a thicker cloud.
!              ---------------------------------------------------
!
GWORK6(:,:) = SPREAD( OTRIG(:), DIM=2, NCOPIES=KLEV )
WHERE ( .NOT. GWORK6(:,:) )
    PUMF(:,:)  = 0.
    PUDR(:,:)  = 0.
    PUER(:,:)  = 0.
    PUTHL(:,:) = PTHL(:,:)
    PURW(:,:)  = PRW(:,:)
    PURC(:,:)  = 0.
    PURI(:,:)  = 0.
END WHERE
!
END SUBROUTINE CONVECT_UPDRAFT_SHAL
!     ######spl
     SUBROUTINE CONVECT_CLOSURE_SHAL( KLON, KLEV,                                 &
                                      PPRES, PDPRES, PZ, PDXDY, PLMASS,           &
                                      PTHL, PTH, PRW, PRC, PRI, OTRIG1,           &
                                      PTHC, PRWC, PRCC, PRIC, PWSUB,              &
                                      KLCL, KDPL, KPBL, KCTL,                     &
                                      PUMF, PUER, PUDR, PUTHL, PURW,              &
                                      PURC, PURI, PCAPE, PTIMEC, KFTSTEPS         )
!    #######################################################################
!
!!**** Uses modified Fritsch-Chappell closure
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted 
!!     (over a time step PTIMEC) environmental values of THETA_l, R_w, R_c, R_i
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PTHC-PTH)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     
!!    CONVECT_CLOSURE_THRVLCL
!!    CONVECT_CLOSURE_ADJUST_SHAL
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV         ! specific heat for dry air and water vapor
!!          XCL, XCI           ! specific heat for liquid water and ice
!!          XTT                ! triple point temperature
!!          XLVTT, XLSTT       ! vaporization, sublimation heat constant
!!
!!      Module MODD_CONVPAR_SHAL
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration 
!!          XSTABC             ! stability factor in CAPE adjustment
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Peter Bechtold 15/11/96 change for enthalpie, r_c + r_i tendencies
!!      Tony Dore   14/10/96 Initialise local variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR_SHAL
USE MODD_CONVPAREXT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                   INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,                   INTENT(IN) :: KLEV   ! vertical dimension
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level 
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
REAL, DIMENSION(KLON),  INTENT(INOUT) :: PTIMEC ! convection time step 
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY  ! grid area (m^2)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH    ! grid scale theta        
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRW    ! grid scale total water  
			                        ! mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRC    ! grid scale r_c 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRI    ! grid scale r_i 
LOGICAL, DIMENSION(KLON),  INTENT(IN) :: OTRIG1 ! logical to keep trace of 
                                                ! convective arrays modified in UPDRAFT
!   
!
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (P)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES ! pressure difference between 
                                                 ! bottom and top of layer (Pa)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PLMASS ! mass of model layer (kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m) 
REAL, DIMENSION(KLON),     INTENT(IN)  :: PCAPE  ! available potent. energy
INTEGER,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                 ! only used for chemical tracers
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
REAL, DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)
!
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i 
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JK, JKP, JKMAX ! vertical loop index
INTEGER :: JI             ! horizontal loop index
INTEGER :: JITER          ! iteration loop index
INTEGER :: JSTEP          ! fractional time loop index
REAL    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
REAL    :: ZCVOCD, ZEPSA  ! C_pv / C_pd, R_v / R_d
!
REAL, DIMENSION(KLON,KLEV) :: ZTHLC       ! convectively adjusted 
                                          ! grid scale enthalpy
REAL, DIMENSION(KLON,KLEV) :: ZOMG        ! conv. environm. subsidence (Pa/s)
REAL, DIMENSION(KLON,KLEV) :: ZUMF        ! non-adjusted updraft mass flux
REAL, DIMENSION(KLON,KLEV) :: ZUER        !   "     updraft  entrainm. rate
REAL, DIMENSION(KLON,KLEV) :: ZUDR        !   "     updraft  detrainm. rate
REAL, DIMENSION(KLON)     :: ZADJ         ! mass adjustment factor
REAL, DIMENSION(KLON)     :: ZADJMAX      ! limit value for ZADJ
REAL, DIMENSION(KLON)     :: ZCAPE        ! new CAPE after adjustment
REAL, DIMENSION(KLON)     :: ZTIMEC       ! fractional convective time step
REAL, DIMENSION(KLON,KLEV):: ZTIMC        ! 2D work array for ZTIMEC
!
REAL, DIMENSION(KLON)     :: ZTHLCL       ! new  theta at LCL
REAL, DIMENSION(KLON)     :: ZRVLCL       ! new  r_v at LCL
REAL, DIMENSION(KLON)     :: ZZLCL        ! height of LCL
REAL, DIMENSION(KLON)     :: ZTLCL        ! temperature at LCL
REAL, DIMENSION(KLON)     :: ZTELCL       ! envir. temper. at LCL
REAL, DIMENSION(KLON)     :: ZTHEUL       ! theta_e for undilute ascent
REAL, DIMENSION(KLON)     :: ZTHES1, ZTHES2! saturation environm. theta_e
REAL, DIMENSION(KLON,KLEV) :: ZTHMFIN, ZTHMFOUT, ZRWMFIN, ZRWMFOUT
REAL, DIMENSION(KLON,KLEV) :: ZRCMFIN, ZRCMFOUT, ZRIMFIN, ZRIMFOUT
                                    ! work arrays for environm. compensat. mass flux
REAL, DIMENSION(KLON)     :: ZPI          ! (P/P00)**R_d/C_pd 
REAL, DIMENSION(KLON)     :: ZLV          ! latent heat of vaporisation
REAL, DIMENSION(KLON)     :: ZLS          ! latent heat of sublimation 
REAL, DIMENSION(KLON)     :: ZCPH         ! specific heat C_ph
INTEGER, DIMENSION(KLON)  :: ITSTEP       ! fractional convective time step
INTEGER, DIMENSION(KLON)  :: ICOUNT       ! timestep counter 
INTEGER, DIMENSION(KLON)  :: ILCL         ! index lifting condens. level
INTEGER, DIMENSION(KLON)  :: IWORK1       ! work array
REAL, DIMENSION(KLON)     :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5
LOGICAL, DIMENSION(KLON)  :: GWORK1, GWORK3! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK4    ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.2    Initialize  local variables
!               ----------------------------
!
!
ZTIMC(:,:) = 0.
ZTHES2(:) = 0.
ZWORK1(:) = 0. 
ZWORK2(:) = 0. 
ZWORK3(:) = 0. 
ZWORK4(:) = 0. 
ZWORK5(:) = 0. 
GWORK1(:) = .FALSE.
GWORK3(:) = .FALSE.  
GWORK4(:,:) = .FALSE.  
ILCL(:)   = KLCL(:)
!
ZCPORD    = XCPD / XRD
ZRDOCP    = XRD / XCPD
ZCVOCD    = XCPV / XCPD 
ZEPSA     = XRV / XRD
!
ZADJ(:)   = 1.
ZWORK5(:) = 1.
WHERE( .NOT. OTRIG1(:) ) ZWORK5(:) = 0. 
!
!
!*       0.3   Compute loop bounds
!              ------------------- 
!
IIE    = KLON
IKB    = 1 + JCVEXB 
IKS    = KLEV
IKE    = KLEV - JCVEXT 
JKMAX  = MAXVAL( KCTL(:) )
!
!
!*       2.     Save initial mass flux values to be used in adjustment procedure
!               ---------------------------------------------------------------
!
ZUMF(:,:)  = PUMF(:,:)
ZUER(:,:)  = PUER(:,:)
ZUDR(:,:)  = PUDR(:,:)
ZOMG(:,:)  = 0.
PWSUB(:,:) = 0. 
!
!
!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger 
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------
!
ZADJMAX(:) = 1000.
IWORK1(:) = ILCL(:)
JKP = MINVAL( KDPL(:) )
DO JK = JKP, IKE
  DO JI = 1, IIE
    IF( JK > KDPL(JI) .AND. JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) / ( ( PUER(JI,JK) + 1.E-5 ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    END IF
  END DO
END DO
!
!
GWORK1(:) = OTRIG1(:)  ! logical array to limit adjustment to not definitively
                       ! adjusted columns
!
DO JK = IKB, IKE
  ZTHLC(:,JK) = PTHL(:,JK) ! initialize adjusted envir. values 
  PRWC(:,JK)  = PRW(:,JK)
  PRCC(:,JK)  = PRC(:,JK)
  PRIC(:,JK)  = PRI(:,JK)
  PTHC(:,JK)  = PTH(:,JK)
END DO
!
!
!
DO JITER = 1, 7  ! Enter adjustment loop to assure that all CAPE is
                 ! removed within the advective time interval TIMEC
!
     ZTIMEC(:) = PTIMEC(:)
     GWORK4(:,:)   = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKS )
     WHERE( GWORK4(:,:) ) PWSUB(:,:) = 0.
     ZOMG(:,:)=0.
!
     DO JK = IKB + 1, JKMAX 
           JKP = MAX( IKB + 1, JK - 1 )
           WHERE ( GWORK1(:) .AND. JK <= KCTL(:) )
!
!
!*       4.     Determine vertical velocity at top and bottom of each layer
!               to satisfy mass continuity.
!               ---------------------------------------------------------------
              ! we compute here Domega/Dp = - g rho Dw/Dz = 1/Dt
!
             ZWORK1(:)   = - ( PUER(:,JKP) - PUDR(:,JKP) ) / PLMASS(:,JKP)
!    
             PWSUB(:,JK) = PWSUB(:,JKP) - PDPRES(:,JK-1) * ZWORK1(:)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer
!
!   
!*       5.     Compute fractional time step. For stability or 
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------
!
             ZWORK1(:) = XSTABT * PDPRES(:,JKP) / ( ABS( PWSUB(:,JK) ) + 1.E-10 )
              ! the factor XSTABT is used for stability reasons
             ZTIMEC(:) = MIN( ZTIMEC(:), ZWORK1(:) ) 
!
              ! transform vertical velocity in mass flux units
             ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / XG 
         END WHERE
     END DO
!
!
     WHERE( GWORK4(:,:) )
           ZTHLC(:,:) = PTHL(:,:) ! reinitialize adjusted envir. values 
           PRWC(:,:)  = PRW(:,:)  ! when iteration criterium not attained
           PRCC(:,:)  = PRC(:,:)
           PRIC(:,:)  = PRI(:,:)
           PTHC(:,:)  = PTH(:,:)
     END WHERE
!
! 
!        6. Check for mass conservation, i.e. ZWORK1 > 1.E-2
!           If mass is not conserved, the convective tendencies
!           automatically become zero.
!           ----------------------------------------------------
!
    DO JI = 1, IIE
       JK=KCTL(JI)
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1 )    &
                                                            - PWSUB(JI,JK)
    END DO
    WHERE( GWORK1(:) .AND. ABS( ZWORK1(:) ) - .01 > 0. )
        GWORK1(:) = .FALSE.
        PTIMEC(:) = 1.E-1
        ZWORK5(:) = 0.
    END WHERE
    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
    END DO
    GWORK4(:,1:IKB) = .FALSE.
    GWORK4(:,IKS)   = .FALSE.
!
    ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1 
    ZTIMEC(:) = PTIMEC(:) / REAL( ITSTEP(:) ) ! adjust  fractional time step
                                           ! to be an integer multiple of PTIMEC
    ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
    ICOUNT(:) = 0
!
!
!
    KFTSTEPS = MAXVAL( ITSTEP(:) )
    DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop here
!
     	 ICOUNT(:) = ICOUNT(:) + 1
!
	     GWORK3(:) =  ITSTEP(:) >= ICOUNT(:) .AND. GWORK1(:) 
!
!
!*       7.     Assign enthalpy and r_w values at the top and bottom of each
!               layer based on the sign of w
!               ------------------------------------------------------------
!
             ZTHMFIN(:,:)   = 0.
             ZRWMFIN(:,:)   = 0.
             ZRCMFIN(:,:)   = 0.
             ZRIMFIN(:,:)   = 0.
             ZTHMFOUT(:,:)  = 0.
             ZRWMFOUT(:,:)  = 0.
             ZRCMFOUT(:,:)  = 0.
             ZRIMFOUT(:,:)  = 0.
!
         DO JK = IKB + 1, JKMAX
         GWORK4(:,JK) = GWORK3(:) .AND. JK <= KCTL(:)
         JKP = MAX( IKB + 1, JK - 1 )
           DO JI = 1, IIE
           IF ( GWORK3(JI) ) THEN
!
               ZWORK1(JI)       = SIGN( 1., ZOMG(JI,JK) )
               ZWORK2(JI)       = 0.5 * ( 1. + ZWORK1(JI) )
               ZWORK1(JI)       = 0.5 * ( 1. - ZWORK1(JI) )
               ZTHMFIN(JI,JK)   = - ZOMG(JI,JK) * ZTHLC(JI,JKP) * ZWORK1(JI)
               ZTHMFOUT(JI,JK)  =   ZOMG(JI,JK) * ZTHLC(JI,JK)  * ZWORK2(JI)
               ZTHMFIN(JI,JKP)  = ZTHMFIN(JI,JKP)  + ZTHMFOUT(JI,JK) * ZWORK2(JI)
               ZTHMFOUT(JI,JKP) = ZTHMFOUT(JI,JKP) + ZTHMFIN(JI,JK)  * ZWORK1(JI)
               ZRWMFIN(JI,JK)   = - ZOMG(JI,JK) * PRWC(JI,JKP) * ZWORK1(JI)
               ZRWMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRWC(JI,JK)  * ZWORK2(JI)
               ZRWMFIN(JI,JKP)  = ZRWMFIN(JI,JKP)  + ZRWMFOUT(JI,JK) * ZWORK2(JI)
               ZRWMFOUT(JI,JKP) = ZRWMFOUT(JI,JKP) + ZRWMFIN(JI,JK)  * ZWORK1(JI)
               ZRCMFIN(JI,JK)   = - ZOMG(JI,JK) * PRCC(JI,JKP) * ZWORK1(JI)
               ZRCMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRCC(JI,JK)  * ZWORK2(JI)
               ZRCMFIN(JI,JKP)  = ZRCMFIN(JI,JKP)  + ZRCMFOUT(JI,JK) * ZWORK2(JI)
               ZRCMFOUT(JI,JKP) = ZRCMFOUT(JI,JKP) + ZRCMFIN(JI,JK)  * ZWORK1(JI)
               ZRIMFIN(JI,JK)   = - ZOMG(JI,JK) * PRIC(JI,JKP) * ZWORK1(JI)
               ZRIMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRIC(JI,JK)  * ZWORK2(JI)
               ZRIMFIN(JI,JKP)  = ZRIMFIN(JI,JKP)  + ZRIMFOUT(JI,JK) * ZWORK2(JI)
               ZRIMFOUT(JI,JKP) = ZRIMFOUT(JI,JKP) + ZRIMFIN(JI,JK)  * ZWORK1(JI)
!
           END IF
           END DO
         END DO
!
         WHERE ( GWORK4(:,:) )
!
!******************************************************************************
!
!*       8.     Update the environmental values of enthalpy and r_w at each level
!               NOTA: These are the MAIN EQUATIONS of the scheme
!               -----------------------------------------------------------------
!
!
           ZTHLC(:,:) = ZTHLC(:,:) + ZTIMC(:,:) / PLMASS(:,:) * (      &
                          ZTHMFIN(:,:) + PUDR(:,:) * PUTHL(:,:)        &
                      - ZTHMFOUT(:,:) - PUER(:,:) * PTHL(:,:)   )
           PRWC(:,:)  = PRWC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
                          ZRWMFIN(:,:) + PUDR(:,:) * PURW(:,:)          &
                      - ZRWMFOUT(:,:) - PUER(:,:) * PRW(:,:)    )    
           PRCC(:,:)  = PRCC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
               ZRCMFIN(:,:) + PUDR(:,:) * PURC(:,:) - ZRCMFOUT(:,:) -  &
                         PUER(:,:) * PRC(:,:)    )    
           PRIC(:,:)  = PRIC(:,:) + ZTIMC(:,:) / PLMASS(:,:) *  (      &
               ZRIMFIN(:,:) + PUDR(:,:) * PURI(:,:) - ZRIMFOUT(:,:) -  & 
                         PUER(:,:) * PRI(:,:)    )    
!
!
!******************************************************************************
!
         END WHERE
!
    END DO ! Exit the fractional time step loop
!
! 
!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------
!
      DO JK = IKB + 1, JKMAX
         DO JI = 1, IIE
         IF( GWORK1(JI) .AND. JK <= KCTL(JI) ) THEN
           ZPI(JI)    = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
           ZCPH(JI)   = XCPD + PRWC(JI,JK) * XCPV
           ZWORK2(JI) = PTH(JI,JK) / ZPI(JI)  ! first temperature estimate
           ZLV(JI)    = XLVTT + ( XCPV - XCL ) * ( ZWORK2(JI) - XTT )
           ZLS(JI)    = XLVTT + ( XCPV - XCI ) * ( ZWORK2(JI) - XTT )
             ! final linearized temperature
           ZWORK2(JI) = ( ZTHLC(JI,JK) + ZLV(JI) * PRCC(JI,JK) + ZLS(JI) * PRIC(JI,JK) &
                       - (1. + PRWC(JI,JK) ) * XG * PZ(JI,JK) ) / ZCPH(JI)
           ZWORK2(JI) = MAX( 180., MIN( 340., ZWORK2(JI) ) )
           PTHC(JI,JK)= ZWORK2(JI) * ZPI(JI) ! final adjusted envir. theta
         END IF
         END DO
      END DO
!
!
!*         11.     Compute new cloud ( properties at new LCL )
!                     NOTA: The computations are very close to
!                           that in routine TRIGGER_FUNCT
!                  ---------------------------------------------
!
      CALL CONVECT_CLOSURE_THRVLCL(  KLON, KLEV,                           &
                                     PPRES, PTHC, PRWC, PZ, GWORK1,        &
                                     ZTHLCL, ZRVLCL, ZZLCL, ZTLCL, ZTELCL, &
                                     ILCL, KDPL, KPBL )
!
!
       ZTLCL(:)  = MAX( 230., MIN( 335., ZTLCL(:) ) )  ! set some overflow bounds
       ZTELCL(:) = MAX( 230., MIN( 335., ZTELCL(:) ) )
       ZTHLCL(:) = MAX( 230., MIN( 345., ZTHLCL(:) ) )
       ZRVLCL(:) = MAX(   0., MIN(   1., ZRVLCL(:) ) )
!
!
!*         12.    Compute adjusted CAPE
!                 ---------------------
!
       ZCAPE(:)  = 0.
       ZPI(:)    = ZTHLCL(:) / ZTLCL(:)
       ZPI(:)    = MAX( 0.95, MIN( 1.5, ZPI(:) ) )
       ZWORK1(:) = XP00 / ZPI(:) ** ZCPORD ! pressure at LCL
!
       CALL CONVECT_SATMIXRATIO( KLON, ZWORK1, ZTELCL, ZWORK3, ZLV, ZLS, ZCPH )
       ZWORK3(:) = MIN(   .1, MAX(   0., ZWORK3(:) ) )
!
	        ! compute theta_e updraft undilute
       ZTHEUL(:) = ZTLCL(:) * ZPI(:) ** ( 1. - 0.28 * ZRVLCL(:) )            &
                                  * EXP( ( 3374.6525 / ZTLCL(:) - 2.5403 )   &
                                  * ZRVLCL(:) * ( 1. + 0.81 * ZRVLCL(:) ) )
!
	        ! compute theta_e saturated environment at LCL
       ZTHES1(:) = ZTELCL(:) * ZPI(:) ** ( 1. - 0.28 * ZWORK3(:) )           &
                                  * EXP( ( 3374.6525 / ZTELCL(:) - 2.5403 )  &
                                  * ZWORK3(:) * ( 1. + 0.81 * ZWORK3(:) ) )
!
      DO JK = MINVAL( ILCL(:) ), JKMAX
        JKP = JK - 1
        DO JI = 1, IIE
          ZWORK4(JI) = 1.
          IF ( JK == ILCL(JI) ) ZWORK4(JI) = 0.
!
           ! compute theta_e saturated environment and adjusted values
           ! of theta
!
          GWORK3(JI)  = JK >= ILCL(JI) .AND. JK <= KCTL(JI) .AND. GWORK1(JI) 
!
          ZPI(JI)     = ( XP00 / PPRES(JI,JK) ) ** ZRDOCP
          ZWORK2(JI)  = PTHC(JI,JK) / ZPI(JI)
        END DO
!
        CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZWORK2, ZWORK3, ZLV, ZLS, ZCPH )
!
!
        DO JI = 1, IIE
          IF ( GWORK3(JI) ) THEN
              ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( 1. - 0.28 * ZWORK3(JI) )   &
                                   * EXP( ( 3374.6525 / ZWORK2(JI) - 2.5403 ) &
                                   * ZWORK3(JI) * ( 1. + 0.81 * ZWORK3(JI) ) )
!
              ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                &
                           ( 1. - ZWORK4(JI) ) * ZZLCL(JI)    ! level thickness
              ZWORK1(JI)  = ( 2. * ZTHEUL(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - 1.
              ZCAPE(JI)   = ZCAPE(JI) + XG * ZWORK3(JI) * MAX( 0., ZWORK1(JI) )
              ZTHES1(JI)  = ZTHES2(JI)
          END IF
        END DO
      END DO
!
!                                                          
!*         13.     Determine mass adjustment factor knowing how much
!                  CAPE has been removed.
!                  -------------------------------------------------
!
       WHERE ( GWORK1(:) )
           ZWORK1(:) = MAX( PCAPE(:) - ZCAPE(:), 0.1 * PCAPE(:) )
           ZWORK2(:) = ZCAPE(:) / ( PCAPE(:) + 1.E-8 )
!       
           GWORK1(:) = ZWORK2(:) > 0.1 .OR. ZCAPE(:) == 0. ! mask for adjustment
       END WHERE
!
       WHERE ( ZCAPE(:) == 0. .AND. GWORK1(:) )  ZADJ(:) = ZADJ(:) * 0.5
       WHERE ( ZCAPE(:) /= 0. .AND. GWORK1(:) )                              &
               ZADJ(:) = ZADJ(:) * XSTABC * PCAPE(:) / ( ZWORK1(:) + 1.E-8 )
       ZADJ(:) = MIN( ZADJ(:), ZADJMAX(:) )  
!
!
!*         13.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization
!                 ----------------------------------------------------
!
       CALL CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, ZADJ,                     &
                                         PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR    )
!
!
      IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached 
                                          ! desired degree of stabilization.
!
END DO  ! end of big adjustment iteration loop
!
!
        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
   PRWC(:,JK) = MAX( 0., PRWC(:,JK) - PRCC(:,JK) - PRIC(:,JK) )
END DO
!
!
END SUBROUTINE CONVECT_CLOSURE_SHAL
!     ######spl
     SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, PADJ,                      &
                                             PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR  )
!    #########################################################################
!
!!**** Uses closure adjustment factor to adjust mass flux and to modify
!!     precipitation efficiency  when necessary. The computations are
!!     similar to routine CONVECT_PRECIP_ADJUST.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust the mass flux using the
!!      factor PADJ computed in CONVECT_CLOSURE
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!
!!    EXTERNAL
!!    --------
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!     
!!    None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    None
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
REAL, DIMENSION(KLON),      INTENT(IN) :: PADJ     ! mass adjustment factor
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IIE, IKB, IKE                 ! horiz. + vert. loop bounds
INTEGER :: JK                            ! vertical loop index
!
!
!-------------------------------------------------------------------------------
!
!*       0.3   Compute loop bounds
!              -------------------
!
IIE  = KLON
IKB  = 1 + JCVEXB 
IKE  = KLEV - JCVEXT
!
!
!*       1.     Adjust mass flux by the factor PADJ to converge to
!               specified degree of stabilization
!               ----------------------------------------------------
!
     DO JK = IKB + 1, IKE
	  PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
     END DO
!
END SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL
