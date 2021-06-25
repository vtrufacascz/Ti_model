      SUBROUTINE IONTIF(CRD,PF107Y,INVDIP,FL,DIMO,B0,
     &                  DIPL,MLT,ALT,DDD,PF107,TI,SIGTI)
C Empirical model of ion temperature (Ti) in the topside ionosphere
C with inclusion of solar activity Ti variation.
C Based on spherical harmonics approximation of measured
C Ti (all available satellites) at altitudes centred on 350km, 430km,
C 600km, and 850km. For intermediate altitudes a Booker
C interpolation is used. Recommended altitude range: 300-1000 km!!!
C Linear extrapolation is used for altitude ranges <300;350)km
C and (2000;2500> km. For days between seasons centred at
C (21.3. = 80; 21.6. = 172; 23.9. 266; 21.12. = 355) Ti is
C interpolated by harmonic functions.
C Inputs: CRD - 0 .. INVDIP
C               1 .. FL, DIMO, B0, DIPL (used for calculation INVDIP inside)
C         PF107Y - 0 .. PF107 dependency NOT included (not recommended!!!)
C                  1 .. PF107 dependency included
C         INVDIP - "mix" coordinate of the dip latitude and of
C                    the invariant latitude;
C                    positive northward, in deg, range <-90.0;90.0>
C         FL, DIMO, BO - McIlwain L parameter, dipole moment in
C                        Gauss, magnetic field strength in Gauss -
C                        parameters needed for invariant latitude
C                        calculation
C         DIPL - dip latitude
C                positive northward, in deg, range <-90.0;90.0>
C         MLT - magnetic local time (central dipole)
C               in hours, range <0;24)
C         ALT - altitude above the Earth's surface;
C               in km, range <500;3000>
C         DDD - day of year; range <0;365>
C         PF107 - Phil Richard's solar radio flux;
C Output: TI - ion temperature in K
C         SIGTI - standard deviation (or model error) of calculated TI in K
C Versions:
C           2021 (FORTRAN)  
C Authors of the model
C                V. Truhlik at al.
C Author of the code:
C         Vladimir Truhlik
C         Institute of Atm. Phys.
C         Bocni II, 1401/1a
C         14100 Praha 4, Sporilov
C         Czech Republic
C         e-mail: vtr@ufa.cas.cz
      REAL INVDIP,FL,DIMO,B0,DIPL,MLT,ALT,PF107,TE,SIGTE
      INTEGER CRD,PF107Y,ISRSAT,DDD,SEZDAY,XDAY
      INTEGER MIRREQ(81)
      REAL D(4,3,81),DERRTI(4,3,81),ASOL(4,3,81),BSOL(4,3,81)
      REAL ASOL2(4,3,81),BSOL2(4,3,81),CSOL2(4,3,81)
      DOUBLE PRECISION B(8),A
      REAL DPI,DTOR,ASA,INVL,RINVL,INVDP,RDIPL,ALFA,BETA
      REAL RMLT,RCOLAT
      REAL C(82)
      INTEGER SEZA,SEZB,DDDA,DDDB,DDDD
      REAL T350,T350A,T350B,T430,T430A,T430B,T600,T600A,T600B,
     &     T850,T850A,T850B
      REAL E350,E350A,E350B,E430,E430A,E430B,E600,E600A,E600B,
     &     E850,E850A,E850B
      REAL TP350A,TP350B,TP430A,TP430B,TP600A,TP600B,
     &     TP850A,TP850B
      REAL ANO(4),AH(4),DNO(2),ST(3)
      INTEGER FUN
      INTEGER I
      COMMON/ARGEXP/ARGMAX
      DATA B/1.259921D0  ,-0.1984259D0 ,-0.04686632D0,-0.01314096D0,
     &      -0.00308824D0, 0.00082777D0,-0.00105877D0, 0.00183142D0/
C////////////////////////////////coefficients - main model part//////////////////////
      DATA (MIRREQ(J),J=1,81)/
     &  1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1,
     &  1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,
     & -1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1,-1, 1, 1,-1,
     &  1,-1, 1,-1, 1,-1, 1,-1, 1, 1,-1, 1, 1,-1, 1,-1, 1, 1/
C      ISRSAT=1
      CALL KOFDTI(MIRREQ,D)
      CALL KERRTI(MIRREQ,DERRTI)
      CALL KOL107(MIRREQ,ASOL,BSOL)
      CALL KOQ107(MIRREQ,ASOL2,BSOL2,CSOL2)
C//////////////////////thresholds of solar activity/////////////////////////////////////
      IF (PF107 .GT. 250) PF107=250
      IF (PF107 .LT. 65) PF107=65
C////////////////////////////////////////////////////////////////////////////////////
      DPI=3.1415926535897
      DTOR=DPI/180.0
      IF (CRD .EQ. 1) THEN
C      calculation of INVDIP from FL, DIMO, BO, and DIPL
C      invariant latitude calculated by highly
C      accurate polynomial expansion
       A=(DIMO/B0)**(1.0D0/3.0D0)/FL
       ASA=A*(B(1)+B(2)*A+B(3)*A**2+B(4)*A**3+B(5)*A**4+
     &        B(6)*A**5+B(7)*A**6+B(8)*A**7)
       IF (ASA .GT. 1.0) ASA=1.0
C      invariant latitude (absolute value)
       RINVL=ACOS(SQRT(ASA))
       INVL=RINVL/DTOR
       RDIPL=DIPL*DTOR
       ALFA=SIN(ABS(RDIPL))**3
       BETA=COS(RINVL)**3
       INVDP=(ALFA*SIGN(1.0,DIPL)*INVL+BETA*DIPL)/(ALFA+BETA)
      ELSE IF	(CRD .EQ. 0) THEN
       INVDP=INVDIP
      ELSE
       RETURN
      END IF
      RMLT=MLT*DTOR*15.0
      RCOLAT=(90.0-INVDP)*DTOR
      CALL SPHARM_IK(C,8,8,RCOLAT,RMLT)
C     21.3. - 20.6.
      IF ((DDD .GE. 80) .AND. (DDD .LT. 172)) THEN
       SEZA=1
       SEZB=2
       DDDA=80
       DDDB=172
       DDDD=DDD
       FUN=0
      END IF
C     21.6. - 22.9.
      IF ((DDD .GE. 172) .AND. (DDD .LT. 266)) THEN
       SEZA=2
       SEZB=1
       DDDA=172
       DDDB=266
       DDDD=DDD
       FUN=1
      END IF
C     23.9. - 20.12.
      IF ((DDD .GE. 266) .AND. (DDD .LT. 355)) THEN
       SEZA=1
       SEZB=3
       DDDA=266
       DDDB=355
       DDDD=DDD
       FUN=0
      END IF
C     21.12. - 20.3.
      IF ((DDD .GE. 355) .OR. (DDD .LT. 80)) THEN
       SEZA=3
       SEZB=1
       DDDA=355
       DDDB=365+80
       DDDD=DDD
        IF (DDD .GE. 355) THEN
         DDDD=DDD
        ELSE
         DDDD=DDD+365
        END IF
       FUN=1 
      END IF
C     model Te
      T350A=0.0
      T350B=0.0
      T430A=0.0
      T430B=0.0
      T600A=0.0
      T600B=0.0
      T850A=0.0
      T850B=0.0          
      DO 30 I=1,81
       T350A=T350A+C(I)*D(1,SEZA,I)
       T350B=T350B+C(I)*D(1,SEZB,I)
       T430A=T430A+C(I)*D(2,SEZA,I)
       T430B=T430B+C(I)*D(2,SEZB,I)
       T600A=T600A+C(I)*D(3,SEZA,I)
       T600B=T600B+C(I)*D(3,SEZB,I)
       T850A=T850A+C(I)*D(4,SEZA,I)
30     T850B=T850B+C(I)*D(4,SEZB,I)
C     model errTe
      E350A=0.0
      E350B=0.0
      E430A=0.0
      E430B=0.0
      E600A=0.0
      E600B=0.0
      E850A=0.0
      E850B=0.0          
      DO 50 I=1,81
       E350A=E350A+C(I)*DERRTI(1,SEZA,I)
       E350B=E350B+C(I)*DERRTI(1,SEZB,I)
       E430A=E430A+C(I)*DERRTI(2,SEZA,I)
       E430B=E430B+C(I)*DERRTI(2,SEZB,I)
       E600A=E600A+C(I)*DERRTI(3,SEZA,I)
       E600B=E600B+C(I)*DERRTI(3,SEZB,I)
       E850A=E850A+C(I)*DERRTI(4,SEZA,I)
50     E850B=E850B+C(I)*DERRTI(4,SEZB,I)
C
      IF (PF107Y .EQ. 1) THEN
       CALL TIF107(PF107,INVDP,C,SEZA,SEZB,
     &                  ASOL,BSOL,ASOL2,BSOL2,CSOL2,  
     &                  TP350A,TP350B,TP430A,TP430B,
     &                  TP600A,TP600B,TP850A,TP850B) 
       T350A=T350A+TP350A
       T350B=T350B+TP350B
       T430A=T430A+TP430A
       T430B=T430B+TP430B
       T600A=T600A+TP600A
       T600B=T600B+TP600B
       T850A=T850A+TP850A
       T850B=T850B+TP850B
      END IF 
C     Ti
      IF (FUN .EQ. 0) THEN
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       SINDAY=SIN(DPI/2.0*XDAY/SEZDAY)
       T350=(T350B-T350A)*SINDAY+T350A
       T430=(T430B-T430A)*SINDAY+T430A
       T600=(T600B-T600A)*SINDAY+T600A
       T850=(T850B-T850A)*SINDAY+T850A
      ELSE
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       COSDAY=COS(DPI/2.0*XDAY/SEZDAY)
       T350=(T350A-T350B)*COSDAY+T350B
       T430=(T430A-T430B)*COSDAY+T430B
       T600=(T600A-T600B)*COSDAY+T600B
       T850=(T850A-T850B)*COSDAY+T850B
      END IF
C     error Ti
      IF (FUN .EQ. 0) THEN
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       SINDAY=SIN(DPI/2.0*XDAY/SEZDAY)
       E350=(E350B-E350A)*SINDAY+E350A
       E430=(E430B-E430A)*SINDAY+E430A
       E600=(E600B-E600A)*SINDAY+E600A
       E850=(E850B-E850A)*SINDAY+E850A
      ELSE
       SEZDAY=(DDDB-DDDA)
       XDAY=DDDD-DDDA
       COSDAY=COS(DPI/2.0*XDAY/SEZDAY)
       E350=(E350A-E350B)*COSDAY+E350B
       E430=(E430A-E430B)*COSDAY+E430B
       E600=(E600A-E600B)*COSDAY+E600B
       E850=(E850A-E850B)*COSDAY+E850B
      END IF      
C ////////////////////////////////////////////////////////
C     Ti linear interpolation in altitude
      IF (ALT .LT. 430.0) THEN
          TI=(T430-T350)/80.0*(ALT-350.0)+T350
       SIGTI=(E430-E350)/80.0*(ALT-350.0)+E350    
      END IF
      IF ((ALT .GE. 430.0) .AND. (ALT .LT. 600.0)) THEN       
          TI=(T600-T430)/170.0*(ALT-430.0)+T430
       SIGTI=(E600-E430)/170.0*(ALT-430.0)+E430    
      END IF
      IF (ALT .GE. 600.0) THEN  
          IF (T850 .LT. T600) THEN
           T850=T600 
          END IF     
          TI=(T850-T600)/250.0*(ALT-600.0)+T600
       SIGTI=(E850-E600)/250.0*(ALT-600.0)+E600    
      END IF
C ////////////////////////////////////////////////////////
      RETURN
      END