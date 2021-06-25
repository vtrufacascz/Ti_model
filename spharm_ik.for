      SUBROUTINE SPHARM_IK(C,L,M,COLAT,AZ)
C CALCULATES THE COEFFICIENTS OF THE SPHERICAL HARMONIC
C FROM IRI 95 MODEL
C NOTE: COEFFICIENTS CORRESPONDING TO COS, SIN SWAPPED!!!
      DIMENSION C(82)
      C(1)=1.
      K=2
      X=COS(COLAT)
      C(K)=X
      K=K+1
      DO 10 I=2,L
      C(K)=((2*I-1)*X*C(K-1)-(I-1)*C(K-2))/I
10    K=K+1
      Y=SIN(COLAT)
      DO 20 MT=1,M
      CAZ=COS(MT*AZ)
      SAZ=SIN(MT*AZ)
      C(K)=Y**MT
      K=K+1
      IF(MT.EQ.L) GOTO 16
      C(K)=C(K-1)*X*(2*MT+1)
      K=K+1
      IF((MT+1).EQ.L) GOTO 16
      DO 15 I=2+MT,L
      C(K)=((2*I-1)*X*C(K-1)-(I+MT-1)*C(K-2))/(I-MT)
15    K=K+1
16    N=L-MT+1
      DO 18 I=1,N
      C(K)=C(K-N)*SAZ
      C(K-N)=C(K-N)*CAZ
18    K=K+1
20    CONTINUE
      RETURN
      END