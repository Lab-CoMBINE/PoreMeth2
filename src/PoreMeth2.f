       SUBROUTINE BIOVITERBII(ETAV,P,EMISSION,T,KTILDE,PATH,PSI)
       
       IMPLICIT NONE
       INTEGER T,KTILDE,IND,PATH(T),PSI(KTILDE,T),I,J,K,COUNTP
       DOUBLE PRECISION NUMMAX
       DOUBLE PRECISION ETAV(KTILDE),EMISSION(KTILDE,T)
       DOUBLE PRECISION P(KTILDE,KTILDE),DELTA(KTILDE,T)

       
       DO 202 I=1,KTILDE
          DELTA(I,1)=ETAV(I)+EMISSION(I,1)
          PSI(I,1)=0
 202   CONTINUE

       COUNTP=0
       DO 203 I=2,T
          DO 213 J=1,KTILDE
             NUMMAX=0.0
             NUMMAX=DELTA(1,I-1)+P(1,(J+COUNTP))
             IND=1
             DO 223 K=2,KTILDE
                IF((DELTA(K,I-1)+P(K,(J+COUNTP))).GT.NUMMAX) THEN
                   NUMMAX=(DELTA(K,I-1)+P(K,(J+COUNTP)))
                   IND=K
                ENDIF
 223         CONTINUE
             PSI(J,I)=IND
             DELTA(J,I)=NUMMAX+EMISSION(J,I)
 213      CONTINUE
          COUNTP=COUNTP+KTILDE
 203   CONTINUE
       
       NUMMAX=0.0
       NUMMAX=DELTA(1,T)
       IND=1
       DO 253 K=2,KTILDE
          IF(DELTA(K,T).GT.NUMMAX) THEN
             NUMMAX=DELTA(K,T)
             IND=K
           ENDIF
 253   CONTINUE
       
       PATH(T)=IND
       DO 263 K=T-1,1,-1
          PATH(K)=PSI(PATH(K+1),K+1)
 263   CONTINUE

       RETURN
       END 


       SUBROUTINE TRANSEMISIHMM(MUK,NCOV,TOTALSEQ,KS,
     c COV,SEPSILON,T,PT,P,EMISSION,TRUNCCOEF)

       IMPLICIT NONE

       INTEGER T,KS,I,J,K,COUNTP,NCOV
       DOUBLE PRECISION PI,ELNSUB,COV(NCOV)
       DOUBLE PRECISION P(KS,KS*NCOV)
       DOUBLE PRECISION EMISSION(KS,T)
       DOUBLE PRECISION MUK(KS),PT(KS)
       DOUBLE PRECISION TOTALSEQ(T),TRUNCCOEF(KS)
       DOUBLE PRECISION SEPSILON(KS),MYZERO
       PARAMETER(PI=3.14159265358979)       

       COUNTP=0
       DO 700 I=1,NCOV
           DO 710 J=1,KS
               DO 720 K=1,KS
                   IF (J.EQ.K) THEN
                   P(J,(K+COUNTP))=ELNSUB(MYZERO,(PT(K)+COV(I)))
                   ELSE
                   P(J,(K+COUNTP))=PT(K)+COV(I)-LOG(KS-1.)
                   ENDIF
 720           CONTINUE
 710       CONTINUE
           COUNTP=COUNTP+KS
 700   CONTINUE



       DO 730 J=1,KS
          DO 740 K=1,T
              EMISSION(J,K)=EMISSION(J,K)+LOG(1/
     c        (SQRT(2*PI)*SEPSILON(J)*TRUNCCOEF(J)))+
     c        (-0.5*((TOTALSEQ(K)-MUK(J))/(SEPSILON(J)))**2)
 740      CONTINUE
 730   CONTINUE
       RETURN
       END 



       SUBROUTINE BIOVITERBI(ETAV,P,EMISSION,T,KTILDE,PATH,PSI)
       
       IMPLICIT NONE
       INTEGER T,KTILDE,IND,PATH(T),PSI(KTILDE,T),I,J,K
       DOUBLE PRECISION NUMMAX
       DOUBLE PRECISION ETAV(KTILDE),EMISSION(KTILDE,T)
       DOUBLE PRECISION P(KTILDE,KTILDE),DELTA(KTILDE,T)
       
       
       DO 302 I=1,KTILDE
          DELTA(I,1)=ETAV(I)+EMISSION(I,1)
          PSI(I,1)=0
 302   CONTINUE

       DO 303 I=2,T
          DO 313 J=1,KTILDE
             NUMMAX=0.0
             NUMMAX=DELTA(1,I-1)+P(1,J)
             IND=1
             DO 323 K=2,KTILDE
                IF((DELTA(K,I-1)+P(K,J)).GT.NUMMAX) THEN
                   NUMMAX=(DELTA(K,I-1)+P(K,J))
                   IND=K
                ENDIF
 323         CONTINUE
             PSI(J,I)=IND
             DELTA(J,I)=NUMMAX+EMISSION(J,I)
 313      CONTINUE

 303   CONTINUE
       
       NUMMAX=0.0
       NUMMAX=DELTA(1,T)
       IND=1
       DO 353 K=2,KTILDE
          IF(DELTA(K,T).GT.NUMMAX) THEN
             NUMMAX=DELTA(K,T)
             IND=K
           ENDIF
 353   CONTINUE
       
       PATH(T)=IND
       DO 363 K=T-1,1,-1
          PATH(K)=PSI(PATH(K+1),K+1)
 363   CONTINUE

       RETURN
       END 



       SUBROUTINE TRANSEMISJSLM(MUK,MI,ETA,TOTALSEQ,KTILDE,
     c NUMSEQ,SMU,SEPSILON,T,G,P,EMISSION,TRUNCCOEF)

       IMPLICIT NONE

       INTEGER T,KTILDE,NUMSEQ,I,J,K
       DOUBLE PRECISION NORM,GSUP,ETA,PI,ELNSUM
       DOUBLE PRECISION MI(NUMSEQ),P(KTILDE,KTILDE)
       DOUBLE PRECISION GVECT(KTILDE),TRUNCCOEF(NUMSEQ,KTILDE)
       DOUBLE PRECISION G(KTILDE,KTILDE)
       DOUBLE PRECISION EMISSION(KTILDE,T)
       DOUBLE PRECISION MUK(NUMSEQ,KTILDE)
       DOUBLE PRECISION TOTALSEQ(NUMSEQ,T)
       DOUBLE PRECISION SMU(NUMSEQ),SEPSILON(NUMSEQ)
       PARAMETER(PI=3.14159265358979)
       
      
       DO 627 I=1,KTILDE
             GSUP=0
             DO 617 J=1,NUMSEQ
                GSUP=GSUP+(-(((MUK(J,I)-MI(J))**2)/(2*SMU(J)**2)))
 617         CONTINUE
             GVECT(I)=GSUP
 627   CONTINUE
       
       NORM=GVECT(1)
       IF (KTILDE.GT.1) THEN
       DO 607 J=2,KTILDE
          NORM=ELNSUM(NORM,GVECT(J))
 607   CONTINUE
       ENDIF

       DO 600 I=1,KTILDE
          DO 602 J=1,KTILDE
             G(I,J)=GVECT(J)-NORM
 602      CONTINUE
 600   CONTINUE
              

       DO 610 J=1,KTILDE
              DO 620 K=1,KTILDE
                 IF (J.EQ.K) THEN
                    P(J,K)=ELNSUM(LOG(1-ETA),(LOG(ETA)+G(J,K)))
                    ELSE
                    P(J,K)=LOG(ETA)+G(J,K)
                 ENDIF
 620          CONTINUE
 610   CONTINUE

       DO 630 J=1,KTILDE
          DO 640 K=1,T
             DO 641 I=1,NUMSEQ
              EMISSION(J,K)=EMISSION(J,K)+LOG(1/
     c        (SQRT(2*PI)*SEPSILON(I)*TRUNCCOEF(I,J)))+
     c        (-0.5*((TOTALSEQ(I,K)-MUK(I,J))/(SEPSILON(I)))**2)
 641         CONTINUE
 640      CONTINUE
 630   CONTINUE
       RETURN
       END





       SUBROUTINE TRANSEMISISLM(MUK,MI,ETA,NCOV,TOTALSEQ,KTILDE,
     c NUMSEQ,SMU,SEPSILON,T,G,P,EMISSION,TRUNCCOEF)

       IMPLICIT NONE

       INTEGER T,KTILDE,NUMSEQ,I,J,K,COUNTP,NCOV
       DOUBLE PRECISION NORM,GSUP,ETA(NCOV),PI,ELNSUM
       DOUBLE PRECISION MI(NUMSEQ),P(KTILDE,KTILDE*NCOV)
       DOUBLE PRECISION GVECT(KTILDE),TRUNCCOEF(NUMSEQ,KTILDE)
       DOUBLE PRECISION G(KTILDE,KTILDE)
       DOUBLE PRECISION EMISSION(KTILDE,T)
       DOUBLE PRECISION MUK(NUMSEQ,KTILDE)
       DOUBLE PRECISION TOTALSEQ(NUMSEQ,T)
       DOUBLE PRECISION SMU(NUMSEQ),SEPSILON(NUMSEQ)
       PARAMETER(PI=3.14159265358979)
       
      
       DO 627 I=1,KTILDE
             GSUP=0
             DO 617 J=1,NUMSEQ
                GSUP=GSUP+(-(((MUK(J,I)-MI(J))**2)/(2*SMU(J)**2)))
 617         CONTINUE
             GVECT(I)=GSUP
 627   CONTINUE
       
       NORM=GVECT(1)
       IF (KTILDE.GT.1) THEN
       DO 607 J=2,KTILDE
          NORM=ELNSUM(NORM,GVECT(J))
 607   CONTINUE
       ENDIF

       DO 700 I=1,KTILDE
          DO 702 J=1,KTILDE
             G(I,J)=GVECT(J)-NORM
 702      CONTINUE
 700   CONTINUE
              
       COUNTP=0
       DO 710 I=1,NCOV
           DO 720 J=1,KTILDE
               DO 750 K=1,KTILDE
                   IF (J.EQ.K) THEN
                      P(J,(K+COUNTP))=ELNSUM(LOG(1-ETA(I)),(LOG(ETA(I))+
     c                G(J,K)))
                      ELSE
                      P(J,(K+COUNTP))=LOG(ETA(I))+G(J,K)
                   ENDIF
 750           CONTINUE
 720       CONTINUE
           COUNTP=COUNTP+KTILDE
 710   CONTINUE


       DO 730 J=1,KTILDE
          DO 740 K=1,T
             DO 741 I=1,NUMSEQ
              EMISSION(J,K)=EMISSION(J,K)+LOG(1/
     c        (SQRT(2*PI)*SEPSILON(I)*TRUNCCOEF(I,J)))+
     c        (-0.5*((TOTALSEQ(I,K)-MUK(I,J))/(SEPSILON(I)))**2)
 741         CONTINUE
 740      CONTINUE
 730   CONTINUE
       RETURN
       END


      DOUBLE PRECISION FUNCTION ELNSUM(X,Y)

      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      IF (X.GT.Y) THEN
         ELNSUM=X+LOG(1+EXP(Y-X))
      ELSE
         ELNSUM=Y+LOG(1+EXP(X-Y))
      ENDIF

      RETURN
      END


      DOUBLE PRECISION FUNCTION ELNSUB(X,Y)

      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      IF (X.GT.Y) THEN
         ELNSUB=X+LOG(1-EXP(Y-X))
      ELSE
         ELNSUB=Y+LOG(1-EXP(X-Y))
      ENDIF

      RETURN
      END


       SUBROUTINE MATRIXCOMPARE3(MCOORD1,MCOORD2,N,M,VECA
     c  ,COORDOVERLAP)
       
       IMPLICIT NONE
       INTEGER N,M,I,J,CONTROL,SUB
       DOUBLE PRECISION OVERLAP,L2,L1
       DOUBLE PRECISION MINOVERLAP,MAXOVERLAP,SCORE
       DOUBLE PRECISION SCOREA,SCOREB,VECA(N),COORDOVERLAP(N,2)
       DOUBLE PRECISION MCOORD1(N,2),MCOORD2(M,2)
       DOUBLE PRECISION COORD1(2),COORD2(2),DIST1,DISTN

       DO 109 I=1,N
        COORD1=MCOORD1(I,:)
        SCORE=0
        DIST1=ABS(COORD1(1)-MCOORD2(1,1))
        DISTN=ABS(COORD1(1)-MCOORD2(M,1))
        IF (DIST1.LE.DISTN) THEN
          J=1
          SUB=1
        ENDIF
        IF (DIST1.GE.DISTN) THEN
          J=M
          SUB= -1
        ENDIF

        CONTROL=1
           DO 209 WHILE (CONTROL.EQ.1)
              COORD2=MCOORD2(J,:)
                 IF (COORD1(2).GT.COORD2(2)) THEN
                    MINOVERLAP=COORD2(2)
                 ELSE
                    MINOVERLAP=COORD1(2)
                 ENDIF
                 IF (COORD1(1).GT.COORD2(1)) THEN
                    MAXOVERLAP=COORD1(1)
                 ELSE
                    MAXOVERLAP=COORD2(1)
                 ENDIF
                 OVERLAP=MINOVERLAP-MAXOVERLAP
                 IF (OVERLAP.GT.0) THEN
                    L2=COORD2(2)-COORD2(1)
                    L1=COORD1(2)-COORD1(1)
                    SCOREA=OVERLAP/L1
                    SCOREB=OVERLAP/L2
                    
                       IF (SCOREA.GT.SCORE) THEN
                        VECA(I)=J
                        COORDOVERLAP(I,1)=MAXOVERLAP
                        COORDOVERLAP(I,2)=MINOVERLAP
                        CONTROL=0
                       ENDIF
                 ENDIF
                IF ((SUB.EQ.1).AND.(COORD1(1).LT.COORD2(2))) THEN
                  CONTROL=0
                ENDIF
                IF ((SUB.EQ.(-1)).AND.(COORD1(2).GT.COORD2(1))) THEN
                  CONTROL=0
                ENDIF
              J=J+SUB
209      CONTINUE
109   CONTINUE
      RETURN
      END
      
      



       SUBROUTINE VECTORCOMPARE2(MCOORD1,MCOORD2,N,M,VECA)
       
       IMPLICIT NONE
       INTEGER N,M,I,J,CONTROL,SUB
       DOUBLE PRECISION VECA(N)
       DOUBLE PRECISION MCOORD1(N),MCOORD2(M,2)
       DOUBLE PRECISION COORD1,COORD2(2),DIST1,DISTN

       DO 407 I=1,N
        COORD1=MCOORD1(I)
        CONTROL=1
        DIST1=ABS(COORD1-MCOORD2(1,1))
        DISTN=ABS(COORD1-MCOORD2(M,1))
        IF (DIST1.LE.DISTN) THEN
          J=1
          SUB=1
        ENDIF
        IF (DIST1.GE.DISTN) THEN
          J=M
          SUB= -1
        ENDIF
           DO 507 WHILE (CONTROL.EQ.1)
              COORD2=MCOORD2(J,:)
              IF ((COORD1.GE.COORD2(1)).AND.(COORD1.LE.COORD2(2))) THEN
                  VECA(I)=J
              ENDIF
              IF ((SUB.EQ.1).AND.(COORD1.LT.COORD2(2))) THEN
                CONTROL=0
              ENDIF
              IF ((SUB.EQ.(-1)).AND.(COORD1.GT.COORD2(1))) THEN
                CONTROL=0
              ENDIF
              J=J+SUB
507      CONTINUE
407   CONTINUE
      RETURN
      END
      
      
      
      
      
      
       SUBROUTINE MATRIXCOMPARE(MCOORD1,MCOORD2,M,COORDOVERLAP,VECA
     c  ,VECB)
       
       IMPLICIT NONE
       INTEGER M,J
       DOUBLE PRECISION OVERLAP,L2,L1
       DOUBLE PRECISION MINOVERLAP,MAXOVERLAP
       DOUBLE PRECISION SCOREA,SCOREB,VECA(M),VECB(M)
       DOUBLE PRECISION COORDOVERLAP(M,2),MCOORD2(M,2)
       DOUBLE PRECISION COORD1(2),COORD2(2),MCOORD1(2)

       
       COORD1=MCOORD1
       DO 208 J=1,M
        COORD2=MCOORD2(J,:)
          IF (COORD1(2).GT.COORD2(2)) THEN
            MINOVERLAP=COORD2(2)
            ELSE
              MINOVERLAP=COORD1(2)
            ENDIF
            IF (COORD1(1).GT.COORD2(1)) THEN
              MAXOVERLAP=COORD1(1)
            ELSE
              MAXOVERLAP=COORD2(1)
            ENDIF
            OVERLAP=MINOVERLAP-MAXOVERLAP
              IF (OVERLAP.GT.0) THEN
                L2=COORD2(2)-COORD2(1)
                L1=COORD1(2)-COORD1(1)
                COORDOVERLAP(J,1)=MAXOVERLAP
                COORDOVERLAP(J,2)=MINOVERLAP
                SCOREA=OVERLAP/L1
                SCOREB=OVERLAP/L2
                VECA(J)=SCOREA
                VECB(J)=SCOREB
              ENDIF
208     CONTINUE
      RETURN
      END


      
      
      
      
       SUBROUTINE MATRIXOVERLAPSUM(MCOORD1,MCOORD2,N,M,OS)
       
       IMPLICIT NONE
       INTEGER N,M,I,J
       DOUBLE PRECISION OVERLAP,OS
       DOUBLE PRECISION MINOVERLAP,MAXOVERLAP
       DOUBLE PRECISION MCOORD1(N,2),MCOORD2(M,2)
       DOUBLE PRECISION COORD1(2),COORD2(2)
       
       OS=0

       DO 107 I=1,N
       COORD1=MCOORD1(I,:)
           DO 207 J=1,M
              COORD2=MCOORD2(J,:)
                 IF (COORD1(2).GT.COORD2(2)) THEN
                    MINOVERLAP=COORD2(2)
                 ELSE
                    MINOVERLAP=COORD1(2)
                 ENDIF
                 IF (COORD1(1).GT.COORD2(1)) THEN
                    MAXOVERLAP=COORD1(1)
                 ELSE
                    MAXOVERLAP=COORD2(1)
                 ENDIF
                 OVERLAP=MINOVERLAP-MAXOVERLAP
                 IF (OVERLAP.GT.0.) THEN
                  OS=OS+OVERLAP
                 ENDIF
207      CONTINUE
107   CONTINUE
      RETURN
      END      
      
      