      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
      
      DIMENSION DE1(6,6),DE2(6,6),DE(6,6)
      DIMENSION DSTRESS1(6),DSTRESS2(6),STRESS1(6),STRESS2(6)
      DIMENSION FG1(6),FG2(6),FG(6)
      DIMENSION DEP1(6,6),DEP2(6,6)
      DIMENSION DDSTRAN(6),ESTRESS(6)
      
      SSTOL=1E-4
          
      FLAMA=PROPS(1)
      FKAPA=PROPS(2)
      GREF=PROPS(3)
      FM=PROPS(4)
      FN0=PROPS(5)
      FGA0=PROPS(6)
      FN=PROPS(7)

      FR=EXP((FN0-FGA0)/(FLAMA-FKAPA))
      patm=101.0
      FTIME=0.0
      FDTIME=1.0   

888   CONTINUE
      FVOID1=STATEV(1)
      FPC1=STATEV(2)
      FVOID0=STATEV(3)     
      
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      
      FP1=SINV1
      FQ1=SINV2
      IF(FQ1.LE.1e-5)FQ1=1e-5
      FSD1=SQRT(FP1**2+FQ1**2)
      FKMOD1=-(1+FVOID1)*FP1/FKAPA
      FGMOD1=GREF*((1+FVOID1)**(-3.0))*SQRT(-FP1/patm)
      
      CALL GETDE(FKMOD1,FGMOD1,DE1)
      
      FPFP=-FN*((-FQ1/FM/FP1)**(FN-1))*(-FQ1/FM)/FP1/FP1+1/FP1/LOG(FR)
      FPFQ=-FN*((-FQ1/FM/FP1)**(FN-1))/FM/FP1
      FATA=FQ1/FP1
      
      FPB1=-FPC1*EXP(-((-FATA/FM)**FN)*LOG(FR))
      FSDB1=-SQRT(1+FATA**2)*FPB1
      
      FKP1=-(1+FVOID1)/(FLAMA-FKAPA)/LOG(FR)*(FM*FM*(FSDB1/FSD1)-FATA*
     & FATA)/2.0/FATA
      FDS1=(FM*FM*(FSD1/FSDB1)-FATA*FATA)/2.0/FATA
      
999   CONTINUE    
      
      DO I=1,6
          DDSTRAN(I)=DSTRAN(I)*FDTIME
      END DO

      CALL GETDEP(FP1,FQ1,STRESS,FDS1,FPFP,FPFQ,DE1,FKP1,DEP1,DDSTRAN,
     & FGS1,FF1,FG1)
            
      DEVP1=FDS1*FGS1
      
      FPC2=FPC1*EXP(-(1+FVOID1)/(FLAMA-FKAPA)*DEVP1)
  
      CALL GETDSTRESS(DEP1,DDSTRAN,DSTRESS1)
 
      DO I=1,6
          STRESS1(I)=STRESS(I)+DSTRESS1(I)
      END DO
      
      DEV=(DDSTRAN(1)+DDSTRAN(2)+DDSTRAN(3))
      FVOID2=EXP(DEV)*(1+FVOID1)-1
 
      CALL SINV(STRESS1,SINV1,SINV2,NDI,NSHR)
      FP2=SINV1
      FQ2=SINV2
      IF(FQ2.LE.1e-5)FQ2=1e-5
      FKMOD2=-(1+FVOID2)*FP2/FKAPA
      FGMOD2=GREF*((1+FVOID2)**(-3.0))*SQRT(-FP2/patm)
      CALL GETDE(FKMOD2,FGMOD2,DE2)
      
      FPFP=-FN*((-FQ2/FM/FP2)**(FN-1))*(-FQ2/FM)/FP2/FP2+1/FP2/LOG(FR)
      FPFQ=-FN*((-FQ2/FM/FP2)**(FN-1))/FM/FP2
      FATA=FQ2/FP2
      
      FSD2=SQRT(FP2**2+FQ2**2) 
      FPB2=-FPC2*EXP(-((-FATA/FM)**FN)*LOG(FR))
      FSDB2=-SQRT(1+FATA**2)*FPB2
     
      FKP2=-(1+FVOID2)/(FLAMA-FKAPA)/LOG(FR)*(FM*FM*(FSDB2/FSD2)-FATA*
     1FATA)/2.0/FATA
      FDS2=(FM*FM*(FSD2/FSDB2)-FATA*FATA)/2.0/FATA
  
      CALL GETDEP(FP2,FQ2,STRESS1,FDS2,FPFP,FPFQ,DE2,FKP2,DEP2,DDSTRAN,
     1FGS2,FF2,FG2)
          
      DEVP2=FDS2*FGS2
                  
      CALL GETDSTRESS(DEP2,DDSTRAN,DSTRESS2)
      
      DO I=1,6
          ESTRESS(I)=0.5*(DSTRESS2(I)-DSTRESS1(I))
          STRESS2(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
      END DO
     
      FEIE=0.0
      FEIS=0.0
      FERR=0.0
      DO I=1,6
          FEIE=FEIE+ESTRESS(I)*ESTRESS(I)
          FEIS=FEIS+STRESS2(I)*STRESS2(I)
      END DO
      FERR=SQRT(FEIE/FEIS)
      IF(FERR.LE.1E-8)FERR=1E-8
      FBETA=0.8*SQRT(SSTOL/FERR)

      IF(FERR.GT.SSTOL)THEN
          
          IF(FBETA.LE.0.1)FBETA=0.1
          FDTIME=FBETA*FDTIME
          GOTO 999
      ELSE
          FTIME=FTIME+FDTIME
 
          IF(FBETA.GE.1.1)FBETA=1.1
          FDTIME=FBETA*FDTIME    
      END IF

      
      DO I=1,6
          STRESS(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
      END DO
      DEVP=0.5*(DEVP1+DEVP2)
      FPC=FPC1*EXP(-(1+FVOID1)/(FLAMA-FKAPA)*DEVP)
      
      STATEV(1)=FVOID2
      STATEV(2)=FPC
      

      IF(FTIME.LT.1)THEN
      IF(FDTIME.GT.(1-FTIME))THEN
          FDTIME=1-FTIME
      END IF

      GOTO 888
      END IF
       
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      FP=SINV1
      FQ=SINV2
      FVOID=FVOID2
      STATEV(1)=FVOID2
      STATEV(2)=FPC
      
      FKMOD=-(1+FVOID)*FP/FKAPA
      FGMOD=GREF*((1+FVOID)**(-3.0))*SQRT(-FP/patm)
      CALL GETDE(FKMOD,FGMOD,DE)
      IF(FQ.LE.1e-5)FQ=1e-5
      
      FPFP=-FN*((-FQ/FM/FP)**(FN-1))*(-FQ/FM)/FP/FP+1/FP/LOG(FR)
      FPFQ=-FN*((-FQ/FM/FP)**(FN-1))/FM/FP
      FATA=FQ/FP
      
      FSD=SQRT(FP**2+FQ**2)     
      FPB=-FPC*EXP(-((-FATA/FM)**FN)*LOG(FR))
      FSDB=-SQRT(1+FATA**2)*FPB
      FKP=-(1+FVOID)/(FLAMA-FKAPA)/LOG(FR)*(FM*FM*(FSDB/FSD)-FATA*
     1FATA)/2.0/FATA
      FDS=(FM*FM*(FSD/FSDB)-FATA*FATA)/2.0/FATA 
     
      CALL GETDEP(FP,FQ,STRESS,FDS,FPFP,FPFQ,DE,FKP,DDSDDE,DSTRAN,FGS
     1,FF,FG)
 
      RETURN
      END
      
      SUBROUTINE GETDE(FKMOD,FGMOD,FDE)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(6,6)
      DO I=1,6
          DO J=1,6
              FDE(I,J)=0.0
          END DO
      END DO
 
      FDE(1,1)=FKMOD+4.0/3.0*FGMOD
      FDE(2,2)=FKMOD+4.0/3.0*FGMOD
      FDE(3,3)=FKMOD+4.0/3.0*FGMOD
      FDE(4,4)=FGMOD
      FDE(5,5)=FGMOD
      FDE(6,6)=FGMOD
      FDE(1,2)=FKMOD-2.0/3.0*FGMOD
      FDE(1,3)=FKMOD-2.0/3.0*FGMOD
      FDE(2,1)=FKMOD-2.0/3.0*FGMOD
      FDE(2,3)=FKMOD-2.0/3.0*FGMOD
      FDE(3,1)=FKMOD-2.0/3.0*FGMOD
      FDE(3,2)=FKMOD-2.0/3.0*FGMOD
      END 
      
      SUBROUTINE GETDSTRESS(FDE,DER,DS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(6,6),DER(6),DS(6)
     
      DO I=1,6
          DS(I)=0.0
      END DO
      DO I=1,6
          DO J=1,6
              DS(I)=DS(I)+FDE(I,J)*DER(J)
          END DO
      END DO
      END
      
      SUBROUTINE GETDEP(FP,FQ,FSTRESS,FDS,FPFP,FPFQ,FDE,FKP,DEP,DER,FGS,
     1FF,FG)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FSTRESS(6),FDE(6,6),DEP(6,6),DER(6)
      DIMENSION FA(6),FB(6),FC(6),FD(6),FE(6),FG(6),FH(6,6),FI(6)
      FA(1)=1.0/3.0
      FA(2)=1.0/3.0
      FA(3)=1.0/3.0
      FA(4)=0.0
      FA(5)=0.0
      FA(6)=0.0

      FB(1)=1.0/2.0/FQ*(2*FSTRESS(1)-FSTRESS(2)-FSTRESS(3))
      FB(2)=1.0/2.0/FQ*(2*FSTRESS(2)-FSTRESS(1)-FSTRESS(3))
      FB(3)=1.0/2.0/FQ*(2*FSTRESS(3)-FSTRESS(1)-FSTRESS(2))
      FB(4)=3.0/2.0/FQ*2*FSTRESS(4)
      FB(5)=3.0/2.0/FQ*2*FSTRESS(5)
      FB(6)=3.0/2.0/FQ*2*FSTRESS(6)


      DO I=1,6
          FC(I)=FDS*FA(I)+FB(I)
          FD(I)=FPFP*FA(I)+FPFQ*FB(I)
      END DO
     
      DO I=1,6
          FE(I)=0.0
          DO J=1,6
              FE(I)=FE(I)+FD(J)*FDE(I,J)
          END DO
      END DO

      FF=0.0
      DO I=1,6
      FF=FF+FE(I)*FC(I)
      END DO

      DO I=1,6
          FG(I)=0.0
          DO J=1,6
              FG(I)=FG(I)+FDE(I,J)*FC(J)
          END DO
      END DO
      DO I=1,6
          DO J=1,6
              FH(I,J)=FG(I)*FE(J)
          END DO
      END DO
      DO I=1,6
          DO J=1,6
              DEP(I,J)=FDE(I,J)-1/(FKP+FF)*FH(I,J)
          END DO
      END DO
      
      DO I=1,6
          FI(I)=0.0 !!!!Stress
          DO J=1,6
              FI(I)=FI(I)+DEP(I,J)*DER(J)
          END DO
      END DO      
      
      FGS=0.0
      DO I=1,6
          FGS=FGS+FD(I)*FI(I)/FKP
      END DO
      FUNLOAD=0.0
      DO I=1,6
          FUNLOAD=FUNLOAD+FD(I)*DER(I)
      END DO
      IF(FUNLOAD.LT.0)THEN
          DO I=1,6
              DO J=1,6
                  DEP(I,J)=FDE(I,J)
              END DO
          END DO
          FGS=0
      END IF
      END
      
       SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
       STATEV(1)=0.654
       !STATEV(1)=0.58! OCR=1 Zhou(2017)
       !STATEV(1)=0.879 ! OCR=1 Zhou(2015)
       !STATEV(1)=0.913 ! OCR=2 Zhou(2015)
       STATEV(2)=150 !OCR=1
       !STATEV(2)=300 !OCR=2
       STATEV(3)=20
       !STATEV(4)=100000
      RETURN
      End