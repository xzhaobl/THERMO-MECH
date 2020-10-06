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
C USER-DEFINED MATRIX     
      DIMENSION DE1(NTENS,NTENS),DE2(NTENS,NTENS),DE(NTENS,NTENS)
      DIMENSION DSTRESS1(NTENS),DSTRESS2(NTENS),STRESS1(NTENS),
     1STRESS2(NTENS)
      DIMENSION FG1(NTENS),FG2(NTENS),FG(NTENS)
      DIMENSION DEP1(NTENS,NTENS),DEP2(NTENS,NTENS)
      DIMENSION DDSTRAN(NTENS),ESTRESS(NTENS)
      
      SSTOL=1E-4
C INPUT MATERIAL PARAMETERS          
      FLAMA=PROPS(1)
      FKAPA=PROPS(2)
      GREF=PROPS(3)
      FM=PROPS(4)
      FN0=PROPS(5)
      FGA0=PROPS(6)
      FN=PROPS(7)

      FR=EXP((FN0-FGA0)/(FLAMA-FKAPA))
      patm=101.0
      FTIME=0.0 !PSEUDO TIME
      FDTIME=1.0  !PSEUDO TIME INCREMENT
C!ABAQUS拉伸为＋，压缩为-，计算先改变主应力主应变的正负号
      DO I=1,NDI
          STRESS(I)=-STRESS(I)
          STRAN(I)=-STRAN(I)
          DSTRAN(I)=-DSTRAN(I)
      END DO
      DO I=1,NSHR
          STRESS(NDI+I)=STRESS(NDI+I)
          STRAN(NDI+I)=STRAN(NDI+I)
          DSTRAN(NDI+I)=DSTRAN(NDI+I)
      END DO     
C LOOP BEGIN
888   CONTINUE
      FVOID1=STATEV(1)
      FPC1=STATEV(2)
      FVOID0=STATEV(3)     
C 第一次应力计算     
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      
      FP1=SINV1 !p
      FQ1=SINV2 !q
      IF(FQ1.LE.1e-5)FQ1=1e-5
      FSD1=SQRT(FP1**2+FQ1**2)
      FKMOD1=(1+FVOID1)*FP1/FKAPA !K
      FGMOD1=GREF*((1+FVOID1)**(-3.0))*SQRT(FP1/patm) !G
      
      CALL GETDE(FKMOD1,FGMOD1,DE1,NDI,NSHR,NTENS)
      
      FPFP=FN*((FQ1/FM/FP1)**(FN-1))*(-FQ1/FM)/FP1/FP1+1/FP1/LOG(FR) !F对p一阶导
      FPFQ=FN*((FQ1/FM/FP1)**(FN-1))/FM/FP1 !F对q一阶导
      FATA=FQ1/FP1 !应力比
      
      FPB1=FPC1*EXP(-((FATA/FM)**FN)*LOG(FR)) !image point 
      FSDB1=SQRT(1+FATA**2)*FPB1

      FDS1=(FM*FM*(FSD1/FSDB1)**2-FATA*FATA)/2.0/FATA
      FKP1=(1+FVOID1)/(FLAMA-FKAPA)/LOG(FR)*(FM*FM*(FSDB1/FSD1)**2-
     1FATA*FATA)/2.0/FATA
        
999   CONTINUE    
      
      DO I=1,NTENS
          DDSTRAN(I)=DSTRAN(I)*FDTIME
      END DO

      CALL GETDEP(FP1,FQ1,STRESS,FDS1,FPFP,FPFQ,DE1,FKP1,DEP1,DDSTRAN,
     1FGS1,FF1,FG1,NDI,NSHR,NTENS)
            
      DEVP1=FDS1*FGS1 !plastic volumetric STRAIN increment & hardening parameter
      
      FPC2=FPC1*EXP((1+FVOID1)/(FLAMA-FKAPA)*DEVP1)
  
      CALL GETDSTRESS(DEP1,DDSTRAN,DSTRESS1,NDI,NSHR,NTENS)
 
      DO I=1,NTENS
          STRESS1(I)=STRESS(I)+DSTRESS1(I)
      END DO
      
      DEV=(DDSTRAN(1)+DDSTRAN(2)+DDSTRAN(3)) !volumetric strain
      FVOID2=EXP(DEV)*(1+FVOID1)-1
C 第二次应力计算 
      CALL SINV(STRESS1,SINV1,SINV2,NDI,NSHR,NDI,NSHR,NTENS)
      FP2=SINV1
      FQ2=SINV2
      IF(FQ2.LE.1e-5)FQ2=1e-5
      FKMOD2=(1+FVOID2)*FP2/FKAPA
      FGMOD2=GREF*((1+FVOID2)**(-3.0))*SQRT(FP2/patm)
      CALL GETDE(FKMOD2,FGMOD2,DE2,NDI,NSHR,NTENS)
      
      FPFP=FN*((FQ2/FM/FP2)**(FN-1))*(-FQ2/FM)/FP2/FP2+1/FP2/LOG(FR)
      FPFQ=FN*((FQ2/FM/FP2)**(FN-1))/FM/FP2
      FATA=FQ2/FP2
      
      FSD2=SQRT(FP2**2+FQ2**2) 
      FPB2=FPC2*EXP(-((FATA/FM)**FN)*LOG(FR))
      FSDB2=SQRT(1+FATA**2)*FPB2
     
      FKP2=(1+FVOID2)/(FLAMA-FKAPA)/LOG(FR)*(FM*FM*(FSDB2/FSD2)**2-
     1FATA*FATA)/2.0/FATA
      FDS2=(FM*FM*(FSD2/FSDB2)**2-FATA*FATA)/2.0/FATA
  
      CALL GETDEP(FP2,FQ2,STRESS1,FDS2,FPFP,FPFQ,DE2,FKP2,DEP2,DDSTRAN,
     1FGS2,FF2,FG2,NDI,NSHR,NTENS)
          
      DEVP2=FDS2*FGS2
     
                  
      CALL GETDSTRESS(DEP2,DDSTRAN,DSTRESS2,NDI,NSHR,NTENS)
      
      DO I=1,NTENS
          ESTRESS(I)=0.5*(DSTRESS2(I)-DSTRESS1(I))
          STRESS2(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
      END DO
      
      error_K=abs((DEVP1-DEVP2)/DEVP2)!error for hardening parameter
           
      FEIE=0.0
      FEIS=0.0
      FERR=0.0
      DO I=1,NTENS
          FEIE=FEIE+ESTRESS(I)*ESTRESS(I)
          FEIS=FEIS+STRESS2(I)*STRESS2(I)
      END DO
      FERR=SQRT(FEIE/FEIS)
      IF(FERR.LE.error_K)FERR=error_K
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

      
      DO I=1,NTENS
          STRESS(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
      END DO
      DEVP=0.5*(DEVP1+DEVP2) !硬化参数更新
      FPC=FPC1*EXP((1+FVOID1)/(FLAMA-FKAPA)*DEVP)
      
      STATEV(1)=FVOID2
      STATEV(2)=FPC
      

      IF(FTIME.LT.1)THEN
      IF(FDTIME.GT.(1-FTIME))THEN
          FDTIME=1-FTIME
      END IF

      GOTO 888
      END IF
C 更新数据       
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      FP=SINV1
      FQ=SINV2
      FVOID=FVOID2
      STATEV(1)=FVOID2
      STATEV(2)=FPC
      
      FKMOD=(1+FVOID)*FP/FKAPA
      FGMOD=GREF*((1+FVOID)**(-3.0))*SQRT(FP/patm)
      CALL GETDE(FKMOD,FGMOD,DE,NDI,NSHR,NTENS)
      IF(FQ.LE.1e-5)FQ=1e-5
      
      FPFP=FN*((FQ/FM/FP)**(FN-1))*(-FQ/FM)/FP/FP+1/FP/LOG(FR)
      FPFQ=FN*((FQ/FM/FP)**(FN-1))/FM/FP
      FATA=FQ/FP
      
      FSD=SQRT(FP**2+FQ**2)     
      FPB=FPC*EXP(-((FATA/FM)**FN)*LOG(FR))
      FSDB=SQRT(1+FATA**2)*FPB
      FKP=(1+FVOID)/(FLAMA-FKAPA)/LOG(FR)*(FM*FM*(FSDB/FSD)**2-
     1FATA*FATA)/2.0/FATA
      FDS=(FM*FM*(FSD/FSDB)**2-FATA*FATA)/2.0/FATA 
     
      CALL GETDEP(FP,FQ,STRESS,FDS,FPFP,FPFQ,DE,FKP,DDSDDE,DSTRAN,FGS
     1,FF,FG,NDI,NSHR,NTENS)
      
C恢复ABAQUS原本的符号
      DO I=1,NDI
          STRESS(I)=-STRESS(I)
          STRAN(I)=-STRAN(I)
          DSTRAN(I)=-DSTRAN(I)
      END DO
      RETURN
      END
      
      SUBROUTINE GETDE(FKMOD,FGMOD,FDE,NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(NTENS,NTENS)
      DO I=1,NTENS
          DO J=1,NTENS
              FDE(I,J)=0.0
          END DO
      END DO
      DO I=1,NDI
          DO J=1,NDI
              FDE(I,J)=FKMOD-2.0/3.0*FGMOD
          END DO
      END DO
      DO I=1,NDI
      FDE(I,I)=FDE(I,I)+2.0*FGMOD
      END DO
      DO I=1,NSHR
      FDE(NDI+I,NDI+I)=FGMOD
      END DO
      END 
      
      SUBROUTINE GETDSTRESS(FDE,DER,DS,NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(NTENS,NTENS),DER(NTENS),DS(NTENS)
     
      DO I=1,NTENS
          DS(I)=0.0
      END DO
      DO I=1,NTENS
          DO J=1,NTENS
              DS(I)=DS(I)+FDE(I,J)*DER(J)
          END DO
      END DO
      END
      
      SUBROUTINE GETDEP(FP,FQ,FSTRESS,FDS,FPFP,FPFQ,FDE,FKP,DEP,DER,FGS,
     1FF,FG,NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FSTRESS(NTENS),FDE(NTENS,NTENS),DEP(NTENS,NTENS),
     1DER(NTENS)
      DIMENSION FA(NTENS),FB(NTENS),FC(NTENS),FD(NTENS),FE(NTENS),
     1FG(NTENS),FH(NTENS,NTENS),FI(NTENS)
      
      DO I=1,NDI
      FA(I)=1.0/3.0
      END DO
      DO I=1,NSHR
      FA(NDI+I)=0.0
      END DO
      
      CALL SINV(FSTRESS,SINV1,SINV2,NDI,NSHR)
      DO I=1,NDI
      FB(I)=1.0/2.0/FQ*(3*FSTRESS(I)-3*SINV1)
      END DO
      DO I=1,NSHR
      FB(NDI+I)=3.0/2.0/FQ*2*FSTRESS(NDI+I)
      END DO


      DO I=1,NTENS
          FC(I)=FDS*FA(I)+FB(I)
          FD(I)=FPFP*FA(I)+FPFQ*FB(I)
      END DO
     
      DO I=1,NTENS
          FE(I)=0.0
          DO J=1,NTENS
              FE(I)=FE(I)+FD(J)*FDE(I,J)
          END DO
      END DO

      FF=0.0
      DO I=1,NTENS
      FF=FF+FE(I)*FC(I)
      END DO

      DO I=1,NTENS
          FG(I)=0.0
          DO J=1,NTENS
              FG(I)=FG(I)+FDE(I,J)*FC(J)
          END DO
      END DO
      DO I=1,NTENS
          DO J=1,NTENS
              FH(I,J)=FG(I)*FE(J)
          END DO
      END DO
      DO I=1,NTENS
          DO J=1,NTENS
              DEP(I,J)=FDE(I,J)-1/(FKP+FF)*FH(I,J)
          END DO
      END DO
      
      DO I=1,NTENS
          FI(I)=0.0 !!!!Stress
          DO J=1,NTENS
              FI(I)=FI(I)+DEP(I,J)*DER(J)
          END DO
      END DO      
      
      FGS=0.0
      DO I=1,NTENS
          FGS=FGS+FD(I)*FI(I)/FKP
      END DO
      FUNLOAD=0.0
      DO I=1,NTENS
          FUNLOAD=FUNLOAD+FD(I)*DER(I)
      END DO
      IF(FUNLOAD.LT.0)THEN
          DO I=1,NTENS
              DO J=1,NTENS
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
       !STATEV(1)=0.634 ! OCR=1 Zhou(2015) itatic Clay
       !STATEV(1)=0.58! OCR=1 Zhou(2017)
       !STATEV(1)=0.879 ! OCR=1 Zhou(2015) Boom Clay
       !STATEV(1)=0.789 ! OCR=2 Zhou(2015)
       !STATEV(2)=150 !OCR=1
       !STATEV(2)=300 !OCR=2
       STATEV(3)=20
       !STATEV(4)=100000
CCCC
       Patm=101
       E1=0.67
       FLAMA=0.09
       FKAPPA=0.01
       M=1.2
       Y=COORDS(2)
       VSTRESS=18.0*(50-Y)+1
       HSTRESS=0.6*VSTRESS
       P=(VSTRESS+HSTRESS*2.0)/3.0
       Q=VSTRESS-HSTRESS
       PC=P*EXP(Q/M/P*LOG(2.718))
       STATEV(2)=PC
       E0=E1-FLAMA*LOG(PC/Patm)+FKAPPA*LOG(PC/P)
       STATEV(1)=E0
      RETURN
      End