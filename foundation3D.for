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
     1 STRESS2(NTENS)
      DIMENSION FG1(NTENS),FG2(NTENS),FG(NTENS)
      DIMENSION DEP1(NTENS,NTENS),DEP2(NTENS,NTENS)
      DIMENSION DDSTRAN(NTENS),ESTRESS(NTENS)
      DIMENSION TH_PLAS(NTENS),TH_PLAS1(NTENS),TH_PLAS2(NTENS)
C DEBUG      
      !Logical :: firstrun=.true.
      !  Integer :: tempvar
      !  If(firstrun)Then
      !      write(*,*)"please input an integer"
      !      read(*,*)tempvar
      !      firstrun=.false.
      !  End If
C INPUT MATERIAL PARAMETERS   
      SSTOL=1E-4
      FLama=PROPS(1)
      FKapa=PROPS(2)
      !Gref=PROPS(3)
      FU=PROPS(3)
      FM=PROPS(4)
      FN0=PROPS(5)
      FGA0=PROPS(6)
      FN=PROPS(7)
      FRN=PROPS(8)
      FRT=PROPS(9)
      alpha=PROPS(10)
      nc=PROPS(11)

      FR=EXP((FN0-FGA0-(FRN-FRT)*DTEMP)/(FLama-FKapa)) !
      BETA_T=(FRN)/(FLama-FKapa)
      patm=101.0
      FTIME=0.0 !PSEUDo TIME
      FDTIME=1.0  !PSEUDo TIME INCREMENT
C!ABAQUS拉伸为＋，压缩为-，计算先改变主应力主应变的正负号
      Do I=1,NDI
          STRESS(I)=-STRESS(I)
          STRAN(I)=-STRAN(I)
          DSTRAN(I)=-DSTRAN(I)
      End Do 
C LOOP BEGIN
888   CONTINUE
      FVOID1=STATEV(1)
      FPC1=STATEV(2)
    
C 第一次应力计算     
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      
      FP1=SINV1 !p
      FQ1=SINV2 !q
      If(FQ1.LE.1e-5)FQ1=1e-5
      FSD1=SQRT(FP1**2+FQ1**2) !distance to the original point
      FKMOD1=(1+FVOID1)*FP1/FKapa !K
      !FGMOD1=Gref*((1+FVOID1)**(-3.0))*SQRT(FP1/patm) !G
      FGMOD1=FKMOD1*3.0*(1-2.0*FU)/2.0/(1+FU)
      
      CALL GETDE(FKMOD1,FGMOD1,DE1,NDI,NSHR,NTENS)
      
      FPFP=FN*((FQ1/FM/FP1)**(FN-1))*(-FQ1/FM)/FP1/FP1+1/FP1/LOG(FR) !first order derivation of F to p
      FPFQ=FN*((FQ1/FM/FP1)**(FN-1))/FM/FP1 !first order derivation of F to q
      FATA=FQ1/FP1 !应力比
      FPFR1=LOG(FP1/FPC1)/LOG(FR)/LOG(FR)*(FRN-FRT)/(FLama-FKapa) ! partial F partial r * partial r partial T
      
      FPB1=FPC1*EXP(-((FATA/FM)**FN)*LOG(FR)) !image point 
      FSDB1=SQRT(1+FATA**2)*FPB1 ! distance of image point to the original point
      Ratio_1=FSDB1/FSD1
      
      FDS1=(FM*FM*(FSD1/FSDB1)**2-FATA*FATA)/2.0/FATA
      FKP1=(1+FVOID1)/(FLama-FKapa)/LOG(FR)*(FM*FM*(FSDB1/FSD1)**(2+nc)-
     & FATA*FATA)/2.0/FATA
        
999   CONTINUE    
      
      Do I=1,NTENS
          DDSTRAN(I)=DSTRAN(I)*FDTIME
          DDTEMP=DTEMP*FDTIME         
      End Do
      
      if (DDTEMP.LT.1e-5)alpha=0.0

      CALL GETDEP(FP1,FQ1,STRESS,FDS1,FPFP,FPFQ,DE1,FKP1,DEP1,DDSTRAN,
     & FGS1,FF1,FG1,NDI,NSHR,NTENS)
            
      DEVP1=FDS1*FGS1 !plastic volumetric STRAIN increment & hardening parameter
      
      !If(DDTEMP.GT.SSTOL)Then
      FPC2=FPC1*EXP((1+FVOID1)/(FLAMA-FKAPA)*DEVP1)*EXP(-BETA_T*DDTEMP)
      !Else
      ! FPC2=FPC1*EXP((1+FVOID1)/(FLAMA-FKAPA)*DEVP1)
      !End if
  
      CALL GETDSTRESS(alpha,DEP1,DDSTRAN,DDTEMP,DSTRESS1,NDI,NSHR,NTENS)
      CALL GETTH_PLAS(FG1,FF1,BETA_T,FR,TH_PLAS1,DDTEMP,FPFR1,
     & NDI,NSHR,NTENS)
      
      
      Do I=1,NTENS
          STRESS1(I)=STRESS(I)+DSTRESS1(I)
     !!&-TH_PLAS1(I)
      End Do
      
      DEV=(DDSTRAN(1)+DDSTRAN(2)+DDSTRAN(3)) !volumetric strain
      FVOID2=EXP(-DEV)*(1+FVOID1)-1 !对数应变&true strain
     & +alpha*DDTEMP
C 第二次应力计算 
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      FP2=SINV1
      FQ2=SINV2
      If(FQ2.LE.1e-5)FQ2=1e-5
      FKMOD2=(1+FVOID2)*FP2/FKapa
      !FGMOD2=Gref*((1+FVOID2)**(-3.0))*SQRT(FP2/patm)
      FGMOD2=FKMOD2*3.0*(1-2.0*FU)/2.0/(1+FU)
      
      CALL GETDE(FKMOD2,FGMOD2,DE2,NDI,NSHR,NTENS)
      
      FPFP=FN*((FQ2/FM/FP2)**(FN-1))*(-FQ2/FM)/FP2/FP2+1/FP2/LOG(FR)
      FPFQ=FN*((FQ2/FM/FP2)**(FN-1))/FM/FP2

      FATA=FQ2/FP2
      FPFR2=LOG(FP2/FPC2)/LOG(FR)/LOG(FR)*(FRN-FRT)/(FLama-FKapa)

      
      FSD2=SQRT(FP2**2+FQ2**2) 
      FPB2=FPC2*EXP(-((FATA/FM)**FN)*LOG(FR))
      FSDB2=SQRT(1+FATA**2)*FPB2
      Ratio_2=FSDB2/FSD2
     
      FKP2=(1+FVOID2)/(FLama-FKapa)/LOG(FR)*(FM*FM*(FSDB2/FSD2)**(nc+2)-
     1 FATA*FATA)/2.0/FATA
      FDS2=(FM*FM*(FSD2/FSDB2)**2-FATA*FATA)/2.0/FATA
  
      CALL GETDEP(FP2,FQ2,STRESS1,FDS2,FPFP,FPFQ,DE2,FKP2,DEP2,DDSTRAN,
     1 FGS2,FF2,FG2,NDI,NSHR,NTENS)
          
      DEVP2=FDS2*FGS2
     
                  
      CALL GETDSTRESS(alpha,DEP2,DDSTRAN,DDTEMP,DSTRESS2,NDI,NSHR,NTENS)
      CALL GETTH_PLAS(FG2,FF2,BETA_T,FR,TH_PLAS2,DDTEMP,
     & FPFR2,NDI,NSHR,NTENS)
      
      Do I=1,NTENS
          ESTRESS(I)=0.5*(DSTRESS2(I)-DSTRESS1(I)) !!!!!!!
          STRESS2(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
     !!&-0.5*TH_PLAS1(I)-0.5*TH_PLAS2(I)
      End Do
              
      FEIE=0.0
      FEIS=0.0
      FERR=0.0
      Do I=1,NTENS
          FEIE=FEIE+ESTRESS(I)*ESTRESS(I)
          FEIS=FEIS+STRESS2(I)*STRESS2(I)
      End Do
      FERR=SQRT(FEIE/FEIS)

      If(FERR.LE.1E-8)FERR=1E-8
      FBETA=0.8*SQRT(SSTOL/FERR)

      If(FERR.GT.SSTOL)Then
          
          If(FBETA.LE.0.1)FBETA=0.1
          FDTIME=FBETA*FDTIME
          GOTO 999
      Else
          FTIME=FTIME+FDTIME
 
          If(FBETA.GE.1.1)FBETA=1.1
          FDTIME=FBETA*FDTIME    
      End If

      
      Do I=1,NTENS       
          STRESS(I)=STRESS(I)+0.5*DSTRESS1(I)+0.5*DSTRESS2(I)
     & -0.5*TH_PLAS1(I)-0.5*TH_PLAS2(I)
      End Do
      DEVP=0.5*(DEVP1+DEVP2) !硬化参数更新
      !If(DDTEMP.GT.SSTOL)Then
      FPC=FPC1*EXP((1+FVOID1)/(FLAMA-FKAPA)*DEVP1)*EXP(-BETA_T*DDTEMP)
      !Else
      !FPC=FPC1*EXP((1+FVOID1)/(FLAMA-FKAPA)*DEVP1)
      !end if
      
      STATEV(1)=FVOID2
      STATEV(2)=FPC
      
      
      If(FTIME.LT.1)Then
      If(FDTIME.GT.(1-FTIME))Then
          FDTIME=1-FTIME
      End If

      GOTO 888
      End If
C 更新数据       
      CALL SINV(STRESS,SINV1,SINV2,NDI,NSHR)
      FP=SINV1
      FQ=SINV2
      FVOID=FVOID2
      
      FKMOD=(1+FVOID)*FP/FKapa
      !FGMOD=Gref*((1+FVOID)**(-3.0))*SQRT(FP/patm)
      FGMOD=FKMOD*3.0*(1-2.0*FU)/2.0/(1+FU)
      CALL GETDE(FKMOD,FGMOD,DE,NDI,NSHR,NTENS)
      If(FQ.LE.1e-5)FQ=1e-5
      
      FPFP=FN*((FQ/FM/FP)**(FN-1))*(-FQ/FM)/FP/FP+1/FP/LOG(FR)
      FPFQ=FN*((FQ/FM/FP)**(FN-1))/FM/FP
      FATA=FQ/FP
      FPFR2=LOG(FP/FPC)/LOG(FR)/LOG(FR)*(FRN-FRT)/(FLama-FKapa)

      
      FSD=SQRT(FP**2+FQ**2)     
      FPB=FPC*EXP(-((FATA/FM)**FN)*LOG(FR))
      FSDB=SQRT(1+FATA**2)*FPB
      FKP=(1+FVOID)/(FLama-FKapa)/LOG(FR)*(FM*FM*(FSDB/FSD)**(nc+2)-
     1 FATA*FATA)/2.0/FATA
      FDS=(FM*FM*(FSD/FSDB)**2-FATA*FATA)/2.0/FATA 
     
      CALL GETDEP(FP,FQ,STRESS,FDS,FPFP,FPFQ,DE,FKP,DDSDDE,DSTRAN,FGS,
     1 FF,FG,NDI,NSHR,NTENS)
      CALL GETTH_PLAS(FG,FF,BETA_T,FR,TH_PLAS,DTEMP,FPFR,NDI,NSHR,NTENS)
      STATEV(1)=FVOID2
      STATEV(2)=FPC
      STATEV(3)=TEMP+DTEMP
      STATEV(4)=FGMOD
     !! if (NOEL==1)then
     !!     if (NPT==1) then
     !!     !Open(101,file='E:\ratio.txt',position='append')
     !!     !Write(101,*),'Ratio_1',Ratio_1,'Ratio_2',Ratio_2
     !!     !Open(101,file='E:\strain.txt',position='append')
     !!     Open(101,file='E:\th_plas.txt',position='append')  
     !!!!     Write(101,*),'time',TIME(1),
     !!!!&'strain',STRAN(1),'Dstrain',DSTRAN(1)
     !!     Write(101,*),'time',TIME(1),'TH_PLAST',TH_PLAS(1)
     !!     Close(101)
     !!     end if
     !! end if 
C恢复ABAQUS原本的符号
      Do I=1,NDI
          STRESS(I)=-STRESS(I)
          STRAN(I)=-STRAN(I)
          DSTRAN(I)=-DSTRAN(I)
      End Do
      RETURN
      End
      
      SUBROUTINE GETDE(FKMOD,FGMOD,FDE,NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(NTENS,NTENS)     
      Do I=1,NTENS
          Do J=1,NTENS
              FDE(I,J)=0.0
          End Do
      End Do
      Do I=1,NDI
          Do J=1,NDI
              FDE(I,J)=FKMOD-2.0/3.0*FGMOD
          End Do
      End Do
      Do I=1,NDI
      FDE(I,I)=FDE(I,I)+2.0*FGMOD
      End Do
      Do I=1,NSHR
      FDE(NDI+I,NDI+I)=FGMOD
      End Do
      End 
      
      SUBROUTINE GETDSTRESS(alpha,FDE,DER,DTEMP,DS,NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FDE(NTENS,NTENS),DER(NTENS),DS(NTENS),alpha_s(NTENS)
     
      Do I=1,NTENS
          DS(I)=0.0
          alpha_s(I)=0.0
      End Do
      Do I=1,NDI
          alpha_s(I)=alpha*DTEMP/3
      END DO
      Do I=1,NTENS
          Do J=1,NTENS
              DS(I)=DS(I)+FDE(I,J)*(DER(J)+alpha_s(J))
          End Do
      End Do
      End
      
      SUBROUTINE GETDEP(FP,FQ,FSTRESS,FDS,FPFP,FPFQ,FDE,FKP,DEP,DER,FGS,
     & FF,FG,NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FSTRESS(NTENS),FDE(NTENS,NTENS),DEP(NTENS,NTENS),
     & DER(NTENS)
      DIMENSION FA(NTENS),FB(NTENS),FC(NTENS),FD(NTENS),FE(NTENS),
     & FG(NTENS),FH(NTENS,NTENS),FI(NTENS)
      
      Do I=1,NDI
      FA(I)=1.0/3.0
      End Do
      Do I=1,NSHR
      FA(NDI+I)=0.0
      End Do
      
      Do I=1,NDI
      FB(I)=3.0/2.0/FQ*(FSTRESS(I)-FP)
      End Do
      Do I=1,NSHR
      FB(NDI+I)=3.0/2.0/FQ*2.0*FSTRESS(NDI+I)
      End Do


      Do I=1,NTENS
          FC(I)=FDS*FA(I)+FB(I)
          FD(I)=FPFP*FA(I)+FPFQ*FB(I)
      End Do
     
      Do I=1,NTENS
          FE(I)=0.0
          Do J=1,NTENS
              FE(I)=FE(I)+FD(J)*FDE(I,J)
          End Do
      End Do

      FF=0.0
      Do I=1,NTENS
      FF=FF+FE(I)*FC(I)
      End Do

      Do I=1,NTENS
          FG(I)=0.0
          Do J=1,NTENS
              FG(I)=FG(I)+FDE(I,J)*FC(J)
          End Do
      End Do
      Do I=1,NTENS
          Do J=1,NTENS
              FH(I,J)=FG(I)*FE(J)
          End Do
      End Do
      Do I=1,NTENS
          Do J=1,NTENS
              DEP(I,J)=FDE(I,J)-1/(FKP+FF)*FH(I,J)
          End Do
      End Do
!unload-reload      
      Do I=1,NTENS
          FI(I)=0.0 !!!!Stress
          Do J=1,NTENS
              FI(I)=FI(I)+DEP(I,J)*DER(J)
          End Do
      End Do      
      
      FGS=0.0
      Do I=1,NTENS
          FGS=FGS+FD(I)*FI(I)/FKP
      End Do
      FUNLOAD=0.0
      Do I=1,NTENS
          FUNLOAD=FUNLOAD+FD(I)*FI(I)
      End Do
      If(FUNLOAD.LT.0)Then
          Do I=1,NTENS
              Do J=1,NTENS
                  DEP(I,J)=FDE(I,J)
              End Do
          End Do
          FGS=0
      End If
      End
   
      SUBROUTINE GETTH_PLAS(FG,FF,BETA_T,FR,TH_PLAS,DTEMP,FPFR,
     & NDI,NSHR,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION FG(NTENS),TH_PLAS(NTENS)
      Do I=1,NTENS
        TH_PLAS(I)=0.0
      End Do
    
      Do I=1,NTENS
       If(DTEMP.GE.1e-5) TH_PLAS(I)=FG(I)/FF*(FPFR+BETA_T/LOG(FR))*DTEMP
      End Do
     
      End
      
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
CCCCCC
       Patm=101
       E1=0.67
       FLAMA=0.25
       FKAPPA=0.05
       M=1.2
       Y=COORDS(3)
       VSTRESS=18.0*Y+1
       HSTRESS=1.2*VSTRESS
       P=(VSTRESS+HSTRESS*2.0)/3.0
       Q=VSTRESS-HSTRESS
       PC=P*EXP(Q/M/P*LOG(2.718))
       STATEV(2)=PC
       E0=E1-FLAMA*LOG(PC/Patm)+FKAPPA*LOG(PC/P)
       STATEV(1)=E0
      RETURN
      End