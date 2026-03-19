!////////////////////////////////////
!  Heat and Water transfer based on water chemical potential
!  Building Envelope Simulation for Multi-Story Model 
!  1-dimensional model made by Hiroaki Saito 
!  dx�̕ύX�A�t���`�����A�ʋC�w�A���I�ʁA�C�ۃf�[�^�A�`�B�����Z
!  For �R������O�Ǒ��w�K���f��  �����l��
!  made by Hiroaki Saito  August 2021
! SI unit
!///////////////////////////////////
IMPLICIT REAL*8(A-H,O-Z)
CHARACTER MOJI*110,OUTFILE*8,FILENAME0*12
!                    ���_���@�w���@�@�@�ǐ��@������
INTEGER,PARAMETER :: NXP=150,KMTLP=10,NWP=50,NRM=50  !,NYP=150,KATLP=70
!
INTEGER,PARAMETER :: MXIT=5000,NMTL=19
REAL,PARAMETER :: EPS1=0.01,EPS2=1.0e+2
INTEGER DAY,DAY2,DAY1
COMMON /MS/TEMPO(25),SJD(25),SJS(25),SJN(25),SAT(NWP,25),SJIN(NWP,25),WSJIN(NWP,25),WSJD(NWP,25)
!//////////COEFFOCIENT//////////////
DIMENSION GMA(NMTL),GCP(NMTL),DX(NWP,KMTLP),NX(NWP)
DIMENSION RMDG(NMTL),RMDL(NMTL),RMDH(NMTL)
DIMENSION DPDU(NWP,NXP),TMPC(NWP,NXP),WGT(NWP,NXP)
DIMENSION DWGX(NWP,NXP+1),DWLX(NWP,NXP+1),RMDX(NWP,NXP+1),DTGX(NWP,NXP+1),DTLX(NWP,NXP+1),DWX(NWP,NXP+1),DTX(NWP,NXP+1)
DIMENSION ADWGX(NWP,NXP+1),ADWLX(NWP,NXP+1),ADTGX(NWP,NXP+1),ADTLX(NWP,NXP+1),ADWX(NWP,NXP+1),ADTX(NWP,NXP+1)
DIMENSION ALF(NWP,KMTLP,2),ALD(NWP,KMTLP,2),ALDP(NWP,KMTLP,2),ALDT(NWP,KMTLP,2),ALFK(NWP,KMTLP,2),ALDK(NWP,KMTLP,2)
DIMENSION NAFX(NWP,KMTLP),NALX(NWP,KMTLP),NBUZAI(NWP,KMTLP),ALP(NWP,KMTLP),ALPK(NWP,KMTLP)
DIMENSION ALP_TOTAL(KMTLP),TEMP_CAVITY(KMTLP),HIGHT_CAVITY(KMTLP)
DIMENSION NOUTD(50,3),NOUTAVD(50,3),TPINT(NWP,KMTLP),RHINT(NWP,KMTLP),KMTL(KMTLP),WOUTAV(50),NOUTMLD(50,3)
!/////////VARIABLE///////////////
DIMENSION WPU(NWP,NXP),WPS(NWP,NXP),WPW(NWP,NXP),TMP(NWP,NXP),RH(NWP,NXP)
DIMENSION XN(NWP,NXP)
DIMENSION AWPU(NWP,NXP),ATMP(NWP,NXP),BWPU(NWP,NXP),BTMP(NWP,NXP)     !,AWPS(NXP),AWPW(NXP)
DIMENSION HWPU(NWP,NXP),HTMP(NWP,NXP)    !,HWPS(NXP),HWPW(NXP)
! DIMENSION STMP(KMTLP,2),SWPS(KMTLP,2),SWPU(KMTLP,2),SWPW(KMTLP,2)&
DIMENSION QS(NWP,KMTLP,2)  !,SRH(NWP,KMTLP,2),ASWPU(NWP,KMTLP,2)
DIMENSION SCWDAY(NWP,25),MCW(NWP),SCW(NWP),CWV(NWP),CWVDAY(NWP,25),NWCN(NWP)
DIMENSION XOD(25),RHOD(25),SATDV(NWP),QQ(NWP,NXP),XM(NWP,NXP),AXN(NWP,NXP),RHDIS(NWP,NXP,12,31),WGTDIS(NWP,NXP,12,31)
DIMENSION RHAVD(NWP,NXP),WGTAVD(NWP,NXP),DGDUQ(NWP,NXP),DGDTQ(NWP,NXP)
DIMENSION ATMPQ(NWP,NXP),SATA(NWP)
!DIMENSION TDES1(24),TDES2(24),THDES1(24),THDES2(24),TSDES1(24),TSDES2(24),THSDES1(24),THSDES2(24)
DIMENSION TPAVD(NWP,NXP),XNAVD(NWP,NXP),TPDIS(NWP,NXP,12,31),XNDIS(NWP,NXP,12,31)
DIMENSION WLOSS(NWP,NXP),WJW(NWP,NXP)  !DAMAGE FUNC
DIMENSION WJRAIN(NWP,NXP),SWJRAIN(NWP,NXP),NRAINPOINT(NWP,3),COFRAIN(NWP,2),VELOUT(24),RAIN(24),SUMRAIN(NWP,NXP,25)
DIMENSION RN(NWP,NXP),HRN(NWP,NXP),DAYSUMRAIN(NWP,NXP),ACOFRAIN(NWP,2),VELDIR(24)
DIMENSION MDAY(12),NORI(NWP,3),WLENG(NWP,4),AW(NWP),VOLR(NRM),QR(NRM)
DIMENSION TMPR(NRM),ATMPR(NRM),RHR(NRM),XR(NRM),WPUR(NRM),HTMPR(NRM),XRM(NRM)!,VPR(NRM)
DIMENSION TICK(NWP,NWP),KDVD(NWP,NWP),KWTYPE(NWP),AKASYA(NWP),RLF(NWP),AS(NWP),EMI(NWP)
DIMENSION DGDUR(NRM),HWPUR(NRM),DGDTR(NRM),AWPUR(NWP),AWSF(NRM),AWS(NRM)
DIMENSION NWORI(NWP,2),AWIN(NWP),WINK(NWP),WINSCW(NWP),WSIN(NWP,25),WSIND(NWP,25)
DIMENSION WSUP(NWP,25),HSUP(NWP,25),IROOMT(NWP),TRAV(NWP),TRDT(NWP),RRAV(NWP),KSUP(20)
DIMENSION IN(NRM,NRM),JN(NRM,NRM),NROPEN(NRM),QRIN(NRM,NRM)!,NWIN(NRM)
!/////////////CONSTANT/////////////////
DATA RG,RW,CPL,ROW/461.5,2.512D+6,4200,998/! [J/kg],CPL[J/kgK],ROW[kg/m3]
!DATA RG,RW,CPL,ROW/461.5,0.,4200,998/! [J/kg],CPL[J/kgK],ROW[kg/m3]
DATA ATP/273.15/
DATA OMG/1.2/
DATA MDAY/31,28,31,30,31,30,31,31,30,31,30,31/
!RW=0.
!
!**** FILE OPEN ****
PRINT *,  '�@�@�`��f�[�^�t�@�C������͂��ĉ�����'
OPEN(UNIT=10,FILE="")
  PRINT *,  '�@�@�C�ۃf�[�^�t�@�C������͂��ĉ�����'
  OPEN(UNIT=7,FILE="")
!PRINT *,  '�@�@�����f�[�^�t�@�C������͂��ĉ�����'
OPEN(UNIT=11,FILE="mcoff1.prn")
  PRINT *,  '�@�@���J�f�[�^�t�@�C������͂��ĉ�����'
OPEN(UNIT=13,FILE="")
READ(13,210)MOJI
READ(13,210)MOJI
PRINT *,  '�@�o�̓t�@�C������͂��ĉ������i8�o�C�g�ȓ��j'
READ(5,*)OUTFILE
PRINT *,  '�@�����v�Z����=1   �����v�Z�Ȃ�=0 '
READ(5,*)NREPT
PRINT *,  '�����ɘa�W������͂��ĉ����� '
READ(5,*)ROTOMG
PRINT *,  '���������̈����@����0�@�l��1 '
READ(5,*)I_HCOFF
PRINT *,  '�ʋC�w���ʏo�̓|�C���g�@����No�y�ёw�ԍ�1,�w�ԍ�2 '
READ(5,*)NQOUT1
READ(5,*)NQOUT2
READ(5,*)NQOUT3
!
READ(10,210)MOJI
!     �ܓx�A�o�x
      READ(10,210)D1
      READ(10,*)HIDO,HKEIDO,HOI
READ(10,210)MOJI
READ(10,*)LYEAR,MON1,MDAY1,LMON,LDAY    ! �N �J�n�� �J�n�� �I���� �I����
READ(10,210)MOJI
READ(10,*)TMPOC,TMPIC,RHO,RHI           !�O�C���x �������x �O�C���x �������x �i�����l�j
READ(10,210)MOJI
READ(10,*)NDVD                          !���ԕ���
READ(10,210)MOJI
READ(10,*)NROOM				!������
READ(10,210)MOJI
DO I=1,NROOM
  READ(10,*)VOLR(I),QR(I),IROOMT(I),TRAV(I),TRDT(I),RRAV(I),NROPEN(I)  !���̐ρ@���C��  ���x�ݒ�L���@���ω��x�@�U���@���ώ��x�A�J����
END DO
READ(10,210)MOJI
!READ(10,*)NOP                           !�J���f�[�^��
READ(10,210)MOJI
DO I=1,NROOM
  K=NROPEN(I)
  IF(K.EQ.0)CYCLE
  DO J=1,K
   READ(10,*)IN(I,J),JN(I,J),QRIN(IN(I,J),JN(I,J))
  END DO
END DO
READ(10,210)MOJI
READ(10,*)NWTYPE			!��TYPE��
DO I=1,NWTYPE
	READ(10,210)MOJI
	READ(10,210)MOJI
	READ(10,*)KMTL(I)		! �ǂ̑w��
        READ(10,210)MOJI
	DO J=1,KMTL(I)
         READ(10,*)NBUZAI(I,J),TPINT(I,J),RHINT(I,J),TICK(I,J),KDVD(I,J)&
	,ALFK(I,J,1),ALFK(I,J,2),ALDK(I,J,1),ALDK(I,J,2),ALPK(I,J)
        END DO
END DO
READ(10,210)MOJI
READ(10,*)NWN				!�Ǖ��ʐ�
READ(10,210)MOJI
READ(10,210)MOJI
DO I=1,NWN  
        !�ǎ��	���� �אڎ� ���� �Z�� ���� �u�����N �X�Ίp �A���x�h ���ˋz�� ���˗�
	READ(10,*)KWTYPE(I),NORI(I,1),NORI(I,2),NORI(I,3),(WLENG(I,J),J=1,4),AKASYA(I),RLF(I),AS(I),EMI(I)
END DO
!AWS=0.
!AWSF=0.
DO LW=1,NWN
    IW=KWTYPE(LW)
    AW(LW)=WLENG(LW,1)*WLENG(LW,2)-WLENG(LW,4)
        DO K=1,KMTL(IW)
		ALP(LW,K)=ALPK(IW,K)*TICK(IW,K)*WLENG(LW,2)         !�����J���ʐσ�A�̎Z�o
	END DO
    !    IF(NORI(LW,1).NE.0)
    IF(NORI(LW,1).EQ.5)THEN
        AWSF(NORI(LW,3))=AWSF(NORI(LW,3))+AW(LW)    ! ���B���˂��󂯂鏰�ʐ�
     ELSE
        AWS(NORI(LW,3))=AWS(NORI(LW,3))+AW(LW)    ! �g�U���˂��󂯂�ǁE�V��ʐ�
    END IF
    IF(NORI(LW,1).EQ.0)THEN
        AWSF(NORI(LW,2))=AWSF(NORI(LW,2))+AW(LW)
        AWS(NORI(LW,2))=AWS(NORI(LW,2))-AW(LW)
    END IF
END DO
!:::::::::::�����J���ʐσ�A�̒��񍇐�:::::::::::
!
DO I=1,NWTYPE     !���  �s�v�i�v�C���j
  DO K=1,KMTL(I)   !�w
    D1=0.
    D2=0.
    IF(ALPK(I,K).GT.1.E-7)THEN
    DO  LW=1,NWN    !�K
        IW=KWTYPE(LW)
	D1=D1+1/(ALP(LW,K)*ALP(LW,K))      
        D2=D2+WLENG(LW,3)                         !�ʋC�w�����̍��v
    END DO
    IF(D1.GT.1E-12)THEN
       ALP_TOTAL(K)=1/SQRT(D1) 
       HIGHT_CAVITY(K)=D2
    END IF
    END IF
  !PRINT *, " I=",I," K=",K,ALP_TOTAL(K),HIGHT_CAVITY(K)
  END DO
END DO
!PAUSE
!:::::::::::::::::::::::::::::::::::::::::::::::
READ(10,210)MOJI
READ(10,*)NWIN				!���Ǖ��ʐ�
READ(10,210)MOJI
READ(10,210)MOJI
DO I=1,NWIN             			
        !����	�אڎ�	�ʐ�	�M�ї����@�Օ��W��
	READ(10,*)NWORI(I,1),NWORI(I,2),AWIN(I),WINK(I),WINSCW(I)
     !   WINSCW(I)=0.
END DO
!:::::::::���I�v�Z���ʓ���:::::::
 READ(10,210)MOJI
 READ(10,*)NSUWC
 READ(10,*)(NWCN(I),I=1,NSUWC)
 ALDC=5./3600.
!:::::�o�͍��W::::::::::::::::::::::::: 
!PAUSE
READ(10,210)MOJI
READ(10,*)NOUT
READ(10,*)(NOUTD(I,1),I=1,NOUT)
READ(10,*)(NOUTD(I,2),I=1,NOUT)
!:::::���ϒl�o�͍ޗ�:::::::::::::::::::: 
READ(10,210)MOJI
READ(10,*)NOUTAV
READ(10,*)(NOUTAVD(I,1),I=1,NOUTAV)
READ(10,*)(NOUTAVD(I,2),I=1,NOUTAV)
!:::::�������ʏo�͍��W:::::::::::::::::::: 
READ(10,210)MOJI
READ(10,*)NOUTML
READ(10,*)(NOUTMLD(I,1),I=1,NOUTML)
READ(10,*)(NOUTMLD(I,2),I=1,NOUTML)
READ(10,*)(NOUTMLD(I,3),I=1,NOUTML)
!***********���M*********
READ(10,210)MOJI
READ(10,*)NSUP
READ(10,*)(KSUP(K),K=1,NSUP)
READ(10,210)MOJI
DO I=1,24
  READ(10,*)(HSUP(KSUP(K),I),WSUP(KSUP(K),I),K=1,NSUP)
END DO
DO K=1,NSUP
DO I=1,24
  HSUP(KSUP(K),I)= HSUP(KSUP(K),I)*1000.
!  HSUP(1,I)=5000.
!  WSUP(1,I)=0.5
END DO
END DO
!***********�Z�����̐ݒ�*********
READ(10,210)MOJI
  READ(10,*)NRAIN
READ(10,210)MOJI
IF(NRAIN.GT.0)THEN
  DO I=1,NRAIN
   READ(10,*)NRAINPOINT(I,1),NRAINPOINT(I,2),NRAINPOINT(I,3),ACOFRAIN(I,1),ACOFRAIN(I,2)
   PRINT *,I,NRAINPOINT(I,1),NRAINPOINT(I,2),NRAINPOINT(I,3),ACOFRAIN(I,1),ACOFRAIN(I,2)
  END DO 
END IF 
!PAUSE
!********************************
CLOSE(10)
210 FORMAT(A)
DTM=1./NDVD*3600.
PXCOF=1./133322.
!PRINT *, DTM,DX
!PAUSE  
!**********************************
      READ(11,210)MOJI
      DO I=1,NMTL
      READ(11,*)K,RMDH(I),RMDG(I),GCP(I),GMA(I)
!      PRINT *,K,RMDH(I),RMDG(I),GCP(I),GMA(I)
      RMDL(I)=0.
      GCP(I)=GCP(I)*GMA(I)
      END DO
      CLOSE(11)
!      PAUSE
!**********************************
!PRINT *, '  �x���o�͂Q�_ ' 
!READ(5,*)NDOS1,NDOS2
!
!
!*********CAL DX,NAFX ***************
DO I=1,NWTYPE
    K=0
    DO J=1,KMTL(I)
	K=K+KDVD(I,J)
	IF(KDVD(I,J).EQ.1)THEN
	  D1=2.
	 ELSE
	  D1=KDVD(I,J)
        END IF
	DX(I,J)=TICK(I,J)/(D1-1.)
	IF(J.EQ.1)THEN
          NAFX(I,J)=1
          NALX(I,J)=KDVD(I,J)
	ELSE
          NAFX(I,J)=NALX(I,J-1)+1
          NALX(I,J)=NAFX(I,J)+KDVD(I,J)-1
	END IF
      END DO
      NX(I)=K
     ! PRINT *, I,K ,(DX(I,J),J=1,KMTL(I))
     ! PAUSE
 END DO
!*******************************
!
!
!
! *******SAT   SINITIAL VALUE *******
  DO I=1,NWN
      SAT(I,1)=TMPOC
      SAT(I,2)=TMPOC
  END DO
      TEMPO(1)=TMPOC
      SJD(1)=0.				
      SJS(1)=0.				
      SJN(1)=0.				
      RHOD(1)=RHINT(1,1)
      D1=TMPOC+ATP				
      RH0=RHINT(1,1)				
      CALL GOFF(D1,FS0,VP)				
!      CALL FUNCX(FS0,RH0,XS)
      XOD(1)=VP*RH0*0.01		
      TEMPO(2)=TMPOC				
!      VELO(2)=0				
      SJD(2)=0.				
      SJS(2)=0.				
      SJN(2)=0.				
      XOD(2)=VP*RH0*0.01				
      RHOD(2)=RHINT(1,1)
!***************************************
!/////////////////////////////////////////
      FILENAME0=OUTFILE//'.DAT'
      OPEN(UNIT=8,FILE=FILENAME0)
      WRITE(8,*)'  FILENAME ',OUTFILE
!      CLOSE(8)
!/////////////////////////////////////////
      FILENAME0=OUTFILE//'.TXT'
      OPEN(UNIT=12,FILE=FILENAME0)
      WRITE(12,*)'  FILENAME ',OUTFILE
      WRITE(12,*)' '
!      WRITE(12,*)' YEAR MON DAY  ','    DAYLY AVERAGE TEMP&VP& SUMMATION OF FLUX  WLOSS(kg/kg) WL(kg/DAY) '
WRITE(12,*)' YEAR MON DAY �Z����(kg/m2DAY)  �ޕ��ϊܐ���(%)  �ܐ������z(%) ���ʌ��������zWloss(kg/kg) ����������WJ(kg/m3DAY) '
WRITE(12,*)'NYEAR MON DAY RAIN1',(K,K=2,NRAIN),'WGTAV1',(K,K=2,NOUTAV),('WGT1',(K,K=2,NOUTMLD(L,3)+1-NOUTMLD(L,2)),L=1,NOUTML),&
           ('WLOSS1',(K,K=2,NOUTMLD(L,3)+1-NOUTMLD(L,2)),L=1,NOUTML),('WJW1',(K,K=2,NOUTMLD(L,3)+1-NOUTMLD(L,2)),L=1,NOUTML)
!      CLOSE(8)
!/////////INITIAL CONDITION//////////////
!PRINT *, RH0,TP1,WP0,WP1,WP2
!PAUSE
DO L=1,NWN
 IW=KWTYPE(L)
 ! WRITE(8,*)AW(L),WLENG(L,1),WLENG(L,2),WLENG(L,3),WLENG(L,4),ALPK(IW,2),ALP(L,2)
  DO K=1,KMTL(IW)
   L1=NAFX(IW,K)
   L2=NALX(IW,K)
   L5=NBUZAI(IW,K)
   TP1=TPINT(IW,K)+ATP
   RH0=RHINT(IW,K)
   CALL SATUWPT(TP1,WP0)
   CALL REWPT(RH0,WP,TP1,RG)
   CALL GOFF(TP1,FS0,VP)				
!  CALL FUNCX(FS0,RH0,XS)
   WP1=WP
   WP2=WP0+WP1
   QQ(L,K)=0.
   DO I=L1,L2
     WPW(L,I)=WP2
     WPU(L,I)=WP1
     WPS(L,I)=WP0
     TMP(L,I)=TP1
     ATMPQ(L,I)=TP1-ATP
     TMPC(L,I)=TP1-ATP
     RH(L,I)=RH0
     XN(L,I)=VP*RH0*0.01
     XM(L,I)=VP*RH0*0.01
     AXN(L,I)=VP*RH0*0.01
!  AWPW(I)=WP2
     AWPU(L,I)=WP1
!  AWPS(I)=WP0
     ATMP(L,I)=TP1
!  HWPW(I)=WP2
     HWPU(L,I)=WP1
!  HWPS(I)=WP0
    HTMP(L,I)=TP1
    RHAVD(L,I)=0.
    WGTAVD(L,I)=0.
    XNAVD(L,I)=0.
    TPAVD(L,I)=0.
    CALL AHGANS(WD,RH0,L5)
    WGT(L,I)=WD
    WLOSS(L,I)=0.
   END DO
 DO I=1,2
  QS(L,K,I)=0.
 END DO
END DO
END DO
DO I=1,NROOM
  TMPR(I)=TMPIC
  ATMPR(I)=TMPIC+ATP
  HTMPR(I)=TMPIC+ATP
  RHR(I)=RHI
  XR(I)=VP*RHI*0.01
  XRM(I)=VP*RHI*0.01
  RH0=RHI
  TP1=TMPIC+ATP
  CALL REWPT(RH0,WP1,TP1,RG)
  WPUR(I)=WP1
  AWPUR(I)=WP1
END DO
DO L=1,NOUTAV
    I=NOUTAVD(L,1)
    J=NOUTAVD(L,2)
    IW=KWTYPE(I)
    L1=NAFX(IW,J)
    WOUTAV(L)=WGT(I,L1)
END DO
 !     SDES1=0.
 !     SDES2=0.
  !    HSDES1=0.
 !     HSDES2=0.
!
DO 444 NYEAR=1,LYEAR
WRITE(8,*)' '
WRITE(8,*) '  YEAR= ',NYEAR
WRITE(8,*) '  DATA  '
WRITE(8,*)'   Year Mon Day IM �~�J��(mm/h) ����(m/s) ���Z����(kg/m2h) �h���w�Z����(kg/m2h) ��C�w������(kg/m2) ' &
          ,' Temp_R1 Rh_R1 Vp_R1  Tmp1',(K,K=2,NOUT),'Vp1',(K,K=2,NOUT),'Rh1',(K,K=2,NOUT),'Mc1',(K,K=2,NOUT),'SAT1',(K,K=2,NWN),'WGTAV1',(K,K=2,NOUTAV) &
          ,' �O�C���x(��) �O�C�����C��(Pa) '
         ! ,' Vp1 2 3 4 5 6 7 8 9   Rh1 2 3 4 5 6 7 8 9  Mc1 2 3 4 5 6 7 8 9  SAT1 2 3 4 5 WGTAV1 2 3 4 5 6(%)'&
!
IF(NYEAR.GT.1) THEN
  READ(7,210)D1
  READ(13,210)D1
  READ(13,210)D1
END IF
      KDAY=0				
      IF(NYEAR.EQ.1)THEN				
        CALL YOMI(MON1,MDAY1,KDAY)				
	MONF=MON1			
	LMONL=LMON			
	IF(LYEAR.GT.1)LMONL=12			
	ELSE			
	MONF=1			
	LMONL=12
        IF(LYEAR.EQ.NYEAR)LMONL=LMON			
      END IF
      IMQ=0
DO 999 MON=MONF,LMONL
      DAY2=MDAY(MON)
      IF(NYEAR.EQ.LYEAR.AND.MON.EQ.LMON)THEN				
        DAY2=LDAY				
        ELSE				
        DAY2=MDAY(MON)				
      END IF				
      IF(NYEAR.EQ.1.AND.MON.EQ.MON1)THEN				
       DAY1=MDAY1				
       ELSE				
       DAY1=1				
      END IF					
DO 888 DAY=DAY1,DAY2
!       IF(DAY.NE.DAY1)THEN
!          NDVD=60
!          DTM=1./NDVD*3600.
!       END IF
	PRINT *, '  NYEAR=',NYEAR,'  MON=',MON,' DAY=',DAY,' IT1-3=',IT2,IT3
	KDAY=KDAY+1
!       CALL ROOM(KDAY,TR,RR,XR)
  100 FORMAT(F4.1,F3.1,3F3.0,2F3.3,2F4.3)				
  101 FORMAT(24F3.0)				
  102 FORMAT(24I3)
      READ(7,*)MOJI
!      PRINT *, MOJI
    DO  I=2,25
      READ(7,100)TEMPO(I),XOD(I),SJD(I),SJS(I),SJN(I),D1,D2,D3,D4
      SJD(I)=SJD(I)*1.16                                           !MKS
      SJS(I)=SJS(I)*1.16                                           !MKS
      SJN(I)=SJN(I)*1.16                                           !MKS
      XOD(I)=XOD(I)*0.001				
      D1=TEMPO(I)+ATP				
      RH0=100.				
      CALL GOFF(D1,FS0,VP)				
      CALL FUNCX(FS0,RH0,XS)				
      RHOD(I)=100.0*(XOD(I)/XS)      				
      XOD(I)=VP*RHOD(I)*0.01                                    !CHANGE
    END DO
!**   �O�������A�J�ʂ̓ǂݍ���***
    DO I=1,24
      READ(13,*)VELOUT(I),RAIN(I),VELDIR(I)
    !  PRINT *,VELOUT(I),RAIN(I)
    END DO
!*******************************
    DO LW=1,NWN
      AKASYA1=AKASYA(LW)
      NORIENT=NORI(LW,1)
      RLF1=RLF(LW)
      AS1=AS(LW)
      EMI1=EMI(LW)
      IF(NORIENT.LE.4.AND.NORIENT.GE.1)THEN
        CALL SOLOCT(KDAY,HIDO,HKEIDO,HOI,AKASYA1,NORIENT,LW)				
        CALL SATCAL(RLF1,AS1,EMI1,LW)
      END IF
    END DO
DO I=1,NWIN
	J=I+NWN
  	NORIENT=NWORI(I,1)
        AKASYA1=90.
        CALL SOLOCT(KDAY,HIDO,HKEIDO,HOI,AKASYA1,NORIENT,J)				
        CALL SATCAL(RLF1,AS1,EMI1,J)
        DO K=2,25
           WSIN(I,K)=WSJIN(J,K)*AWIN(I)*WINSCW(I)    !�g�U����
           WSIND(I,K)=WSJD(J,K)*AWIN(I)*WINSCW(I)     ! ���B����
	END DO
END DO
DO LW=1,NWN
 IW=KWTYPE(LW)
 DO I=1,NX(IW)
   WGTDIS(LW,I,MON,DAY)=WGT(LW,I)
 END DO
END DO
!*************���ԃ��[�v*****************
DO 777 IM=1,24
!PRINT *,'  IM=',IM
!/////KETURO/////////////
      DO LW=1,NWIN
	  CWVDAY(LW,IM)=0.
      END DO
!/////�Z�����̌v�Z�i���C�u�����z���l�����������̊֐��j////////
    IF(NRAIN.GT.0)THEN
      DO K=1,NRAIN
           I=NRAINPOINT(K,1)
           J=NRAINPOINT(K,2)
           SUMRAIN(I,J,IM)=0.
           HIGHT_SER=(I-0.5)*WLENG(I,3)+0.7  !���������ǒ����ɐݒ� ��b��0.7m
           VEL_MOD=VELOUT(IM)*(HIGHT_SER/6.5)**0.25    !�����̍����␳�i�����h=6.5m�j0.25��
           D1=0.
           D2=0.
           D3=0.
           D4=0.
          !IF(VELOUT(IM).GE.1.)THEN
          IF(VELOUT(IM).GE.0.1)THEN
           DO I=1,30
             D1=I*0.1*VEL_MOD 
             !                                                       !v      ���̐Z����������
             !IF(D1.GE.8.5)THEN
             !  D2=0.0029*D1*D1 - 0.0284*D1 + 0.0346                                     !g(v)
             !  IF(D2.GT.1.)D2=1.
             ! ELSE
             !  D2=0.
             ! END IF                                                                   !CIM=9.2  �ΐ�_��
             !  CIM=9.2                                                                  !�O�ǂ�Cim
             ! D2=ACOFRAIN(K,1)*D1*COS(HOI*3.1415/180)/CIM                              !�Z����=ACOFRAIN(K,1)
               D2=ACOFRAIN(K,1)*D1*COS((VELDIR(IM)-(180+90*(NORI(NRAINPOINT(K,1),1)-1)))/180*3.1415)*1.5*0.2   !�����̍l�� ASHRAE�̖\�I�W��(SEVERE)�𗘗p
               IF(D2.LT.0.)D2=0.
             !
             D1=3.14*0.5*D1/VEL_MOD/VEL_MOD*EXP(3.14*0.25*(D1/VEL_MOD)**2)      !f(v)
             D3=D3+D1*D2                                                                ! sig(f(v)*g(v))
             D4=D4+D1                                                                   ! sig(f(v))
            END DO
            COFRAIN(K,1)=D3/D4
            !COFRAIN(K,1)=0.1
           ELSE 
            COFRAIN(K,1)=0
           END IF
      END DO
    END IF
!//////////////////////////////////////////////////////////
      IMM=IM+1 
      IMQ=IMQ+1
!   �������ǋy�я��ʂ̓��ߓ��˗ʂ̌v�Z(������)
 WSINSUM=0.
 WSINDSUM=0.
 DO  I=1,NWIN
     WSINSUM=WSINSUM+WSIN(I,IM)     !�g�U�����i�e�ʊg�U�j
     WSINDSUM=WSINDSUM+WSIND(I,IM)  ! ���B�����i���ʓ��ˁj
 END DO
DO LW=1,NWN
  IW=KWTYPE(LW)
  K=KMTL(IW)
 ! IF(NORI(LW,1).EQ.5) QS(LW,K,2)=WSINDSUM*AW(LW)/AWSF(NORI(LW,3))+WSINSUM*AW(LW)/AWS(NORI(LW,3))  !���K��
 ! IF(NORI(LW,1).EQ.0) QS(LW,K,1)=WSINDSUM*AW(LW)/AWSF(NORI(LW,2))+WSINSUM*AW(LW)/AWS(NORI(LW,2))  !���ԏ��A�V��i����������j
 ! IF(NORI(LW,1).GE.1.AND.NORI(LW,1).LE.4) QS(LW,K,2)=WSINSUM*AW(LW)/AWS(NORI(LW,3))  !�O��
END DO
!      NCONT=1
!
!        TDES1(IM)=DES1
!        TDES2(IM)=DES2
!        THDES1(IM)=HDES1
!        THDES2(IM)=HDES2
!        TSDES1(IM)=SDES1
!        TSDES2(IM)=SDES2
!        THSDES1(IM)=HSDES1
!        THSDES2(IM)=HSDES2
!
IF(NYEAR.GE.1)THEN
 ! IF(MON.EQ.1.OR.MON.EQ.8)THEN
      !D1=SJD(IM)+SJS(IM)+SJN(IM)
       IF(IM.EQ.1)THEN
         I=24
        ELSE
         I=IM-1
       END IF
       D2=SUMRAIN(NRAINPOINT(1,1),NRAINPOINT(1,2),I)
       D1=SWJRAIN(NRAINPOINT(1,1),NRAINPOINT(1,2))*3600.
       D3=RN(NRAINPOINT(1,1),NRAINPOINT(1,2))
       WRITE(8,*)NYEAR,MON,DAY,IM-1,RAIN(I),VELOUT(I),D1,D2,D3& !   Year Mon Day IM �~�J��(mm/h) ����(m/s) ���Z����(kg/m2h) �h���w�Z����(kg/m2h) ��C�w������(kg/m2) 
               ,(TMPR(K),RHR(K),XR(K),K=1,NROOM),(TMPC(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(XN(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(RH(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT)&
               ,(WGT(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(SATDV(K),K=1,NWN),(WOUTAV(K),K=1,NOUTAV),TEMPO(I),XOD(I),QQ(NQOUT1,NQOUT2),QQ(NQOUT1,NQOUT3)
  !     WRITE(8,*)NYEAR,MON,DAY,IM-1,RAIN(I),COFRAIN(1,1),SUMRAIN(NRAINPOINT(2,1),NRAINPOINT(2,2),I),SWJRAIN(NRAINPOINT(2,1),NRAINPOINT(2,2)),WJRAIN(2,6)*3600,WJRAIN(2,7)*3600&
  !             ,WJRAIN(4,6)*3600,WJRAIN(4,7)*3600,(TMPR(K),RHR(K),XR(K),K=1,NROOM),(TMPC(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(XN(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(RH(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT)&
  !             ,(WGT(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(SATDV(K),K=1,NWN),TEMPO(I),XOD(I)
  !     WRITE(8,*)NYEAR,MON,DAY,IM,SUMRAIN(NRAINPOINT(1,1),NRAINPOINT(1,2),IM)&
  !             ,TMPR(1),(TMPC(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),XR(1),(XN(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),RHR(1),(RH(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT)&
  !             ,(WPU(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(CWVDAY(NWCN(K),I),K=1,NSUWC),(SATDV(K),K=1,NWN),TEMPO(IM),XOD(IM),(WSIN(K,IM),K=1,NWIN),WSJIN(11,IM+1)
  ! END IF
END IF
DO 666 NF=1,NDVD
     !    WRITE(8,*)NYEAR,MON,DAY,IM-1,SWJRAIN(2,6)/AW(2),SUMRAIN(2,6,IM-1)/AW(2),WJRAIN(2,6)*3600/NDVD/AW(2),&
     !   SWJRAIN(2,7)/AW(2),SUMRAIN(2,7,IM-1)/AW(2),WJRAIN(2,7)*3600/NDVD/AW(2)
     !   WRITE(8,*)NYEAR,MON,DAY,IM,NF,SJIN(NORI(1,1),IM)&
     !          ,TMPR(1),(TMPC(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),XR(1),(XN(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),RHR(1),(RH(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT)&
     !          ,(WPU(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(SCWDAY(NWCN(K),IM),K=1,NSUWC),(SAT(K,IM),K=1,4),TEMPO(IM),XOD(IM),SJD(IM),SJS(IM),SJN(IM)
!     ******�C�ۃf�[�^�̓ǂݍ���**********
      DO I=1,NWN
	IF(NORI(I,1).GE.1.AND.NORI(I,1).LE.4)THEN
          D1=(SAT(I,IMM)-SAT(I,IM))*DTM/3600
	  SATDV(I)=SAT(I,IM)+(NF-1)*D1
	END IF
        IF(NORI(I,1).EQ.5)THEN
	  D5=(TEMPO(IMM)-TEMPO(IM))*DTM/3600     !  �������x�i�ቺ��0.7�j
	 ! SATDV(I)=TEMPO(IM)+(NF-1)*D5
	  D5=TEMPO(IM)+(NF-1)*D5
          SATDV(I)=TMPR(NORI(I,3))-0.7*(TMPR(NORI(I,3))-D5)
	  SAT(I,IM)=SATDV(I)
	END IF
!SATDV(I)=5.0
      SATA(I)=SATDV(I)+ATP
      END DO
      D5=(TEMPO(IMM)-TEMPO(IM))*DTM/3600				
      D6=(XOD(IMM)-XOD(IM))*DTM/3600			
      D7=(RHOD(IMM)-RHOD(IM))*DTM/3600				
      TPO=TEMPO(IM)+(NF-1)*D5				
      XO=XOD(IM)+(NF-1)*D6				
      RHO=RHOD(IM)+(NF-1)*D7	
!XO=436.5749692
!TPO=5.0
!RHO=50.
!     ***************************
!   �O�C�̉��w�|�e���V�����̌v�Z
TMPOC=TPO
TP1=TMPOC+ATP
CALL REWPT(RHO,WP1,TP1,RG)
WPUO=WP1
TMPO=TMPOC+ATP
CALL CALDGDU(RHO,TMPO,RG,D1,D2)
DGDUO=1.2*PXCOF*D1
DGDTO=1.2*PXCOF*D2
!PRINT *, '  BEFORE 393'
!PAUSE
!****�Z���ʂ̐ݒ�*******
!�~�J�ʂ̊��Z 1mm/h=1L/m2h=1Kg/m2h
IF(NRAIN.GT.0)THEN
DO K=1,NRAIN
   I=NRAINPOINT(K,1)
   J=NRAINPOINT(K,2)
   IF(COFRAIN(K,1).GE.0.)THEN
      SWJRAIN(I,J)=0.01*RAIN(IM)*COFRAIN(K,1)/3600.     !�ޕ\�ʂւ̐Z����kg/m2S 
     !SWJRAIN(I,J)=SWJRAIN(I,J)+AW(I)*0.01*RAIN(IM)*COFRAIN(K,1)/NDVD   !  ��C�w�̕ێ�������
   END IF
END DO
END IF

!
!////////PUT IN PREBIOUS VALUE///////////
DO I=1,NROOM
  HWPUR(I)=WPUR(I)
  HTMPR(I)=ATMPR(I)
  XRM(I)=XR(I)
END DO
DO LW=1,NWN
 IW=KWTYPE(LW)
 DO K=1,KMTL(IW)
   L1=NAFX(IW,K)
   L2=NALX(IW,K)
   L5=NBUZAI(IW,K)
   DO I=L1,L2
!  HWPW(I)=WPW(I)
     HWPU(LW,I)=WPU(LW,I)
!  HWPS(I)=WPS(I)
     HTMP(LW,I)=TMP(LW,I)
     XM(LW,I)=XN(LW,I)
     HRN(LW,I)=RN(LW,I)
   END DO
  END DO
END DO
!
!
!///////////////COMBAIN METHIOD//////////////////////////
IT3=0
810 AMAX3=0.
    AMAX4=0.
!**** CULCURATION CONDUCTANCE****
DO LW=1,NWN
  IW=KWTYPE(LW)
  DO K=1,KMTL(IW)
    L1=NAFX(IW,K)
    L2=NALX(IW,K)
    L5=NBUZAI(IW,K)
!   GCP(L5)=0.
    DO I=L1,L2
     IF(RH(LW,I).GT.99.98)CYCLE
  ! IF(RH(LW,I).GT.98.0.AND.L5.EQ.3)CYCLE
     IF(RH(LW,I).GT.90.0.AND.L5.EQ.5)THEN
     ! RMDL(L5)=5.0E-9
     !  RMDL(L5)=0. 
      D1=RH(LW,I)
      D2=TMP(LW,I)
      CALL DIFF(L5,D1,D2,RG,D3)     !Liquid conductivity (kg/msJ/kg)
      RMDL(L5)=D3
     ELSE 
      RMDL(L5)=0.
    END IF
!*****************************
   RH1=RH(LW,I)
   TP1=TMP(LW,I)
   CALL CALDGDU(RH1,TP1,RG,DGDU,DGDT)
!*****************************
!
   D3=AW(LW)
!   D1=4.78D-13*EXP(0.104*RH(LW,I))    !�σV�[�g�̃R���_�N�^���X�i�o�T�F�c��j
   IF(I.EQ.L1)THEN
!     IF(L5.EQ.6)THEN                          !�σV�[�g�̃R���_�N�^���X�i�{�[�h���ʁj
!       ALD(LW,K,1)=D1*DGDU*D3
!       ALDT(LW,K,1)=D1*DGDT*D3
!       ALDP(LW,K,1)=D1*D3
!      ELSE 
       ALD(LW,K,1)=ALDK(IW,K,1)*DGDU*D3
       ALDT(LW,K,1)=ALDK(IW,K,1)*DGDT*D3
       ALDP(LW,K,1)=ALDK(IW,K,1)*D3
 !    END IF
     ALF(LW,K,1)=ALFK(IW,K,1)*D3
   END IF
   IF(I.EQ.L2)THEN
!     IF(L5.EQ.3)THEN                           !�σV�[�g�̃R���_�N�^���X�iGW�h���w�j
!       ALD(LW,K,2)=D1*DGDU*D3
!       ALDT(LW,K,2)=D1*DGDT*D3
!       ALDP(LW,K,2)=D1*D3
!     ELSE
       ALD(LW,K,2)=ALDK(IW,K,2)*DGDU*D3
       ALDT(LW,K,2)=ALDK(IW,K,2)*DGDT*D3
       ALDP(LW,K,2)=ALDK(IW,K,2)*D3
 !    END IF
     ALF(LW,K,2)=ALFK(IW,K,2)*D3
   END IF
   IF(L5.EQ.2)THEN   !��C�w
    DGDUQ(LW,K)=1.2*DGDU*PXCOF
    DGDTQ(LW,K)=1.2*DGDT*PXCOF 
    CYCLE
   END IF
!
   RMDX(LW,I)=RMDH(L5)*D3/DX(IW,K)
     WP=WPU(LW,I)
     TP=TMP(LW,I)
     D1=GMA(L5)
  IF(L5.NE.2)THEN
   IF(L5.NE.3)THEN
     CALL CALDPDU(WP,TP,D2,D1,L5,RG,ROW)
     DPDU(LW,I)=D2
    ELSE
     IF(RH(LW,I).LE.99.98)THEN
      CALL CALDPDU(WP,TP,D2,D1,L5,RG,ROW)
      DPDU(LW,I)=D2
     END IF
   END IF
   END IF
   D2=RMDG(L5)
   D1=RH(LW,I)*0.01
!   IF(D1.GE.0.92)D1=0.92
!   IF(D1.LE.0.25)D1=0.25
!    DPDU(I)=0.
!   ������2004.4
!   IF(L5.EQ.3)D2=(0.070066+0.0090852*D1*0.01+0.036*D1*D1*1.E-4)*2.778E-4
!   IF(L5.EQ.5)D2=(0.00068874-0.0012938*D1*0.01+0.004379*D1*D1*1.E-4)*2.778E-4
!   IF(L5.EQ.6)D2=(0.014025-0.0086462*D1*0.01+0.018157*D1*D1*1.E-4)*2.778E-4
!   IF(L5.EQ.8)D2=(0.010511-0.00097486*D1*0.01+0.0030587*D1*D1*1.E-4)*2.778E-4
!   IF(L5.EQ.11)D2=(0.0039139-0.0048556*D1*0.01+0.011291*D1*D1*1.E-4)*2.778E-4
!   IF(L5.EQ.12)D2=(0.0064129-0.005005*D1*0.01+0.012101*D1*D1*1.E-4)*2.778E-4
! ������2005
!   IF(L5.EQ.5)D2=2.778E-4*(0.0126*D1**1.1685*EXP(-1.4542*(1-D1**26.511)))
!   IF(L5.EQ.6)D2=2.778E-4*(0.049055*D1**0.08379*EXP(-1.2994*(1-D1**20.307)))
!   IF(L5.EQ.8)D2=2.778E-4*(0.020624*D1**0.067194*EXP(-0.554*(1-D1**27.37)))
!   IF(L5.EQ.11)D2=2.778E-4*(0.09277*D1**0.16041*EXP(-1.6124*(1-D1**15.565)))
!   IF(L5.EQ.12)D2=2.778E-4*(0.040355*D1**0.0066503*EXP(-1.2173*(1-D1**13.177)))
!
    IF(L5.EQ.5)D2=1.87E-11*D1**2.4019*EXP(-0.78864*(1-D1**1.1471))  !*3.45   !  Moisture conductivity (kg/msPa)  ASHRAE
   ADWGX(LW,I)=D2*DGDU/DX(IW,K)*D3
   ADWLX(LW,I)=RMDL(L5)*D3/DX(IW,K)
  ! ADWLX(LW,I)=RMDL(L5)*DPDU(LW,I)*D3/DX(IW,K)
   ADWX(LW,I)=ADWGX(LW,I)+ADWLX(LW,I)
   ADTGX(LW,I)=D2*DGDT*D3/DX(IW,K)
   ADTLX(LW,I)=0.*RMDL(L5)*DPDU(LW,I)*D3/DX(IW,K)
   ADTX(LW,I)=ADTGX(LW,I)+ADTLX(LW,I)
  END DO
  DO I=L1+1,L2
   DWGX(LW,I)=(ADWGX(LW,I)+ADWGX(LW,I-1))*0.5
   DWLX(LW,I)=(ADWLX(LW,I)+ADWLX(LW,I-1))*0.5
   DWX(LW,I)=(ADWX(LW,I)+ADWX(LW,I-1))*0.5
   DTGX(LW,I)=(ADTGX(LW,I)+ADTGX(LW,I-1))*0.5
   DTLX(LW,I)=(ADTLX(LW,I)+ADTLX(LW,I-1))*0.5
   DTX(LW,I)=(ADTX(LW,I)+ADTX(LW,I-1))*0.5
  END DO
 END DO
END DO
!
DO I=1,NROOM
   RH1=RHR(I)
   TP1=ATMPR(I)
   CALL CALDGDU(RH1,TP1,RG,DGDU,DGDT)
   DGDUR(I)=1.2*DGDU*PXCOF
   DGDTR(I)=1.2*DGDT*PXCOF
END DO


!
!
!////////////CAL_TEMP BY OVER RELAXATION METHOD///////////
IT1=0                                                                     
ITQ=0
IT4=IT4+1
!***** conbained  temp and ventiration****
  217 CONTINUE
      ITQ=ITQ+1
      IF(ITQ.GT.100)THEN
      PRINT *, '  no convagence   ventiration'
      GO TO 218
      END IF
!*****************************************
510 AMAX1=0. 
!******ROOM TEMP*****
DO I=1,NROOM
  UHEN=1300.*VOLR(I)/DTM*HTMPR(I)
  SAHEN=1300.*VOLR(I)/DTM
  DO LW=1,NWN     !   ���ǂ���̎擾�M��
   IW=KWTYPE(LW)
   IF(NORI(LW,3).EQ.I)THEN
      UHEN=UHEN+ALF(LW,KMTL(IW),2)*TMP(LW,NX(IW))!+RW*ALDP(LW,KMTL(IW),2)*(XN(LW,NX(IW))-XR(I)) ! SENSIBLE
      SAHEN=SAHEN+ALF(LW,KMTL(IW),2)
   END IF
   IF(NORI(LW,2).EQ.I)THEN
      UHEN=UHEN+ALF(LW,1,1)*TMP(LW,1)!+RW*ALDP(LW,1,1)*(XN(LW,1)-XR(I)) ! SENSIBLE
      SAHEN=SAHEN+ALF(LW,1,1)
   END IF
  END DO
  DO LW=1,NWIN     !   ���ǂ���̎擾�M��
    IF(NWORI(LW,2).EQ.I)THEN
      UHEN=UHEN+WINK(LW)*AWIN(LW)*TMPO      !+WSIN(LW,IM)*AWIN(LW)*WINSCW(LW)
      SAHEN=SAHEN+WINK(LW)*AWIN(LW)
    END IF
  END DO
  !QR(I)=0.
  IF(NROPEN(I).NE.0)THEN
    DO J=1,NROPEN(I)
      IF(JN(I,J).EQ.0)THEN
        UHEN=UHEN+1300.*QRIN(I,JN(I,J))*TMPO/3600
       ELSE 
        UHEN=UHEN+1300.*QRIN(I,JN(I,J))*ATMPR(JN(I,J))/3600      !  ���C���ׂ̐ώZ
      END IF
      SAHEN=SAHEN+1300.*QRIN(I,JN(I,J))/3600.
    END DO
  END IF
  !UHEN=UHEN+1300.*QR(I)*VOLR(I)*TMPO/3600.+HSUP(I,IM)   !   ���C����
  !SAHEN=SAHEN+1300.*QR(I)*VOLR(I)/3600.
  UHEN=UHEN+HSUP(I,IM)
  ATMPR(I)=UHEN/SAHEN
 ! ATMPR(I)=295.5
  TMPR(I)=ATMPR(I)-ATP
!****Constant temp********************
  IF(IROOMT(I).EQ.1)THEN
   D1=TRAV(I)
   D2=TRDT(I)
   D3=RRAV(I)
   CALL ROOM(KDAY,TRR,RR,XRR,D1,D2,D3)
   TMPR(I)=TRR
   ATMPR(I)=TMPR(I)+ATP
  END IF
END DO
!*************************************
!
!
DO LW=1,NWN
  IW=KWTYPE(LW)
  DO K=1,KMTL(IW)
   L1=NAFX(IW,K)
   L2=NALX(IW,K)
   L5=NBUZAI(IW,K)
!*******  VENTED CAVITY  ******
 IF(L5.EQ.2)THEN
   I=L1
   D1=1300.*DX(IW,K)*AW(LW)/DTM
   D2=ALF(LW,K,1)+ALF(LW,K,2)+RW*(ALDT(LW,K,1)+ALDT(LW,K,2))
   D4=(ALF(LW,K,1)+RW*ALDT(LW,K,1))*TMP(LW,I-1)+(ALF(LW,K,2)+RW*ALDT(LW,K,2))*TMP(LW,I+1)
   D5=RW*(ALD(LW,K,1)*(WPU(LW,I-1)-WPU(LW,I))+ALD(LW,K,2)*(WPU(LW,I+1)-WPU(LW,I)))
   IF(LW.EQ.1)THEN
     UHEN=D1*HTMP(LW,I)+D4+D5+1300.*QQ(LW,K)*TMPO
    ELSE
     UHEN=D1*HTMP(LW,I)+D4+D5+1300.*QQ(LW,K)*TMP(LW-1,I)
   END IF
   SAHEN=D1+D2+1300.*QQ(LW,K)
   TMP(LW,I)=UHEN/SAHEN
   CYCLE
  END IF
!*******************************
 DO I=L1,L2
  IF(I.EQ.L1)THEN               !*********�[���O
   D1=0.5*GCP(L5)*DX(IW,K)*AW(LW)/DTM
   IF(I.EQ.1)THEN                                               !�ύX
    IF(NORI(LW,2).EQ.0)THEN                                         !             �׎��O�C
      D2=ALF(LW,K,1)+RMDX(LW,I+1)+RW*(ALDT(LW,K,1)+DTGX(LW,I+1))
      D3=CPL*(DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
!      D4=(ALF(LW,K,1)+RW*ALDT(LW,K,1))*SATA(LW)+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
      D4=ALF(LW,K,1)*SATA(LW)+RW*ALDT(LW,K,1)*TMPO+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
      D5=RW*(ALD(LW,K,1)*(WPUO-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
      UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I+1)+QS(LW,K,1)
      SAHEN=D1+D2+D3
     ELSE                                                              !             �׎�����
      D2=ALF(LW,K,1)+RMDX(LW,I+1)+RW*(ALDT(LW,K,1)+DTGX(LW,I+1))
      D3=CPL*(DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
!      D4=(ALF(LW,K,1)+RW*ALDT(LW,K,1))*SATA(LW)+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
      D4=ALF(LW,K,1)*ATMPR(NORI(LW,2))+RW*ALDT(LW,K,1)*ATMPR(NORI(LW,2))+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
      D5=RW*(ALD(LW,K,1)*(WPUR(NORI(LW,2))-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
      UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I+1)+QS(LW,K,1)
      SAHEN=D1+D2+D3
     END IF
    ELSE   
!!      IF(NBUZAI(LW,K-1).EQ.2)QS(LW,K,1)=0.818*5.66*((ATMP(LW,I-2)/100.)**4.0-(ATMP(LW,I)/100.)**4.0)
      D2=ALF(LW,K,1)+RMDX(LW,I+1)+RW*(ALDT(LW,K,1)+DTGX(LW,I+1))
      D3=CPL*(DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
      D4=(ALF(LW,K,1)+RW*ALDT(LW,K,1))*TMP(LW,I-1)+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
      D5=RW*(ALD(LW,K,1)*(WPU(LW,I-1)-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
      UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I+1)+QS(LW,K,1)
      SAHEN=D1+D2+D3
   END IF
   TMP(LW,I)=UHEN/SAHEN
  END IF
  IF(I.GT.L1.AND.I.LT.L2)THEN   !*********������
   D1=GCP(L5)*DX(IW,K)*AW(LW)/DTM
   D2=CPL*(DWLX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))&
       +DTLX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
   D3=RMDX(LW,I)+RMDX(LW,I+1)+RW*(DTGX(LW,I)+DTGX(LW,I+1))
   D4=(RMDX(LW,I)+RW*DTGX(LW,I))*TMP(LW,I-1)+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
   D5=RW*(DWGX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
   D6=D2*(TMP(LW,I-1)+TMP(LW,I+1))
   D7=D2*2
   UHEN=D1*HTMP(LW,I)+D4+D5+D6
   SAHEN=D1+D3+D7
   TMP(LW,I)=UHEN/SAHEN
  END IF
  IF(I.EQ.L2)THEN               !*********�[����
   D1=0.5*GCP(L5)*DX(IW,K)*AW(LW)/DTM
   IF(I.EQ.NX(IW))THEN
      D2=ALF(LW,K,2)+RMDX(LW,I)+RW*(ALDT(LW,K,2)+DTGX(LW,I))
      D3=CPL*(DWLX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DTLX(LW,I)*(TMP(LW,I-1)-TMP(LW,I)))
      D4=(ALF(LW,K,2)+RW*ALDT(LW,K,2))*ATMPR(NORI(LW,3))+(RMDX(LW,I)+RW*DTGX(LW,I))*TMP(LW,I-1)
      D5=RW*(ALD(LW,K,2)*(WPUR(NORI(LW,3))-WPU(LW,I))+DWGX(LW,I)*(WPU(LW,I-1)-WPU(LW,I)))
      UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I-1)+QS(LW,K,2)
 !     IF(NORI(LW,1).EQ.5)THEN
 !       UHEN=UHEN+WSIND(LW,IM)*AWIN(LW)*WINSCW(LW)  !CHANGE
      SAHEN=D1+D2+D3
     ELSE
 !      IF(NBUZAI(K+1).EQ.2)QS(K,2)=0.818*5.66*((ATMP(I+2)/100.)**4.0-(ATMP(I)/100.)**4.0)
      D2=ALF(LW,K,2)+RMDX(LW,I)+RW*(ALDT(LW,K,2)+DTGX(LW,I))
      D3=CPL*(DWLX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DTLX(LW,I)*(TMP(LW,I-1)-TMP(LW,I)))
      D4=(ALF(LW,K,2)+RW*ALDT(LW,K,2))*TMP(LW,I+1)+(RMDX(LW,I)+RW*DTGX(LW,I))*TMP(LW,I-1)
      D5=RW*(ALD(LW,K,2)*(WPU(LW,I+1)-WPU(LW,I))+DWGX(LW,I)*(WPU(LW,I-1)-WPU(LW,I)))
      UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I-1)+QS(LW,K,2)
      SAHEN=D1+D2+D3
    END IF
    TMP(LW,I)=UHEN/SAHEN
  END IF
  ! WRITE(8,*)LW,IW,I,TMP(LW,I)
 END DO
!
!PRINT *, '  JISSITU SYUURYOU'  
!PRINT *, TMP(1),TMP(L2)
!PAUSE

!
END DO   !K
END DO   !LW
!/////////OVER RELAXATION ////// 
DO LW=1,NWN
  IW=KWTYPE(LW)
  DO K=1,KMTL(IW)
   L1=NAFX(IW,K)
   L2=NALX(IW,K)
   L5=NBUZAI(IW,K)
   DO I=L1,L2
    TMP(LW,I)=(TMP(LW,I)-ATMP(LW,I))*OMG+ATMP(LW,I)
    D8=TMP(LW,I)-ATMP(LW,I)
    IF(ABS(D8).GT.AMAX1)AMAX1=ABS(D8)
    ATMP(LW,I)=TMP(LW,I)
 END DO
END DO
END DO
!/////////JUDGEMENT CONVERGENCE//////
 IT1=IT1+1
 IF(IT1.GT.MXIT)THEN
!************************************************
!     OPEN(UNIT=8,FILE=FILENAME0,POSITION='APPEND')
      WRITE(8,*) ' HASSAN TEMP IT1= ',IT1
      WRITE(8,*)(TMP(1,I),I=1,NX(1))
      WRITE(8,*)' ������ '
    !PRINT *, ' HASSAN TEMP IT1= ',IT1
  201 FORMAT(A)
!      CLOSE(8)
!****************************************
      GO TO 555
    !   GO TO 218
END IF
IF(AMAX1.GT.EPS1) GO TO 510
!*************���C�v�Z��������********************
    IF(ITQ.GT.1)THEN
      D1=0
      D2=0
      DO LW=1,NWN
        IW=KWTYPE(LW)
        DO K=1,KMTL(IW)
         L1=NAFX(IW,K)
         L2=NALX(IW,K)
         L5=NBUZAI(IW,K)
          IF(L5.EQ.2.AND.ALPK(IW,K).GT.0.)THEN
           D1=ABS(TMP(LW,L1)-ATP-ATMPQ(LW,L1))
           IF(D1.GT.D2)D2=D1
          END IF
          ATMPQ(LW,L1)=TMP(LW,L1)-ATP
        END DO
      END DO
       	IF(D2.LT.EPS1)GO TO 218
       END IF
!*************���C�ʂ̎Z�o(Q=m3/s)****************
 DO I=1,NWTYPE
   DO K=1,KMTL(I)
    D1=0.
    L1=NAFX(I,K) !!!
    L2=NALX(I,K)
    L5=NBUZAI(I,K)
    DO LW=1,NWN
      IW=KWTYPE(LW)
      IF(L5.EQ.2.AND.ALPK(IW,K).GT.0)THEN
       D1=D1+TMP(LW,L1)*WLENG(LW,3)
      END IF
    END DO
!
    IF(HIGHT_CAVITY(K).GT.0.1)THEN
      TEMP_CAVITY(K)=D1/HIGHT_CAVITY(K)
      GMAQ1=353.25/(TEMP_CAVITY(K))
      GMAO=353.25/(TMPO)
      PQ1=GMAO-GMAQ1                       
      D1=0.5*ABS(PQ1)*HIGHT_CAVITY(K)*9.8                !�����с@�������� 0.5, �d�͉����x g(��-��)
    END IF
!
    DO LW=1,NWN
      IW=KWTYPE(LW)
      IF(L5.EQ.2.AND.ALPK(IW,K).GT.0)QQ(LW,K)=DSQRT(D1)*4.*ALP_TOTAL(K)
    END DO
  END DO
END DO
!PRINT *, '  ���C SYUURYOU'  
!PAUSE
 GOTO 217
218 CONTINUE
! ************END VENT****************************
!
!
 DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
    L1=NAFX(IW,K)
    L2=NALX(IW,K)
    L5=NBUZAI(IW,K)
     DO I=L1,L2
      TP=TMP(LW,I)
      CALL SATUWPT(TP,SW)
      WPS(LW,I)=SW
     END DO
    END DO
 END DO
!
!
!PRINT *, '  TMP SYUURYOU'  
!PAUSE
!////////////CAL_WATER POTENTIAL BY OVER RELAXATION METHOD/////////
  IT2=0
  CWMAX=2.8E-2
  !CWMAX=2.8
  ITC=0
710 AMAX2=0.D0 
!PXCOF=1./133322.
!******ROOM MOISTURE*****
!�������x�v�Z
DO I=1,NROOM
!  IF(IROOMT(I).EQ.1)CYCLE
  UHEN=DGDUR(I)*VOLR(I)/DTM*HWPUR(I)
  SAHEN=DGDUR(I)*VOLR(I)/DTM
  DO LW=1,NWN     !   ���ǂ���̕���
   IW=KWTYPE(LW)
   IF(NORI(LW,3).EQ.I)THEN
      UHEN=UHEN+ALD(LW,KMTL(IW),2)*WPU(LW,NX(IW))+ALDT(LW,KMTL(IW),2)*(TMP(LW,NX(IW))-ATMPR(I))
      SAHEN=SAHEN+ALD(LW,KMTL(IW),2)
   END IF
   IF(NORI(LW,2).EQ.I)THEN
      UHEN=UHEN+ALD(LW,1,1)*WPU(LW,1)+ALDT(LW,1,1)*(TMP(LW,1)-ATMPR(I))
      SAHEN=SAHEN+ALD(LW,1,1)
   END IF
  END DO
  DO LW=1,NWIN   !���K���X�̌��I
    IF(NWORI(LW,2).EQ.I)THEN
           IF(MCW(LW).EQ.1)THEN
           T1=ATMPR(I)
           WP=WPUR(I)
           CALL WPTRE(RH0,WP,T1,RG)
           CALL GOFF(T1,FS,VP)
           CALL FUNCX(FS,RH0,X1)
	   KMC=MCW(LW)
           T1=ATMPR(I)-WINK(LW)/9.3*(ATMPR(I)-TMPO)
	   CALL CWIF1(T1,X1,XSUT,KMC)
	   MCW(LW)=KMC
	   UHEN=UHEN+ALDC*(XSUT-X1)*AW(LW)
	  END IF
      END IF
  END DO
  IF(NROPEN(I).NE.0)THEN
    DO J=1,NROPEN(I)
      D1=QRIN(I,JN(I,J))/3600.
      IF(JN(I,J).EQ.0)THEN
        UHEN=UHEN+DGDUO*WPUO*D1+D1*(DGDTO*TMPO-DGDTR(I)*ATMPR(I)) 
       ELSE 
        UHEN=UHEN+DGDUR(JN(I,J))*WPUR(JN(I,J))*D1+D1*(DGDTR(JN(I,J))*ATMPR(JN(I,J))-DGDTR(I)*ATMPR(I))      !  ���C���ׂ̐ώZ
      END IF
      SAHEN=SAHEN+DGDUR(I)*D1
    END DO
  END IF

!  D1=QR(I)*VOLR(I)/3600.
!  UHEN=UHEN+DGDUO*WPUO*D1+D1*(DGDTO*TMPO-DGDTR(I)*ATMPR(I))+WSUP(I,IM)/3600.   !   ���C����
!  SAHEN=SAHEN+DGDUR(I)*QR(I)*VOLR(I)/3600.
  UHEN=UHEN+WSUP(I,IM)/3600.
  WPUR(I)=UHEN/SAHEN
  D1=(WPUR(I)-AWPUR(I))*OMG+AWPUR(I)
!  D8=(D1-AWPU(LW,I))/(D1+AWPU(LW,I))*2.
   D8=D1-AWPUR(I)
  AWPUR(I)=D1
!  WRITE(8,*)NF,LW,I,WPUR(1),UHEN,SAHEN
!****Constant humidity********************
  IF(IROOMT(I).EQ.1)THEN
   D1=TRAV(I)
   D2=TRDT(I)
   D3=RRAV(I)
   CALL ROOM(KDAY,TRR,RR,XRR,D1,D2,D3)
   !TRR=TMPR(I)
   D1=RR
   D2=TRR+ATP
   CALL REWPT(D1,WP,ATP,RG)  !  ???? D2 ?
   AWPUR(I)=WP
   WPUR(I)=WP
   RHR(I)=RR
   D8=0.
  END IF
  IF(ABS(D8).GT.AMAX2)AMAX2=ABS(D8)
END DO
!*************************************
!************�����̐����ێ��ʋy�ыz���ʂ̌v�Z**********
IF(NRAIN.GT.0)THEN
DO K=1,NRAIN
   I=NRAINPOINT(K,1)
   J=NRAINPOINT(K,2)
   IW=KWTYPE(I)
   L1=NBUZAI(IW,NRAINPOINT(K,3))
   L5=NBUZAI(IW,L1)
   IF(J.EQ.NALX(IW,NRAINPOINT(K,3)))THEN     !  �����Ƃ̗אڎ��_�̑I��
      D2=WPU(I,J-1)
      D3=XM(I,J+1)
     ELSE
      D2=WPU(I,J+1)
      IF(J.EQ.1)THEN
         D3=XO                               !�O���ޕ\�ʁA�O�C�Ƃ̎��x
       ELSE
         D3=XM(I,J-1)
      END IF
   END IF
   IF(SWJRAIN(I,J).GE.0.OR.RN(I,J).GE.0)THEN
     T1=TMP(I,J)
     CALL GOFF(T1,FS,VP)
     RH1=RH(I,J)
     CALL CALDGDU(RH1,T1,RG,DGDU,DGDT)
     D4=3.43E-08*(D3-VP)*0.3                                       !****�G��ʗ�0.3  ��������̏�����
     IF(L5.GE.10.AND.L5.LE.12)THEN                                 !�o�b�N�V�[���[������R 2.4e+5 m2sPa/kg by�@����
        D1=1/(DX(IW,J)/3.73E-6+2.4E+5/DGDU)                        !�O�a���̐����`���� 3.73e-6 kg/ms(J/kg)�@���Ԃ��� by �ɒ�@D�_
       ELSE
        D1=1/(DX(IW,J)/(3.73E-6*0.05)+2.4E+5/DGDU)                 !�ؐ��i��z��
     END IF
     RN(I,J)=(D1*D2+D4+SWJRAIN(I,J))*3600/NDVD+HRN(I,J)            !******�����̐�����kg/m2
     WJRAIN(I,J)=-D1*D2*AW(I)
     ! RN(I,J)=(D1*WPU(I,J+1)/DX(IW,J)+SWJRAIN(I,J))*3600/NDVD+HRN(I,J)  !******�����̐�����kg/m2
     ! WJRAIN(I,J)=-D1*WPU(I,J+1)*AW(I)/DX(IW,J)                         
     IF(RN(I,J).LE.0)THEN
         RN(I,J)=0.
         WJRAIN(I,J)=SWJRAIN(I,J)*AW(I)+HRN(I,J)*AW(I)                       !******�Z����WJRAIN kg/s
     END IF
!
     !D1=3.43E-08*(XM(I,J)-VP)*0.3                               !****�G��ʗ�0.3
     !IF(HRN(I,J).GE.0.005)THEN                                  !******����o����臒l5g/m2
     !   RN(I,J)=(D1-2.68E-06+SWJRAIN(I,J))*3600/NDVD+HRN(I,J)   !******���[�t�B���O�Z��
     !   WJRAIN(I,J+1)=2.68E-06*AW(I)                            !******�Z����WJRAIN kg/s
     !  ELSE                                                     !****C*GMA*SQRT(2gh)= 2.68E-06 kg/m2s
     ! RN(I,J)=(D1+SWJRAIN(I,J))*3600/NDVD+HRN(I,J)
     !  WJRAIN(I,J+1)=0.
     !END IF
     !   RN(I,J)=(-SWJRAIN(I,J))*3600/NDVD+HRN(I,J)   !******���[�t�B���O�Z���@�Z���͑S�č���
     !   WJRAIN(I,J+1)=SWJRAIN(I,J)*AW(I)             !******�Z����WJRAIN kg/s
      !IF(RN(I,J).LE.0)THEN
      !   RN(I,J)=0.
      !   WJRAIN(I,J)=0.
      !ELSE
      !   WJRAIN(I,J)=D1*AW(I)
      !END IF
  END IF
END DO
END IF
!****************************************
!
!
DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
    L1=NAFX(IW,K)
    L2=NALX(IW,K)
    L5=NBUZAI(IW,K)
!*******  VENTED CAVITY  ******
 IF(L5.EQ.2)THEN
   I=L1
   D4=0.
   D5=0.
   D1=ALD(LW,K,1)+ALD(LW,K,2)+DGDUQ(LW,K)*QQ(LW,K)
   D2=ALD(LW,K,1)*WPU(LW,I-1)+ALD(LW,K,2)*WPU(LW,I+1)     !+DGDUQ(LW,K)*QQ(LW,K)*WPUO
   D3=ALDT(LW,K,1)*(TMP(LW,I-1)-TMP(LW,I))+ALDT(LW,K,2)*(TMP(LW,I+1)-TMP(LW,I))   !+DGDTQ(LW,K)*QQ(LW,K)*(TMPO-TMP(LW,I))
   IF(LW.EQ.1)THEN
     D2=D2+DGDUQ(LW,K)*QQ(LW,K)*WPUO                              !�ŉ��w
     D3=D3+DGDTQ(LW,K)*QQ(LW,K)*(TMPO-TMP(LW,I))
     UHEN=D1*HTMP(LW,I)+D4+D5+1300.*QQ(LW,K)*TMPO                 !D4,D5 �s�v?
    ELSE
     D2=D2+DGDUQ(LW,K)*QQ(LW,K)*WPU(LW-1,I)                        !2�K�ȏ�
     D3=D3+DGDTQ(LW,K)*QQ(LW,K)*(TMP(LW-1,I)-TMP(LW,I))
   END IF
   IF(RN(LW,I+1).GT.0.)THEN                                       ! D4 �����̐��������� �G��ʗ�0.3
      T1=TMP(LW,I+1)
      CALL GOFF(T1,FS,VP)
      D4=3.43E-08*(VP-XM(LW,L1))*AW(LW)*0.3
   END IF
   IF(RN(LW,I-1).GT.0.)THEN                                       ! D5 �O���̐��������� �G��ʗ�0.3
      T1=TMP(LW,I-1)
      CALL GOFF(T1,FS,VP)
      D5=3.43E-08*(VP-XM(LW,L1))*AW(LW)*0.3
   END IF
   UHEN=D2+D3+D4+D5
   SAHEN=D1
   WPU(LW,I)=UHEN/SAHEN
   IF(WPU(LW,I).GE.0.00)WPU(LW,I)=-1.3E-4
!   WRITE(8,*) IT2,RH0,XN(LW,I),WP
   CYCLE
 END IF
! **************
 DO I=L1,L2
 D7=WJRAIN(LW,I)
 IF(I.EQ.L1)THEN               !***********�[���O
  D1=ROW*DPDU(LW,I)*DX(IW,K)*AW(LW)/DTM*0.5
   IF(I.EQ.1)THEN
    IF(NORI(LW,2).EQ.0)THEN
     D2=DWX(LW,I+1)+ALD(LW,K,1)
     D3=DWX(LW,I+1)*WPU(LW,I+1)+ALD(LW,K,1)*WPUO
!     D4=DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))+ALDT(LW,K,1)*(SATA(LW)-TMP(LW,I))
     D4=DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))+ALDT(LW,K,1)*(TMPO-TMP(LW,I))
     SAHEN=D1+D2
     UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)*DX(IW,K)*AW(LW)
    ELSE
     D2=DWX(LW,I+1)+ALD(LW,K,1)
     D3=DWX(LW,I+1)*WPU(LW,I+1)+ALD(LW,K,1)*WPUR(NORI(LW,2))
     D4=DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))+ALDT(LW,K,1)*(ATMPR(NORI(LW,2))-TMP(LW,I))
     SAHEN=D1+D2
     UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)*DX(IW,K)*AW(LW)
    END IF
   ELSE
        D2=DWX(LW,I+1)+ALD(LW,K,1)
        D3=DWX(LW,I+1)*WPU(LW,I+1)+ALD(LW,K,1)*WPU(LW,I-1)
        D4=DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))+ALDT(LW,K,1)*(TMP(LW,I-1)-TMP(LW,I))
        SAHEN=D1+D2
        UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)*DX(IW,K)*AW(LW)  !
   END IF
!    WPU(LW,I)=UHEN/SAHEN
 END IF
 IF(I.GT.L1.AND.I.LT.L2)THEN   !*********������
     D1=ROW*DPDU(LW,I)*DX(IW,K)*AW(LW)/DTM
     D2=DWX(LW,I)+DWX(LW,I+1)
     D3=DWX(LW,I)*WPU(LW,I-1)+DWX(LW,I+1)*WPU(LW,I+1)
     D4=DTX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))
     SAHEN=D1+D2
     UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)*DX(IW,K)*AW(LW)  !WJRAIN(LW,I)!
!     WPU(LW,I)=UHEN/SAHEN
!  UHEN=D3*D5+D4*D6+D2+D1*HWPU(LW,I)
 END IF
 IF(I.EQ.L2)THEN                !************�[����
  D1=ROW*DPDU(LW,I)*DX(IW,K)*AW(LW)/DTM*0.5
   IF(I.EQ.NX(IW))THEN
     D2=DWX(LW,I)+ALD(LW,K,2)
     D3=DWX(LW,I)*WPU(LW,I-1)+ALD(LW,K,2)*WPUR(NORI(LW,3))
     D4=DTX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+ALDT(LW,K,2)*(ATMPR(NORI(LW,3))-TMP(LW,I))
     SAHEN=D1+D2
     UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)*DX(IW,K)*AW(LW)  !WJRAIN(LW,I)
    ELSE
       D2=DWX(LW,I)+ALD(LW,K,2)
       D3=DWX(LW,I)*WPU(LW,I-1)+ALD(LW,K,2)*WPU(LW,I+1)
!     IF(NBUZAI(K+1).EQ.2)D3=DWX(I)*WPU(I-1)+ALD(K,2)*HWPU(I+1)
       D4=DTX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+ALDT(LW,K,2)*(TMP(LW,I+1)-TMP(LW,I))
       SAHEN=D1+D2
       UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)*DX(IW,K)*AW(LW)  !WJRAIN(LW,I)
  END IF
!  WPU(LW,I)=UHEN/SAHEN
 END IF
  WPU(LW,I)=UHEN/SAHEN
   IF(WPU(LW,I).GE.0.00)WPU(LW,I)=-1.3E-3      !-100.
 END DO
END DO   !K
END DO   !LW
!
!
!/////////OVER RELAXATION ////// 
DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
      L1=NAFX(IW,K)
      L2=NALX(IW,K)
      L5=NBUZAI(IW,K)
      DO I=L1,L2
       D1=(WPU(LW,I)-AWPU(LW,I))*OMG+AWPU(LW,I)
!  D8=(D1-AWPU(LW,I))/(D1+AWPU(LW,I))*2.
       D8=D1-AWPU(LW,I)
       IF(ABS(D8).GT.AMAX2)AMAX2=ABS(D8)
       AWPU(LW,I)=D1
      END DO
    END DO
END DO
DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
      L1=NAFX(IW,K)
      L2=NALX(IW,K)
      L5=NBUZAI(IW,K)
      DO I=L1,L2
        WPW(LW,I)=WPU(LW,I)+WPS(LW,I)
        WP0=WPU(LW,I)
        TP=TMP(LW,I)
        CALL WPTRE(RH0,WP0,TP,RG)
        RH(LW,I)=RH0
        CALL GOFF(TP,FS,VP)
!  CALL FUNCX(FS,RH0,D1)
        XN(LW,I)=RH0*VP*0.01
!  IF(RH0.GE.100.)RH0=100.  ! :::���U�h�~
        CALL AHGANS(WD,RH0,L5)
        WGT(LW,I)=WD
        TMPC(LW,I)=TMP(LW,I)-ATP
      END DO
    END DO
END DO

!/////////JUDGEMENT CONVERGENCE//////
 IT2=IT2+1
 IF(IT2.GT.MXIT)THEN
!************************************************
!     OPEN(UNIT=8,FILE=FILENAME0,POSITION='APPEND')
   !   WRITE(8,*) ' HASSAN WPU IT2= ',IT2
   !   WRITE(8,*)(WPU(I),I=1,NX)
   !   WRITE(8,*)' ������ '
!     CLOSE(8)
!****************************************
!   PRINT *, ' HASSAN WPU AMAX2= ',AMAX2
   AMAX2=0.
!   PAUSE
!      GO TO 555
END IF
IF(AMAX2.GT.EPS2) GO TO 710
!****************************************
!
!
!IF(DAY.EQ.1.AND.IM.EQ.1)&
!WRITE(8,*)NF,RHR(1),XR(1),WPUR(1),RH(1,NX(1)),XN(1,NX(1)),WPU(1,NX(1))

DO I=1,NROOM
IF(IROOMT(I).EQ.1)CYCLE
     WP0=WPUR(I)
     TP=ATMPR(I)
     CALL WPTRE(RH0,WP0,TP,RG)
     RHR(I)=RH0
     CALL GOFF(TP,FS,VP)
!  CALL FUNCX(FS,RH0,D1)
     XR(I)=RH0*VP*0.01
     TMPR(I)=ATMPR(I)-ATP
END DO
!::::::::���K���X�̌��I���蕔�A�������蕔�ɒǉ�
          IF(ITC.GE.30)THEN
          PRINT *, ' HASSAN KETURO '
          PAUSE
          GO TO 555
          END IF
	  K=0
	  ITC=ITC+1
  DO I=1,NROOM
	DO LW=1,NWIN
        IF(NWORI(LW,2).EQ.I)THEN
	  T1=ATMPR(I)
          RH0=RHR(I)
          CALL GOFF(T1,FS,VP)
          CALL FUNCX(FS,RH0,X1)
	  KMC=MCW(LW)
          T1=ATMPR(I)-WINK(LW)/9.3*(ATMPR(I)-TMPO)
	  CALL CWIF1(T1,X1,XSUT,KMC)
	  MCW(LW)=KMC
	  IF(MCW(LW).EQ.1)THEN
	!    T1=TMP(I)
        !    RH0=RH(I)
        !    CALL GOFF(T1,FS,VP)
        !    CALL FUNCX(FS,RH0,X1)
	    CW=(X1-XSUT)*DTM*ALDC*AW(LW)     !X1:ROOM, XSUT:WIN
            CWV(LW)=CW
	    SCW(LW)=CW+SCW(LW)
	    IF(SCW(LW).GT.CWMAX)SCW(LW)=CWMAX
	    IF(SCW(LW).LE.0)THEN
	     MCW(LW)=0
	     SCW(LW)=0.
             CWV(LW)=0.
	    END IF
	  ELSE
	!    T1=TMP(I)
        !    RH0=RH(I)
         !   CALL GOFF(T1,FS,VP)
         !   CALL FUNCX(FS,RH0,X1)
	    KMC=MCW(LW)
	    CALL CWIF1(T1,X1,XSUT,KMC)
	    MCW(LW)=KMC
	  IF(MCW(LW).EQ.1)THEN
	    K=1
	  END IF
	  END IF
	  SCWDAY(LW,IM)=SCW(LW)
        END IF
        END DO
     END DO
         IF(K.EQ.1)GO TO 710
!::::::::���I�ʓ���
   DO LW=1,NWIN
 	  IF(CWV(LW).GE.0.)CWVDAY(LW,IM)=CWVDAY(LW,IM)+CWV(LW)
   END DO
!
!:::::::::���I�I��
!
! ///////JUDGEMENT CONVERGENCE OF HEAT AND MOISTURE/////
IF(NREPT.EQ.1)THEN  !////COMBINE
DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
      L1=NAFX(IW,K)
      L2=NALX(IW,K)
      L5=NBUZAI(IW,K)
      DO I=L1,L2
!   D8=(WPU(LW,I)-BWPU(LW,I))/(WPU(LW,I)+BWPU(LW,I))*2.
        D8=WPU(LW,I)-BWPU(LW,I)
        D9=TMP(LW,I)-BTMP(LW,I)
   	BWPU(LW,I)=WPU(LW,I)
   	BTMP(LW,I)=TMP(LW,I)
   	IF(ABS(D8).GT.AMAX3)AMAX3=ABS(D8)
   	IF(ABS(D9).GT.AMAX4)AMAX4=ABS(D9)
 	END DO
    END DO
 END DO
!/////////JUDGEMENT CONVERGENCE//////
 IT3=IT3+1
 IF(IT3.GT.MXIT)THEN
 PRINT *, "   NO CONVERGENCE   IM=",IM
!************************************************
!     OPEN(UNIT=8,FILE=FILENAME0,POSITION='APPEND')
!       WRITE(8,*) ' HASSAN TMP&WPU IT3= ',IT3
!      WRITE(8,*)(WPU(I),I=1,NX)
!      WRITE(8,*)' ������ '
!     CLOSE(8)
!****************************************
!      GO TO 555
       GO TO 556
END IF
!IF(AMAX2.GT.EPS2.OR.AMAX1.GT.EPS1) GO TO 810
IF(AMAX3.GT.EPS2.OR.AMAX4.GT.EPS1) GO TO 810
END IF  !/////COMBINE
556 CONTINUE
!PRINT *,  SWPU(1,1),SWPU(1,2), ' ROOPOUT '
!PAUSE
!
!////////////HEAT FLUX J/m2s  MOISTURE FLUX  g/m2s//////////////////
 !     HDES1=ALF(1,1)*(TMPO-TMP(1,1))   !+RW*ALDP(1,1)*(XO-XN(1))  SENSITIVE 
 !     HDES2=ALF(KMTL,2)*(TMP(NX)-ATMPR(1,1))   !+RW*ALDP(KMTL,2)*(XN(NX)-XR)  SENSITIVE
 !     DES1=(ALD(1,1)*(WPUO-WPU(1))+ALDT(1,1)*(TMPO-TMP(1)))*1000.
 !     DES2=(ALD(KMTL,2)*(WPU(NX)-WPUI)+ALDT(KMTL,2)*(TMP(NX)-TMPI))*1000.
    !  DES2=(ALDP(KMTL,2)*(XN(NX)-XR))*1000.
 !     SDES1=SDES1+DES1*DTM
 !     SDES2=SDES2+DES2*DTM
 !     HSDES1=HSDES1+HDES1/NDVD
 !     HSDES2=HSDES2+HDES2/NDVD
!//////////CAL AVERAGE VALUE//////////////
DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
      L1=NAFX(IW,K)
      L2=NALX(IW,K)
      L5=NBUZAI(IW,K)
      DO I=L1,L2
       RHAVD(LW,I)=RH(LW,I)+RHAVD(LW,I)
       WGTAVD(LW,I)=WGT(LW,I)+WGTAVD(LW,I)
       XNAVD(LW,I)=XN(LW,I)+XNAVD(LW,I)
       TPAVD(LW,I)=TMPC(LW,I)+TPAVD(LW,I)
      END DO
    END DO
END DO
IF(NRAIN.GT.0)THEN      !�Z���ʂ̐ώZ SUMRAIN (kg/m2h)
DO K=1,NRAIN
   I=NRAINPOINT(K,1)
   J=NRAINPOINT(K,2)
       !SWJRAIN(I,J)=SWJRAIN(I,J)+(WJRAIN(I,J)-WJRAIN(I,J+1))*3600/NDVD
       !IF(SWJRAIN(I,J).LT.0)SWJRAIN(I,J)=0.
       SUMRAIN(I,J,IM)=SUMRAIN(I,J,IM)+WJRAIN(I,J)*3600/NDVD/AW(I)
END DO
END IF
666 CONTINUE
IF(NOUTAV.GT.0)THEN                           !�ޗ��̕��ϊܐ����̌v�Z
DO L=1,NOUTAV
  LW=NOUTAVD(L,1)
  J=NOUTAVD(L,2)
  IW=KWTYPE(LW)
  L1=NAFX(IW,J)
  L2=NALX(IW,J)
  D1=DX(IW,J)
  D2=WGT(LW,L1)*0.5   !  *DX(IW,J)
  D2=D2+WGT(LW,L2)*0.5  !*DX(IW,J)
  DO K=L1+1,L2-1
   ! D1=D1+DX(IW,J)
   ! D2=D2+WGT(LW,K)*DX(IW,J)
    D2=D2+WGT(LW,K)
  END DO
 ! WOUTAV(L)=D2/D1/(L2-L1+1)
  WOUTAV(L)=D2/(L2-L1+1)
END DO
END IF
  !  I=1
  !  J=2
  !  CALL DOSU(NYEAR,MON,I,RH(NDOS1)) !I=1 OR 2
  !  CALL DOSU(NYEAR,MON,J,RH(NDOS2))
!IF(NYEAR.GE.1.AND.NBOUT.EQ.0)THEN
!    IF(NYEAR.EQ.1.OR.NYEAR.EQ.5.OR.NYEAR.EQ.10)THEN
!      IF(MON.EQ.1.OR.MON.EQ.7)WRITE(8,*)MON,DAY,IM,SAT(NORI,IMM),(TMPC(NOUTD(K)),K=1,NOUT),(RH(NOUTD(K)),K=1,NOUT),(WGT(NOUTD(K)),K=1,NOUT),&
!              (XN(NOUTD(K)),K=1,NOUT),SCWDAY(1,IM),CWVDAY(1,IM)  !,QQ(2)*3600
 !     WRITE(8,*)NYEAR,MON,DAY,IM,SAT(NORI,IMM),(TMPC(NOUTD(K)),K=1,NOUT),(RH(NOUTD(K)),K=1,NOUT),(WGT(NOUTD(K)),K=1,NOUT)&
!               ,(XN(NOUTD(K)),K=1,NOUT),(SCWDAY(1,IM),K=1,NSUWC),TDES2(IM),THDES2(IM),TSDES2(IM),QQ(2)/DX(2)
!    END IF
777 CONTINUE
! *******SAT   IREKAE *******
      DO I=1,NWN
	SAT(I,1)=SAT(I,25)
        SJIN(I,1)=SJIN(I,25)
      END DO
      TEMPO(1)=TEMPO(25)				
!      VELO(1)=VELO(25)				
      SJD(1)=SJD(25)				
      SJS(1)=SJS(25)				
      SJN(1)=SJN(25)				
      XOD(1)=XOD(25)				
      RHOD(1)=RHOD(25)				
!***************************************
DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
      L1=NAFX(IW,K)
      L2=NALX(IW,K)
      L5=NBUZAI(IW,K)
      DO I=L1,L2
        RHDIS(LW,I,MON,DAY)=RHAVD(LW,I)/24./NDVD
!       WGTDIS(LW,I,MON,DAY)=WGTAVD(LW,I)/24./NDVD
        XNDIS(LW,I,MON,DAY)=XNAVD(LW,I)/24./NDVD
        TPDIS(LW,I,MON,DAY)= TPAVD(LW,I)/24./NDVD
        RHAVD(LW,I)=0.
        WGTAVD(LW,I)=0.
        XNAVD(LW,I)=0.
        TPAVD(LW,I)=0.
      END DO
    END DO 
END DO
!//////////�J���Z���ʐώZ�lDAYSUMRAIN(kg/m2DAY)/////////
DO K=1,NRAIN
   I=NRAINPOINT(K,1)
   J=NRAINPOINT(K,2)
   DAYSUMRAIN(I,J)=0.
   D1=0.
   DO L=1,24
       D1=SUMRAIN(I,J,L)+D1
   END DO
   DAYSUMRAIN(I,J)=D1
END DO
!////////Calculation Mass Loss�i�����ϒl�ɂ��Z�o�j//////////
!IF(D1.NE.10000)GO TO 1221
DO LW=1,NWN
    IW=KWTYPE(LW)
    DO K=1,KMTL(IW)
	L1=NAFX(IW,K)
 	L2=NALX(IW,K)
 	L5=NBUZAI(IW,K)
	WLOSSMAX=0.6
	IF(L5.EQ.4.OR.L5.EQ.5)THEN
	IF(I_HCOFF.EQ.1)THEN
		HCOFF=0.319  			  !����������
	   ELSE
		HCOFF=0.
	END IF
	DO I=L1,L2
	  D1=TPDIS(LW,I,MON,DAY)
	  D2=RHDIS(LW,I,MON,DAY)
	  CALL wood rot(LW,I,D1,D2,DLOSS,ROTOMG)
	  WLOSS(LW,I)=WLOSS(LW,I)+DLOSS
	  IF(WLOSS(LW,I) > WLOSSMAX)DLOSS=0.
	  WJW(LW,I)=HCOFF*DLOSS*GMA(L5)/86400.   ! Time unit:h=24,s=86400
	END DO 
	END IF 
    END DO 
END DO
1221 CONTINUE
!
!/////////////End Mass Loss//////////////

!
!WRITE(12,*)NYEAR,MON,DAY,(TPDIS(NOUTD(K,1),NOUTD(K,2),MON,DAY),K=1,NOUT),(XNDIS(NOUTD(K,1),NOUTD(K,2),MON,DAY),K=1,NOUT),&
!	   (RHDIS(NOUTD(K,1),NOUTD(K,2),MON,DAY),K=1,NOUT)
WRITE(12,*)NYEAR,MON,DAY,(DAYSUMRAIN(NRAINPOINT(K,1),NRAINPOINT(K,2)),K=1,NRAIN),(WOUTAV(K),K=1,NOUTAV),((WGTDIS(NOUTMLD(L,1),K,MON,DAY),K=NOUTMLD(L,2),NOUTMLD(L,3)),L=1,NOUTML),&
           ((WLOSS(NOUTMLD(L,1),K),K=NOUTMLD(L,2),NOUTMLD(L,3)),L=1,NOUTML),((WJW(NOUTMLD(L,1),K)*86400,K=NOUTMLD(L,2),NOUTMLD(L,3)),L=1,NOUTML)
!WRITE(12,*)NYEAR,MON,DAY,(DAYSUMRAIN(NRAINPOINT(K,1),NRAINPOINT(K,2)),K=1,NRAIN),(WOUTAV(K),K=1,NOUTAV),(WGTDIS(1,K,MON,DAY),K=7,12),(WGTDIS(5,K,MON,DAY),K=7,12),&
 !          (WLOSS(1,K),K=7,12),(WLOSS(5,K),K=7,12),(WJW(1,K)*86400,K=7,12),(WJW(5,K)*86400,K=7,12)
888 CONTINUE
999 CONTINUE
REWIND 7
REWIND 13
!WRITE(12,*)'NDOS1=',NDOS1,'NDOS1=',NDOS2
!PRINT *,'NDOS1=',NDOS1,'NDOS1=',NDOS2
!CALL SIGMA(NYEAR)
!DO L=1,LYEAR
WRITE(12,*)' '
WRITE(12,*)'  MOISTURE CONTENT DISTRIBUTION (mass%) '
WRITE(12,*)'  YEAR MON I '
 DO K=1,12
  J=1
  WRITE(12,*)NYEAR,K,(WGTDIS(1,I,K,J),I=1,NX(2))
 END DO
WRITE(12,*)' '
!END DO
444 CONTINUE
WRITE(12,*)' '
      WRITE(12,*)'  '
      PRINT *,'  POINT,TMPC(I),RH(I),WGT(I),XN(I)'
      DO I=1,NX(1)
       PRINT *,I,TMPC(1,I),RH(1,I),WGT(1,I),XN(1,I)
      END DO

!WRITE(8,*)'  WALL VALUE'
!DO LW=1,NWN
!   WRITE(8,*)' '
!  WRITE(8,*)'  WALL VALUE  No=',LW
!  IW=KWTYPE(LW)
!   WRITE(8,*)LW,IW,(TMP(LW,I),I=1,NX(IW))
!   WRITE(8,*)LW,IW,(RH(LW,I),I=1,NX(IW))
!   WRITE(8,*)LW,IW,(WPU(LW,I),I=1,NX(IW))
!   WRITE(8,*)LW,IW,(XN(LW,I),I=1,NX(IW))
!END DO

!
!*******************************************
PRINT *, MOJI,WPW(1,1),WGT(1,1),XN(1,1),XO,SCWDAY(1,1)
PRINT *,XO,AXN(1,1),XM(1,1),ATMPQ(1,1),HWPUR(1),DGDUR(1)
555 CONTINUE
CLOSE(8)
CLOSE(7)
CLOSE(12)
CLOSE(13)
PAUSE
END


SUBROUTINE SATUWPT(TP,SW)  ! cpw[J/kgK],sw[J/kg]
IMPLICIT REAL*8(A-H,O-Z)
 CPW=(30.36+0.009615*TP+0.00000118*TP*TP)/(0.018016)
 CALL GOFF(TP,FS,VP)
 SW=(644243)+CPW*(TP-273.15)-TP*CPW*DLOG(TP/273.15)+461.5*TP*DLOG(VP/101325)
 RETURN
END
!
SUBROUTINE GOFF(TMP,FS,VP)  ! TP(K),FS(mmHg),VP(Pa)
 IMPLICIT REAL*8(A-H,O-Z)
 EW=-6096.938/TMP+21.2409642-TMP*2.711193E-2+TMP*TMP*1.673952E-5+2.433502*LOG(TMP) 
 VP=EXP(EW)  
 FS=VP/133.322 
 RETURN 
END         
!
SUBROUTINE FUNCX(FS,RH,X)
IMPLICIT REAL*8(A-H,O-Z)
 FD=FS*RH*.01
 X=FD*.6217/(760.-FD)
 RETURN      
END
!
SUBROUTINE REWPT(RH,WP,TP,RG)
IMPLICIT REAL*8(A-H,O-Z)
 RH0=RH*0.01
 WP=RG*TP*DLOG(RH0)
 RETURN
END
!
SUBROUTINE WPTRE(RH,WP,TP,RG)
IMPLICIT REAL*8(A-H,O-Z)
 D1=WP/RG/TP
 D2=DEXP(D1)
 RH=D2*100.
 RETURN
END
!
SUBROUTINE CALDPDU(WPT,TP,DPDU,GMA,ML0,RG,ROW)
IMPLICIT REAL*8(A-H,O-Z)
!  DW=5.D+3
   D1=WPT
   IF(D1>0.) D1=-100.  !-1.D-3
    DW=ABS(D1*0.01)
    W1=D1+DW
    W2=D1-DW
!   ELSE
!    DW=50
!    W1=-50
!    W2=-150
!   END IF
!  D1=W1-W2
!   IF(W1.GT.-100.)W1=-100.
  CALL WPTRE(RH1,W1,TP,RG)
  CALL WPTRE(RH2,W2,TP,RG)
  CALL AHGANS(WD1,RH1,ML0)
  CALL AHGANS(WD2,RH2,ML0)
  VGT1=0.01*WD1*GMA/ROW
  VGT2=0.01*WD2*GMA/ROW
!   DPDU=(VGT1-VGT2)/(W1-W2)
  DPDU=0.5*(VGT1-VGT2)/DW
!  D1=(VGT1-VGT2)*1000.
!  D1=D1/(W1-W2)
!  DPDU=D1*0.001
  RETURN
END
!
!
SUBROUTINE CALDGDU(RH,TMP,RG,DGDU,DGDT)
IMPLICIT REAL*8(A-H,O-Z)
   RH1=RH+0.01
   RH2=RH-0.01
   CALL REWPT(RH1,WP1,TMP,RG)
   CALL REWPT(RH2,WP2,TMP,RG)
   CALL GOFF(TMP,FS1,VP)
   X1=VP*RH1*0.01
   X2=VP*RH2*0.01
   DGDU=ABS((X1-X2)/(WP1-WP2))
!
   TP1=TMP+0.1
   TP2=TMP-0.1
   CALL GOFF(TP1,FS1,VP1)
   CALL GOFF(TP2,FS2,VP2)
   X1=VP1*RH*0.01
   X2=VP2*RH*0.01
   DGDT=ABS((X1-X2)/(TP1-TP2))
  RETURN
END


!
     SUBROUTINE AHGANS(WD,RHM,ML0)
     IMPLICIT REAL*8(A-H,O-Z)
      GO TO (11,10,13,14,15,16,17,18,19,20,21,22,13,13,13,13,13,28,29),ML0
!******SEASING BOARD******
   11 IF(RHM.LE.89)WD=-LOG(1-.01*RHM)/.1612944
      IF(RHM.GT.89.AND.RHM.LT.95)WD=639.22-14.988*RHM+.08943*RHM*RHM
      IF(RHM.GE.95)WD=2.446*RHM-209.9
      GO TO 10
!*****GLASS WOOL*****
   13 D1=RHM*.01
      D2=EXP(-1.3016*(1-D1**101.667))
      D3=0.16*D2*D1**7.6448
      IF(D1.GT.0.98)  D3=47.79595*D1-46.79595      !�O�a�ܐ���=100wt%
      !IF(D1.GT.0.98)  D3=4347.796*D1-4260.796      !�O�a�ܐ���=8700wt%
      WD=D3*100.
     !  IF(D1.LT.0.98)THEN
     !   WD= 0.022449*D1*100.
     !  ELSE
     !   WD= (4348.9*D1 - 4261.9)*100.
     ! END IF
      GO TO 10
!*****WOOD*****
   14 D1=RHM*.01
      !D2=EXP(-0.479*(1-D1**15.7))
      !D3=0.36*D2*D1**0.861
      !
     ! D2=EXP(-1.9854*(1-D1**91.36))
     ! D3=1.2725*D2*D1**2.3836
      D2=EXP(-2.1708*(1-D1**90.988))
      D3=1.58*D2*D1**1.3664
      WD=D3*100.
!�����ߒ��Fw=0.36���^0.861exp[-0.479(1 - ��^15.7)]
!�z���ߒ��Fw=0.36���^0.949exp[-0.657(1 - ��^15.2)]
      GO TO 10
!*****PLY WOOD*****
   15 D1=RHM*.01
      !D2=EXP(-0.94044*(1-D1**12.064))
      !D3=0.39932*D2*D1**0.87105
      D2=EXP(-1.6471*(1-D1**32.051))    !ASHRAE PLYWOOD 1
      D3=1.0625*D2*D1**1.5373
      WD=D3*100.
      GO TO 10
!  COEFF={0.39932,0.87105,-0.94044,12.064}
!*****PLASTER BOARD*****
   16 D1=RHM*.01
      D2=EXP(-1.6445*(1-D1**75.811))    !ASHRAE GYPSUM
      D3=0.56978*D2*D1**0.18552
!      D2=EXP(-1.5212*(1-D1**16.289))
!      D3=0.3*D2*D1**0.6741
      WD=D3*100.
      GO TO 10!
!*****ALC***********
   17 D1=RHM*.01
      D2=EXP(-0.66866*(1-D1**17.316))
      D3=0.071101*D2*D1**1.0146
      WD=D3*100.
      GO TO 10
!  COEFF={0.071101,1.0146,-0.66866,17.316}
!*****NANSHITSUSENIBAN*****
   18 D1=RHM*.01
      D2=EXP(-0.93059*(1-D1**4.6145))
      D3=0.35093*D2*D1**0.48739
      WD=D3*100.
      GO TO 10
!*****THERMOPLY*****
   19 D1=RHM*.01
      D2=EXP(-1.2566*(1-D1**5.483))
      D3=0.4367*D2*D1**0.3513
      WD=D3*100.
      GO TO 10
!*****SAIDHING*****
   20 D1=RHM*.01
      D2=EXP(-1.05*(1-D1**11.1))
      D3=0.34*D2*D1**0.786
!     ASHRAE Siding No.38    0.635	2.5972	-1.1225	15.878
      D2=EXP(-1.1225*(1-D1**15.878))
      D3=0.635*D2*D1**2.5972
      WD=D3*100.
      GO TO 10
!*****tutikabe*****
   21 D1=RHM*.01
      D2=EXP(-0.69*(1-D1**25.06))
      D3=0.0645*D2*D1**0.658
      WD=D3*100.
      GO TO 10
!*****�y�ʃ����^��*****
   22 D1=RHM*.01
      D2=EXP(-1.75*(1-D1**2.015))
      D3=0.203*D2*D1**(-0.120)
      WD=D3*100.
      GO TO 10
!*****�W���ށiOM���-�j*****
   28 D1=RHM*.01
      D2=EXP(-0.849*(1-D1**9.527))
      D3=0.381*D2*D1**(0.738)
      WD=D3*100.
      GO TO 10
!*****�\���p���iOM���-�j*****
   29 D1=RHM*.01
      D2=EXP(-0.907*(1-D1**9.412))
      D3=0.382*D2*D1**(0.670)
      WD=D3*100.
      GO TO 10
!
   10 CONTINUE
      RETURN      
      END 

!
!::::::::�T�u���[�`���ǉ�
      SUBROUTINE CWIF1(T1,X1,X0,MCW)
      IMPLICIT REAL*8(A-H,O-Z)
      RH0=100.
      CALL GOFF(T1,FS0,VP)
      CALL FUNCX(FS0,RH0,X0)
      IF(X1.GE.X0.AND.MCW.EQ.0)MCW=1
      RETURN
      END

     SUBROUTINE DOSU(NYEAR,MON,I,D1)
     IMPLICIT REAL*8(A-H,O-Z)				
      COMMON /DA/DSU(2,60,20),D(2,60),NSU(2,60),DG(2,60),DN(2,60),DNG(2,60)				
      J=MON+12*(NYEAR-1)				
      D(I,J)=D1				
	IF(D(I,J).GT.0.AND.D(I,J).LE.5)  DSU(I,J,1)=DSU(I,J,1)+1			
	IF(D(I,J).GT.5.AND.D(I,J).LE.10)  DSU(I,J,2)=DSU(I,J,2)+1			
	IF(D(I,J).GT.10.AND.D(I,J).LE.15)  DSU(I,J,3)=DSU(I,J,3)+1			
	IF(D(I,J).GT.15.AND.D(I,J).LE.20)  DSU(I,J,4)=DSU(I,J,4)+1			
	IF(D(I,J).GT.20.AND.D(I,J).LE.25)  DSU(I,J,5)=DSU(I,J,5)+1			
	IF(D(I,J).GT.25.AND.D(I,J).LE.30)  DSU(I,J,6)=DSU(I,J,6)+1			
	IF(D(I,J).GT.30.AND.D(I,J).LE.35)  DSU(I,J,7)=DSU(I,J,7)+1			
	IF(D(I,J).GT.35.AND.D(I,J).LE.40)  DSU(I,J,8)=DSU(I,J,8)+1			
	IF(D(I,J).GT.40.AND.D(I,J).LE.45)  DSU(I,J,9)=DSU(I,J,9)+1			
	IF(D(I,J).GT.45.AND.D(I,J).LE.50)  DSU(I,J,10)=DSU(I,J,10)+1			
	IF(D(I,J).GT.50.AND.D(I,J).LE.55)  DSU(I,J,11)=DSU(I,J,11)+1			
	IF(D(I,J).GT.55.AND.D(I,J).LE.60)  DSU(I,J,12)=DSU(I,J,12)+1			
	IF(D(I,J).GT.60.AND.D(I,J).LE.65)  DSU(I,J,13)=DSU(I,J,13)+1			
	IF(D(I,J).GT.65.AND.D(I,J).LE.70)  DSU(I,J,14)=DSU(I,J,14)+1			
	IF(D(I,J).GT.70.AND.D(I,J).LE.75)  DSU(I,J,15)=DSU(I,J,15)+1			
	IF(D(I,J).GT.75.AND.D(I,J).LE.80)  DSU(I,J,16)=DSU(I,J,16)+1			
	IF(D(I,J).GT.80.AND.D(I,J).LE.85)  DSU(I,J,17)=DSU(I,J,17)+1			
	IF(D(I,J).GT.85.AND.D(I,J).LE.90)  DSU(I,J,18)=DSU(I,J,18)+1			
	IF(D(I,J).GT.90.AND.D(I,J).LE.95)  DSU(I,J,19)=DSU(I,J,19)+1			
	IF(D(I,J).GT.95.AND.D(I,J).LE.500)  DSU(I,J,20)=DSU(I,J,20)+1			
	NSU(I,J)=NSU(I,J)+1			
	DG(I,J)=DG(I,J)+D(I,J)			
	DN(I,J)=D(I,J)**2			
	DNG(I,J)=DNG(I,J)+DN(I,J)
!        PRINT *, I,J,NSU(I,J),DNG(I,J)
!        PAUSE
	RETURN			
	END			
!				
!				
       SUBROUTINE SIGMA(NYEAR)				
       IMPLICIT REAL*8(A-H,O-Z)
       COMMON /DA/DSU(2,60,20),D(2,60),NSU(2,60),DG(2,60),DN(2,60),DNG(2,60)  				
       DIMENSION DGN(2,60),AS(2,60),V(2,60),SS(2,60)  				
       L1=(NYEAR-1)*12+1				
       L2=(NYEAR-1)*12+12			
       DO 300 I=1,2				
       DO 300 K=L1,L2
       IF(NSU(I,K).EQ.0)GO TO 300
!       PRINT *, I,K,NSU(I,K),DGN(I,K),DNG(I,K),AS(I,K),V(I,K)
!       PAUSE
       DGN(I,K)=DG(I,K)*DG(I,K)
       AS(I,K)=DNG(I,K)-(DGN(I,K)/NSU(I,K))
       V(I,K)=AS(I,K)/(NSU(I,K)-1)
       SS(I,K)=SQRT(V(I,K))
  300 CONTINUE				
      WRITE(12,*)' '				
      DO 350 I=1,2				
      WRITE(12,*)' I=',I,'  NSU  SUM  HEIHOUWA   BUNSAN   HENSA  '
      PRINT *,' I=',I,'  NSU  SUM  HEIHOUWA   BUNSAN   HENSA  '				
      DO 350 J=L1,L2				
       WRITE(12,*)I,NSU(I,J),DG(I,J),AS(I,J),V(I,J),SS(I,J)
       PRINT *,I,NSU(I,J),DG(I,J),AS(I,J),V(I,J),SS(I,J)				
  350 CONTINUE				
      WRITE(12,*)' '				
      WRITE(12,*)'  DOSUU BUNNPU'				
      DO 360 I=1,2				
      WRITE(12,*)' '				
      WRITE(12,*)' I=',I,'  MON '				
      WRITE(12,*)' HNAI 1  2  3  4  5  6  7  8  9  10  11  12'				
      DO 360 K=1,20				
        WRITE(12,*)K,(DSU(I,J,K),J=L1,L2)
         PRINT *,K,(DSU(I,J,K),J=L1,L2)				
  360 CONTINUE				
      RETURN				
      END		
 !
!********YOMITOBASI******************				
      SUBROUTINE YOMI(MON1,MDAY1,KDAY)
      CHARACTER MOJI*50				
      DIMENSION MDAY(12)				
      DATA MDAY/31,28,31,30,31,30,31,31,30,31,30,31/				
      L1=0				
      KDAY=0				
      DO 60 I = 1,MON1-1				
        L1=MDAY(I)+L1				
   60 CONTINUE      				
      L2=MDAY1+L1-1				
      KDAY=L2				
      L3=L2*25				
!      PRINT *,' YOMITOBASI=',L3				
      DO 70 I = 1,L3+1				
        READ(7,*)MOJI
   70 CONTINUE
      L4=L2*24
      DO I=1,L4
        READ(13,*)D1,D2
      END DO
      PRINT *, '  MOJI ',MOJI,D1,D2				
!  103 FORMAT(A)				
      RETURN				
      END				
!********************************				
!      								
!      				
      SUBROUTINE SOLOCT(KW,SIDO,SKEIDO,HOI,AKASYA,NORIENT,L)
     IMPLICIT REAL*8(A-H,O-Z)
     INTEGER,PARAMETER::NWP=50				
!     REAL KEIDO,IDO			! COUTION!!!	
      COMMON /SS/W(30),DV4(NWP,30),BKAKU				
      DIMENSION V2(30),V3(30)				
!     DIMENSION DW(4,30),DV2(4,30),DV3(4,30)
      DATA C1,C2,C3,C4,C5,C6,C7/0.006322,0.405748,0.153231,&				
     0.00588,0.207099,0.003233,0.620129/				
      DATA SD1,SD2,SD3,SD4,SD5,SD6,SD7/0.000279,0.122772,1.498311,&				
     0.165458,1.261546,0.005354,1.1571/				
!      IF(IWEEK.LE.1) THEN				
!       GO TO 10				
!        ELSE				
!         GO TO 20				
!      END IF				
!  10 CONTINUE				
! *****calculate solar locasion subroutine***** 				
      ENSYU=3.1415				
      HOIM=HOI*ENSYU/180.				
      BKAKU=AKASYA*ENSYU/180.				
!     FAI=HKEIDO*ENSYU/180.				
      HIDO=SIDO*ENSYU/180.				
      OMEGA=2*ENSYU*KW/366.				
      SEKI=C1-C2*COS(OMEGA+C3)-C4*COS(2.*OMEGA+C5)				
      SEKI=SEKI-C6*COS(3.*OMEGA+C7)				
      EKNJ=-SD1+SD2*COS(OMEGA+SD3)-SD4*COS(2.*OMEGA-SD5)				
      EKNJ=EKNJ-SD6*COS(3.*OMEGA-SD7)				
      DO 30 I = 4, 20				
        K=I+1				
        TKU=ENSYU*((I+EKNJ-12)+(SKEIDO-135.)/15.)/12.				
        W(K)=SIN(HIDO)*SIN(SEKI)+COS(HIDO)*COS(SEKI)*COS(TKU)				
        V2(K)=COS(SEKI)*SIN(TKU)				
        V3(K)=-SIN(SEKI)*COS(HIDO)+COS(SEKI)*SIN(HIDO)*COS(TKU)				
        IF(V3(K).LT.0) V3(K)=0				
	IF(W(K).LT.0)W(K)=0.			
        HKODO=ASIN(W(K))				
!       DO 30 J=1,4				
!	E1=(J-1)*90.			
!        IF(NORIENT.LE.4)THEN
          J=NORIENT
!        ELSE
!          DV4(L,K)=0.
!         GOTO 30
!        END IF
        E1=(J-1)*90
	E2=E1*ENSYU/180.			
	HOIM2=HOIM+E2			
        S1=V2(K)/COS(HKODO)				
        S3=V3(K)/COS(HKODO)				
!        S2=S1*COS(HOIM2)-S3*SIN(HOIM2)				
        S4=S3*COS(HOIM2)+S1*SIN(HOIM2)				
!        DW(J,K)=W(K)*COS(BKAKU)+COS(HKODO)*SIN(BKAKU)*S4				
!        DV2(J,K)=COS(HKODO)*S2				
!        DV3(J,K)=-(W(K))*SIN(BKAKU)+COS(HKODO)*COS(BKAKU)*S4				
	DV4(L,K)=COS(BKAKU)*SIN(HKODO)+SIN(BKAKU)*COS(HKODO)*S4			
	IF(DV4(L,K).LE.0.)DV4(L,K)=0.			
   30 CONTINUE				
   20 CONTINUE				
!      IWEEK=IWEEK+1				
!      IF(IWEEK.GE.7) IWEEK=1				
!      WRITE(8,*)SEKI,EKNJ,TKU,SIDO,SKEIDO,OMEGA,HOI,AKASYA,HOIM,BKAKU,HIDO
      RETURN				
      END				
      				
      SUBROUTINE SATCAL(RLF,AS,EMI,L)				
!      				
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER,PARAMETER::NWP=50
      COMMON /MS/TEMPO(25),SJD(25),SJS(25),SJN(25),SAT(NWP,25),SJIN(NWP,25),WSJIN(NWP,25),WSJD(NWP,25)
      COMMON /SS/W(30),DV4(NWP,30),BKAKU				
      DATA AFO/22.4/				
      DO 11 I=2,25				
  !     DO 10 J=1,NSAT				
        D1=(1+COS(BKAKU))/2				
        DSRS=SJD(I)*DV4(L,I)				
	SRS=D1*SJS(I)			
        DSRH=SJD(I)*W(I)				
        HNISYA=SJS(I)+DSRH				
        RSRS=(1-(1+COS(BKAKU))/2)*RLF*HNISYA				
        SNISYA=DSRS+SRS+RSRS				
        SAT(L,I)=TEMPO(I)+(1/AFO)*(AS*SNISYA-D1*SJN(I)*EMI)
        SJIN(L,I)=AS*SNISYA-D1*SJN(I)*EMI
!	WSJIN(L,I)=SNISYA
        WSJIN(L,I)=AS*(SRS+RSRS)
	WSJD(L,I)=AS*DSRS
 !  10  CONTINUE
      ! WRITE(8,*)I,(SAT(J,I),J=1,4),(DV4(J,I),J=1,4),W(I),DSRS,SRS,DSRH,RSRS,SNISYA
   11 CONTINUE		
      RETURN	
      END				
      				
!//////////////////////////
SUBROUTINE DIFF(K,RH,TMP,RG,RML)
! calculation liquid conductance
!Y=e^(a+bX+cX2+dX3+eX4+fX5+gX6)
IMPLICIT REAL*8(A-H,O-Z)
!
 CALL AHGANS(WD,RH,K)
 IF(WD>45.AND.WD<110)THEN
   D1=WD*0.01
   DW=EXP(-30.91+4.2967*D1-0.22017*D1*D1)
   RH1=RH+0.005
   RH2=RH-0.005
   CALL REWPT(RH1,WP1,TMP,RG)
   CALL REWPT(RH2,WP2,TMP,RG)
   CALL AHGANS(WD1,RH1,K)
   CALL AHGANS(WD2,RH2,K)
   RML=998.*DW*ABS((WD1-WD2)*0.01/(WP1-WP2))
  ELSE
   RML=0.
 END IF
RETURN 
END 
 !    
      SUBROUTINE ROOM(KDAY,TR,RR,XR,TRAV,TRDT,RRAV)
      IMPLICIT REAL*8(A-H,O-Z)
      D1=2*3.1415*(KDAY-212)*24./8760
      TR=TRDT*COS(D1)+TRAV
      RR=RRAV
      T0=TR+273.15				
      CALL GOFF(T0,FS0,VP)				
      !CALL FUNCX(FS0,RR,X0)				
      XR=VP*RR*0.01				
      RETURN				
      END
!
!!//////////////////////////////////////////////
! Biological damage function for wood rot decay 
! File name:hdamage.f95
! Reaction model
! Hiroaki Saito
! ////////////////////////////////////////////
SUBROUTINE wood rot(K,I,TMP,RH,DLOSS,ROTOMG)
!
!
IMPLICIT REAL*8(A-H,O-Z)
!INTEGER,PARAMETER :: NXP=150,NYP=150
INTEGER,PARAMETER :: NXP=150,NWP=50  !,NYP=150,KATLP=70
DIMENSION TIME_S(NWP,NXP),L_STAGE(NWP,NXP)
!
DT_DAMAGE=60*60*24  !Unit:sec
W=-1.   ! OSB
RH_GROWTH=98.
DLOSS=0.
!
!
! �����  Initial response time

IF(TMP <= 0.)TIME_S(K,I)=0.
IF(TMP > 0.) THEN
  CALL  RH_CRITICAL(TMP,RHC)
  IF(RH<RHC)TIME_S(K,I)=0.
  IF(RH>=RHC)THEN
    CALL TIN(TMP,RH,W,GC,FC)
     D1=-GC/FC
    !D1=2.0
     IF(D1>0)THEN
	TIME_INT=D1
	ELSE
        IF(D1>-0.59)THEN
	   TIME_INT=0.
	  ELSE
	   TIME_INT=100.
	END IF
     END IF
     IF(TIME_INT<=0.5)THEN
	TIME_INT=0.5*DT_DAMAGE*30.
	ELSE
	TIME_INT= TIME_INT*DT_DAMAGE*30.
     END IF
     TIME_S(K,I)=TIME_S(K,I)+DT_DAMAGE
     IF(TIME_S(K,I) > TIME_INT)L_STAGE(K,I)=1
  END IF
END IF

!������  Growth Stage
IF(TMP>0.AND.TMP<=40)THEN
  IF(L_STAGE(K,I)==1.OR.L_STAGE(K-1,I)==1.OR.L_STAGE(K+1,I)==1)THEN
     IF(RH.GE.RH_GROWTH)THEN
	REACTION_K=(2.77-3.23*TMP+0.865*TMP*TMP -0.0189*TMP**3)*1.D-10*ROTOMG
	DLOSS=REACTION_K*DT_DAMAGE
     END IF
  END IF
END IF
!WRITE(12,*)I,TMP,RH,RHC,TIME_S(K,I),TIME_INT,L_STAGE(K,I),REACTION_K,DLOSS
RETURN 
END
	
!////////////////////////
!////////////////////////

SUBROUTINE RH_CRITICAL(TMP,RHC)
 IMPLICIT REAL*8(A-H,O-Z)
 !IMPLICIT NONE
 ! REAL, INTENT(in) :: T;
 ! REAL :: rot_critical_relative_humidity;
  IF (TMP<=15.)RHC = -0.5*TMP+100.
  IF (TMP>15.)RHC= 92.5
  IF( RHC >= 100.) RHC=100.
RETURN
END

SUBROUTINE TIN(T,RH,W,GC,FC)
 IMPLICIT REAL*8(A-H,O-Z)
 FC= 0.1384*T+0.4370*RH-42.9450+W*(0.034*T-0.021*RH+1.721)
 GC= -2.2270*T-0.0347*RH+0.0244*T*RH + W*(-0.504*T+0.0096*RH+0.0047*T*RH)
RETURN
END 



