import numpy as np
from dataclasses import dataclass

from modules.room import Room, InputRoom
from modules.thermo_dynamics import ATP, get_wp, RG, GOFF, FUNCX, get_dgdu, ROW_CP, WPTRE, DIFF, AHGANS, CPL, ROW, get_rho
from modules.config import TMPOC, HTC_TW, HOI
from modules.wall import Wall, WallType, Layer
from modules.solar_position import get_solar_position
from modules.weather import Weather
from modules.date_operation import Period
from modules.nrain import NRAINPOINT, ACOFRAIN, set_default_nrains
from scipy.stats import rayleigh
from modules.date_operation import get_step_d_t, get_step_n


INTEGER,PARAMETER :: NXP=150,KMTLP=10,NWP=50,NRM=50

# サブルーチンで内部保持しているパラメータ（time_s, l_stage）があったため、クラス化した。
wdm = WoodDecayModel(nwp=NWP, nxp=NXP)


INTEGER,PARAMETER :: MXIT=5000
REAL,PARAMETER :: EPS1=0.01,EPS2=1.0e+2
INTEGER DAY,DAY2,DAY1


!//////////COEFFOCIENT//////////////
DIMENSION NX(NWP)
DIMENSION WGT(NWP,NXP)
DIMENSION KMTL(KMTLP),WOUTAV(50)

!/////////VARIABLE///////////////
DIMENSION WPU(NWP,NXP),WPS(NWP,NXP),WPW(NWP,NXP),TMP(NWP,NXP)
DIMENSION XN(NWP,NXP)
DIMENSION ATMP(NWP,NXP),BWPU(NWP,NXP),BTMP(NWP,NXP)
DIMENSION MCW(NWP)
DIMENSION QQ(NWP,NXP),XM(NWP,NXP),RHDIS(NWP,NXP,12,31)
DIMENSION RHAVD(NWP,NXP),WGTAVD(NWP,NXP)
DIMENSION ATMPQ(NWP,NXP)
DIMENSION TPAVD(NWP,NXP),XNAVD(NWP,NXP),TPDIS(NWP,NXP,12,31)
DIMENSION WLOSS(NWP,NXP),WJW(NWP,NXP)
DIMENSION WJRAIN(NWP,NXP),RAIN(24)
DIMENSION RN(NWP,NXP),HRN(NWP,NXP)

DATA OMG/1.2/

CHARACTER MOJI*110, FILENAME0*12

PRINT *,  '　収束計算あり=1   収束計算なし=0 '
READ(5,*)NREPT
PRINT *,  '腐朽緩和係数を入力して下さい '
READ(5,*)ROTOMG
PRINT *,  '水分生成の扱い　無視0　考慮1 '
READ(5,*)I_HCOFF


# 水分化学ポテンシャルによる建物の温湿度計算

# 緯度
HIDO = 36.7

# 経度
HKEIDO = 137.21


# 時間分割（1/DT {h}）
NDVD = 20


ipt_room = InputRoom.read(d={})
room = Room.init(ipt_room=ipt_room, n_div=NDVD)



# 壁体部位数 NWTYPE 壁TYPE数
NWTYPE = 2

# 部材数
KMTL = [7, 7]

nrains = set_default_nrains()


walls = Wall.read_default()


# :::::::::結露計算部位入力:::::::

# 結露計算質点
# 出力質点数
NSUWC = 4
# 窓部位
NWCN = [1, 2, 3, 4]

ALDC = 5./3600.

# :::::出力座標::::::::::::::::::::::::: 
# 出力質点
NOUT = 12

NOUTD = [[1, 6], [1, 12], [1, 17], [1, 22], [4, 6], [4, 12], [4, 17], [4, 22], [7, 6], [7, 12], [7, 17], [7, 22]]

# :::::平均値出力材料::::::::::::::::::::
# 平均値出力材料
NOUTAV = 9
NOUTAVD = [[1, 1], [1, 3], [1, 4], [4, 1], [4, 3], [4, 4], [7, 1], [7, 3], [7, 4]]

# :::::腐朽部位出力座標::::::::::::::::::::
# 腐朽部位出力質点
NOUTML = 3
NOUTMLD = [[1, 12, 17], [4, 12, 17], [7, 12, 17]]

dt = 1. / NDVD * 3600.
PXCOF = 1. / 133322.


K = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]


# *********CAL DX,NAFX ***************
for I, walltype in enumerate(walltypes):
    K=0

    for J, layer in enumerate(walltype.layers):

        if J == 1:
            NALX(I, J) = layer.n_div
        else:
            NALX(I, J) = NALX(I, J-1) + layer.n_div

    for J, layer in enumerate(walltype.layers):

        if J == 1:
            NAFX(I, J) = 1
        else:
            NAFX(I, J) = NALX(I, J-1)+1

    NX(I) = walltype.n_div_total

# *******SAT   SINITIAL VALUE *******

RH0 = walltypes[0][0].initial_humidity


# /////////INITIAL CONDITION//////////////

for L, wall in enumerate(walls):

    IW = wall.kwtype

    for K, layer in enumerate(wall.layers):

        TP1 = layer.initial_temp + ATP
        RH0 = layer.initial_humidity
        WP0 = SATUWPT(TP1)
        WP = get_wp(RH0, TP1)
   
        FS0, VP = GOFF(TP1)
		
        WP1=WP
        WP2=WP0 + WP1
        QQ(L,K)=0.

        WPW(L,I) = WP2
        WPU(L,I) = WP1
        WPS(L,I) = WP0
        # TMP(L,I) = TP1
        ATMPQ(L,I) = TP1 - ATP
        XN(L,I) = VP * RH0 * 0.01
        XM(L,I) = VP * RH0 * 0.01
        ATMP(L,I) = TP1
        RHAVD(L,I) = 0.
        WGTAVD(L,I) = 0.
        XNAVD(L,I) = 0.
        TPAVD(L,I) = 0.
        WD = AHGANS(rhm=RH0, ml0=L5)
        WGT(L,I) = WD
        WLOSS(L,I) = 0.

for L in range(NOUTAV):
    I = NOUTAVD(L,1)
    J = NOUTAVD(L,2)
    IW = walls[I].kwtype
    L1 = NAFX(IW,J)
    WOUTAV(L) = WGT(I,L1)


### プログラム ここから ###

period = Period()

weather = Weather.load_data(file_name='weather.csv', longitude=HIDO, latitude=HKEIDO)


t_ws_is = [np.zeros(wall.n_mesh_total, dtype=float) for wall in walls]
wp_ws_is = [np.zeros(wall.n_mesh_total, dtype=float) for wall in walls]

#DO 444
for NYEAR in range(1, period.LYEAR + 1):

    KDAY=0				
    
    MONF = period.get_start_month(year=NYEAR)
    LMONL = period.get_last_month(year=NYEAR)

    IMQ=0

    # 999
    for MON in range(MONF, LMONL + 1):
    
        DAY1 = period.get_start_day(year=NYEAR)
        DAY2 = period.get_end_day(year=NYEAR, month=MON)

        # 三浦コメント：ここから１日のループ計算
        # 888
        for DAY in range(DAY1, DAY2):
    
            KDAY=KDAY+1

            # 777
            for IM in range(24):

                # //////////////////////////////////////////////////////////
                IMM=IM+1 
                IMQ=IMQ+1

                # 三浦コメント：おそらくここから1時間の分割す数に応じたループがはじまる。
                # 666
                for n in range(NDVD):
                    NF = n + 1

                    frac = n / NDVD

                    oc = weather.get_condition(month=MON, day=DAY, hour=IM, frac=frac)


                    # ////////PUT IN PREBIOUS VALUE///////////

                    for n_LW, wall in enumerate(walls):
                        LW = n_LW
                        IW=walls[LW].kwtype

                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2

                                wall.t_n_is = wall.t_n_pls
                                # 1ステップ前の値を入れ替えている。
                                XM(LW,I)=XN(LW,I)
                                HRN(LW,I)=RN(LW,I)
                            END DO
                        END DO
                    END DO

                    !///////////////COMBAIN METHIOD//////////////////////////
                    IT3=0
                    810 AMAX3=0.
                    AMAX4=0.

                    !////////////CAL_TEMP BY OVER RELAXATION METHOD///////////
                    ITQ=0
                    IT4=IT4+1

                    !***** conbained  temp and ventiration****
    
                    # GOTO 文でここに戻ってくる
                    217 CONTINUE

                    ITQ=ITQ+1
                    IF(ITQ.GT.100)THEN
                        PRINT *, '  no convagence   ventiration'
                        GO TO 218
                    END IF
                    !*****************************************
                    
                    n = get_step_n(month=month, day=day, hour=hour, n_hour=n, n_div=NDVD)

                    # GOTO 文でここに戻ってくる
                    510 AMAX1=0. 

                    t_ws_is_next = []

                    delta_ws = []

                    for w, wall in enumerate(walls):

                        t_is = t_ws_is[w]

                        t_is_next = wall.get_t_n_pls(t_is=t_is, dt=dt, oc=oc, theta_r_n=room.theta_n(n=n), wp_r_n=room.wp_n(n=n), QQ=QQ, t_upstream=t_upstream)

                        t_ws_is_next.append((t_is_next - t_is) * OMG + t_is)

                        delta_is = t_is_next - t_is

                        delta_ws.append(np.abs(delta_is).max())

                    if max(delta_ws) > EPS1:
                        GO TO 510

                    # *************換気計算収束判定********************
                    # 最初はこの計算を飛ばして、2回目からこの判定を行う。
                    # ITQ：換気と温度の収束計算の繰り返し回数
                    if ITQ > 1:

                        D1=0

                        D2=0

                        for LW, wall in enumerate(walls):

                            IW = wall.kwtype

                            for i in range(wall.n_mesh_total):

                                if wall.is_air_layer_is[i] and wall[K].alpha > 0.0:
                                    
                                    # TMP 絶対温度
                                    D1 = abs(TMP(LW, i)-ATP-ATMPQ(LW, i))
                                    if D1 > D2:
                                        D2 = D1

                                ATMPQ(LW, i)=TMP(LW, i)-ATP

                        IF(D2.LT.EPS1)GO TO 218


                    # *************換気量の算出(Q=m3/s)****************
                    rho_o = get_rho(t=oc.t_k)
                    
                    # 換気量, m3/s
                    QQ = []
                    for wall in walls:
                        qq = 0.0
                        for i in range(wall.n_mesh_total):
                            if wall.is_air_layer_is[i]:
                                rho_i = get_rho(t=wall.t_n_is[i])
                                qq += wall.alpha_a_ls * (2 / rho_o * abs(rho_o - rho_i) * 9.8 * wall.height) ** 0.5
                        QQ.append(qq)


                    # 三浦コメント：ここで強制的に前の方に戻される。
                    # 戻された部分からここまでの間に、何かの条件を満たせばGOTO218で飛び、ループを抜けられる。
                    GOTO 217

                    218 CONTINUE

                    # ************END VENT****************************

                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                TP=TMP(LW,I)
                                SW = CALL SATUWPT(TP)
                                WPS(LW,I)=SW
                            END DO
                        END DO
                    END DO

                    # ////////////CAL_WATER POTENTIAL BY OVER RELAXATION METHOD/////////
                    IT2=0
                    CWMAX=2.8E-2

                    ITC=0
                    710 AMAX2=0.c_liquid_i_mns 

                    # *************************************
                    # ************水膜の水分保持量及び吸水量の計算**********
                    IF(NRAINS.GT.0)THEN
                        DO K=1,NRAINS
                            # 壁の番号
                            I=NRAINPOINT(K,1)
                            # 質点の番号
                            J=NRAINPOINT(K,2)
                            IW=walls[I].kwtype
                            L1=walltypes[IW][NRAINPOINT(k,3)].num
                            # 材料番号
                            L5=walltypes[IW][L1].num

                            # NALX:室内側の接点番号
                            IF(J.EQ.NALX(IW,NRAINPOINT(K,3)))THEN     !  水膜との隣接質点の選択
                                D2=WPU(I,J-1)
                                D3=XM(I,J+1)
                            ELSE
                                D2=WPU(I,J+1)
                                IF(J.EQ.1)THEN
                                    D3=oc.xod                              !外装材表面、外気との収支
                                ELSE
                                    D3=XM(I,J-1)
                                END IF
                            END IF

                            IF(walls[I].get_swjrain(oc=oc, i=J) >= 0.0 or RN(I,J).GE.0)THEN
                                t_srf=TMP(I,J)
                                FS, VP = GOFF(t_srf)
                                RH1 = walls[I].rh[J]

                                DGDU = get_dgdu(rh=RH1, t=t_srf)
                                DGDT = get_dgdt(rh=RH1, t=t_srf)

                                # 湿気伝達率 kg/(m2 s Pa))
                                D4=3.43E-08*(D3-VP)*0.3 # ****濡れ面率0.3  水膜からの蒸発量

                                IF(L5.GE.10.AND.L5.LE.12)THEN                                 !バックシーラー透水抵抗 2.4e+5 m2sPa/kg by　長村
                                    # 3.73：水分伝導率（材料番号が10～12）（セメント系材料：サイディングとか）　kg/ms(J/kg)
                                    # 2.4：バックシーラー　塗膜が塗ってある　塗膜の抵抗　m2sPa/kg
                                    # 抵抗値にして逆数にしてコンダクタンスになおしている。
                                    c_liquid_i_pls=1/( walltypes[IW].layers[J].dx /3.73E-6+2.4E+5/DGDU)                        !飽和時の水分伝導率 3.73e-6 kg/ms(J/kg)　いぶし瓦 by 伊庭　D論
                                ELSE
                                    # 木製品を想定　木製品の水分伝導率は、3.73e-6 kg/ms(J/kg)　いぶし瓦 by 伊庭　D論
                                    # 瓦の5%になっている。（根拠はいまのところ不明）
                                    c_liquid_i_pls=1/( walltypes[IW].layers[J].dx /(3.73E-6*0.05)+2.4E+5/DGDU)                 !木製品を想定
                                END IF

                                # 水幕の水分量（保持している水の量）, kg/m2
                                # D2 はポテンシャル
                                # D4 は蒸発量
                                # HRN は前のステップの水分量
                                RN(I,J)=(c_liquid_i_pls*D2+D4+walls[I].get_swjrain(oc=oc, i=J))*3600/NDVD + HRN(I,J)            !******水膜の水分量kg/m2
                                # 水幕からの吸水量（表面に水幕があって材料に吸われる分）
                                # 水幕が残っている場合は飽和水蒸気圧（水分伝導率で計算）

                                WJRAIN(I,J)=-c_liquid_i_pls*D2*walls[LW].area

                                IF(RN(I,J).LE.0)THEN
                                    RN(I,J)=0
                                    # RNがゼロの場合は水幕がないので、雨水が直接材料に吸われることになる。
                                    WJRAIN(I,J) = walls[I].get_swjrain(oc=oc, i=J) *walls[LW].area+HRN(I,J)*walls[LW].area                      !******浸水量WJRAIN kg/s
                                END IF
                            END IF
                        END DO
                    END IF

                    # ****************************************

                    theta_r_n = room.theta_n(n=n)
                    wp_r_n = room.wp_n(n=n)

                    wp_ws_is_next = []

                    dlta_ws = []

                    for w, wall in enumerate(walls):

                        wp_is = wp_ws_is[w]
                        
                        wp_is_next = wall.get_wp_n_pls(
                            dt=dt,
                            oc=oc,
                            theta_r_n=theta_r_n,
                            wp_r_n=wp_r_n,
                            wp_is=wp_is,
                            QQ=QQ[w],
                            RN=RN[w],
                            XM=XM[w],
                            WJRAIN=WJRAIN[w],
                            WJW=WJW[w]
                        )

                        wp_ws_is_next.append((wp_is_next - wp_is) * OMG + wp_is)

                        delta_is = wp_is_next - wp_is

                        delta_ws.append(np.abs(delta_is).max())

                    AMAX2 = max(delta_ws)

                    # wp_ws_is_next が次のWPUのこと。

                    for n_LW, wall in enumerate(walls):
                        LW = n_LW

                        IW=walls[LW].kwtype

                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num

                            DO I=L1,L2
                                WPW(LW,I)=WPU(LW,I)+WPS(LW,I)
                                WP0=WPU(LW,I)
                                TP=TMP(LW,I)
                                RH0 = WPTRE(WP0,TP)
                                wall.rh[I] =RH0
                                FS, VP = GOFF(TP)

                                XN(LW,I)=RH0*VP*0.01

                                WD = AHGANS(rhm=RH0, ml0=L5)
                                WGT(LW,I)=WD
                                wall.set_TMPC(I, TMPC=TMP(LW, I) - ATP)
                            END DO
                        END DO
                    END DO

                    # /////////JUDGEMENT CONVERGENCE//////
                    IT2=IT2+1

                    IF(IT2.GT.MXIT)THEN
                        AMAX2=0.
                    END IF

                    IF(AMAX2.GT.EPS2) GO TO 710

                    # ::::::::窓ガラスの結露判定部、収束判定部に追加
                    IF(ITC.GE.30)THEN
                        PRINT *, ' HASSAN KETURO '
                        PAUSE
                        GO TO 555
                    END IF

                    K=0
                    ITC=ITC+1
    
                    IF(K.EQ.1)GO TO 710
    
                    # ///////JUDGEMENT CONVERGENCE OF HEAT AND MOISTURE/////

                    # 収束計算ありの場合　NREPT=1, 収束計算なしの場合　NREPT=0
                    IF(NREPT.EQ.1)THEN  !////COMBINE

                        DO LW=1,len(walls)
                            IW=walls[LW].kwtype
                            DO K=1,KMTL(IW)

                                L1=NAFX(IW,K)
                                L2=NALX(IW,K)
                                L5=walltypes[IW][K].num

                                DO I=L1,L2
                                    D8=WPU(LW,I)-BWPU(LW,I)
                                    D9=TMP(LW,I)-BTMP(LW,I)
                                    BWPU(LW,I)=WPU(LW,I)
                                    BTMP(LW,I)=TMP(LW,I)
                                    IF(ABS(D8).GT.AMAX3)AMAX3=ABS(D8)
                                    IF(ABS(D9).GT.AMAX4)AMAX4=ABS(D9)
                                END DO
                            END DO
                        END DO

                        # /////////JUDGEMENT CONVERGENCE//////
                        IT3=IT3+1

                        IF(IT3.GT.MXIT)THEN
                            PRINT *, "   NO CONVERGENCE   IM=",IM
                            GO TO 556
    
                        END IF

                        IF(AMAX3.GT.EPS2.OR.AMAX4.GT.EPS1) GO TO 810

                    END IF  !/////COMBINE

                    # 556 CONTINUE

                    # ////////////HEAT FLUX J/m2s  MOISTURE FLUX  g/m2s//////////////////
                    # //////////CAL AVERAGE VALUE//////////////
                    for n_LW, wall in enumerate(walls):
                        LW = n_LW
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                RHAVD(LW,I) = wall.rh[I]  +RHAVD(LW,I)
                                WGTAVD(LW,I)=WGT(LW,I)+WGTAVD(LW,I)
                                XNAVD(LW,I)=XN(LW,I)+XNAVD(LW,I)
                                TPAVD(LW,I) = wall.TMPC[I] + TPAVD(LW,I)
                            END DO
                        END DO
                    END DO

                # 666 CONTINUE

                IF(NOUTAV.GT.0)THEN                           !材料の平均含水率の計算
                    DO L=1,NOUTAV
                        LW=NOUTAVD(L,1)
                        J=NOUTAVD(L,2)
                        IW=walls[LW].kwtype
                        L1=NAFX(IW,J)
                        L2=NALX(IW,J)
                        c_liquid_i_pls=walltypes[IW].layers[J].dx
                        D2=WGT(LW,L1)*0.5
                        D2=D2+WGT(LW,L2)*0.5

                        DO K=L1+1,L2-1
                            D2=D2+WGT(LW,K)
                        END DO

                        WOUTAV(L)=D2/(L2-L1+1)
                    END DO
                END IF

            # 777 CONTINUE

            # *******SAT   IREKAE *******
            DO I=1,len(walls)
                SJIN(I,1)=SJIN(I,25)
            END DO

            # ***************************************
            DO LW=1,len(walls)
                IW=walls[LW].kwtype
                DO K=1,KMTL(IW)
                    L1=NAFX(IW,K)
                    L2=NALX(IW,K)
                    L5=walltypes[IW][K].num
                    DO I=L1,L2
                        RHDIS(LW,I,MON,DAY)=RHAVD(LW,I)/24./NDVD
                        TPDIS(LW,I,MON,DAY)= TPAVD(LW,I)/24./NDVD
                        RHAVD(LW,I)=0.
                        WGTAVD(LW,I)=0.
                        XNAVD(LW,I)=0.
                        TPAVD(LW,I)=0.
                    END DO
                END DO 
            END DO

            # ////////Calculation Mass Loss（日平均値による算出）//////////

            for n_LW, wall in enumerate(walls):
                LW = n_LW + 1

                IW=walls[LW].kwtype

                DO K=1,KMTL(IW)

                    L1=NAFX(IW,K)
                    L2=NALX(IW,K)
                    L5=walltypes[IW][K].num
                    WLOSSMAX=0.6

                    IF(L5.EQ.4.OR.L5.EQ.5)THEN

                        IF(I_HCOFF.EQ.1)THEN
                            HCOFF=0.319  			  !水分生成量
                        ELSE
                            HCOFF=0.
                        END IF
    
                        DO I=L1,L2
                            c_liquid_i_pls=TPDIS(LW,I,MON,DAY)
                            D2=RHDIS(LW,I,MON,DAY)
                            DLOSS = wdm.WOOD_ROT(LW, I, c_liquid_i_pls, D2, ROTOMG)
                            WLOSS(LW,I)=WLOSS(LW,I)+DLOSS
                            IF(WLOSS(LW,I) > WLOSSMAX)DLOSS=0.
                            WJW(LW,I)=HCOFF*DLOSS * wall.gma_is[i] /86400.   ! Time unit:h=24,s=86400
                        END DO 
                    END IF 

                END DO 
            END DO

        # 888 CONTINUE

    # 999 CONTINUE

    REWIND 7
    REWIND 13

# 444 CONTINUE


# *******************************************

555 CONTINUE

CLOSE(8)
CLOSE(7)
CLOSE(12)
CLOSE(13)

PAUSE

END


def SATUWPT(TP)

    # J/(kg K)
    CPW = ( 30.36 + 0.009615 * TP + 0.00000118 * TP * TP ) / ( 0.018016 )

    FS, VP = GOFF(TP)

    # J/ kg
    SW=(644243)+CPW*(TP-273.15)-TP*CPW*DLOG(TP/273.15)+461.5*TP*DLOG(VP/101325)

    return SW


#::::::::サブルーチン追加
def CWIF1(t: float, x_air: float, mcw: int):
    """_summary_

    Args:
        t: 表面温度, K
        x_air: 空気の混合比, kg/kg(DA)
        mcw (_type_): _description_

    Returns:
        _type_: _description_
    """
    
    # 表面の飽和水蒸気圧, Pa
    _, vp = GOFF(t=t)
    
    # 表面の混合比, kg/kg(DA)
    x_srf = FUNCX(rh=100.0, vp=vp)

    # 空気の混合比の方が表面の混合比よりも大きい場合
    if x_air >= x_srf and mcw == 0:
        mcw = 1

    return mcw, x_srf


class WoodDecayModel:
    """
        Biological damage function for wood rot decay 
        File name:hdamage.f95
        Reaction model
        Hiroaki Saito
    """

    def __init__(self, nwp=50, nxp=150):
        # AIコメント：状態を保持するための配列 (FortranのDIMENSION相当)
        self.time_s = np.zeros((nwp + 2, nxp + 2))   # インデックス余裕を持たせる
        self.l_stage = np.zeros((nwp + 2, nxp + 2), dtype=int)

    def WOOD_ROT(self, k, i, tmp, rh, rotomg):
        """
        木材の腐朽被害関数
        k, i: 現在のグリッド等のインデックス
        tmp: 温度
        rh: 相対湿度
        rotomg: 腐朽速度係数
        """
        w = -1.0  # OSB (TIN関数に渡すパラメータと推測)
        rh_growth = 98.0
        dloss = 0.0
        reaction_k = 0.0
        dt_damage = 60 * 60 * 24  # 1日 (秒)

        # --- 発芽期 (Initial response time) ---
        if tmp <= 0.0:
            self.time_s[k, i] = 0.0
        else:
            # 以前定義した RH_CRITICAL を呼び出し
            rhc = RH_CRITICAL(tmp)
            
            if rh < rhc:
                self.time_s[k, i] = 0.0
            else:
                # TIN関数の呼び出し (引数wが必要)
                gc, fc = TIN(tmp, rh, w)
                
                # ゼロ除算回避
                if fc != 0:
                    d1 = -gc / fc
                else:
                    d1 = 100.0 # 仮の大きな値
                
                # 発芽時間の計算 (TIME_INT)
                if d1 > 0:
                    time_int = d1
                elif d1 > -0.59:
                    time_int = 0.0
                else:
                    time_int = 100.0
                
                # 係数調整
                if time_int <= 0.5:
                    time_int = 0.5 * dt_damage * 30.0
                else:
                    time_int = time_int * dt_damage * 30.0
                
                # 時間経過の蓄積
                self.time_s[k, i] += dt_damage
                if self.time_s[k, i] > time_int:
                    self.l_stage[k, i] = 1

        # --- 成長期 (Growth Stage) ---
        if 0 < tmp <= 40:
            # 自身または隣接するノードが発芽しているかチェック
            # (k-1, k+1 の範囲エラーを防ぐため境界チェックが必要)
            if (self.l_stage[k, i] == 1 or 
                self.l_stage[k-1, i] == 1 or 
                self.l_stage[k+1, i] == 1):
                
                if rh >= rh_growth:
                    # 腐朽反応速度の計算
                    reaction_k = (2.77 - 3.23 * tmp + 0.865 * (tmp**2) - 0.0189 * (tmp**3)) * 1e-10 * rotomg
                    dloss = reaction_k * self.dt_damage

        return dloss
	

def RH_CRITICAL(tmp):
    """
    AIのコメント：温度(tmp)に基づき、臨界相対湿度(rhc)を計算する。
    """

    if tmp <= 15.0:
        rhc = -0.5 * tmp + 100.0
    else:
        rhc = 92.5
        
    # 上限を100に制限
    if rhc >= 100.0:
        rhc = 100.0
        
    return rhc


def TIN(t, rh, w):
    """
    AIのコメント：温度(t), 相対湿度(rh), 含水率(?)w に基づき、FCとGCを計算する。
    """
    
    fc = (0.1384 * t + 0.4370 * rh - 42.9450 + w * (0.034 * t - 0.021 * rh + 1.721))
    
    gc = (-2.2270 * t - 0.0347 * rh + 0.0244 * t * rh + w * (-0.504 * t + 0.0096 * rh + 0.0047 * t * rh))
          
    return gc, fc


