import numpy as np
from dataclasses import dataclass

from modules.room import Room, InputRoom
from modules.thermo_dynamics import ATP, get_wp, RG, GOFF, FUNCX, get_dgdu, ROW_CP, WPTRE, DIFF, AHGANS, CALDPDU, CPL, ROW
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
DIMENSION DTX(NWP,NXP+1)
DIMENSION ADTX(NWP,NXP+1)
DIMENSION ALDT(NWP,KMTLP,2)
DIMENSION ALP_TOTAL(KMTLP),TEMP_CAVITY(KMTLP),HIGHT_CAVITY(KMTLP)
DIMENSION KMTL(KMTLP),WOUTAV(50)

!/////////VARIABLE///////////////
DIMENSION WPU(NWP,NXP),WPS(NWP,NXP),WPW(NWP,NXP),TMP(NWP,NXP)
DIMENSION XN(NWP,NXP)
DIMENSION AWPU(NWP,NXP),ATMP(NWP,NXP),BWPU(NWP,NXP),BTMP(NWP,NXP)
DIMENSION HWPU(NWP,NXP),HTMP(NWP,NXP)    !,HWPS(NXP),HWPW(NXP)
! DIMENSION STMP(KMTLP,2),SWPS(KMTLP,2),SWPU(KMTLP,2),SWPW(KMTLP,2)&
DIMENSION QS(NWP,KMTLP,2)
DIMENSION SCWDAY(NWP,25),MCW(NWP),SCW(NWP),CWV(NWP)
DIMENSION SATDV(NWP),QQ(NWP,NXP),XM(NWP,NXP),AXN(NWP,NXP),RHDIS(NWP,NXP,12,31)
DIMENSION RHAVD(NWP,NXP),WGTAVD(NWP,NXP),DGDUQ(NWP,NXP),DGDTQ(NWP,NXP)
DIMENSION ATMPQ(NWP,NXP)
!DIMENSION TDES1(24),TDES2(24),THDES1(24),THDES2(24),TSDES1(24),TSDES2(24),THSDES1(24),THSDES2(24)
DIMENSION TPAVD(NWP,NXP),XNAVD(NWP,NXP),TPDIS(NWP,NXP,12,31),XNDIS(NWP,NXP,12,31)
DIMENSION WLOSS(NWP,NXP),WJW(NWP,NXP)  !DAMAGE FUNC
DIMENSION WJRAIN(NWP,NXP),RAIN(24)
DIMENSION RN(NWP,NXP),HRN(NWP,NXP)
DIMENSION MDAY(12)
DIMENSION WSIND(NWP,25)


DATA OMG/1.2/
DATA MDAY/31,28,31,30,31,30,31,31,30,31,30,31/

!**** FILE OPEN ****

  PRINT *,  '　　気象データファイルを入力して下さい'
  OPEN(UNIT=7,FILE="")

!PRINT *,  '　　物性データファイルを入力して下さい'
OPEN(UNIT=11,FILE="mcoff1.prn")

  PRINT *,  '　　風雨データファイルを入力して下さい'
OPEN(UNIT=13,FILE="")

CHARACTER MOJI*110,OUTFILE*8,FILENAME0*12
READ(13,210)MOJI
READ(13,210)MOJI
PRINT *,  '　出力ファイルを入力して下さい（8バイト以内）'
READ(5,*)OUTFILE
PRINT *,  '　収束計算あり=1   収束計算なし=0 '
READ(5,*)NREPT
PRINT *,  '腐朽緩和係数を入力して下さい '
READ(5,*)ROTOMG
PRINT *,  '水分生成の扱い　無視0　考慮1 '
READ(5,*)I_HCOFF
PRINT *,  '通気層流量出力ポイント　厚壁No及び層番号1,層番号2 '
READ(5,*)NQOUT1
READ(5,*)NQOUT2
READ(5,*)NQOUT3

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


# :::::::::::相当開口面積αAの直列合成:::::::::::

#種類  不要（要修正）
for I in range(len(walltypes)):

    # 層
    for K in range(KMTL(I)):
        c_liquid_i_pls = 0.0
        D2 = 0.0
        if walls[I][K].alpha > 1.E-7:
            # 階
            for LW, wall in enumerate(walls):
                IW = wall.kwtype
                c_liquid_i_pls = c_liquid_i_pls + 1 / (wall.alpha_a_ls[K] * wall.alpha_a_ls[K])
                # 通気層高さの合計
                D2 = D2 + walls[LW].height
            if c_liquid_i_pls > 1E-12:
                # αAの2乗の合計値
                # この場合のαAは下端と上端の合成したαA
                ALP_TOTAL(K)=1/SQRT(c_liquid_i_pls) 
                HIGHT_CAVITY(K)=D2
# 薄壁部位数
NWIN = 4



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
        L1 = wall.first_mesh_indices[k]
        L2 = wall.last_mesh_indices[k]
        L5 = layer.num
        TP1 = layer.initial_temp + ATP
        RH0 = layer.initial_humidity
        WP0 = SATUWPT(TP1)
        WP = get_wp(RH0, TP1)
   
        FS0, VP = GOFF(TP1)
		
        WP1=WP
        WP2=WP0 + WP1
        QQ(L,K)=0.

        for I in range(L1,L2):
            WPW(L,I) = WP2
            WPU(L,I) = WP1
            WPS(L,I) = WP0
            # TMP(L,I) = TP1
            ATMPQ(L,I) = TP1 - ATP
            XN(L,I) = VP * RH0 * 0.01
            XM(L,I) = VP * RH0 * 0.01
            AXN(L,I) = VP * RH0 * 0.01
            AWPU(L,I) = WP1
            ATMP(L,I) = TP1
            RHAVD(L,I) = 0.
            WGTAVD(L,I) = 0.
            XNAVD(L,I) = 0.
            TPAVD(L,I) = 0.
            WD = AHGANS(rhm=RH0, ml0=L5)
            WGT(L,I) = WD
            WLOSS(L,I) = 0.

        for I in range(1, 2):
            QS(L,K,I) = 0.

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

                    # Pa / (J/K)
                    c_liquid_i_pls = get_dgdu(rh=oc.rh, t=oc.t_k)

                    # Pa / K
                    D2 = get_dgdt(rh=oc.rh, t=oc.t_k)
                    
                    # 外気空気の水蒸気圧を変換するための
                    # 1.2: 空気の密度, kg(DA) / m3
                    # 絶対湿度
                    # 通気層の時に使う
                    # 換気によって入ってくる水分量を化学ポテンシャルで処理するための換算
                    # kg(DA) / m3 * kg / kg(DA) * m3 / s = kg(水蒸気) / s
                    # PXCOF: 絶対湿度を圧力にかえる係数
                    # PXCOF = 1. / 133322. 
                    # J = rho * V * dX = (dg/du) * V * dmu + 
                    # 133322: エクセルで絶対湿度と水蒸気圧
                    DGDUO = 1.2 * PXCOF * c_liquid_i_pls

                    DGDTO = 1.2 * PXCOF * D2

                    # ////////PUT IN PREBIOUS VALUE///////////

                    for n_LW, wall in enumerate(walls):
                        LW = n_LW
                        IW=walls[LW].kwtype

                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                HWPU(LW,I) = WPU(LW,I)

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
                    !**** CULCURATION CONDUCTANCE****

                    for n_lw, wall in enumerate(walls):
                        LW = n_lw + 1

                        IW = wall.kwtype

                        for n_k, layer in enumerate(wall.layers):
                            K = n_k + 1

                            L1 = wall.first_mesh_indices[K]
                            L2 = wall.last_mesh_indices[K]
                            L5 = layer.num

                            DO I=L1,L2

                                IF(I.EQ.L1)THEN
                                    # (kg / s) / K = kg/(m2 s Pa) * Pa / K * m2
                                    ALDT(LW,K,1) = wall.layers[k].cond_m_o* wall.dgdt(i) * wall.area
                                END IF

                                IF(I.EQ.L2)THEN
                                    ALDT(LW,K,2) = wall.layers[k].cond_m_i* wall.dgdt(i) * wall.area
                                END IF

                                IF(L5.EQ.2)THEN   !空気層

                                    # 絶対湿度の場合 kg(DA)/m3 m3/s kg/kg(DA) = kg/s
                                    # kg/m3 * 1/pa *  Pa / (J/K) = kg / m3 (J/K) 
                                    DGDUQ(LW,K)=1.2* wall.dgdu(i) *PXCOF
                                    # kg/m3 * 1/pa *  Pa / K = kg / m3 K 
                                    DGDTQ(LW,K)=1.2* wall.dgdt(i) *PXCOF 
                                    CYCLE
                                END IF


                                ADTX(LW,I) = wall.ADTGX(i) + wall.ADTLX(i=I, wpt=WPU(LW,I), tp=TMP(LW,I))
                            END DO

                            DO I=L1+1,L2
                                # 計算した値を平均する方が良い。逆数同士足した値の逆数がベター
                                # 本来であれば、逆数の和の逆数にすべきだが、片方がゼロになる場合はゼロ割の可能性があるので注意が必要。
                                DTX(LW,I)=(ADTX(LW,I)+ADTX(LW,I-1))*0.5
                            END DO
                        END DO
                    END DO

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
                    
                    # GOTO 文でここに戻ってくる
                    510 AMAX1=0. 

                    n = get_step_n(month=month, day=day, hour=hour, n_hour=n, n_div=NDVD)

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

                            for K in range(wall.layers):
                                L1 = NAFX(IW,K)
                                L2=NALX(IW,K)
                                L5=walltypes[IW][K].num
                                if wall.materials_is[i] == 2 and wall[K].alpha > 0.0:
                                    
                                    # TMP 絶対温度
                                    D1 = abs(TMP(LW,L1)-ATP-ATMPQ(LW,L1))
                                    if D1 > D2:
                                        D2 = D1

                                ATMPQ(LW,L1)=TMP(LW,L1)-ATP
                            END DO
                        END DO
                        IF(D2.LT.EPS1)GO TO 218
                    END IF

                    # *************換気量の算出(Q=m3/s)****************
                    DO I=1,NWTYPE
                        DO K=1,KMTL(I)
                            D1=0.
                            L1=NAFX(I,K) !!!
                            L2=NALX(I,K)
                            L5=walltypes[I][K].num

                            for w, wall in enumerate(walls):
                                LW = w + 1
                                IW=walls[LW].kwtype

                                # alpha はαAのこと。αは密閉空気層の場合もありうる。
                                # まず公開するのは複雑ではない計算。
                                # まずは単純な一次元計算にする。
                                IF( L5 == 2 and wall[K].alpha > 0.0)THEN
                                    # 温度×高さの積算
                                    D1=D1+TMP(LW,L1)*walls[LW].height
                                END IF
                            END DO

                            IF(HIGHT_CAVITY(K).GT.0.1)THEN
                                # 加重平均温度を計算
                                TEMP_CAVITY(K) = D1 / HIGHT_CAVITY(K)
                                # 353.25: kg K/m3
                                # P V = n R T   P: N/m2, V: m3, n: kg, R: J/(kg K), T: K
                                # n R / V = P / T * 何か定数
                                # 353.25 は、大気圧(pa)101325わる乾燥空気の気体常数(J/(kg K) 287.05で計算した値。
                                # kg/m3 が求まる。
                                GMAQ1=353.25/(TEMP_CAVITY(K))
                                GMAO=353.25/(oc.t_k)
                                PQ1=GMAO-GMAQ1
                                # 0.5 Δρ g h                       
                                D1=0.5*ABS(PQ1)*HIGHT_CAVITY(K)*9.8                !中性帯　高さ中央 0.5, 重力加速度 g(ρ-ρ)
                            END IF

                            for i, wall in enumerate(walls):
                                LW = i + 1
                                IW = wall.kwtype
                                IF(L5.EQ.2.AND.walls[IW][K].alpha.GT.0)
                                    QQ(LW,K)=DSQRT(D1)*4.*ALP_TOTAL(K)
                            END DO
                        END DO
                    END DO

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
                            I=NRAINPOINT(K,1)
                            J=NRAINPOINT(K,2)
                            IW=walls[I].kwtype
                            L1=walltypes[IW][NRAINPOINT(k,3)].num
                            L5=walltypes[IW][L1].num

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
                                D4=3.43E-08*(D3-VP)*0.3                                       !****濡れ面率0.3  水膜からの蒸発量

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

                    for i, wall in enumerate(walls):
                        
                        LW = i + 1
                        
                        IW = wall.kwtype

                        for i in range(wall.n_mesh_total):
                            
                            # (kg/s)/(J/kg)
                            SAHEN = 0.0
                            # kg/s
                            UHEN = 0.0

                            if wall.get_layer(i).num == 2:
                                m_cap = 0.0
                            else:
                                m_cap = ROW * wall.DPDU(i=i, wpt=WPU(LW,i), tp=TMP(LW,i)) * wall.v_is[i] / dt

                            SAHEN += wall.m_cap(i=i, wp_is=WPU(LW,i))
                            UHEN += wall.m_cap(i=i, wp_is=WPU(LW,i)) * HWPU(LW,I)

                            t_i_mns = oc.t_k if wall.is_outside_surface(i) else TMP(LW,I-1)
                            t_i = TMP(LW,I)
                            t_i_pls = room_t_n(n=n) if wall.is_inside_surface(i) else TMP(LW,I+1)

                            # 水分化学ポテンシャルによる蒸気と液水の移動量に関するコンダクタンス, kg/s / (J / kg)
                            if wall.is_outside_end_point_is[i]:
                                # (kg / s) / (J / kg) = kg/(m2 s Pa) * Pa / (J / kg) * m2
                                D21 = wall.layers[k].cond_m_o * wall.dgdu(i) * wall.area
                            else:
                                D21 = wall.DWX(i=i)
                            if wall.is_inside_end_point_is[i]:
                                D22 = wall.layers[k].cond_m_i * wall.dgdu(i) * wall.area
                            else:
                                D22 = wall.DWX(i=i+1)
                            SAHEN += D21 + D22

                            D31 = D21 * (oc.wp if wall.is_outside_surface(i) else WPU(i - 1))
                            D32 = D22 * (room.wp_n(n=n) if wall.is_inside_surface(i) else WPU(i+1))
                            UHEN += D31 + D32

                            # 温度差による移動コンダクタンス (kg/s / K) に温度差をかけた値,  kg/s
                            D41 = (ALDT(LW,K,1) if wall.is_outside_end_point_is[i] else DTX(LW,I)) * (t_i_mns - t_i)
                            D42 = (ALDT(LW,K,2) if wall.is_inside_end_point_is[i] else DTX(LW,I+1)) * (t_i_pls - t_i)
                            UHEN += D41 + D42

                            if wall.get_layer(i).num == 2:

                                D4=0.
                                D5=0.

                                IF(LW.EQ.1)THEN
                                    D2 = DGDUQ(LW,K) * QQ(LW,K) * oc.wp                              !最下層
                                    D3 = DGDTQ(LW,K) * QQ(LW,K) * ( oc.t_k -TMP(LW,I))
                                ELSE
                                    # DGDUQ：水分化学ポテンシャルによる蒸気と液水の移動量に関するコンダクタンス, kg/s / (J / kg)
                                    # DGDTQ：温度差による移動コンダクタンス (kg/s / K) に温度差をかけた値,  kg/s
                                    D2 = DGDUQ(LW,K) * QQ(LW,K) * WPU(LW-1,I)                        !2階以上
                                    D3 = DGDTQ(LW,K) * QQ(LW,K) * (TMP(LW-1,I)-TMP(LW,I))
                                END IF

                                IF(RN(LW,I+1).GT.0.)THEN                                       ! D4 内側の水膜蒸発量 濡れ面率0.3
                                    t_srf=TMP(LW,I+1)
                                    FS, VP = GOFF(t_srf)
                                    # 3.43E-08: 湿気伝達率　kg/(m2 s Pa)
                                    D4=3.43E-08*(VP-XM(LW,i))*walls[LW].area*0.3
                                END IF

                                IF(RN(LW,I-1).GT.0.)THEN                                       ! D5 外側の水膜蒸発量 濡れ面率0.3
                                    t_srf=TMP(LW,I-1)
                                    FS, VP = GOFF(t_srf)
                                    D5=3.43E-08*(VP-XM(LW,i))*walls[LW].area*0.3
                                END IF

                                UHEN += D2 + D3 + D4 + D5
                                SAHEN += DGDUQ(LW,K)*QQ(LW,K)
                                WPU(LW,I)=UHEN/SAHEN
                                IF(WPU(LW,I).GE.0.00)WPU(LW,I)=-1.3E-4
                                CYCLE
                            END IF

                            # WJRAIN 雨水由来の浸入量
                            # WJW：木材が分解した場合にセルロースが分解された場合に発生する水分量
                            D5 = WJRAIN(LW,I) + WJW(LW,I) * wall.dx_is[i] * wall.area # WJRAIN(LW,I)!

                            UHEN += D5

                            WPU(LW,I)=UHEN/SAHEN

                            IF(WPU(LW,I).GE.0.00)WPU(LW,I)=-1.3E-3      !-100.

                        END DO   !K
                    END DO   !LW

                    # /////////OVER RELAXATION ////// 
                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                c_liquid_i_pls=(WPU(LW,I)-AWPU(LW,I))*OMG+AWPU(LW,I)
                                D8=c_liquid_i_pls-AWPU(LW,I)
                                IF(ABS(D8).GT.AMAX2)AMAX2=ABS(D8)
                                AWPU(LW,I)=c_liquid_i_pls
                            END DO
                        END DO
                    END DO

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
                        XNDIS(LW,I,MON,DAY)=XNAVD(LW,I)/24./NDVD
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
                            WJW(LW,I)=HCOFF*DLOSS * wall.GMA(I) /86400.   ! Time unit:h=24,s=86400
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
PRINT *,' YOMITOBASI=',L3				
DO 70 I = 1,L3+1				
    READ(7,*)MOJI
70 CONTINUE

L4=L2*24

DO I=1,L4
    READ(13,*)D1,D2
END DO

PRINT *, '  MOJI ',MOJI,D1,D2				
RETURN				
END				
!********************************				


def ROOM(kday, trav, trdt, rrav):
    """
    室内の温湿度計算
    kday: 通算日(または経過時間)
    trav: 平均気温
    trdt: 気温の振幅
    rrav: 平均相対湿度
    """
    # D1 = 2 * PI * (KDAY - 212) * 24 / 8760
    # 212日はおそらく8月1日付近を基準（位相）にしている
    d1 = 2 * math.pi * (kday - 212) * 24.0 / 8760.0
    
    # 室温の計算 (余弦曲線での近似)
    tr = trdt * math.cos(d1) + trav
    rr = rrav
    
    t0 = tr + 273.15  # 摂氏 -> ケルビン
    
    # GOFF関数から飽和水蒸気圧等を取得
    fs0, vp = GOFF(t0)
    
    # 水蒸気圧の計算 (VP * 相対湿度[%] * 0.01)
    xr = vp * rr * 0.01
    
    return tr, rr, xr


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


