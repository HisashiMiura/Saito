import numpy as np
from dataclasses import dataclass

from modules.room import Room, InputRoom
from modules.thermo_dynamics import ATP, get_wp, RG, GOFF, FUNCX, get_dgdu, ROW_CP, WPTRE, CPL, ROW, get_rho, get_x
from modules.config import NREPT, HIDO, HKEIDO, NDVD, I_HCOFF, OMG, EPS1, EPS2
from modules.wall import Wall
from modules.solar_position import get_solar_position
from modules.weather import Weather
from modules.date_operation import Period
from modules.nrain import NRAINPOINT, ACOFRAIN, set_default_nrains
from scipy.stats import rayleigh
from modules.date_operation import get_step_d_t, get_step_n
from modules.input_data import InputRoom
from modules.wood_decay import rh_threshold, tin
from modules import wood_decay


INTEGER,PARAMETER :: NXP=150,KMTLP=10,NWP=50,NRM=50


INTEGER,PARAMETER :: MXIT=5000

INTEGER DAY,DAY2,DAY1


!//////////COEFFOCIENT//////////////
DIMENSION NX(NWP)
DIMENSION WGT(NWP,NXP)
DIMENSION KMTL(KMTLP)

!/////////VARIABLE///////////////
DIMENSION WPU(NWP,NXP),WPS(NWP,NXP),TMP(NWP,NXP)

DIMENSION ATMP(NWP,NXP),BWPU(NWP,NXP),BTMP(NWP,NXP)
DIMENSION MCW(NWP)
DIMENSION QQ(NWP,NXP),XM(NWP,NXP)
DIMENSION ATMPQ(NWP,NXP)
DIMENSION WLOSS(NWP,NXP),WJW(NWP,NXP)
DIMENSION WJRAIN(NWP,NXP),RAIN(24)
DIMENSION RN(NWP,NXP),HRN(NWP,NXP)

CHARACTER MOJI*110, FILENAME0*12



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

        WPU(L,I) = WP1
        WPS(L,I) = WP0
        # TMP(L,I) = TP1
        ATMPQ(L,I) = TP1 - ATP
        XN(L,I) = VP * RH0 * 0.01
        XM(L,I) = VP * RH0 * 0.01
        ATMP(L,I) = TP1
        WD = wall.material_is[i].get_psi(rh=RH0)
        WGT(L,I) = WD
        WLOSS(L,I) = 0.

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

                    # 外界気象条件
                    oc = weather.get_condition(month=MON, day=DAY, hour=IM, frac=frac)


                    # ////////PUT IN PREBIOUS VALUE///////////

                    for LW, wall in enumerate(walls):

                        for i in range(wall.nrains):

                                wall.t_n_is = wall.t_n_pls

                                # 1ステップ前の値を入れ替えている。
                                XM(LW,i) = XN(LW,i)

                    for w, wall in enumerate(walls):
                        
                        wall.rn = rn_ws_is[w]

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

                    # 温度に関する収束計算が終了した時点での温度, [W, I] (ジャグ配列)
                    t_ws_is = []
                    
                    n = get_step_n(month=month, day=day, hour=hour, n_hour=n, n_div=NDVD)


                    # *************換気量の算出(Q=m3/s)****************
                    rho_o = oc.rho
                    
                    for w, wall in enumerate(walls):
                        
                        # 通気層の換気量, m3/s, [I]
                        v_air = wall.get_v_air(oc=oc)

                        t_is = wall.t_n_is

                        # 収束判定を全ての壁について一気にやっていたところを、壁１枚１枚で判定するように変更した。
                        # ある壁が他の壁に与える影響が無い（あるいは少ない）場合はこの方法の方が良いと思われる。
                        # 仮にある壁が他の壁に与える影響がある場合（例えば、ある壁の壁体内通気層温度が連通している他の壁の通気層の上流側温度になる場合など）であっても、
                        # その影響度合いが小さいのであれば、陽解法的にステップnの値を用いることによって簡易化する方が、計算速度の観点から良いと思われる。
                        for i in range(500):

                            t_is_next = wall.get_t_n_pls(t_is=t_is, dt=dt, oc=oc, theta_r_n=room.theta_n(n=n), wp_r_n=room.wp_n(n=n), v_air=v_air, t_upstream=oc.t_k)

                            t_is = (t_is_next - t_is) * OMG + t_is

                            delta_is = t_is_next - t_is

                            # ステップn+1とステップnの値の差が許容誤差範囲内に収まったらループを抜ける。
                            if EPS1 >= np.abs(delta_is).max():
                                break

                        # for文が最後までまわりきってしまった（収束しなかった）場合の措置。
                        else:
                            raise Exception('温度計算が収束しませんでした。')
                
                        t_ws_is.append(t_is)


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
                    rn_ws_is = []
                    wjrain_ws_is = []

                    for w, wall in enumerate(walls):

                        rn_is, wjrsin_is = wall.get_rn_n_pls(wp_is=wp_is, oc=oc, p_v_rm=room.p_v_n(n), dt=dt)

                        rn_ws_is.append(rn_is)
                        wjrain_ws_is.append(wjrsin_is)


                    theta_r_n = room.theta_n(n=n)
                    wp_r_n = room.wp_n(n=n)

                    wp_ws_is_next = []

                    delta_ws = []

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
                            WJRAIN=wjrain[w],
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
                                WP0=WPU(LW,I)
                                TP=TMP(LW,I)
                                RH0 = WPTRE(WP0,TP)
                                wall.rh[I] =RH0
                                FS, VP = GOFF(TP)

                                XN(LW,I)=RH0*VP*0.01

                                WD = wall.material_is[i].get_u(rh=RH0)
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

                    K=0
                    ITC=ITC+1
    
                    IF(K.EQ.1)GO TO 710
    
                    # ///////JUDGEMENT CONVERGENCE OF HEAT AND MOISTURE/////

                    # 収束計算ありの場合　NREPT=1, 収束計算なしの場合　NREPT=0
                    # NREPT=1のときは温度と水分の収束計算をする。
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

                        if IT3 > MXIT:
                            print('NO CONVERGENCE')
                            GO TO 556
    

                        if AMAX3 > EPS2 or AMAX4 > EPS1:
                            GO TO 810

                    END IF  !/////COMBINE

                    # 556 CONTINUE

                # 666 CONTINUE

            # 777 CONTINUE

            for wall in walls:

                for i in range(wall.n_mesh_total):

                    time_s = wood_decay.update_time_s(theta=wall.theta_n_is(i), rh=wall.rh[i], time_s=wall.time_s_is[i])

                    l_stage = wood_decay.update_stage(theta=wall.theta_n_is(i), rh=wall.rh[i], time_s=wall.time_s_is[i], l_stage=wall.l_stage_is[i])

                    wall.time_s_is[i] = time_s
                    wall.l_stage_is[i] = l_stage
                
                # 質量減少分, kg/(kg d)
                mass_loss_is = wood_decay.get_mass_loss_is(l_stage_is=wall.l_stage_is, theta_is=wall.theta_n_is, rh_is=wall.rh)

                wall.m_loss += mass_loss_is

                WLOSS_MAX = 0.6

                # I_HCOFF：水分生成を計算するかどうかの判定パラメータ
                if I_HCOFF == 1:
                    HCOFF=0.319 # 水分生成量（実験値）
                else:
                    HCOFF=0.

                for i in range(wall.n_mesh_total):

                    if wall.m_loss > WLOSS_MAX:
                        wall.wjw[i] = 0.0
                    
                    else:
                        # HCOFF: 不朽した分は菌の代謝に使われ、残りは水分になる。質量現象に対する水分量, 無次元
                        # gma: 密度, kg/m3
                        # 86400: s/d
                        # WJW：木材が分解した場合にセルロースが分解された場合に発生する水分量, kg/(m3 s)
                        # kg（水）/(m3（木材） s) = HCOFF（m3（水）/m3（材料））/(kg（材料）/kg（材料）) * DLOSS(kg（木材）/kg（木材）/d) * gma（kg（木材）/m3（木材））
                        # TODO: wall.gma_is は水の密度が正しい                            
                        wall.wjw[i] = HCOFF * mass_loss_is[i] * wall.gma_is[i] / 86400

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


def SATUWPT(TP):

    # J/(kg K)
    CPW = ( 30.36 + 0.009615 * TP + 0.00000118 * TP * TP ) / ( 0.018016 )

    FS, VP = GOFF(TP)

    # J/ kg
    SW=(644243)+CPW*(TP-273.15)-TP*CPW*DLOG(TP/273.15)+461.5*TP*DLOG(VP/101325)

    return SW

