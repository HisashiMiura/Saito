import numpy as np

from modules.room import Room, InputRoom
from modules.thermo_dynamics import GOFF
from modules.config import HIDO, HKEIDO, NDVD, I_HCOFF, OMG, EPS1, EPS2
from modules.wall import Wall
from modules.weather import Weather
from modules.date_operation import Period
from modules.nrain import set_default_nrains
from modules.date_operation import get_step_n
from modules.input_data import InputRoom
from modules import wood_decay


ipt_room = InputRoom.read(d={})
room = Room.init(ipt_room=ipt_room, n_div=NDVD)


nrains = set_default_nrains()


walls = Wall.read_default()

dt = 1. / NDVD * 3600.


period = Period()

weather = Weather.load_data(file_name='weather.csv', longitude=HIDO, latitude=HKEIDO)


t_ws_is = [np.zeros(wall.n_mesh_total, dtype=float) for wall in walls]
wp_ws_is = [np.zeros(wall.n_mesh_total, dtype=float) for wall in walls]

for NYEAR in range(1, period.LYEAR + 1):

    KDAY=0				
    
    MONF = period.get_start_month(year=NYEAR)
    LMONL = period.get_last_month(year=NYEAR)

    IMQ=0

    for MON in range(MONF, LMONL + 1):
    
        DAY1 = period.get_start_day(year=NYEAR)
        DAY2 = period.get_end_day(year=NYEAR, month=MON)

        # 三浦コメント：ここから１日のループ計算
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

                    for w, wall in enumerate(walls):
                        
                        wall.rn = rn_ws_is[w]

                    n = get_step_n(month=month, day=day, hour=hour, n_hour=n, n_div=NDVD)

                    # *************換気量の算出(Q=m3/s)****************
                    rho_o = oc.rho

                    v_air_ws_is = []
                    
                    # 温度に関する収束計算が終了した時点での温度, [W, I] (ジャグ配列)
                    t_ws_is = []

                    for w, wall in enumerate(walls):
                        
                        # 通気層の換気量, m3/s, [I]
                        v_air = wall.get_v_air(oc=oc)

                        v_air_ws_is.append(v_air)

                        t_is = wall.t_n_is

                        # 収束判定を全ての壁について一気にやっていたところを、壁１枚１枚で判定するように変更した。
                        # ある壁が他の壁に与える影響が無い（あるいは少ない）場合はこの方法の方が良いと思われる。
                        # 仮にある壁が他の壁に与える影響がある場合（例えば、ある壁の壁体内通気層温度が連通している他の壁の通気層の上流側温度になる場合など）であっても、
                        # その影響度合いが小さいのであれば、陽解法的にステップnの値を用いることによって簡易化する方が、計算速度の観点から良いと思われる。
                        for i in range(500):

                            t_is_next = wall.get_t_n_pls(t_is=t_is, dt=dt, oc=oc, theta_r_n=room.theta_n(n=n), wp_r_n=room.wp_n(n=n), v_air=v_air, t_upstream=oc.t_k)

                            delta_is = t_is_next - t_is

                            t_is = (t_is_next - t_is) * OMG + t_is

                            # ステップn+1とステップnの値の差が許容誤差範囲内に収まったらループを抜ける。
                            if EPS1 >= np.abs(delta_is).max():
                                break

                        # for文が最後までまわりきってしまった（収束しなかった）場合の措置。
                        else:
                            raise Exception('温度計算が収束しませんでした。')
                
                        t_ws_is.append(t_is)

                    # ////////////CAL_WATER POTENTIAL BY OVER RELAXATION METHOD/////////
                    IT2=0
                    CWMAX=2.8E-2

                    theta_r_n = room.theta_n(n=n)
                    wp_r_n = room.wp_n(n=n)

                    rn_ws_is = []
                    wjrain_ws_is = []

                    theta_r_n = room.theta_n(n=n)
                    wp_r_n = room.wp_n(n=n)

                    for w, wall in enumerate(walls):

                        rn_is, wjrsin_is = wall.get_rn_n_pls(wp_is=wp_is, oc=oc, p_v_rm=room.p_v_n(n), dt=dt)

                        rn_ws_is.append(rn_is)
                        wjrain_ws_is.append(wjrsin_is)

                    wp_ws_is = []

                    for w, wall in enumerate(walls):

                        v_air_is = v_air_ws_is[w]
                        wjrain_is = wjrain_ws_is[w]

                        wp_is = wp_ws_is[w]

                        for i in range(500):
                        
                            wp_is_next = wall.get_wp_n_pls(dt=dt, oc=oc, theta_r_n=theta_r_n, wp_r_n=wp_r_n, wp_is=wp_is, v_air_is=v_air_is, RN=rn_ws_is[w], WJRAIN=wjrain_is)

                            delta_is = wp_is_next - wp_is

                            wp_is = (wp_is_next - wp_is) * OMG + wp_is

                            if EPS2 >= np.abs(delta_is).max():
                                break
                        
                        else:
                            raise Exception('湿度計算が収束しませんでした。')
                        
                        wp_ws_is.append(wp_is)

                    # 一旦、熱と水分の同時収束計算の評価式は削除した。


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


def SATUWPT(TP):

    # J/(kg K)
    CPW = ( 30.36 + 0.009615 * TP + 0.00000118 * TP * TP ) / ( 0.018016 )

    FS, VP = GOFF(TP)

    # J/ kg
    SW=(644243)+CPW*(TP-273.15)-TP*CPW*DLOG(TP/273.15)+461.5*TP*DLOG(VP/101325)

    return SW

