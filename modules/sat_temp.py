import numpy as np


def get_theta_sat_d_t(
        sin_h_d_t: float,
        cos_h_d_t: float,
        sin_a_d_t: float,
        cos_a_d_t: float,
        alpha_k: float,
        beta_k: float,
        i_dn_d_t: float,
        i_sky_d_t: float,
        rho_g: float,
        theta_o_d_t: float,
        alpha_s_k: float,
        epsilon_k: float,
        r_n_d_t: float,
        h_o_k: float
    ) -> float:
    """相当外気温度を計算する。

    Args:
        sin_h_d_t: 日付dの時刻tにおける太陽高度の正弦, -
        cos_h_d_t: 日付dの時刻tにおける太陽高度の余弦, -
        sin_a_d_t: 日付dの時刻tにおける太陽方位角の正弦, -
        cos_a_d_t: 日付dの時刻tにおける太陽方位角の余弦, -
        alpha_k: 部位kの方位角, rad
        beta_k: 部位kの傾斜角, rad
        i_dn_d_t: 日付dの時刻tにおける法線面直達日射量, W/m2K
        i_sky_d_t: 日付dの時刻tにおける水平面天空日射量, W/m2K
        rho_g: 地面反射率, -
        theta_o_d_t: 日付dの時刻tにおける外気温度, ℃
        alpha_s_k: 部位kの日射吸収率, -
        beta_k: 部位kの長波長放射率, -
        r_n_d_t: 日付dの時刻tにおける夜間放射量, W/m2K
        h_o_k: 部位kの室外側総合熱伝達率, W/m2K
    
    Returns:
        日付dの時刻tにおける部位kの相当外気温度, ℃
    """


    # 日付dの時刻tにおける部位kの入射角の余弦, -, eq.1
    cos_i_d_t_k = np.cos(beta_k) * sin_h_d_t + np.sin(beta_k) * cos_h_d_t * (cos_a_d_t * np.cos(alpha_k) + sin_a_d_t * np.sin(alpha_k))

    # 日付dの時刻tにおける部位kに入射する日射量の直達成分, W/m2K, eq.2
    i_d_d_t_k = max(i_dn_d_t * cos_i_d_t_k, 0.0)

    # 部位kの天空に対する形態係数, -, eq.3
    f_sky_k = (1 + np.cos(beta_k)) / 2

    # 日付dの時刻tにおける部位kに入射する日射量の天空成分, W/m2K, eq.4
    i_s_d_t_k = i_sky_d_t * f_sky_k

    # 日付dの時刻tにおける水平面全天日射量, W/m2, eq.5
    i_t_d_t = i_dn_d_t * sin_h_d_t + i_sky_d_t

    # 日付dの時刻tにおける部位kに入射する日射量の地面反射成分, W/m2, eq.6
    i_r_d_t_k = i_t_d_t * rho_g * (1 - f_sky_k)

    # 日付dの時刻tにおける部位kに入射する日射量, W/m2, eq.7
    i_d_t_k = i_d_d_t_k + i_s_d_t_k + i_r_d_t_k

    # 日付dの時刻tにおける部位kの相当外気温度, deg.C, eq.8
    theta_sat_d_t_k = theta_o_d_t + (alpha_s_k * i_d_t_k - f_sky_k * epsilon_k * r_n_d_t) / h_o_k

    net_gain = alpha_s_k * i_d_t_k - f_sky_k * epsilon_k * r_n_d_t

    d_gain = alpha_s_k * i_d_d_t_k

    sky_gain = alpha_s_k * (i_s_d_t_k + i_r_d_t_k)

    return theta_sat_d_t_k, net_gain, d_gain, sky_gain