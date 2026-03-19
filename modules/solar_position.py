import numpy as np

def get_solar_position(d: int, t: float, phi_lon: float, phi_lat: float) -> tuple[float, float, float, float]:
    """太陽位置を計算する。

    Args:
        d: 日付
        t: 時刻
        phi_lon: 経度, rad
        phi_lat: 緯度, rad
    
    Returns:
        日付dの時刻tにおける太陽高度の正弦, -
        日付dの時刻tにおける太陽高度の余弦, -
        日付dの時刻tにおける太陽方位角の正弦, -
        日付dの時刻tにおける太陽方位角の余弦, -
    """

    # 日付dにおける日角, rad, eq.1
    omega_d = 2.0 * np.pi * d / 366.0

    # 日付dにおける日赤緯, rad, eq.2
    delta_d = (0.006322 - 0.405748 * np.cos(omega_d + 0.153231)
            - 0.00588 * np.cos(2 * omega_d + 0.207099)
            - 0.003233 * np.cos(3 * omega_d + 0.620129))

    # 日付dにおける均時差, h, eq.3
    e_d = (-0.000279 + 0.122772 * np.cos(omega_d + 1.498311)
            - 0.165458 * np.cos(2 * omega_d - 1.261546)
            - 0.005354 * np.cos(3 * omega_d - 1.1571))

    # 日付dの時刻tにおける時角, rad, eq.4
    t_d_t = np.pi / 12.0 * ((t + e_d - 12.0) + (phi_lon * 180 / np.pi - 135.0) / 15.0)

    # 日付dの時刻tにおける太陽高度の正弦, -, eq.5
    sin_h_d_t = np.sin(phi_lat) * np.sin(delta_d) + np.cos(phi_lat) * np.cos(delta_d) * np.cos(t_d_t)

    # 日付dの時刻tの太陽高度の余弦, eq.6
    cos_h_d_t = np.cos(np.arcsin(sin_h_d_t))

    # 日付dの時刻tの太陽方位角の正弦, -, eq.7
    sin_a_d_t = np.cos(delta_d) * np.sin(t_d_t) / cos_h_d_t
    
    # 日付dの時刻tの入射角の余弦, eq.8
    cos_a_d_t = (
        - np.sin(delta_d) * np.cos(phi_lat) + np.cos(delta_d) * np.sin(phi_lat) * np.cos(t_d_t)
        ) / cos_h_d_t

    return sin_h_d_t, cos_h_d_t, sin_a_d_t, cos_a_d_t