from dataclasses import dataclass
import numpy as np
import math

from .config import HTC_TW
from .sat_temp import get_theta_sat_d_t
from .weather import OutdoorCondition
from .config import HOI

@dataclass
class ThinWall:

    # 方位
    direction: int

    # 隣接室
    room_sideB: int

    # 面積, m2
    area: float

    # 熱貫流率, W/(m2 K)
    u_value: float

    # 遮蔽係数
    shading_factor: float

    @property
    def cond(self):
        """室外から室内へのコンダクタンス, W/K"""
        return self.u_value * self.area
    
    def get_surf_temp(self, t_rm: float, t_out: float) -> float:
        """室内側の表面温度を計算する。

        Args:
            t_rm: 室内温度, K
            t_out: 室外温度, K

        Returns:
            室内側表面温度, K
        """

        return t_rm - self.u_value / HTC_TW * (t_rm - t_out)
    
    def get_alpha_w(self):
    
        e = (self.direction - 1) * 90.0
        return np.radians(e + HOI)

    
    def get_k_sat(self, oc: OutdoorCondition):

        angle = 90.0
        albedo = 0.3
        absorption = 0.9
        emissivity = 0.9

        sin_h_d_t = oc.sin_h
        cos_h_d_t = oc.cos_h
        sin_a_d_t = oc.sin_a
        cos_a_d_t = oc.cos_a
        i_dn_d_t = oc.i_dn
        i_sky_d_t = oc.i_sky
        theta_o_d_t = oc.t
        r_n_d_t = oc.r_n

        alpha_k = self.get_alpha_w
        beta_k = np.radians(angle)

        theta_sat_d_t_k, net_gain, d_gain, sky_gain = get_theta_sat_d_t(
            sin_h_d_t=sin_h_d_t,
            cos_h_d_t=cos_h_d_t,
            sin_a_d_t=sin_a_d_t,
            cos_a_d_t=cos_a_d_t,
            alpha_k=alpha_k,
            beta_k=beta_k,
            i_dn_d_t=i_dn_d_t,
            i_sky_d_t=i_sky_d_t,
            rho_g=albedo,
            theta_o_d_t=theta_o_d_t,
            alpha_s_k=absorption,
            epsilon_k=emissivity,
            r_n_d_t=r_n_d_t,
            h_o_k=22.4
        )

        return theta_sat_d_t_k, net_gain, d_gain, sky_gain

    @classmethod
    def read(cls, d:dict):

        return cls(
            direction=d['direction'],
            room_sideB=d['room_sideB'],
            area=d['area'],
            u_value=d['u_value'],
            shading_factor=d['shading_factor']
        )

    @classmethod
    def read_default(cls):

        # 南壁, 西壁, 北壁, 東壁
        ds = [
            {
                'direction': 1,
                'room_sideB': 4,
                'area': 0.09,
                'u_value': 6.5,
                'shading_factor': 0.01
            },
            {
                'direction': 2,
                'room_sideB': 4,
                'area': 0.09,
                'u_value': 6.5,
                'shading_factor': 0.01
            },
            {
                'direction': 3,
                'room_sideB': 4,
                'area': 0.09,
                'u_value': 6.5,
                'shading_factor': 0.01
            },
            {
                'direction': 4,
                'room_sideB': 4,
                'area': 0.09,
                'u_value': 6.5,
                'shading_factor': 0.01
            },
        ]

        return [cls.read(d=d) for d in ds]
