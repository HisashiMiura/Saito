from dataclasses import dataclass
import numpy as np


from .input_data import InputWall
from .direction import Direction
from .weather import OutdoorCondition
from .config import HOI
from .surface_solar import get_surface_solar_d_t


@dataclass
class WallSurface:

    # 方位
    direction: Direction

    # 傾斜角
    angle: float

    # アルベド
    albedo: float

    @classmethod
    def read(cls, iw: InputWall):

        # 方位
        direction = iw.direction

        # 傾斜角
        angle = iw.angle

        # アルベド
        albedo = iw.albedo

        return WallSurface(
            direction=direction,
            angle=angle,
            albedo=albedo
        )

    def get_q_sol_d_t_k(self, oc: OutdoorCondition) -> float:
        """壁面の外側に入る日射・放射由来の熱を計算する。

        Args:
            oc: 外気のコンディションを保持するクラス

        Returns:
            float: 壁面の外側に入る日射・放射由来の熱量, W/m2
        """

        drct = self.direction

        # 壁の向きが下向きの場合は入射する日射量は計算するまでもなく0とする。
        if drct == Direction.BOTTOM:

            return 0.0
        
        else:

            # 壁が上向きの場合は方位角を定義できない。
            if drct == Direction.TOP:
                alpha = None
            else:
                alpha = np.radians(drct.alpha + HOI)

            beta = np.radians(self.angle)

            # 壁面にあたる短波長日射量と夜間放射量を計算する。
            # 短波長日射量, W/m2
            # 夜間放射量, W/m2
            i_d_t_k, r_d_t_k = get_surface_solar_d_t(
                sin_h_d_t=oc.sin_h,
                cos_h_d_t=oc.cos_h,
                sin_a_d_t=oc.sin_a,
                cos_a_d_t=oc.cos_a,
                alpha_k=alpha,
                beta_k=beta,
                i_dn_d_t=oc.i_dn,
                i_sky_d_t=oc.i_sky,
                rho_g=self.albedo,
                r_n_d_t=oc.r_n
            )

            q_sol_d_t_k = self.absorption * i_d_t_k - r_d_t_k * self.emissivity

            return q_sol_d_t_k

    def get_confrains