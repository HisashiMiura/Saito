from dataclasses import dataclass
import numpy as np


from modules.input_wall import InputWall
from modules.direction import Direction
from modules.weather import OutdoorCondition
from modules.config import HOI
from modules.surface_solar import get_surface_solar_d_t


@dataclass
class WallSurface:

    # 方位
    direction: Direction

    # 傾斜角
    angle: float

    # アルベド
    albedo: float

    # 評価高さ, m
    eva_heigt: float

    @classmethod
    def read(cls, iw: InputWall):

        # 方位
        direction = iw.direction

        # 傾斜角
        angle = iw.angle

        # アルベド
        albedo = iw.albedo

        # 評価高さ, m
        eva_height = iw.eva_height

        return WallSurface(
            direction=direction,
            angle=angle,
            albedo=albedo,
            eva_heigt=eva_height
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

    def get_v_wind_eva_k(self, v_wind: float) -> float:
        """評価高さにおける風速を求める。

        Args:
            v_wind: 風速, m/s

        Returns:
            評価高さにおける風速, m/s
        """

        # 基準風速は6.5m高さ
        return v_wind * (self.eva_height / 6.5)**0.25
    
    def get_wind_angle(self, wind_direction: float) -> float:
        """風向と壁の法線のなす角度を求める。

        Args:
            wind_direction: 風向き, 度

        Returns:
            風向と壁の法線のなす角度, 度
            
        Notes:
            TODO: 壁が傾斜している場合も考慮するべきではないか。現行の方法だと垂直壁しか考えていないように思われる。
            TODO: 風向きの定義が、北側から時計回りで定義されているのでは？
        """

        return wind_direction - (180.0 + self.direction.alpha)


