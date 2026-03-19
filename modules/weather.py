import numpy as np
from dataclasses import dataclass
import pandas as pd
import os

from config import ATP
from .thermo_dynamics import GOFF, get_wp
import solar_position
from date_operation import get_step_d, get_step_d_t


@dataclass
class OutdoorCondition:
    """瞬時の外気コンディションを保持するためのデータクラス"""

    # 太陽高度の正弦
    sin_h: float
    
    # 太陽高度の余弦
    cos_h: float
    
    # 太陽方位角の正弦
    sin_a: float
    
    # 太陽方位角の余弦
    cos_a: float
    
    # 温度, deg.C
    t: float

    # 水蒸気圧, Pa
    xod: float

    # 直達日射
    i_dn: float

    # 天空日射
    i_sky: float

    # 夜間放射
    r_n: float

    # 降雨量, mm/h
    rainfall: float

    # 風速, m/s
    v_wind: float

    # 風向, deg.
    wind_direction: float

    @property
    def t_k(self):
        """絶対温度, K"""
        return self.t + ATP
    
    @property
    def wp(self):
        """水分化学ポテンシャル, J/kg"""
        return get_wp(rh=self.rh, t=self.t_k)
    
    @property
    def vp_sat(self):
        """飽和水蒸気圧, Pa"""
        return GOFF(t=self.t_k)[1]
    
    @property
    def rh(self):
        """相対湿度, %"""
        return self.xod / self.vp_sat * 100.0


@dataclass
class Weather:

    # 経度, 度
    longitude: float

    # 緯度, 度
    latitude: float

    # 気温, deg.C, [8760]
    t_ns: np.array

    # 絶対湿度, g/kg(DA), [8760]
    x_ns: np.array

    # 南北風（南風：正　北風：負）, [8760]
    v_sn_ns: np.array

    # 東西風（西風：正　東風：負）, [8760]
    v_we_ns: np.array

    # 大気圧, Pa, [8760]
    p_atm_ns: np.array

    # (1時間前から当該時刻までの)降水量の積算値, mm/h, [8760]
    rf_ns: np.array

    # 風速（瞬時値）, m/s, [8760]
    v_wind_ns: np.array

    # (1時間前から当該時刻までの)法線面直達日射量の平均値, W/m2, [8760]
    i_dn_ns: np.array

    # (1時間前から当該時刻までの)水平面天空日射量の平均値, W/m2, [8760]
    i_sky_ns: np.array

    # (1時間前から当該時刻までの)夜間放射量の平均値, W/m2, [8760]
    r_n_ns: np.array

    @property
    def p_ns(self) -> np.ndarray:
        """水蒸気圧(Pa)を求める。"""

        # 水蒸気圧, Pa, [8760]
        # TODO: 外部モジュールにまとめて書くのが望ましい。
        p_ns = self.p_atm_ns * self.x_ns / (622.0 + self.x_ns)

        return p_ns
    
    @property
    def d_wind_ns(self) -> np.ndarray:
        """風向を求める。"""

        # 風向, °, [8760]
        # おそらく南から時計回りの度。
        # TODO: 仮置き。このプログラムがどこ基準、どちら周りで定義しているのかを調べて守成すること。
        d_wind_ns = np.rad2deg(np.atan2(self.v_sn_ns, self.v_we_ns))

        return d_wind_ns

    def get_at_step(self, month: int, day: int, hour: int, frac: float = 0.0):

        return self.get_condition(month=month, day=day, hour=hour, frac=frac)

    def get_condition(self, month: int, day: int, hour: int, frac: float):

        t = float(hour) + frac

        d = get_step_d(month=month, day=day)

        phi_lon = np.radians(self.longitude)

        phi_lat = np.radians(self.latitude)

        sin_h_d_t, cos_h_d_t, sin_a_d_t, cos_a_d_t = solar_position.get_solar_position(
            d=d, t=t, phi_lon=phi_lon, phi_lat=phi_lat
        )

        n = get_step_d_t(month=month, day=day, hour=hour)

        n_next = 0 if n == 8760 else n + 1

        def f(v: np.ndarray) -> float:

            return v[n] + frac * (v[n_next] - v[n])

        t = f(self.t_ns)

        xod = f(self.p_ns)

        i_dn = f(self.i_dn_ns)

        i_sky = f(self.i_sky_ns)

        r_n = f(self.r_n_ns)

        rf = f(self.rf_ns)

        v_wind = f(self.v_wind_ns)

        d_wind = f(self.d_wind_ns)

        return OutdoorCondition(
            sin_h=sin_h_d_t,
            cos_h=cos_h_d_t,
            sin_a=sin_a_d_t,
            cos_a=cos_a_d_t,
            t=t,
            xod=xod,
            i_dn=i_dn,
            i_sky=i_sky,
            r_n=r_n,
            rainfall=rf,
            v_wind=v_wind,
            wind_direction=d_wind
        )

    @classmethod
    def load_data(cls, file_name: str, longitude: float, latitude: float):
        """_summary_

        Args:
            file_name: 気象データファイル名
            longitude: 緯度, 度
            latitude: 経度, 度
        """

        file_path = os.path.join(os.path.dirname(__file__), file_name)

        d = pd.read_csv(file_path, encoding='shift_jis')
        
        # 気温, deg.C, [8760]
        t_ns = np.array(d['気温'])

        # 絶対湿度, g/kg(DA), [8760]
        x_ns = np.array(d['絶対湿度'])

        # 南北風（南風：正　北風：負）, [8760]
        v_sn_ns = np.array(d['南北風'])

        # 東西風（西風：正　東風：負）, [8760]
        v_we_ns = np.array(d['東西風'])

        # 大気圧, Pa, [8760]
        p_atm_ns = np.array(d['気圧'])

        # (1時間前から当該時刻までの)降水量の積算値, mm/h, [8760]
        rf_ns = np.array(d['降水量'])

        # 風速, m/s, [8760]
        v_wind_ns = np.array(d['風速'])

        # 法線面直達日射量, MJ/m2, [8760]
        i_dn_ns = np.array(d['法線面直達日射量'])

        # 水平面天空日射量, MJ/m2, [8760]
        i_sky_ns = np.array(d['水平面天空日射量'])

        # 夜間放射量, MU/m2, [8760]
        r_n_ns = np.array(d['夜間放射'])

        return Weather(
            longitude=longitude,
            latitude=latitude,
            t_ns=t_ns,
            x_ns=x_ns,
            v_sn_ns=v_sn_ns,
            v_we_ns=v_we_ns,
            p_atm_ns=p_atm_ns,
            rf_ns=rf_ns,
            v_wind_ns=v_wind_ns,
            i_dn_ns=convert_MJ_to_W(i_dn_ns),
            i_sky_ns=convert_MJ_to_W(i_sky_ns),
            r_n_ns=convert_MJ_to_W(r_n_ns)
        )


def convert_MJ_to_W(v):
    """MJからWへの変換"""

    return v * 1000000 / 3600

