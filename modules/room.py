from dataclasses import dataclass
import numpy as np
import math

from .thermo_dynamics import ATP, GOFF, get_wp, ROW_CP, WPTRE, FUNCX, get_dgdt, get_dgdu
from .config import TMPIC, RHI
from .date_operation import get_step_total


@dataclass
class InputRoom:

    # 平均温度, ℃
    average_temperature: float

    # 振幅, ℃
    amplitude_temperature: float

    # 平均湿度
    average_relative_humidity: float

    @classmethod
    def read(cls, d: dict):

        try:
            average_temperature = d.get('average temperature', 22.5)
        except:
            raise ValueError('\'average temperature\' に不正な値が指定されました。')
        
        try:
            amplitude_temperature = d.get('amplitude temperature', 4.5)
        except:
            raise ValueError('\'amplitude temperature\' に不正な値が指定されました。')

        try:
            average_relative_humidity = d.get('average relative humidity', 60.0)
        except:
            raise ValueError('\'average relative humidity\' に不正な値が指定されました。')

        InputRoom(
            average_temperature=average_temperature,
            amplitude_temperature=amplitude_temperature,
            average_relative_humidity=average_relative_humidity
        )


@dataclass
class Room:

    # 温度, ℃, [N]
    theta_ns: np.array

    # 相対湿度, %, [N]
    rh_ns: np.array

    @classmethod
    def init(cls, ipt_room: InputRoom, n_div: int):
        """Room クラスを作成

        Args:
            ipt_room: InputRoom クラス
            n_div: 1時間の分割数
        """

        theta_ave = ipt_room.average_temperature
        theta_amp = ipt_room.amplitude_temperature
        rh_ave = ipt_room.average_relative_humidity

        # 総ステップ数 (=N)
        # 例えば1時間の分割数が60(つまり1秒間隔)の場合は、365*24*60
        n_total = get_step_total(n_div=n_div)

        # 0 から n_total までの連番を作成する。
        ns = np.arange(n_total)

        # 位相(-2 pi ~ 2 pi), [N]
        phase_ns = 2 * np.pi * (ns - 212 * 24 * n_div) / n_total 

        # 温度, ℃, [N]
        theta_ns = theta_ave + theta_amp * np.cos(phase_ns)

        # 相対湿度, %, [N]
        rh_ns = np.full_like(a=theta_ns, fill_value=rh_ave, dtype=float)

        return Room(
                theta_ns=theta_ns,
                rh_ns=rh_ns
        )

    def rh_n(self, n: int) -> float:
        """ステップnの相対湿度を求める。

        Args:
            n: ステップ
        
        Returns:
            ステップnの相対湿度, %
        """

        return self.rh_ns[n]
    
    def theta_n(self, n: int) -> float:
        """ステップnの室内の温度を求める。
        
        Args:
            n: ステップ
        
        Returns:
            ステップnの温度, ℃
        """

        return self.theta_ns[n]
    
    def t_n(self, n: int) -> float:
        """ステップnの絶対温度を求める。

        Args:
            n: ステップ

        Returns:
            ステップnの絶対温度, K
        """

        return self.theta_n(n=n) + ATP
    
    def p_vsat_n(self, n: int) -> float:
        """ステップnの飽和水蒸気圧を求める。

        Args:
            n: ステップ

        Returns:
            ステップnの飽和水蒸気圧, Pa
        """

        return GOFF(t=self.t_n(n=n))[1]

    def p_v_n(self, n: int) -> float:
        """ステップnの水蒸気圧を求める。

        Args:
            n: ステップ

        Returns:
            ステップnの水蒸気圧, Pa
        """

        return self.p_vsat_n(n=n) * self.rh_n(n=n) * 0.01
    
    def wp_n(self, n: int) -> float:
        """ステップnの水分化学ポテンシャルを計算する。

        Args:
            n: ステップ

        Returns:
            ステップnの水分化学ポテンシャル, J/K
        """

        return get_wp(rh=self.rh_n(n=n), t=self.t_n(n=n))

    def x_v_n(self, n: int) -> float:
        """ステップnの混合比を求める。

        Args:
            n: ステップ

        Returns:
            ステップnの混合比, kg/kg(DA)
        """

        return FUNCX(rh=self.rh_n(n=n), vp=self.p_vsat_n(n=n))


            

            
