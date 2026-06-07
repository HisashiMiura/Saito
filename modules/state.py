from dataclasses import dataclass
import numpy as np

from modules.thermo_dynamics import ATP
from modules.thermo_dynamics import get_wp, GOFF, get_rh


@dataclass
class State:

    # 絶対温度, K
    t: float

    # 水分化学ポテンシャル, J/kg
    mu: float

    @classmethod
    def init(cls, t: float, rh: float):
        """Stateクラスの初期化を行うクラスメソッド

        Args:
            t: 絶対温度, K
            rh: 相対湿度, %
        
        Returns:
            Stateクラスのインスタンス
        """

        # 水分化学ポテンシャルの初期値を計算
        wp = get_wp(rh=rh, t=t)

        return State(t=t, mu=wp)

    @property
    def theta(self) -> float:
        return self.t - ATP
    
    @theta.setter
    def theta(self, value):
        self.t = value + ATP
    
    @property
    def p_v_sat(self) -> float:
        """飽和水蒸気圧, Pa"""

        return GOFF(self.t)[1]
    
    @property
    def rh(self) -> float:
        """相対湿度, %"""

        return get_rh(mu=self.mu, t=self.t)
    
    @property
    def p_v(self) -> float:
        """水蒸気圧, Pa"""

        return self.rh * self.p_v_sat * 0.01
    
