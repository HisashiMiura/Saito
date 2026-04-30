from dataclasses import dataclass
import numpy as np
import math
from scipy.stats import rayleigh
from numpy.typing import NDArray 

from .surface_solar import get_surface_solar_d_t
from .weather import OutdoorCondition
from config import HOI
from nrain import NRAIN, get_nrains_of_walls
from .thermo_dynamics import ATP, DIFF, get_dgdt, get_dgdu, get_wp, RW, CALDPDU, ROW, CPL, GOFF, WPTRE
from .direction import Direction
from .materials import Materials, Material
from .input_data import InputWall
from .wall_surface import WallSurface

@dataclass
class Layer:

    # 部材No.
    num: int

    # 名称
    name: str

    # 初期温度, ℃
    initial_temp: float

    # 初期湿度, %
    initial_humidity: float

    # 幅, m
    thick: float

    # 分割数
    n_div: int

    # 熱コンダクタンス, W/(m2 K)（外面）
    cond_h_o: float

    # 熱コンダクタンス, W/(m2 K)（内面）
    cond_h_i: float

    # 湿気コンダクタンス, kg/(m2 s Pa)（外面）
    cond_m_o: float

    # 湿気コンダクタンス, kg/(m2 s Pa)（内面）
    cond_m_i: float

    # 斎藤先生コメント
    # 伝達率に関しては正確にはコンダクタンスとして扱っており、防湿層がある場合は透湿抵抗値の逆数（ここでは2.00E-11）を入れます。
    # 材料同士が接触する場合のコンダクタンスは、上の表では熱コンダクタンスは50（W/m2K）、湿気コンダクタンスは7.75E-07（kg/msPa）を暫定的に入れています。
    # 塗膜やクロスなどの透湿抵抗も、ここで調整します。

    # 通気層の場合の流量係数α（高さ1mあたりの値）　三浦コメント：これは高さ1mあたりではなく、幅1mあたりの間違いでは？
    alpha: float

    @property
    def dx(self):
        if self.n_div == 1:
            d1 = 1
        else:
            d1 = self.n_div - 1
        
        return self.thick / d1


@dataclass
class WallType:

    layers: list[Layer]

    @property
    def n_div_total(self):
        return sum([layer.n_div for layer in self.layers])


@dataclass
class WallState:

    # 相対湿度, %
    rh: np.ndarray

    # 温度, deg.C
    TMPC: np.ndarray

    # ステップ n+1 における絶対温度, K
    t_n_pls: np.ndarray

    # ステップ n における絶対温度, K
    HTMP: np.ndarray

    @classmethod
    def init(cls, layers: list[Layer], lookup_table: list[int]):

        ls = [layers[layer_index] for layer_index in lookup_table]

        rh = np.array([l.initial_humidity for l in ls])
        theta = np.array([l.initial_temp for l in ls])
        t = theta + ATP

        return WallState(
            rh=rh,
            TMPC=theta,
            t_n_pls=t,
            HTMP=t
        )









@dataclass
class Wall:

    # 壁表面
    wsurf: WallSurface

    # 壁種類
    kwtype: int

    # 方位
    direction: Direction

    # 相当開口面積（αA）, m2
    alpha_a_ls: list[float]

    # 高さ, m
    height: float

    # 面積, m2
    area: float

    # layers
    layers: np.ndarray

    # メッシュ数
    n_mesh_total: int

    # 状態量
    state: WallState

    # 評価高さ
    eva_height: float

    # 浸水ポイント
    # 厚壁No.
    # wall_no: int
    # 座標
    # pos: int
    # 層No
    # layer_no: int
    # 浸水率
    # ratio: float
    # 閾値風速, m/s
    # v: float

    nrains: list[NRAIN]

    # 雨水浸入ポイントかどうか
    is_rainpoint: list[bool]

    # 雨水浸入の率
    wall_fall_ratio: list[float]

    # 閾値風速, m/s
    wall_fall_wind_threshold: list[float]

    # メッシュ番号に対するレイヤーインデックス
    lookup_table: np.ndarray

    # 質点iの質点間距離, m, [I]
    dx_is: np.ndarray

    # 室外側の端点かどうか, [I]
    is_outside_end_point_is: np.ndarray

    # 室内側の端点かどうか, [I]
    is_inside_end_point_is: np.ndarray

    # 通気層かどうか, [I]
    is_air_layer_is: np.ndarray

    # 質点iの体積, m3, [I]
    v_is: np.ndarray

    # 熱容量, J/(m3 K), [I]
    gcp_is: np.ndarray

    # 熱伝導率, W/(m K), [I]
    lambda_is: np.ndarray

    # 密度, kg/m3, [I]
    gma_is: np.ndarray

    # 材料, [I]
    material_is: list[Material]

    # 室外側の端点における熱コンダクタンス, W/(m2 K), [I]
    cond_h_o_is: np.ndarray

    # 室内側の端点における熱コンダクタンス, W/(m2 K), [I]
    cond_h_i_is: np.ndarray

    # 室外側の端点における湿気コンダクタンス, kg/(m2 s Pa), [I]
    cond_m_o_is: np.ndarray

    # 室内側の端点における湿気コンダクタンス, kg/(m2 s Pa), [I]
    cond_m_i_is: np.ndarray

    # ステップnの絶対温度, K, [I]
    t_n_is: np.ndarray

    # ステップnの水分化学ポテンシャル, J/kg, [I]
    wp_n_is: np.ndarray

    # 水分保有量, kg/s
    rn: np.ndarray

    def p_sv(self, i: int) -> float:
        """飽和水蒸気圧, Pa"""
        return GOFF(t=self.t_n_is[i])[1]

    @property    
    def rh(self):
        """相対湿度, %"""
        return WPTRE(wp=self.wp_n_is, t=self.t_n_is)
    
    @property
    def TMPC(self):
        return self.state.TMPC
    
    def set_TMPC(self, i: int, TMPC: float):
        self.state.TMPC[i] = TMPC

    @property
    def t_n_pls(self):
        return self.state.t_n_pls
    
    def set_t_n_pls(self, i: int, t_n_pls: float):
        self.state.t_n_pls[i] = t_n_pls

    def RMDL(self, i: int):
        """水分伝導率, kg/(m s (J/kg))"""

        if self.rh[i] > 90.0 and self.material_is[i].id == 5:
            return DIFF(rh=self.rh[i], k=self.t_n_pls[i], material=self.material_is[i])
        else:
            return 0.0

    def get_layer(self, i) -> Layer:
        return self.layers[self.lookup_table[i]]
    
    def RMDG(self, i: int) -> float:
        """湿気伝導率, kg /(m s Pa) """
        # 後退差分において湿気伝導率の値は相対湿度に依存性があるが前の時刻の相対湿度を用いることにより線形化している。
        if self.material_is[i].id == 5:
            d1 = self.rh[i] * 0.01
            return 1.87E-11 * d1**2.4019 * math.exp(-0.78864 * ( 1 - d1**1.1471))  #*3.45   !  Moisture conductivity (kg/msPa)  ASHRAE
        else:
            return self.material_is[i].lambda_m
    
    def DPDU(self, i: int):
        """(m3/m3)/(J/kg)"""
        return CALDPDU(wpt=self.wp_n_is[i], tp=self.t_n_is[i], gma=self.gma_is[i], ml0=self.get_layer(i).num, material=self.material_is[i])
    
    def is_outside_surface(self, i: int) -> bool:
        """質点iが室外側表面"""

        return i == 0
    
    def is_inside_surface(self, i: int) -> bool:
        """質点iが室内側表面"""

        return i == self.n_mesh_total - 1

    def dgdu(self, i: int) -> float:
        """水分化学ポテンシャルに対する水蒸気圧の微分, Pa / (J / kg)"""
        return get_dgdu(rh=self.rh[i], t=self.t_n_pls[i])

    def dgdt(self, i: int) -> float:
        """絶対温度に対する水蒸気圧の微分, Pa / K"""
        return get_dgdt(rh=self.rh[i], t=self.t_n_pls[i])
    
    def cap(self, i: int) -> float:
        """熱容量, J/K"""

        # 質点iが空気層の場合
        if self.is_air_layer_is[i]:
            #1300は空気の体積熱容量, J/(kg K)
            return 1300.0 * self.dx_is[i] * self.area
        else:
            return self.v_is[i] * self.gcp_is[i]

    def cap_m(self, i: int) -> float:
        """水分移動に伴う容量項, kg/(J/kg)"""

        # 質点iが空気層の場合
        if self.is_air_layer_is[i]:
            # 空気層の場合の熱容量はなし
            # ステップnにおける水分量は移流の計算項で表現されている。
            return 0.0
        else:
            # ROW: 水の密度, kg/m3
            # DPDU: 容積含水率を水分化学ポテンシャルで微分した値, (m3/m3)/(J/kg)
            # v_is: 容積, m3
            # kg/(J/kg) = kg/m3 * (m3/m3)/(J/kg) * m3
            return ROW * self.DPDU(i=i) * self.v_is[i]

    def _c_h_t_i_mns(self, i: int) -> float:
        """温度差を駆動力とする熱移動に関する係数（室外側）, W/(K m2)"""
    
        # 各レイヤの室外側の端点の場合
        if self.is_outside_end_point_is[i]:
                                    
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            return self.cond_h_o_is[i]

        else:
            return self.get_series_combination(
                x=(self.lambda_is[i-1]) / self.dx_is[i-1],
                y=(self.lambda_is[i]) / self.dx_is[i]
            )
    
    def _c_h_t_i_pls(self, i: int) -> float:
        """温度差を駆動力とする熱移動に関する係数（室内側）, W/(K m2)"""

        # 各レイヤの室内側の端点の場合
        if self.is_inside_end_point_is[i]:

            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            return self.cond_h_i_is[i]

        else:
            return self.get_series_combination(
                x=(self.lambda_is[i+1]) / self.dx_is[i+1], 
                y=(self.lambda_is[i]) / self.dx_is[i]
            )

    def _c_vap_t_i_mns(self, i: int) -> float:
        """温度差を駆動力とする水蒸気移動に関する係数（室外側）, (kg/(s m2))/K"""

        # 各レイヤの室外側の端点の場合
        if self.is_outside_end_point_is[i]:
                                    
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            return self.cond_m_o_is[i] * self.dgdt(i)

        else:
            # RMDG：湿気伝導率, kg / (m s Pa)
            # dgdt：水蒸気圧の温度勾配, Pa/K
            return self.get_series_combination(
                x=(self.RMDG(i-1) * self.dgdt(i-1)) / self.dx_is[i-1],
                y=(self.RMDG(i) * self.dgdt(i)) / self.dx_is[i]
            )

    def _c_vap_t_i_pls(self, i: int) -> float:
        """温度差を駆動力とする水蒸気移動に関する係数（室内側）, (kg/(s m2))/K"""

        # 各レイヤの室内側の端点の場合
        if self.is_inside_end_point_is[i]:

            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            return self.cond_m_i_is[i] * self.dgdt(i)

        else:
            return self.get_series_combination(
                x=(self.RMDG(i+1) * self.dgdt(i+1)) / self.dx_is[i+1],
                y=(self.RMDG(i) * self.dgdt(i)) / self.dx_is[i]
            )

    def _c_vap_wp_i_mns(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする水蒸気移動に関する係数（室外側）, (kg/(s m2))/(J/kg)"""

        # 各レイヤの室外側の端点の場合
        if self.is_outside_end_point_is[i]:
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            # 湿気コンダクタンス, kg/(m2 s Pa) * Pa/(J/kg) * m2 = (kg/s)/(J/kg) 
            return self.cond_m_o_is[i] * self.dgdu(i)
        else:
            # RMDG: 湿気伝導率, kg/(m s Pa)
            # dgdu: 水分化学ポテンシャルに対する水蒸気圧の微分, Pa / (J/kg)
            return self.get_series_combination(
                x=self.RMDG(i-1) * self.dgdu(i-1) / self.dx_is[i-1],
                y=self.RMDG(i) * self.dgdu(i) / self.dx_is[i]
            )
    
    def _c_vap_wp_i_pls(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする水蒸気移動に関する係数（室内側）, (kg/(s m2))/(J/kg)"""

        # 各レイヤの室内側の端点の場合
        if self.is_inside_end_point_is[i]:
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            # 湿気コンダクタンス, kg/(m2 s Pa) * Pa/(J/kg) * m2 = W/(J/kg) 
            return self.cond_m_i_is[i] * self.dgdu(i)
        else:
            return self.get_series_combination(
                x=self.RMDG(i+1) * self.dgdu[i+1] / self.dx_is[i+1],
                y=self.RMDG(i) * self.dgdu[i] / self.dx_is[i]
            )
    
    def _c_liq_wp_i_mns(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする液水移動に関する係数（室外側）, (kg/(s m2))/(J/kg)"""

        if self.is_outside_end_point_is[i]:
            return 0.0
        else:
            return self.get_series_combination(
                x=self.RMDL(i-1) / self.dx_is[i-1],
                y=self.RMDL(i) / self.dx_is[i]
            )
    
    def _c_liq_wp_i_pls(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする液水移動に関する係数（室内側）, (kg/(s m2))/(J/kg)"""

        if self.is_inside_end_point_is[i]:
            return 0.0
        else:
            return self.get_series_combination(
                x=self.RMDL(i+1) / self.dx_is[i+1],
                y=self.RMDL(i) / self.dx_is[i]
            )

    def get_t_n_pls(self, t_is: np.ndarray, dt: float, oc: OutdoorCondition, theta_r_n: float, wp_r_n: float, QQ: float, t_upstream: float):

        t_is_next = np.zeros_like(t_is, dtype=float)

        for i in range(self.n_mesh_total):

            ### 質点の温度
            # 室外側
            t_i_mns = oc.t if self.is_outside_surface(i=i) else t_is[i - 1]
            # 中央（後退差分計算において温度については繰り返し計算に用いる温度を採用する。）
            t_i = t_is[i]
            # 室内側
            t_i_pls = theta_r_n if self.is_inside_surface(i=i) else t_is[i + 1]

            ### 質点の水分化学ポテンシャル
            # 室外側
            wp_i_mns = oc.wp if self.is_outside_surface(i=i) else self.wp_n_is[i - 1]
            # 中央（後退差分計算において水分化学ポテンシャルについては前のステップの値を用いる。）
            wp_i = self.wp_n_is[i]
            # 室内側
            wp_i_pls = wp_r_n if self.is_inside_surface(i=i) else self.wp_n_is[i + 1]

            # 日射による吸収熱量, W
            q_sol_d_t_k = self.wsurf.get_q_sol_d_t_k(oc=oc) if self.is_outside_surface(i=i) else 0.0

            # 温度差を駆動力とする熱移動に関する係数, W/K
            c_h_t_i_mns = (self._c_h_t_i_mns(i) + RW * self._c_vap_t_i_mns(i)) * self.area
            c_h_t_i_pls = (self._c_h_t_i_pls(i) + RW * self._c_vap_t_i_pls(i)) * self.area

            # 水分化学ポテンシャル差を駆動力とする水蒸気移動に伴う熱移動に関する係数, W/(J/kg)
            c_wp_i_mns = RW * self._c_vap_wp_i_mns(i) * self.area
            c_wp_i_pls = RW * self._c_vap_wp_i_pls(i) * self.area

            # 液水移動量, kg/s
            # 温度差駆動の液水移動量は十分小さいため無視する。
            j_liq_i_mns = self._c_liq_wp_i_mns(i) * (wp_i_mns - wp_i) * self.area
            j_liq_i_pls = self._c_liq_wp_i_pls(i) * (wp_i_pls - wp_i) * self.area

            cap = self.cap(i) / dt

            # W
            UHEN = (
                cap * self.t_n_is[i]
                + q_sol_d_t_k
                + c_h_t_i_mns * t_i_mns
                + c_h_t_i_pls * t_i_pls
                + c_wp_i_mns * (wp_i_mns - wp_i)
                + c_wp_i_pls * (wp_i_pls - wp_i)
                + CPL * j_liq_i_mns * t_i_mns
                + CPL * j_liq_i_pls * t_i_pls
            )

            # W/K
            SAHEN = (
                cap
                + c_h_t_i_mns
                + c_h_t_i_pls
                + CPL * j_liq_i_mns
                + CPL * j_liq_i_pls
            )

            if self.is_air_layer_is:
                # 空気層の場合に移流分を考慮する。
                # 空気の容積比熱, J/(m3 K)                             
                UHEN =+ 1300.0 * QQ * t_upstream
                SAHEN =+ 1300.0 * QQ
            
            t_is_next[i] = UHEN / SAHEN

        return t_is_next
          
    def get_wp_n_pls(
            self, dt: float, oc: OutdoorCondition, theta_r_n: float, wp_r_n: float, wp_is: np.ndarray,
            QQ: float, RN: np.ndarray, XM: np.ndarray, WJRAIN: np.ndarray, WJW: np.ndarray):
        """ステップn+1における水分化学ポテンシャルを求める。

        Args:
            dt (float): _description_
            oc (OutdoorCondition): _description_
            theta_r_n (float): _description_
            wp_r_n (float): _description_
            wp_is (np.ndarray): _description_
            QQ: 換気量, m3/s
            RN (np.ndarray): _description_
            XM (np.ndarray): _description_
            WJRAIN: 雨水浸入量, kg/s
            WJW (np.ndarray): _description_

        Returns:
            _type_: _description_
        """

        wp_n_pls = np.zeros_like(self.wp_n_is, dtype=float)

        for i in range(self.n_mesh_total):
            
            # 質点の温度
            # 室外側
            t_i_mns = oc.t_k if self.is_outside_surface(i=i) else self.t_n_is[i - 1]
            # 中央（後退差分計算において温度については前のステップの値を用いる。）
            t_i = self.t_n_is[i]
            # 室内側
            t_i_pls = theta_r_n + ATP if self.is_inside_surface(i=i) else self.t_n_is[i + 1]

            # 質点の水分化学ポテンシャル
            # 室外側
            wp_i_mns = oc.wp if self.is_outside_surface(i=i) else wp_is[i - 1]
            # 中央
            wp_i = wp_is[i]
            # 室内側
            wp_i_pls = wp_r_n if self.is_inside_surface(i=i) else wp_is[i + 1]

            # 水分化学ポテンシャル差を駆動力とする水分（水蒸気＋液水）移動に関する係数, (kg/s) / (J/kg)
            c_vap_liq_wp_i_mns = (self._c_vap_wp_i_mns(i=i) + self._c_liq_wp_i_mns(i=i)) * self.area
            c_vap_liq_wp_i_pls = (self._c_vap_wp_i_pls(i=i) + self._c_liq_wp_i_pls(i=i)) * self.area

            # 温度差を駆動力とする水蒸気移動量, kg/s
            # 温度差を駆動力とする液水移動量は非常に小さいため無視する。
            j_vap_i_mns = self._c_vap_t_i_mns(i=i) * self.area * (t_i_mns - t_i)
            j_vap_i_pls = self._c_vap_t_i_pls(i=i) * self.area * (t_i_pls - t_i)

            # (kg/s)/(J/kg)
            # m_cap: kg/(J/kg)
            SAHEN = (
                self.cap_m(i=i) / dt
                + c_vap_liq_wp_i_mns
                + c_vap_liq_wp_i_pls
            )

            # kg/s
            UHEN = (
                self.cap_m(i=i) / dt * self.wp_n_is[i]
                + c_vap_liq_wp_i_mns * wp_i_mns
                + c_vap_liq_wp_i_pls * wp_i_pls
                + j_vap_i_mns
                + j_vap_i_pls
            )

            PXCOF = 1.0 / 133322.0

            # 湿気伝達率, kg / m2 s Pa
            alpha_dsh_m = 3.43e-8

            # 水膜の濡れ面積率
            r_wet = 0.3

            if self.is_air_layer_is[i]:

                # PXCOF: 絶対湿度を水蒸気圧にかえる係数
                # PXCOF = 1. / 133322.  (kg/m3)/Pa
                # 133322: エクセルで絶対湿度と水蒸気圧 

                # (kg/s)/(J/kg) = ??? * Pa/(J/kg) * ((kg/m3)/Pa) * m3/s
                # もとのプログラムから1.2を消した
                j_wtr_vent_wp = self.dgdu(i) * PXCOF * QQ
                # (kg/s)/K = kg/m3 * Pa/K * (1/Pa) * m3/s
                j_wtr_vent_k = self.dgdt(i) * PXCOF * QQ

                # RN: 水膜の保持水分量, kg/m2

                # 水膜からの水分移動量, kg/s
                j_dsh_vap_i_pls = alpha_dsh_m * (self.p_sv(i+1) - XM(i)) * self.area * r_wet if RN[i+1] > 0.0 else 0.0

                j_dsh_vap_i_mns = alpha_dsh_m * (self.p_sv(i-1) - XM(i)) * self.area * r_wet if RN[i-1] > 0.0 else 0.0

                # 風上側は外気にした。（もともとは上流側の通気層の状態量が入ることになっていた。）
                UHEN += (
                    j_wtr_vent_wp * oc.wp
                    + j_wtr_vent_k * (oc.t_k - t_i)
                    + j_dsh_vap_i_pls
                    + j_dsh_vap_i_mns
                )

                SAHEN += j_wtr_vent_wp

            # WJRAIN 雨水由来の浸入量
            D7 = WJRAIN[i]

            # WJW：木材が分解した場合にセルロースが分解された場合に発生する水分量, kg/(m3 s)
            D5 = WJW(i) * self.dx_is[i] * self.area

            UHEN += D5 + D7

            wp_next = UHEN / SAHEN

            # ゼロを超えない措置。
            if self.is_air_layer_is[i]:
                if wp_next >= 0.0:
                    wp_next = -1.3E-4
            else:
                if wp_next >= 0.0:
                    wp_next = -1.3E-3
            
            wp_n_pls[i] = wp_next

        return wp_n_pls

    def get_series_combination(x: float, y: float) -> float:

        if x + y != 0.0:
            return x * y / (x + y)
        else:
            return 0.0

    def get_v_wind_eva_k(self, v_wind: float):
        """評価高さにおける風速を求める。

        Args:
            v_wind: 風速, m/s

        Returns:
            評価高さにおける風速, m/s
        """

        # 基準風速は6.5m高さ
        return v_wind * (self.eva_height / 6.5)**0.25
    
    def get_rn_n_pls(self, wp_is: np.ndarray, oc: OutdoorCondition, p_v_rm: float, dt: float):
        
        rn = np.zeros(shape=self.n_mesh_total, dtype=float)
        wjrain = np.zeros(shape=self.n_mesh_total, dtype=float)

        for i in range(self.n_mesh_total):
            
            if self.is_rainpoint:

                if self.is_outside_surface:
                    mns = oc.xod
                elif self.is_air_layer_is[i-1]:
                    rh = WPTRE(wp=wp_is[i-1], t=self.t_n_is[i-1])
                    _, vsp = GOFF(t=self.t_n_is[i-1])
                    mns = vsp * rh * 0.01
                else:
                    mns = wp_is[i-1]

                if self.is_inside_surface:
                    pls = p_v_rm
                elif self.is_air_layer_is[i+1]:
                    rh = WPTRE(wp=wp_is[i+1], t=self.t_n_is[i+1])
                    _, vsp = GOFF(t=self.t_n_is[i+1])
                    pls = vsp * rh * 0.01
                else:
                    pls = wp_is[i+1]

                # （材料番号が10～12）（セメント系材料：サイディングとか）
                # 水分伝導率, kg/ms(J/kg)
                if self.material_is[i].id >= 10 and self.material_is[i].id <= 12:
                    # 飽和時の水分伝導率 3.73e-6 kg/ms(J/kg)　いぶし瓦 by 伊庭　D論
                    rmdl = 3.73e-6
                else:
                    # 木製品を想定
                    rmdl = 3.73e-6 * 0.05
                # バックシーラー透水抵抗 2.4e+5 m2sPa/kg by 長村
                r = 2.4e5
                # コンダクタンス (kg/s) / m2 (J/kg)
                c = 1 / (self.dx_is[i] / rmdl) + r / self.dgdu(i)

                    
                # 当該質点が飽和している前提で計算する。
                # 3.43E-08：湿気伝達率 kg/(m2 s Pa))
                # 水分流の計算, kg/(m2 s)                            
                # 0.3 = 濡れ面率
                if self.is_outside_surface or self.is_air_layer_is[i-1]:
                    p_sv = GOFF(t=self.t_n_is[i])
                    f_mns = 3.43e-8 * (mns - p_sv) * 0.3
                else:
                    f_mns = c * (mns - wp_is[i])                                
                if self.is_inside_surface or self.is_air_layer_is[i+1]:
                    p_sv = GOFF(t=self.t_n_is[i])
                    f_pls = 3.43e-8 * (pls - p_sv) * 0.3
                else:
                    f_pls = c * (pls - wp_is[i])

                # 水幕の水分量（保持している水の量）, kg/m2
                # D2 はポテンシャル
                # D4 は蒸発量
                # HRN は前のステップの水分量
                
                # 水膜の水分量の計算, kg/m2
                rn = (f_mns + f_pls + self.get_swjrain(oc=oc, i=i)) * dt + self.rn[i]

                # 水幕からの吸水量（表面に水幕があって材料に吸われる分）
                # 水幕が残っている場合は飽和水蒸気圧（水分伝導率で計算）
                if rn < 0.0:
                    rn_n_pls = 0.0
    
                    # RNがゼロの場合は水幕がないので、雨水が直接材料に吸われることになる。
                    w = (self.get_swjrain(oc=oc, i=i) + self.rn[i]) * self.area
                else:
                    rn_n_pls = rn
                    w = c * (0 - self.wp_n_is) * self.area
                
                rn[i] = rn_n_pls
                wjrain[i] = w
        
        return rn, wjrain




    def get_confrains(self, v_wind: float, wind_direction: float):
        """_summary_

        Args:
            v_wind: 風速, m/s
        """

        v_mod = self.get_v_wind_eva_k(v_wind=v_wind)

        x = np.linspace(0.1 * v_mod, 3.0 * v_mod, 30)

        confrains = np.zeros_like(self.nrains, dtype=float)

        for i, nrain in enumerate(self.nrains):

            if v_wind > 0.1:

                # 1.5 = 

                # ASHRAE 160-2009
                # 雨水暴露係数
                FE = 1.5
                # 雨水付着係数
                FD = 1.0
                # 経験的な定数, kg s / (m3 mm)
                FL = 0.2

                # TODO: Wind_direction の定義をきちんと確認しないといけない。この式のままだと、北側から時計回りか？
                # 水平方向などにも対応させないといけないのではないか？
                # この式は垂直壁にしか対応していないので、3次元的にcosを計算する必要あり。
                d2 = np.maximum(
                    nrain.ratio * x * np.cos(np.radians(wind_direction - (180.0 + self.direction.alpha))) * FE * FD * FL,
                    0.0
                )

                sigma =  v_mod / np.sqrt(np.pi / 2)
                weight = rayleigh.pdf(x=x, scale=sigma)

                confrain = np.sum(weight * d2) / np.sum(weight)

            else:

                confrain = 0.0
            
            confrains[i] = confrain

        return confrains
    
    def get_swjrain(self, oc: OutdoorCondition, i: int):
        """メッシュ番号iで指定されたセルへの浸水量(kg/m2s)を計算する。

        Args:
            oc (OutdoorCondition): _description_
            i: メッシュ番号

        Returns:
            _type_: _description_
        """

        confrains = self.get_confrains(v_wind=oc.v_wind, wind_direction=oc.wind_direction)

        # 材表面への浸水量 kg/(m2 s)
        # rainfall mm/h
        # confrains: 
        # 浸水率をパーセントでいれているので、ここで単位換算している。
        swjrains = 0.01 * oc.rainfall * confrains / 3600.0

        for (nrain, swjrain) in (self.nrains, swjrains):
            if nrain.pos - 1 == i:
                return swjrain
        
        return 0.0

    @classmethod
    def read(cls, d: dict, i: int):

        ipt_wall = InputWall.read(d=d)

        layers = [
            Layer(num=10, name='サイディング', initial_temp=15.0, initial_humidity=50.0, thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
            Layer(num=2, name='通気層', initial_temp=15.0, initial_humidity=50.0, thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
            Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
            Layer(num=5, name='構造用合板1', initial_temp=10.0, initial_humidity=80.0, thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
            Layer(num=3, name='グラスウール1', initial_temp=10.0, initial_humidity=50.0, thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
            Layer(num=5, name='構造用合板1', initial_temp=10.0, initial_humidity=80.0, thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
            Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
        ]

        # 面積, m2
        area = ipt_wall.len_long * ipt_wall.len_short

        # Layerそれぞれにおける室外側と室内側のメッシュ番号, [L], [L]
        outside_end_point_mesh_indices, inside_end_point_mesh_indices = _get_first_and_last_mesh_indices(layers=layers)

        # 相当開口面積（αA）, m2
        alpha_a_ls = [layer.alpha * layer.thick * ipt_wall.len_short for layer in layers]

        # ある質点がどのレイヤーに対応するかを保持するリスト, [I]
        lookup_table = create_lookup_table(layers)

        # 質点の総数
        n_mesh_total = sum([layer.n_div for layer in layers])

        # 室外側の端点かどうか, [I]
        is_outside_end_point_is = np.array([outside_end_point_mesh_indices[layer_index] == i for (i, layer_index) in enumerate(lookup_table)])

        # 室内側の端点かどうか, [I]
        is_inside_end_point_is = np.array([inside_end_point_mesh_indices[layer_index] == i for (i, layer_index) in enumerate(lookup_table)])

        # 熱コンダクタンス（室外側）, W/(m2 K), [I]
        cond_h_o_is = np.array([layers[layer_index].cond_h_o if is_outside_end_point_is[i] else 0.0 for (i, layer_index) in enumerate(lookup_table)])

        # 熱コンダクタンス（室内側）, W/(m2 K), [I]
        cond_h_i_is = np.array([layers[layer_index].cond_h_i if is_inside_end_point_is[i] else 0.0 for (i, layer_index) in enumerate(lookup_table)])

        # 湿気コンダクタンス（室外側）, kg/(m2 s Pa), [I]
        cond_m_o_is = np.array([layers[layer_index].cond_m_o if is_outside_end_point_is[i] else 0.0 for (i, layer_index) in enumerate(lookup_table)])

        # 湿気コンダクタンス（室内側）, kg/(m2 s Pa), [I]
        cond_m_i_is = np.array([layers[layer_index].cond_m_i if is_inside_end_point_is[i] else 0.0 for (i, layer_index) in enumerate(lookup_table)])

        # 通気層かどうか, [I]
        is_air_layer_is = np.array([layers[layer_index].name == '通気層' for layer_index in lookup_table])

        # 質点iの質点間距離, m, [I]
        # レイヤー表面の質点において、SideA側の質点の場合はSideB側のみが、SideB側の質点の場合はSideA側のみが定義される。
        dx_is = np.array([layers[layer_index].dx for layer_index in lookup_table])

        # 材料, Material クラス, [I]
        ms = Materials()
        material_is = [ms.get_material(name=layers[layer_index].name) for layer_index in lookup_table]

        # 質点iの容積比熱, J/(m3 K), [I]
        gcp_is = np.array([material_i.c * material_i.rho for material_i in material_is])

        # 質点iの熱伝導率, W/(m K), [I]
        lambda_is = np.array([material_i.lambda_h for material_i in material_is])

        # 質点iの密度, kg/m3, [I]
        gma_is = np.array([material_i.rho for material_i in material_is])

        # 質点iの体積, m3, [I]
        # 端点の場合は体積が端点以外の部分の半分になる。（空気層は除く。）
        v_is = np.where(is_outside_end_point_is | is_inside_end_point_is, 0.5, 1.0) * dx_is * area
        
        state = WallState.init(layers, lookup_table)

        wsurf = WallSurface.read(iw=ipt_wall)

        # 方位
        direction = ipt_wall.direction

        # 壁の下端と上端の高さの差（換気計算に用いられる）, m
        height = ipt_wall.height

        # 評価高さ, m
        eva_height = ipt_wall.eva_height

        # NRAIN Class
        # 厚壁No. 
        # wall_no: int
        # 座標
        # pos: int
        # 層No
        # layer_no: int
        # 浸水率
        # ratio: float
        # 閾値風速, m/s
        # v: float

        nrains_list: list[NRAIN] = get_nrains_of_walls(wall_index=i)

        # 雨水浸入ポイントかどうか
        is_rainpoint = np.full(n_mesh_total, False)

        # 雨水浸入の率
        wall_fall_ratio = np.zeros(n_mesh_total)

        # 閾値風速, m/s
        wall_fall_wind_threshold = np.zeros(n_mesh_total)

        for nrain in nrains_list:
            is_rainpoint[nrain.pos] = True
            wall_fall_ratio[nrain.pos] = nrain.ratio
            wall_fall_wind_threshold[nrain.pos] = nrain.v

        # 初期温度, K, [I]
        t_init_is = [layers[layer_index].initial_temp + ATP for layer_index in lookup_table]

        rn = np.zeros(n_mesh_total)

        return Wall(
            kwtype=d['kwtype'],
            direction=direction,
            alpha_a_ls=alpha_a_ls,
            height=height,
            area=area,
            layers=layers,
            lookup_table=lookup_table,
            n_mesh_total=n_mesh_total,
            state=state,
            eva_height=eva_height,
            nrains=nrains_list,
            dx_is=dx_is,
            is_outside_end_point_is=is_outside_end_point_is,
            is_inside_end_point_is=is_inside_end_point_is,
            is_air_layer_is=is_air_layer_is,
            v_is=v_is,
            gcp_is=gcp_is,
            lambda_is=lambda_is,
            gma_is=gma_is,
            material_is=material_is,
            cond_h_o_is=cond_h_o_is,
            cond_h_i_is=cond_h_i_is,
            cond_m_o_is=cond_m_o_is,
            cond_m_i_is=cond_m_i_is,
            t_n_is=t_init_is,
            wsurf=wsurf,
            is_rainpoint=is_rainpoint,
            wall_fall_ratio=wall_fall_ratio,
            wall_fall_wind_threshold=wall_fall_wind_threshold,
            rn=rn
        )

    @classmethod
    def read_default(cls):

        ds = get_ds()

        walls = [Wall.read(d=d, i=i) for i, d in enumerate(ds)]

        return walls


def create_lookup_table(layers: list[Layer]) -> np.ndarray:

    lookup_table = []

    for i, layer in enumerate(layers):
        lookup_table.extend([i] * layer.n_div)
    
    return np.array(lookup_table)


def _get_first_and_last_mesh_indices(layers: list[Layer]) -> tuple[list[int], list[int]]:
    
    first_mesh_indices = [0] * len(layers)
    last_mesh_indices = [0] * len(layers)

    for i, layer in enumerate(layers):

        if i == 0:
            last_mesh_indices[i] = 0
        else:
            last_mesh_indices[i] = first_mesh_indices[i - 1] + 1

        first_mesh_indices[i] = last_mesh_indices[i] + layer.n_div - 1
   
    return first_mesh_indices, last_mesh_indices


def get_ds():

    # 換気計算の高さと壁の高さが傾斜がある場合は一致しない。
    return [
        {
            'kwtype': 1,
            'direction': 'e',
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'angle': 90.0,
            'emissivity': 0.9,
            'walltypes': 1,
            'eva_height': 4.2,
        },
    ]

