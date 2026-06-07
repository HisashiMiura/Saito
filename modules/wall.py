from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray 

from .weather import OutdoorCondition
from config import HOI
from nrain import NRAIN, get_nrains_of_walls
from .thermo_dynamics import ATP, DIFF, get_dpv_dt, get_dpv_dmu, RW, ROW, CPL, GOFF, get_rh, get_rho
from .direction import Direction
from .materials import Materials, Material
from .input_wall import InputWall
from .wall_surface import WallSurface
from modules.state import State
from modules.config import INITIAL_WALL_TEMPERATURE, INITIAL_WALL_RERATIVE_HUMIDITY
from modules import nrain
from modules.cell import (
    Cell, 
    CellOutsideEndPoint, CellInsideEndPoint, CellOutsideEndPointAirLayer, CellInsideEndPointAirLayer,
    CellInterior, CellOutsideSurface, CellInsideSurface, CellAirLayer
)
from modules import wood_decay


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

    # 絶対温度, K
    t: np.ndarray

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
            t=t,
            HTMP=t
        )
    
    @property
    def theta(self):
        return self.t - ATP
    
    @theta.setter
    def theta(self, value):
        self.t = value + ATP







@dataclass
class Wall:

    # 壁表面
    wsurf: WallSurface

    # 相当開口面積（αA）, m2
    alpha_a_ls: list[float]

    # 高さ, m
    height: float

    # 面積, m2
    area: float

    # メッシュ数
    n_mesh_total: int

    # 質点iの質点間距離, m, [I]
    dx_is: np.ndarray

    # 通気層かどうか, [I]
    is_air_layer_is: np.ndarray

    # 密度, kg/m3, [I]
    gma_is: np.ndarray

    # 材料, [I]
    material_is: list[Material]

    # ステップnの絶対温度, K, [I]
    t_n_is: np.ndarray

    # ステップnの水分化学ポテンシャル, J/kg, [I]
    wp_n_is: np.ndarray

    # 水分保有量, kg/s
    rn: np.ndarray

    # 発芽までの蓄積時間, [I]
    time_s_is: np.ndarray

    # 発芽の有無, [I]
    l_stage_is: np.ndarray

    # 腐朽により減少した質量の割合, [I]
    m_loss: np.ndarray

    wjw: np.ndarray

    # ステップnの状態量, [I]
    state_n_is: list[State]

    rain_leakage_is: list[nrain.RainLeakage]

    cells_is: list[Cell]

    def dpv_dmu(self, i: int) -> float:
        """水蒸気圧を水分化学ポテンシャルで偏微分（絶対温度一定）した値, Pa / (J / kg)"""
        return get_dpv_dmu(rh=self.state_n_is[i].rh, t=self.state_n_is[i].t)

    def dpv_dt(self, i: int) -> float:
        """水蒸気圧を絶対温度で偏微分（水分化学ポテンシャル一定）した値, Pa / K"""
        return get_dpv_dt(rh=self.state_n_is[i].rh, t=self.state_n_is[i].t)
    
    def _c_h_t_i_mns(self, i: int) -> float:
        """温度差を駆動力とする熱移動に関する係数（室外側）, W/(K m2)"""
        return 1 / (self.cells_is[i-1].r_h_pls + self.cells_is[i].r_h_mns)
    
    def _c_h_t_i_pls(self, i: int) -> float:
        """温度差を駆動力とする熱移動に関する係数（室内側）, W/(K m2)"""
        return 1 / (self.cells_is[i+1].r_h_mns + self.cells_is[i].r_h_pls)

    def _c_vap_t_i_mns(self, i: int) -> float:
        """温度差を駆動力とする水蒸気移動に関する係数（室外側）, (kg/(s m2))/K"""
        # dpv_dt：水蒸気圧の温度勾配, Pa/K
        # r_m：透湿抵抗, (m2 s Pa)/kg
        return 1 / (self.cells_is[i-1].r_m_pls / self.dpv_dt(i-1) + self.cells_is[i].r_m_mns / self.dpv_dt(i))

    def _c_vap_t_i_pls(self, i: int) -> float:
        """温度差を駆動力とする水蒸気移動に関する係数（室内側）, (kg/(s m2))/K"""
        return 1 / (self.cells_is[i+1].r_m_mns / self.dpv_dt(i+1) + self.cells_is[i].r_m_pls / self.dpv_dt(i))

    def _c_vap_wp_i_mns(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする水蒸気移動に関する係数（室外側）, (kg/(s m2))/(J/kg)"""
        # RMDG: 湿気伝導率, kg/(m s Pa)
        # dgdu: 水分化学ポテンシャルに対する水蒸気圧の微分, Pa / (J/kg)
        return 1 / (self.cells_is[i-1].r_h_pls / self.dpv_dmu(i-1) + self.cells_is[i].r_h_mns / self.dpv_dmu(i))
    
    def _c_vap_wp_i_pls(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする水蒸気移動に関する係数（室内側）, (kg/(s m2))/(J/kg)"""
        return 1 / (self.cells_is[i+1].r_h_mns / self.dpv_dmu[i+1] + self.cells_is[i].r_h_pls / self.dpv_dmu[i])
    
    def _c_liq_wp_i_mns(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする液水移動に関する係数（室外側）, (kg/(s m2))/(J/kg)"""
        return 1 / (self.cells_is[i-1].r_liq_wp_pls + self.cells_is[i].r_liq_wp_mns)

    def _c_liq_wp_i_pls(self, i: int) -> float:
        """水分化学ポテンシャル差を駆動力とする液水移動に関する係数（室内側）, (kg/(s m2))/(J/kg)"""
        return 1 / (self.cells_is[i+1].r_liq_wp_mns + self.cells_is[i].r_liq_wp_pls)

    def get_t_next_is(self, t_is: np.ndarray, dt: float, oc_n_pls: OutdoorCondition, t_r_n_pls: float, wp_r_n_pls: float, v_air_n: float, t_upstream_n_pls: float):
        """反復法における次の計算の温度を求める。

        Args:
            t_is: 温度（反復法における一時的な温度）, K, [I]
            dt: 時間刻み幅, s
            oc_n_pls: ステップn+1における外気条件
            t_r_n_pls: ステップn+1における室温, K
            wp_r_n_pls: ステップn+1における室内の水分化学ポテンシャル, J/kg
            v_air_n: ステップnからステップn+1における通気層内の空気の平均流速, m/s
            t_upstream_n_pls: ステップn+1における通気層内に流入する空気の温度, K

        Returns:
            収束計算における次の計算の温度, K
        """

        # 反復法における次の計算の温度（の入れ物）, K, [I]
        t_next_is = np.zeros_like(t_is, dtype=float)

        for i in range(self.n_mesh_total):

            ### 質点の温度, K
            # マイナス側の温度
            # 質点iが室外側の端点の場合は外気温を用いる。
            t_i_mns = oc_n_pls.t_k if isinstance(self.cells_is[i], CellOutsideSurface) else t_is[i - 1]
            # 中央（後退差分計算において温度については繰り返し計算に用いる温度を採用する。）
            t_i = t_is[i]
            # プラス側の温度
            # 質点iが室内側の端点の場合は室内温度を用いる。
            t_i_pls = t_r_n_pls if isinstance(self.cells_is[i], CellInsideSurface) else t_is[i + 1]

            ### 質点の水分化学ポテンシャル, J/kg, [I]
            # 室外側
            # 質点iが室外側の端点の場合は外気の水分化学ポテンシャルを用いる。
            wp_i_mns = oc_n_pls.wp if isinstance(self.cells_is[i], CellOutsideSurface) else self.wp_n_is[i - 1]
            # 中央（後退差分計算において水分化学ポテンシャルについては前のステップの値を用いる。）
            wp_i = self.wp_n_is[i]
            # 室内側
            # 質点iが室内側の端点の場合は室内の水分化学ポテンシャルを用いる。
            wp_i_pls = wp_r_n_pls if isinstance(self.cells_is[i], CellInsideSurface) else self.wp_n_is[i + 1]

            # 日射による吸収熱量, W
            # 質点iが室外側の端点の場合は日射による吸収熱量を考慮する。
            q_sol_d_t_k = self.wsurf.get_q_sol_d_t_k(oc=oc_n_pls) if isinstance(self.cells_is[i], CellOutsideSurface) else 0.0

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

            # 熱容量を時間刻みで除した値, W/K
            cap = self.cells_is[i].cap * self.area / dt

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
                UHEN =+ 1300.0 * v_air_n[i] * t_upstream_n_pls
                SAHEN =+ 1300.0 * v_air_n[i]
            
            t_next_is[i] = UHEN / SAHEN

        return t_next_is

    def get_wp_n_pls(
            self, wp_is: np.ndarray, dt: float, oc: OutdoorCondition, theta_r_n: float, wp_r_n: float,
            v_air_is: np.ndarray, RN: np.ndarray, WJRAIN: np.ndarray):
        """ステップn+1における水分化学ポテンシャルを求める。

        Args:
            dt (float): _description_
            oc (OutdoorCondition): _description_
            theta_r_n (float): _description_
            wp_r_n (float): _description_
            wp_is (np.ndarray): _description_
            QQ: 換気量, m3/s
            RN (np.ndarray): _description_
            WJRAIN: 雨水浸入量, kg/s

        Returns:
            _type_: _description_
        """

        wp_n_pls = np.zeros_like(self.wp_n_is, dtype=float)

        for i in range(self.n_mesh_total):
            
            # 質点の温度
            # 室外側
            t_i_mns = oc.t_k if isinstance(self.cells_is[i], CellOutsideSurface) else self.t_n_is[i - 1]
            # 中央（後退差分計算において温度については前のステップの値を用いる。）
            t_i = self.t_n_is[i]
            # 室内側
            t_i_pls = theta_r_n + ATP if isinstance(self.cells_is[i], CellInsideSurface) else self.t_n_is[i + 1]

            # 質点の水分化学ポテンシャル
            # 室外側
            wp_i_mns = oc.wp if isinstance(self.cells_is[i], CellOutsideSurface) else wp_is[i - 1]
            # 中央
            wp_i = wp_is[i]
            # 室内側
            wp_i_pls = wp_r_n if isinstance(self.cells_is[i], CellInsideSurface) else wp_is[i + 1]

            # 水分化学ポテンシャル差を駆動力とする水分（水蒸気＋液水）移動に関する係数, (kg/s) / (J/kg)
            c_vap_liq_wp_i_mns = (self._c_vap_wp_i_mns(i=i) + self._c_liq_wp_i_mns(i=i)) * self.area
            c_vap_liq_wp_i_pls = (self._c_vap_wp_i_pls(i=i) + self._c_liq_wp_i_pls(i=i)) * self.area

            # 温度差を駆動力とする水蒸気移動量, kg/s
            # 温度差を駆動力とする液水移動量は非常に小さいため無視する。
            j_vap_i_mns = self._c_vap_t_i_mns(i=i) * self.area * (t_i_mns - t_i)
            j_vap_i_pls = self._c_vap_t_i_pls(i=i) * self.area * (t_i_pls - t_i)

            cap_m = self.cells_is[i].cap_m(state=self.state_n_is)

            # (kg/s)/(J/kg)
            # m_cap: kg/(J/kg)
            SAHEN = (
                cap_m / dt
                + c_vap_liq_wp_i_mns
                + c_vap_liq_wp_i_pls
            )

            # kg/s
            UHEN = (
                cap_m / dt * self.wp_n_is[i]
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
                j_wtr_vent_wp = self.dpv_dmu(i) * PXCOF * v_air_is[i]
                # (kg/s)/K = kg/m3 * Pa/K * (1/Pa) * m3/s
                j_wtr_vent_k = self.dpv_dt(i) * PXCOF * v_air_is[i]

                # RN: 水膜の保持水分量, kg/m2

                # 水膜からの水分移動量, kg/s

                j_dsh_vap_i_pls = alpha_dsh_m * (self.state_n_is[i+1].p_v_sat - self.state_n_is[i].p_v) * self.area * r_wet if RN[i+1] > 0.0 else 0.0

                j_dsh_vap_i_mns = alpha_dsh_m * (self.state_n_is[i-1].p_v_sat - self.state_n_is[i].p_v) * self.area * r_wet if RN[i-1] > 0.0 else 0.0

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
            D5 = self.wjw[i] * self.dx_is[i] * self.area

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

    def get_rn_n_pls(self, wp_is: np.ndarray, oc: OutdoorCondition, p_v_rm: float, dt: float):

        # 評価風速, m/s
        v_mod = self.wsurf.get_v_wind_eva_k(v_wind=oc.v_wind)

        # 風向と壁の法線のなす角度, 度
        angle = self.wsurf.get_wind_angle(wind_direction=oc.wind_direction)

        # 降水量, mm/s
        rf = oc.rainfall / 3600.0

        # 湿気伝達率, kg/(m s Pa)
        alpha_dsh_m = 3.43e-8
        
        rn = np.zeros(shape=self.n_mesh_total, dtype=float)
        wjrain = np.zeros(shape=self.n_mesh_total, dtype=float)

        for i in range(self.n_mesh_total):
            
            if self.rain_leakage_is[i].is_rainpoint:

                # 材料面への浸水量, kg?(m2 s)
                swjrain = self.rain_leakage_is[i].get_confrain(v_mod=v_mod, angle=angle, rf=rf) * rf

                if isinstance(self.cells_is[i], CellOutsideSurface):
                    mns = oc.xod
                elif self.is_air_layer_is[i-1]:
                    rh = get_rh(mu=wp_is[i-1], t=self.t_n_is[i-1])
                    _, vsp = GOFF(t=self.t_n_is[i-1])
                    mns = vsp * rh * 0.01
                else:
                    mns = wp_is[i-1]

                if isinstance(self.cells_is[i], CellInsideSurface):
                    pls = p_v_rm
                elif self.is_air_layer_is[i+1]:
                    rh = get_rh(mu=wp_is[i+1], t=self.t_n_is[i+1])
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
                c = 1 / (self.dx_is[i] / rmdl + r / self.dpv_dmu(i))

                    
                # 当該質点が飽和している前提で計算する。
                # 3.43E-08：湿気伝達率 kg/(m2 s Pa))
                # 水分流の計算, kg/(m2 s)                            
                # 0.3 = 濡れ面率
                if isinstance(self.cells_is[i], CellOutsideSurface) or self.is_air_layer_is[i-1]:
                    p_sv = GOFF(t=self.t_n_is[i])
                    f_mns = alpha_dsh_m * (mns - p_sv) * 0.3
                else:
                    f_mns = c * (mns - wp_is[i])                                
                if isinstance(self.cells_is[i], CellInsideSurface) or self.is_air_layer_is[i+1]:
                    p_sv = GOFF(t=self.t_n_is[i])
                    f_pls = alpha_dsh_m * (pls - p_sv) * 0.3
                else:
                    f_pls = c * (pls - wp_is[i])

                # 水幕の水分量（保持している水の量）, kg/m2
                # D2 はポテンシャル
                # D4 は蒸発量
                # HRN は前のステップの水分量
                
                # 水膜の水分量の計算, kg/m2
                rn = (f_mns + f_pls + swjrain) * dt + self.rn[i]

                # 水幕からの吸水量（表面に水幕があって材料に吸われる分）
                # 水幕が残っている場合は飽和水蒸気圧（水分伝導率で計算）
                if rn < 0.0:
                    rn_n_pls = 0.0
    
                    # RNがゼロの場合は水幕がないので、雨水が直接材料に吸われることになる。
                    w = (swjrain + self.rn[i]) * self.area
                # 水膜がある。
                else:
                    rn_n_pls = rn
                    w = c * (0 - self.wp_n_is) * self.area
                
                rn[i] = rn_n_pls
                wjrain[i] = w
        
        return rn, wjrain

    def get_v_air(self, oc: OutdoorCondition):
        """通気層の換気量を求める。
        """

        v_air = np.zeros(self.n_mesh_total, dtype=float)

        for i in range(self.n_mesh_total):

            if self.is_air_layer_is[i]:

                rho_i = get_rho(t=self.t_n_is)
                v_air[i] = self.alpha_a_ls * (2 / oc.rho * abs(oc.rho - rho_i) * 9.8 * self.height) ** 0.5
        
        return v_air

    def update_t_wd_gen(self):

        for i in range(self.n_mesh_total):

            state = self.state_n_is[i]

            time_s = wood_decay.update_time_s(theta=state.theta, rh=state.rh, time_s=self.time_s_is[i])

            l_stage = wood_decay.update_stage(theta=state.theta, rh=state.rh, time_s=self.time_s_is[i])

            self.time_s_is[i] = time_s

            self.l_stage_is[i] = l_stage

    @classmethod
    def read(cls, ipt_wall: InputWall):

        i = 0

        layers = [
            Layer(num=10, name='サイディング', thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
            Layer(num=2, name='通気層', thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
            Layer(num=6, name='石膏ボード', thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
            Layer(num=5, name='構造用合板1', thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
            Layer(num=3, name='グラスウール1', thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
            Layer(num=5, name='構造用合板1', thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
            Layer(num=6, name='石膏ボード', thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
        ]

        # 面積, m2
        area = ipt_wall.len_vertical * ipt_wall.len_horizontal

        # レイヤー分割数, [L]
        n_div_ls = [layer.n_div for layer in layers]

        # Layerそれぞれにおける室外側と室内側のメッシュ番号, [L], [L]
        outside_end_point_mesh_index_ls, inside_end_point_mesh_index_ls = _get_outside_and_inside_end_point_mesh_index_ls(n_div_ls=n_div_ls)

        # 相当開口面積（αA）, m2
        alpha_a_ls = [layer.alpha * layer.thick * ipt_wall.len_horizontal for layer in layers]

        # ある質点がどのレイヤーに対応するかを保持するリスト, [I]
        lookup_table = create_lookup_table(layers)

        # 質点の総数
        n_mesh_total = sum([layer.n_div for layer in layers])

        # 室外側の端点かどうか, [I]
        is_outside_end_point_is = np.array([outside_end_point_mesh_index_ls[layer_index] == i for (i, layer_index) in enumerate(lookup_table)])

        # 室内側の端点かどうか, [I]
        is_inside_end_point_is = np.array([inside_end_point_mesh_index_ls[layer_index] == i for (i, layer_index) in enumerate(lookup_table)])

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

        # 質点iの密度, kg/m3, [I]
        gma_is = np.array([material_i.rho for material_i in material_is])

        # 質点iの体積, m3, [I]
        # 端点の場合は体積が端点以外の部分の半分になる。（空気層は除く。）
        v_is = np.where(is_outside_end_point_is | is_inside_end_point_is, 0.5, 1.0) * dx_is * area
        
        cells_is = []

        for (i, layer_index) in enumerate(lookup_table):

            # 質点が室外側表面の場合
            if i == 0:
                cell = CellOutsideSurface(
                    x=dx_is[i]/2,
                    material=material_is[i]
                )

            # 質点が室内側表面の場合
            elif i == n_mesh_total - 1:
                cell = CellInsideSurface(
                    x=dx_is[i]/2,
                    material=material_is[i]
                )

            # 質点が空気層の場合
            elif is_air_layer_is[i]:
                cell = CellAirLayer(x=dx_is[i])

            # 質点が室外側の端点の場合
            elif is_outside_end_point_is[i]:

                # 質点が室外側の端点でかつ空気層に面する場合
                if is_air_layer_is[i-1]:
                    cell = CellOutsideEndPointAirLayer(
                        x=dx_is[i]/2,
                        material=material_is[i]
                    )

                else:
                    cell = CellOutsideEndPoint(
                        x=dx_is[i]/2,
                        material=material_is[i],
                        cond_h_o=cond_h_o_is[i],
                        cond_m_o=cond_m_o_is[i]
                    )

            # 質点が室内側の端点の場合
            elif is_inside_end_point_is[i]:

                # 質点が室内側の端点でかつ空気層に面する場合
                if is_air_layer_is[i+1]:
                    cell = CellInsideEndPointAirLayer(
                        x=dx_is[i]/2,
                        material=material_is[i]
                    )

                else:
                    cell = CellInsideEndPoint(
                        x=dx_is[i]/2,
                        material=material_is[i],
                        cond_h_i=cond_h_i_is[i],
                        cond_m_i=cond_m_i_is[i]
                    )

            else:
                cell = CellInterior(
                    x=dx_is[i],
                    material=material_is[i]
                )

            cells_is.append(cell)

        states = [
            State.init(t=INITIAL_WALL_TEMPERATURE + ATP, rh=INITIAL_WALL_RERATIVE_HUMIDITY)
            for _ in range(n_mesh_total)
        ]

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
        is_rainpoint_is = np.full(n_mesh_total, False)

        # 雨水浸入の率
        wall_fall_ratio_is = np.zeros(n_mesh_total)

        # 閾値風速, m/s
        wall_fall_wind_threshold_is = np.zeros(n_mesh_total)

        for nrain in nrains_list:
            is_rainpoint_is[nrain.pos] = True
            wall_fall_ratio_is[nrain.pos] = nrain.ratio
            wall_fall_wind_threshold_is[nrain.pos] = nrain.v

        # 初期温度, K, [I]
        t_init_is = [layers[layer_index].initial_temp + ATP for layer_index in lookup_table]

        rn = np.zeros(n_mesh_total)

        time_s_is = np.zeros(n_mesh_total)

        l_stage_is = np.full(shape=n_mesh_total, fill_value=False)

        m_loss = np.zeros(n_mesh_total)

        wjw = np.zeros(n_mesh_total)

        rain_leakage_is = [
            nrain.RainLeakage(is_rainpoint=is_rainpoint, ratio=wall_fall_ratio, wind_threshold=wall_fall_wind_threshold)
            for (is_rainpoint, wall_fall_ratio, wall_fall_wind_threshold) in zip(is_rainpoint_is, wall_fall_ratio_is, wall_fall_wind_threshold_is)]

        return Wall(
            alpha_a_ls=alpha_a_ls,
            height=height,
            area=area,
            n_mesh_total=n_mesh_total,
            dx_is=dx_is,
            is_air_layer_is=is_air_layer_is,
            gma_is=gma_is,
            material_is=material_is,
            t_n_is=t_init_is,
            wsurf=wsurf,
            rn=rn,
            time_s_is=time_s_is,
            l_stage_is=l_stage_is,
            m_loss=m_loss,
            wjw=wjw,
            state_n_is=states,
            rain_leakage_is=rain_leakage_is,
            cells_is=cells_is
        )


def create_lookup_table(layers: list[Layer]) -> np.ndarray:

    lookup_table = []

    for i, layer in enumerate(layers):
        lookup_table.extend([i] * layer.n_div)
    
    return np.array(lookup_table)


def _get_outside_and_inside_end_point_mesh_index_ls(n_div_ls: list[int]) -> tuple[list[int], list[int]]:
    """各レイヤーについて室外側の質点番号と室内側の質点番号を計算する。

    Args:
        n_div_ls: レイヤーの分割数, [L]

    Returns:
        室外側の質点番号, [L]
        室内側の質点番号, [L]
    """
    
    outside_end_point_mesh_index_ls = [0] * len(n_div_ls)
    inside_end_point_mesh_index_ls = [0] * len(n_div_ls)

    for i, n_div_l in enumerate(n_div_ls):

        if i == 0:
            outside_end_point_mesh_index_ls[i] = 0
        else:
            outside_end_point_mesh_index_ls[i] = inside_end_point_mesh_index_ls[i - 1] + 1

        inside_end_point_mesh_index_ls[i] = outside_end_point_mesh_index_ls[i] + n_div_l - 1
   
    return outside_end_point_mesh_index_ls, inside_end_point_mesh_index_ls


