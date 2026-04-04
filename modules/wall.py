from dataclasses import dataclass
import numpy as np
import math
from scipy.stats import rayleigh
from numpy.typing import NDArray 

from .sat_temp import get_theta_sat_d_t
from .weather import OutdoorCondition
from config import HOI
from nrain import NRAIN, get_nrains_of_walls
from .thermo_dynamics import ATP, DIFF, get_dgdt, get_dgdu, get_wp, RW, CALDPDU, ROW, CPL
from .direction import Direction
from .materials import Materials, Material

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

    # 傾斜角, 度
    angle: float

    # アルベド
    albedo: float

    # 日射吸収率
    absorption: float

    # 長波長放射率
    emissivity: float

    # layers
    layers: np.ndarray

    # first mesh indices
    first_mesh_indices: list[int]

    # last mesh indices
    last_mesh_indices: list[int]

    # メッシュ数
    n_mesh_total: int

    # 状態量
    state: WallState

    # 評価高さ
    eva_height: float

    # 浸水ポイント
    nrains: list[NRAIN]

    # メッシュ番号に対するレイヤーインデックス
    lookup_table: np.ndarray

    # 質点iの質点間距離, m, [I]
    dx_is: np.ndarray

    # 室外側の端点かどうか, [I]
    is_outside_end_point_is: np.ndarray

    # 室内側の端点かどうか, [I]
    is_inside_end_point_is: np.ndarray

    # 質点iの体積, m3, [I]
    v_is: np.ndarray

    # 熱容量, J/(m3 K), [I]
    gcp_is: np.ndarray

    # 熱伝導率, W/(m K), [I]
    lambda_is: np.ndarray

    # 材料, [I]
    materials_is: list[Material]

    # ステップnの絶対温度, K, [I]
    t_n_is: np.ndarray

    # ステップnの水分化学ポテンシャル, J/kg, [I]
    wp_n_is: np.ndarray

    @property    
    def rh(self):
        """相対湿度, %"""
        return self.state.rh
    
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

        if self.rh[i] > 90.0 and self.materials_is[i].id == 5:
            return DIFF(n_mat=self.materials_is[i].id, rh=self.rh[i], k=self.t_n_pls[i])
        else:
            return 0.0

    def get_layer(self, i) -> Layer:
        return self.layers[self.lookup_table[i]]
    
    def RMDG(self, i: int) -> float:
        """湿気伝導率, kg /(m s Pa) """
        if self.materials_is[i].id == 5:
            d1 = self.rh[i] * 0.01
            return 1.87E-11 * d1**2.4019 * math.exp(-0.78864 * ( 1 - d1**1.1471))  #*3.45   !  Moisture conductivity (kg/msPa)  ASHRAE
        else:
            return self.materials_is[i].rmdd
    
    def ADWLX(self, i: int) -> float:
        """水分化学ポテンシャル駆動による液水移動量（コンダクタンス）, kg/s / (J/kg)"""
        # RMDL: 水分伝導率, kg/(m s (J/kg))
        # kg/(m s (J/kg)) * m2 / m = kg/s / (J/kg)
        return self.RMDL(i) * self.area / self.dx_is[i]
    
    def ADTGX(self, i: int) -> float:
        """"""
        # RMDG: 湿気伝導率, kg / (m s Pa)
        # dgdt: Pa / K
        # (kg / s) / K
        return self.RMDG(i) * self.dgdt(i) * self.area / self.dx_is[i]
    
    def ADWGX(self, i: int) -> float:
        """水分化学ポテンシャル駆動による水蒸気移動量（コンダクタンス）, kg/s / (J/kg)"""
        # RMDG：湿気伝導率, kg /(m s Pa)
        # kg / (m s Pa) * Pa / (J/kg) * m2 / m = kg/s / (J/kg)
        return self.RMDG(i) * self.dgdu(i) * self.area / self.dx_is[i]

    def ADWX(self, i: int) -> float:
        # kg/s / (J/kg)
        return self.ADWGX(i) + self.ADWLX(i)
    
    def DWX(self, i: int) -> float:
        # 計算した値を平均する方が良い。逆数同士足した値の逆数がベター
        # 本来であれば、逆数の和の逆数にすべきだが、片方がゼロになる場合はゼロ割の可能性があるので注意が必要。
        return (self.ADWX(i) + self.ADWX(i - 1)) / 2

    def DPDU(self, i: int, wpt: float, tp: float):
        """(J/kg)/K"""
        # 1/(J/kg)
        return CALDPDU(wpt=wpt, tp=tp, gma=self.GMA[i], ml0=self.get_layer(i).num, row=ROW)
    
    def ADTLX(self, i: int, wpt: float, tp: float):
        return self.RMDL(i) * self.DPDU(i=i, wpt=wpt, tp=tp) * self.area / self.dx_is[i]
    
    def DTLX(self, i: int, wpts: list[float], tps: list[float]):
        return (self.ADTXL(i=i, wpt=wpts[i], tp=tps[i]) + self.ADTXL(i=i-1, wpt=wpts[i-1], tp=tps[i-1]) )*0.5

    def GMA(self, i: int) -> float:
        """密度, kg/m3"""
        return self.materials_is[i].GMA
    
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
    
    def get_theta_surf_out(self, oc: OutdoorCondition, theta_r: float) -> float:
        """室外側の表面温度を求める。

        Args:
            oc: 外気条件
            theta_r: 室内温度, deg.C

        Returns:
            室外側表面温度, deg.C
        """

        if self.direction == Direction.BOTTOM:
            return theta_r * 0.3 + oc.t * 0.7
        else:
            return self.get_k_sat(oc=oc)
    
    def get_t_surf_out(self, oc: OutdoorCondition, theta_r: float) -> float:
        """室外側の表面温度を求める。

        Args:
            oc: 外気条件
            theta_r: 室内温度, deg.C

        Returns:
            室外側表面温度, K
        """

        return self.get_theta_surf_out(oc=oc, theta_r=theta_r) + ATP

    def cap(self, i: int) -> float:
        """熱容量, J/K"""

        # 質点iが空気層の場合
        if self.get_layer(i).num == 2:
            #1300は空気の体積熱容量, J/(kg K)
            return 1300.0 * self.dx_is[i] * self.area
        else:
            return self.v_is[i] * self.gcp_is[i]

    def m_cap(self, i: int, wp_is: np.ndarray) -> float:
        """水分移動に伴う熱容量, kg/(J/kg)"""

        # 質点iが空気層の場合
        if self.get_layer(i).num == 2:
            # 空気層の場合の熱容量はなし？
            return 0.0
        else:
            # kg/m3 * 1/(J/kg) * m3 = kg/(J/kg)
            # ここの式の単位変換はあっているのか？
            return ROW * self.DPDU(i=i, wpt=wp_is[i], tp=self.t_n_is[i]) * self.v_is[i]

    def get_c_t_i_mns(self, i: int) -> float:
        """温度差を駆動力とする熱移動に関する係数（室外側）, W/K"""
    
        # 各レイヤの室外側の端点の場合
        if self.is_outside_end_point_is[i]:
                                    
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            return (self.get_layer(i).cond_h_o + RW * self.get_layer(i).cond_m_o * self.dgdt(i)) * self.area

        else:
            # ( 熱伝導率 W/(m K) + 水の蒸発潜熱 J/kg * 温度勾配に対する水分伝導率 (kg/(m s K)) ) / 質点間距離 m * 面積 m2
            # 水の蒸発潜熱 2.512 * 10^6 J/kg
            # 温度勾配に対する水分伝導率 lambda_(T,g)
            # RMDG：湿気伝導率, kg / (m s Pa)
            return self.get_series_combination(
                x=(self.lambda_is[i-1] + RW * self.RMDG(i-1) * self.dgdt(i-1)) / self.dx_is[i-1],
                y=(self.lambda_is[i] + RW * self.RMDG(i) * self.dgdt(i)) / self.dx_is[i]
            ) * self.area
    
    def get_c_t_i_pls(self, i: int) -> float:
        """温度差を駆動力とする熱移動に関する係数（室内側）, W/K"""

        # 各レイヤの室内側の端点の場合
        if self.is_inside_end_point_is[i]:

            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            return (self.get_layer(i).cond_h_i + RW * self.get_layer(i).cond_m_i * self.dgdt(i)) * self.area
        
        else:
            return self.get_series_combination(
                x=(self.lambda_is[i+1] + RW * self.RMDG(i+1) * self.dgdt(i+1)) / self.dx_is[i+1],
                y=(self.lambda_is[i] + RW * self.RMDG(i) * self.dgdt(i)) / self.dx_is[i]
            ) * self.area
        
    def get_c_wp_i_mns(self, i: int) -> float:
        """水分ポテンシャル差を駆動力とする水蒸気移動に伴う熱移動に関する係数（室外側）, W/(J/kg)"""

        # 各レイヤの室外側の端点の場合
        if self.is_outside_end_point_is[i]:
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            # 水の蒸発潜熱 J/kg * 湿気コンダクタンス, kg/(m2 s Pa) * Pa/(J/kg) * m2 = W/(J/kg) 
            return RW * self.get_layer(i).cond_m_o * self.dgdu(i) * self.area
        else:
            # RW: 蒸発潜熱, J/kg    2.512 * 10^6
            # DWGX: 蒸気成分の伝導率, kg/s / (J/kg)
            # J/kg * (kg/s) / (J/kg) = W / (J/kg)
            # RMDG: 湿気伝導率, kg/(m s Pa)
            # J/kg * kg/(m s Pa) * Pa/(J/kg) / m * m2
            return self.get_series_combination(
                x=RW * self.RMDG(i-1) * self.dgdu[i-1] / self.dx_is[i-1],
                y=RW * self.RMDG(i) * self.dgdu[i] / self.dx_is[i]
            ) * self.area
    
    def get_c_wp_i_pls(self, i: int) -> float:
        """水分ポテンシャル差を駆動力とする水蒸気移動に伴う熱移動に関する係数（室内側）, W/(J/kg)"""

        # 各レイヤの室内側の端点の場合
        if self.is_inside_end_point_is[i]:
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            # 水の蒸発潜熱 J/kg * 湿気コンダクタンス, kg/(m2 s Pa) * Pa/(J/kg) * m2 = W/(J/kg) 
            return RW * self.get_layer(i).cond_m_i * self.dgdu(i) * self.area
        else:
            return self.get_series_combination(
                x=RW * self.RMDG(i+1) * self.dgdu[i+1] / self.dx_is[i+1],
                y=RW * self.RMDG(i) * self.dgdu[i] / self.dx_is[i]
            ) * self.area

    def get_c_liquid_i_mns(self, i: int, wp_i_mns: float, wp_i: float, t_i_mns: float, t_i: float):
        """液水移動による熱の移動"""

        if self.is_outside_end_point_is[i]:
            return 0.0
        else:
            # 水分化学ポテンシャル駆動の液水移動の係数, (kg/s)/(m2 J/kg)
            dwlx = self.get_series_combination(
                x=self.RMDL(i-1) / self.dx_is[i-1],
                y=self.RMDL(i) / self.dx_is[i]
            )

            # 温度差駆動の液水移動の係数, (kg/s)/(m2 K) = kg/(m s (J/kg)) * (J/kg)/K / m
            dtlx = self.get_series_combination(
                x=self.RMDL(i-1) * self.DPDU(i=i-1, wpt=wp_i_mns, tp=t_i_mns) / self.dx_is[i-1],
                y=self.RMDL(i) * self.DPDU(i=i, wpt=wp_i, tp=t_i) / self.dx_is[i]
            )

            # CPL: 水の比熱 = 4200.0 J/(kg K)
            return CPL * (dwlx * (wp_i_mns - wp_i) + dtlx * (t_i_mns - t_i)) * self.area

    def get_c_liquid_i_pls(self, i: int, wp_i_pls: float, wp_i: float, t_i_pls: float, t_i: float):
        """液水移動による熱の移動"""

        if self.is_inside_end_point_is[i]:
            return 0.0
        else:
            # 水分化学ポテンシャル駆動の液水移動の係数, (kg/s)/(m2 J/kg)
            dwlx = self.get_series_combination(
                x=self.RMDL(i+1) / self.dx_is[i+1],
                y=self.RMDL(i) / self.dx_is[i]
            )

            # 温度差駆動の液水移動の係数, (kg/s)/(m2 K) = kg/(m s (J/kg)) * (J/kg)/K / m
            dtlx = self.get_series_combination(
                x=self.RMDL(i+1) * self.DPDU(i=i+1, wpt=wp_i_pls, tp=t_i_pls) / self.dx_is[i+1],
                y=self.RMDL(i) * self.DPDU(i=i, wpt=wp_i, tp=t_i) / self.dx_is[i]
            )

            # CPL: 水の比熱 = 4200.0 J/(kg K)
            return CPL * (dwlx * (wp_i_pls - wp_i) + dtlx * (t_i_pls - t_i)) * self.area

    def get_t_n_pls(self, t_is: np.ndarray, dt: float, oc: OutdoorCondition, theta_r_n: float, wp_r_n: float, QQ: float, t_upstream: float):

        t_is_next = np.zeros(t_is, dtype=float)

        for i in range(self.n_mesh_total):

            t_i = t_is(i)

            # 室外側の質点の温度
            # 室外側の質点の水分ポテンシャル

            # 質点が室外側表面の場合
            if self.is_outside_surface(i=i):
                t_i_mns = self.get_t_surf_out(oc=oc, theta_r=theta_r_n)
                wp_i_mns = oc.wp
            else:
                t_i_mns = t_is[i - 1]
                wp_i_mns = self.wp_n_is[i - 1]

            wp_i = self.wp_n_is[i]

            # 室内側の質点の温度
            # 質点が室内側表面の場合
            if self.is_inside_surface(i=i):
                t_i_pls = theta_r_n
                wp_i_pls = wp_r_n
            else:
                t_i_pls = t_is[i + 1]
                wp_i_pls = self.wp_n_is[i + 1]


            # 温度差を駆動力とする熱移動に関する係数（室外側）, W/K
            c_t_i_mns = self.get_c_t_i_mns(i)

            # 温度差を駆動力とする熱移動に関する係数（室内側）, W/K
            c_t_i_pls = self.get_c_t_i_pls(i)

            # 水分ポテンシャル差を駆動力とする水蒸気移動に伴う熱移動に関する係数（室外側）, W/(J/kg)
            c_wp_i_mns = self.get_c_wp_i_mns(i)

            # 水分ポテンシャル差を駆動力とする水蒸気移動に伴う熱移動に関する係数（室内側）, W/(J/kg)
            c_wp_i_pls = self.get_c_wp_i_pls(i)

            # 液水移動に伴う熱移動（室外側）, W/K
            c_liquid_i_mns = self.get_c_liquid_i_mns(i=i, wp_i_mns=wp_i_mns, wp_i=wp_i, t_i_mns=t_i_mns, t_i=t_i)

            # 液水移動に伴う熱移動（室内側）, W/K
            c_liquid_i_pls = self.get_c_liquid_i_pls(i=i, wp_i_pls=wp_i_pls, wp_i=wp_i, t_i_pls=t_i_pls, t_i=t_i)

            cap = self.cap(i) / dt

            # W
            UHEN = (
                cap * self.t_n_is[i]
                + c_t_i_mns * t_i_mns
                + c_t_i_pls * t_i_pls
                + c_wp_i_mns * (wp_i_mns - wp_i)
                + c_wp_i_pls * (wp_i_pls - wp_i)
                + c_liquid_i_mns * t_i_mns
                + c_liquid_i_pls * t_i_pls
            )

            # W/K
            SAHEN = (
                cap
                + c_t_i_mns
                + c_t_i_pls
                + c_liquid_i_mns
                + c_liquid_i_pls
            )

            if self.get_layer(i).num == 2:
                # 空気層の場合に移流分を考慮する。
                # 空気の容積比熱, J/(m3 K)                             
                UHEN =+ 1300.0 * QQ * t_upstream
                SAHEN =+ 1300.0 * QQ
            
            t_is_next[i] = UHEN / SAHEN

        return t_is_next
          




    def get_series_combination(x: float, y: float) -> float:

        if x + y != 0.0:
            return x * y / (x + y)
        else:
            return 0.0


    def get_k_sat(self, oc: OutdoorCondition):

        sin_h_d_t = oc.sin_h
        cos_h_d_t = oc.cos_h
        sin_a_d_t = oc.sin_a
        cos_a_d_t = oc.cos_a
        i_dn_d_t = oc.i_dn
        i_sky_d_t = oc.i_sky
        theta_o_d_t = oc.t
        r_n_d_t = oc.r_n

        drct = self.direction

        if drct in [Direction.TOP, Direction.BOTTOM]:
            alpha_k = None
        else:
            alpha_k = np.radians(drct.alpha + HOI)
            
        beta_k = np.radians(self.angle)

        theta_sat_d_t_k, net_gain, d_gain, sky_gain = get_theta_sat_d_t(
            sin_h_d_t=sin_h_d_t,
            cos_h_d_t=cos_h_d_t,
            sin_a_d_t=sin_a_d_t,
            cos_a_d_t=cos_a_d_t,
            alpha_k=alpha_k,
            beta_k=beta_k,
            i_dn_d_t=i_dn_d_t,
            i_sky_d_t=i_sky_d_t,
            rho_g=self.albedo,
            theta_o_d_t=theta_o_d_t,
            alpha_s_k=self.absorption,
            epsilon_k=self.emissivity,
            r_n_d_t=r_n_d_t,
            h_o_k=22.4
        )

        return theta_sat_d_t_k, net_gain, d_gain, sky_gain

    
    def get_v_wind_eva_k(self, v_wind: float):
        """評価高さにおける風速を求める。

        Args:
            v_wind: 風速, m/s

        Returns:
            評価高さにおける風速, m/s
        """

        # 基準風速は6.5m高さ
        return v_wind * (self.eva_height / 6.5)**0.25

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
        """_summary_

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
    def read(cls, d: dict, walltypes: list[WallType], i: int):

        # 長辺, m
        len_long = d['len_long']

        # 長辺, m
        len_short = d['len_short']

        # 面積, m2
        area = len_long * len_short

        # 壁タイプのID
        wall_type_id = int(d['walltypes'])

        # （IDで指定された）壁タイプ
        wall_type: WallType = walltypes[wall_type_id]

        # レイヤー
        layers = np.array(wall_type.layers)

        # Layerそれぞれにおける室外側と室内側のメッシュ番号, [L], [L]
        outside_end_point_mesh_indices, inside_end_point_mesh_indices = _get_first_and_last_mesh_indices(layers=wall_type.layers)

        # 相当開口面積（αA）, m2
        alpha_a_ls = [layer.alpha * layer.thick * len_short for layer in wall_type.layers]
            
        # ある質点がどのレイヤーに対応するかを保持するリスト, [I]
        lookup_table = create_lookup_table(wall_type.layers)

        # 質点の総数
        n_mesh_total = sum([layer.n_div for layer in layers])

        # 室外側の端点かどうか, [I]
        is_outside_end_point_is = np.array([outside_end_point_mesh_indices[layer_index] == i for (i, layer_index) in enumerate(lookup_table)])

        # 室内側の端点かどうか, [I]
        is_inside_end_point_is = np.array([inside_end_point_mesh_indices[layer_index] == i for (i, layer_index) in enumerate(lookup_table)])

        # 質点iの質点間距離, m, [I]
        # レイヤー表面の質点において、SideA側の質点の場合はSideB側のみが、SideB側の質点の場合はSideA側のみが定義される。
        dx_is = np.array([wall_type.layers[layer_index].dx for layer_index in lookup_table])

        # 材料, Material クラス, [I]
        ms = Materials()
        material_is = [ms.get_material(name=wall_type.layers[layer_index].name) for layer_index in lookup_table]

        # 質点iの熱容量, J/(m3 K), [I]
        gcp_is = np.array([material_i.cgg * material_i.gma for material_i in material_is])

        # 質点iの熱伝導率, W/(m K), [I]
        lambda_is = np.array([material_i.rmd for material_i in material_is])

        # 質点iの体積, m3, [I]
        # 端点の場合は体積が端点以外の部分の半分になる。（空気層は除く。）
        v_is = np.where(is_outside_end_point_is | is_inside_end_point_is, 0.5, 1.0) * dx_is * area
        
        state = WallState.init(layers, lookup_table)

        # 方位, 0, 1, 2, 3, 4
        direction = Direction(d['direction'])

        # 傾斜角
        if 'angle' not in d:
            if direction == Direction.TOP:
                angle = 0.0
            elif direction == Direction.BOTTOM:
                angle = 180.0
            else:
                angle = 90.0
        else:
            angle = d['angle']

        # アルベド
        albedo = d.get('albedo', 0.1)

        # 室外側表面の日射吸収率
        absorption = d.get('absoption', 0.9)

        # 室外側表面の長波長放射率
        emissivity = d.get('emissivity', 0.9)

        nrains: list[NRAIN] = get_nrains_of_walls(wall_index=i)

        # 初期温度, K, [I]
        t_init_is = [wall_type.layers[layer_index].initial_temp + ATP for layer_index in lookup_table]

        return Wall(
            kwtype=d['kwtype'],
            direction=direction,
            alpha_a_ls=alpha_a_ls,
            height=d['height'],
            area=area,
            angle=angle,
            albedo=albedo,
            absorption=absorption,
            emissivity=emissivity,
            layers=walltypes[d['walltypes']].layers,
            first_mesh_indices=outside_end_point_mesh_indices,
            last_mesh_indices=inside_end_point_mesh_indices,
            lookup_table=lookup_table,
            n_mesh_total=n_mesh_total,
            state=state,
            eva_height=d['eva_height'],
            nrains=nrains,
            dx_is=dx_is,
            is_outside_end_point_is=is_outside_end_point_is,
            is_inside_end_point_is=is_inside_end_point_is,
            v_is=v_is,
            gcp_is=gcp_is,
            lambda_is=lambda_is,
            material_is=material_is,
            t_n_is=t_init_is
        )

    @classmethod
    def read_default(cls):

        ds = get_ds()

        walltypes = get_wall_types()

        walls = [Wall.read(d=d, walltypes=walltypes, i=i) for i, d in enumerate(ds)]

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


def get_wall_types():

    layers1 = [
        Layer(num=10, name='サイディング', initial_temp=15.0, initial_humidity=50.0, thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
        Layer(num=2, name='通気層', initial_temp=15.0, initial_humidity=50.0, thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=5, name='構造用合板1', initial_temp=10.0, initial_humidity=80.0, thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=3, name='グラスウール1', initial_temp=10.0, initial_humidity=50.0, thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
        Layer(num=5, name='構造用合板1', initial_temp=10.0, initial_humidity=80.0, thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
    ]

    layers2 = [
        Layer(num=10, name='サイディング', initial_temp=15.0, initial_humidity=50.0, thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
        Layer(num=2, name='通気層', initial_temp=15.0, initial_humidity=50.0, thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=5, name='構造用合板1', initial_temp=10.0, initial_humidity=80.0, thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=3, name='グラスウール1', initial_temp=10.0, initial_humidity=50.0, thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
        Layer(num=5, name='構造用合板1', initial_temp=10.0, initial_humidity=80.0, thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
    ]

    return [
        WallType(layers=layers1),
        WallType(layers=layers2)
    ]


