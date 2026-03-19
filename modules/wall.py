from dataclasses import dataclass
import numpy as np
import math
from scipy.stats import rayleigh

from .sat_temp import get_theta_sat_d_t
from .weather import OutdoorCondition
from config import HOI
from nrain import NRAIN, get_nrains_of_walls
from .thermo_dynamics import ATP, DIFF, get_dgdt, get_dgdu, get_wp, RW


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


class Materials:

    def __init__(self):

        self._K = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
        self._RMDH = [0.0058, 0.0058, 0.045, 0.0824, 0.16, 0.2192, 1.6, 0.0586, 0.2111, 0.963, 0.522, 0.522, 0.021, 0.0339, 0.028, 0.036, 0.0379, 0.1397, 0.1596]
        self._RMDG = [8.69E-12, 8.69E-12, 0.00E+00, 1.58E-10, 1.86E-12, 4.58E-12, 2.93E-11, 4.58E-12, 2.61E-11, 2.29E-11, 4.17E-10, 2.03E-11, 2.79E-11, 5.00E-14, 3.48E-12, 2.61E-12, 1.30E-10, 1.30E-10, 8.94E-12, 1.36E-12]
        self._GCP = [0, 0, 8.40E+02, 1.62E+03, 1.86E+03, 2.84E+03, 9.00E+02, 1.50E+03, 1.31E+03, 8.79E+02, 1.05E+03, 1.05E+03, 8.78E+02, 1.26E+03, 1.26E+03, 8.37E+02, 8.38E+02, 1.26E+03, 1.30E+03]
        self._GMA = [0, 16, 484, 600, 318, 2300, 306, 715, 1095, 1334, 1091, 38.6, 30, 40, 32, 16, 413, 644]
        self._RMDL = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self._GCP = [x * y for (x, y) in zip (self._GCP, self._GMA)]
    
    def GCP(self, i: int):
        return self._GCP[i]
    
    def RMDH(self, i: int):
        """熱伝導率, W/(m K)"""
        return self._RMDH[i]
    
    def RMDG(self, i: int):
        return self._RMDG[i]
    
    def GMA(self, i: int):
        return self._GMA[i]


materials = Materials()


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
    direction: int

    # 隣接室1
    room_sideA: int

    # 隣接室2
    room_sideB: int

    # 長辺
    len_long: float

    # 短辺
    len_short: float

    # 高さ
    height: float

    # ブランク
    blank: float

    # 面積, m2
    area: float

    # 傾斜角, 度
    angle: float

    # アルベド
    albedo: float

    # 日射吸収率
    absorption: float

    # 放射率
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
        
    # 相当開口面積（αA）, m2
    @property 
    def alpha_a(self) -> list[float]:
        return [layer.alpha * layer.thick * self.len_short for layer in self.layers]

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

    def HTMP(self, i: int):
        return self.HTMP[i]

    def set_HTMP(self, i: int, HTMP: float):
        self.state.HTMP = HTMP

    @property
    def t_side_A_n_pls(self):
        return self.t_n_pls[0]
    
    @property
    def t_side_B_n_pls(self):
        return self.t_n_pls[self.n_mesh_total - 1]

    def RMDL(self, i: int):
        """液体伝導率, kg/(m s (J/kg))"""

        if self.rh[i] > 90.0 and self.material_num[i] == 5:
            return DIFF(n_mat=self.material_num[i], rh=self.rh[i], k=self.t_n_pls[i])
        else:
            return 0.0
        
    @property
    def layer_indices(self) -> np.ndarray:
        return self.lookup_table
    
    def get_layer(self, i) -> Layer:
        return self.layers[self.layer_indices[i]]
    
    def is_first_mesh(self, i):
        layer_index = self.layer_indices[i]
        return i == self.first_mesh_indices[layer_index]
    
    def is_last_mesh(self, i):
        layer_index = self.layer_indices[i]
        return i == self.last_mesh_indices[layer_index]

    @property
    def material_num(self):    
        
        return np.vectorize(lambda layer: layer.num)(self.layers[self.layer_indices])
    
    def RMDX(self, i: int) -> float:
        """温度を駆動力とする熱コンダクタンス, W/K"""
        return self.lambda_is[i] * self.area / self.dx[i]
    
    def RMDG(self, i: int) -> float:
        material_num =self.material_num
        if material_num == 5:
            d1 = self.rh[i] * 0.01
            return 1.87E-11 * d1**2.4019 * math.exp(-0.78864 * ( 1 - d1**1.1471))  #*3.45   !  Moisture conductivity (kg/msPa)  ASHRAE
        else:
            return materials.RMDG(material_num)
    
    def ADWGX(self, i: int) -> float:
        return self.RMDG(i) * self.dgdu[i] / self.dx[i] * self.area
    
    def ADWLX(self, i: int) -> float:
        return self.RMDL(i) * self.area / self.dx[i]
    
    def ADWX(self, i: int) -> float:
        return self.ADMGX(i) + self.ADWLX(i)
    
    def ADTGX(self, i: int) -> float:
        """"""
        # RMDG: 湿気伝導率, kg / (m s Pa)
        # dgdt: Pa / K
        # (kg / s) / K
        return self.RMDG(i) * self.dgdt(i) * self.area / self.dx[i]
    
    @property
    def GMA(self) -> np.ndarray:
        return np.vectorize(materials.GMA)(self.material_num)

    @property
    def dx(self) -> np.ndarray:
        """メッシュ間距離, m"""
        
        return np.vectorize(lambda layer: layer.dx)(self.layers[self.layer_indices])
    
    def is_outside_surface(self, i: int) -> bool:

        return i == 0
    
    def is_inside_surface(self, i: int) -> bool:

        return i == self.n_mesh_total - 1

    def dgdu(self, i: int) -> float:
        """水分化学ポテンシャルに対する水蒸気圧の微分, Pa / (J / kg)"""
        return get_dgdu(rh=self.rh[i], t=self.t_n_pls[i])

    def dgdt(self, i: int) -> float:
        """絶対温度に対する水蒸気圧の微分, Pa / K"""
        return get_dgdt(rh=self.rh[i], t=self.t_n_pls[i])
    
    def get_ALF_side_A(self, i: int):
        """Side A側の熱コンダクタンス, W/K"""

        if self.is_first_mesh(i):
            return self.get_layer(i).cond_h_o * self.area
        else:
            return 0.0
    
    def get_ALF_side_B(self, i: int):
        """Side B側の熱コンダクタンス, W/K"""

        if self.is_last_mesh(i):
            return self.get_layer(i).cond_h_i * self.area
        else:
            return 0.0
    
    def alpha_o(self, i: int):
        """熱コンダクタンス（室外側）, W/K"""

        if self.is_outside_end_point_is[i]:
            return self.get_layer(i).cond_h_o
        else:
            return 0.0
    
    def alpha_i(self, i: int):
        """熱コンダクタンス（室内側）, W/K"""

        if self.is_inside_end_point_is[i]:
            return self.get_layer(i).cond_h_i
        else:
            return 0.0

    def alpha_dsh_m_o(self, i: int):
        """湿気コンダクタンス（室外側）, kg/(m2 s Pa)"""
        if self.is_outside_end_point_is[i]:
            return self.get_layer(i).cond_m_o
        else:
            return 0.0
    
    def alpha_dsh_m_i(self, i: int):
        """湿気コンダクタンス（室内側）, kg/(m2 s Pa)"""
        if self.is_inside_end_point_is[i]:
            return self.get_layer(i).cond_m_i
        else:
            return 0.0

    def get_theta_surf_sideB(self, oc: OutdoorCondition, theta_r: float) -> float:
        """室外側の表面温度を求める。

        Args:
            oc: 外気条件
            theta_r: 室内温度, ℃

        Returns:
            室外側表面温度, ℃
        """

        if self.direction == 5:
            return theta_r * 0.3 + oc.t * 0.7
        else:
            return self.get_k_sat(oc=oc)
    
    def get_t_surf_out(self, oc: OutdoorCondition, theta_r: float) -> float:
        """室外側の表面温度を求める。

        Args:
            oc: 外気条件
            theta_r: 室内温度, ℃

        Returns:
            室外側表面温度, K
        """

        return self.get_theta_surf_sideB(oc=oc, theta_r=theta_r) + ATP

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
            c_o = (self.lambda_is[i-1] + RW * self.RMDG(i-1) * self.dgdt(i-1)) / self.dx[i-1]
            c = (self.lambda_is[i] + RW * self.RMDG(i) * self.dgdt(i)) / self.dx[i]
            return 1 / (1 / c_o + 1 / c) * self.area
    
    def get_c_t_i_pls(self, i: int) -> float:
        """温度差を駆動力とする熱移動に関する係数（室内側）, W/K"""

        # 各レイヤの室内側の端点の場合
        if self.is_inside_end_point_is[i]:

            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            return (self.get_layer(i).cond_h_i + RW * self.get_layer(i).cond_m_i * self.dgdt(i)) * self.area
        
        else:

            c_i = (self.lambda_is[i+1] + RW * self.RMDG(i+1) * self.dgdt(i+1)) / self.dx[i+1]
            c = (self.lambda_is[i] + RW * self.RMDG(i) * self.dgdt(i)) / self.dx[i]
            return 1 / (1 / c_i + 1 / c) * self.area
    
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
            c_o = RW * self.RMDG(i-1) * self.dgdu[i-1] / self.dx[i-1] * self.area
            c = RW * self.RMDG(i) * self.dgdu[i] / self.dx[i] * self.area
            return 1 / (1 / c_o + 1 / c) * self.area
    
    def get_c_wp_i_pls(self, i: int) -> float:
        """水分ポテンシャル差を駆動力とする水蒸気移動に伴う熱移動に関する係数（室内側）, W/(J/kg)"""

        # 各レイヤの室内側の端点の場合
        if self.is_inside_end_point_is[i]:
            # TODO: 本来であれば、隣との接触抵抗は無限大の値の入力を必要としない抵抗としてもっておき、隣側の抵抗も加算した上での逆数とすべきではないか。
            # 水の蒸発潜熱 J/kg * 湿気コンダクタンス, kg/(m2 s Pa) * Pa/(J/kg) * m2 = W/(J/kg) 
            return RW * self.get_layer(i).cond_m_i * self.dgdu(i) * self.area
        else:
            c_i = RW * self.RMDG(i+1) * self.dgdu[i+1] / self.dx[i+1] * self.area
            c = RW * self.RMDG(i) * self.dgdu[i] / self.dx[i] * self.area
            return 1 / (1 / c_i + 1 / c) * self.area

    def get_k_sat(self, oc: OutdoorCondition):

        sin_h_d_t = oc.sin_h
        cos_h_d_t = oc.cos_h
        sin_a_d_t = oc.sin_a
        cos_a_d_t = oc.cos_a
        i_dn_d_t = oc.i_dn
        i_sky_d_t = oc.i_sky
        theta_o_d_t = oc.t
        r_n_d_t = oc.r_n

        alpha_k = self.get_alpha_w()
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

    def get_alpha_w(self):

        if self.direction == 5:

            return None

        else:
            
            return np.radians((self.direction - 1) * 90.0 + HOI)
    
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
                # 0.2 = 
                # ASHRAE の係数
                d2 = np.maximum(
                    nrain.ratio * x * np.cos(np.radians(wind_direction - (180.0 + 90.0 * (self.direction - 1)))) * 1.5 * 0.2,
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
        swjrains = 0.01 * oc.rainfall * confrains / 3600.0

        for (nrain, swjrain) in (self.nrains, swjrains):
            if nrain.pos - 1 == i:
                return swjrain
        
        return 0.0

    @classmethod
    def read(cls, d: dict, walltypes: list[WallType], i: int):

        len_long = d['len_long']

        len_short = d['len_short']

        blank = d['blank']

        # 面積, m2
        area = len_long * len_short - blank

        # 壁タイプのID
        wall_type_id = int(d['walltypes'])

        # （IDで指定された）壁タイプ
        wall_type: WallType = walltypes[wall_type_id]

        # レイヤー
        layers = np.array(wall_type.layers)

        # Layerそれぞれにおける室外側と室内側のメッシュ番号, [L], [L]
        outside_end_point_mesh_indices, inside_end_point_mesh_indices = _get_first_and_last_mesh_indices(layers=wall_type.layers)
            
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

        # 質点iの熱容量, J/(m3 K), [I]
        gcp_is = np.array([materials.GCP(wall_type.layers[layer_index].num) for layer_index in lookup_table])

        # 質点iの熱伝導率, W/(m K), [I]
        lambda_is = np.array([materials.RMDH(wall_type.layers[layer_index].num) for layer_index in lookup_table])

        # 質点iの体積, m3, [I]
        # 端点の場合は体積が端点以外の部分の半分になる。（空気層は除く。）
        v_is = np.where(is_outside_end_point_is | is_inside_end_point_is, 0.5, 1.0) * dx_is * area
        
        state = WallState.init(layers, lookup_table)

        # 方位, 0, 1, 2, 3, 4
        direction = d['direction']

        # 傾斜角
        angle = d['angle']

        nrains: list[NRAIN] = get_nrains_of_walls(wall_index=i)
        
        return Wall(
            kwtype=d['kwtype'],
            direction=direction,
            room_sideA=d['room_sideA'],
            room_sideB=d['room_sideB'],
            len_long=len_long,
            len_short=len_short,
            height=d['height'],
            blank=blank,
            area=area,
            angle=angle,
            albedo=d['albedo'],
            absorption=d['absorption'],
            emissivity=d['emissivity'],
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
            lambda_is=lambda_is
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

    return [
        {
            'kwtype': 1,
            'direction': 4,
            'room_sideA': 0,
            'room_sideB': 1,
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'blank': 0,
            'angle': 90.0,
            'albedo': 0.3,
            'absorption': 0.9,
            'emissivity': 0.9,
            'walltypes': 1,
            'eva_height': 4.2,
        },
        {
            'kwtype': 2,
            'direction': 4,
            'room_sideA': 0,
            'room_sideB': 1,
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'blank': 0,
            'angle': 90.0,
            'albedo': 0.3,
            'absorption': 0.9,
            'emissivity': 0.9,
            'walltypes': 2,
            'eva_height': 11.2,
        },
        {
            'kwtype': 2,
            'direction': 4,
            'room_sideA': 0,
            'room_sideB': 1,
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'blank': 0,
            'angle': 90.0,
            'albedo': 0.3,
            'absorption': 0.9,
            'emissivity': 0.9,
            'walltypes': 2,
            'eva_height': 18.2,
        },
        {
            'kwtype': 2,
            'direction': 4,
            'room_sideA': 0,
            'room_sideB': 1,
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'blank': 0,
            'angle': 90.0,
            'albedo': 0.3,
            'absorption': 0.9,
            'emissivity': 0.9,
            'walltypes': 2,
            'eva_height': 25.2,
        },
        {
            'kwtype': 2,
            'direction': 4,
            'room_sideA': 0,
            'room_sideB': 1,
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'blank': 0,
            'angle': 90.0,
            'albedo': 0.3,
            'absorption': 0.9,
            'emissivity': 0.9,
            'walltypes': 2,
            'eva_height': 32.2,
        },
        {
            'kwtype': 2,
            'direction': 4,
            'room_sideA': 0,
            'room_sideB': 1,
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'blank': 0,
            'angle': 90.0,
            'albedo': 0.3,
            'absorption': 0.9,
            'emissivity': 0.9,
            'walltypes': 2,
            'eva_height': 39.2,
        },
        {
            'kwtype': 1,
            'direction': 4,
            'room_sideA': 0,
            'room_sideB': 1,
            'len_long': 7.0,
            'len_short': 0.42,
            'height': 7.0,
            'blank': 0,
            'angle': 90.0,
            'albedo': 0.3,
            'absorption': 0.9,
            'emissivity': 0.9,
            'walltypes': 2,
            'eva_height': 46.2,
        },
    ]


def get_wall_types():

    layers1 = [
        Layer(num=10, name='サイディング', initial_temp=15.0, initial_humidity=50.0, thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
        Layer(num=2, name='通気層', initial_temp=15.0, initial_humidity=50.0, thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=3, name='GW', initial_temp=10.0, initial_humidity=50.0, thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
        Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
    ]

    layers2 = [
        Layer(num=10, name='サイディング', initial_temp=15.0, initial_humidity=50.0, thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
        Layer(num=2, name='通気層', initial_temp=15.0, initial_humidity=50.0, thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=3, name='GW', initial_temp=10.0, initial_humidity=50.0, thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
        Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
        Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
    ]

    return [
        WallType(layers=layers1),
        WallType(layers=layers2)
    ]


