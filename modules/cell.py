from dataclasses import dataclass
from abc import ABC, abstractmethod


from modules.materials import Material
from modules.thermo_dynamics import ROW_CP, ROW, COND_H_I, COND_H_O, COND_H_AIR, COND_M_I, COND_M_O, COND_M_AIR
from modules.state import State


@dataclass
class Cell(ABC):

    # 幅, m
    x: float

    @property
    @abstractmethod
    def cap(self) -> float:
        """熱容量, J/(m2 K)"""
        pass

    @abstractmethod
    def cap_m(self, state: State) -> float:
        """水分移動に伴う容積項, (kg/m2)/(J/kg)"""
        pass

    @property
    @abstractmethod
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        pass

    @property
    @abstractmethod
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        pass

    @property
    @abstractmethod
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        pass

    @property
    @abstractmethod
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        pass

    @abstractmethod
    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        pass

    @abstractmethod
    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        pass


@dataclass
class CellAirLayer(Cell):

    @property
    def cap(self) -> float:
        """熱容量, J/(m2 K)"""
        return self.x * ROW_CP
    
    def cap_m(self, state: State) -> float:
        """水分移動に伴う容積項, (kg/m2)/(J/kg)"""
        return 0.0

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（－）, (m2 K)/W"""
        return 0.0

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return 0.0

    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return 0.0

    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return 0.0

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return float('inf')

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return float('inf')


@dataclass
class CellMaterial(Cell):

    material: Material

    @property
    def cap(self) -> float:
        """熱容量, J/(m2 K)"""
        return self.x * self.material.c_rho

    def cap_m(self, state: State) -> float:
        """水分移動に伴う容積項, (kg/m2)/(J/kg)"""

        # ROW: 水の密度, kg/m3
        # dpsi_dmu: 容積含水率を水分化学ポテンシャルで微分した値, (m3/m3)/(J/kg)
        # x: 幅, m
        return ROW * self.material.get_dpsi_dmu(mu=state.mu, t=state.t) * self.x

    @property
    @abstractmethod
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        pass

    @property
    @abstractmethod
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        pass

    @property
    @abstractmethod
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        pass

    @property
    @abstractmethod
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        pass

    @abstractmethod
    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        pass

    @abstractmethod
    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        pass


@dataclass
class CellOutsideSurface(CellMaterial):

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        return 1 / COND_H_O

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return self.x / self.material.lambda_h
    
    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return 1 / COND_M_O
    
    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return float('inf')

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)


@dataclass
class CellInsideSurface(CellMaterial):

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        return self.x / self.material.lambda_h

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return 1 / COND_H_I

    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m
    
    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return 1 / COND_M_I

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return float('inf')


# 室外側端点（空気層に面しない・外気に面しない）
@dataclass
class CellOutsideEndPoint(CellMaterial):

    # 熱コンダクタンス（＋）, W/(m2 K)
    cond_h_o: float

    # 湿気コンダクタンス（室外側）, (kg/s) / (m2 Pa)
    cond_m_o: float

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        return 1 / self.cond_h_o

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return self.x / self.material.lambda_h

    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return 1 / self.cond_m_o
    
    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return float('inf')

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)


# 室内側端点（空気層に面しない・室内に面しない）
@dataclass
class CellInsideEndPoint(CellMaterial):

    # 熱コンダクタンス（＋）, W/(m2 K)
    cond_h_i: float

    # 湿気コンダクタンス（室内側）, (kg/s) / (m2 Pa)
    cond_m_i: float

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        return self.x / self.material.lambda_h

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return 1 / self.cond_h_i

    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m
    
    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return 1 / self.cond_m_i

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return float('inf')


@dataclass
class CellOutsideEndPointAirLayer(CellMaterial):

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        return 1 / COND_H_AIR

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return self.x / self.material.lambda_h

    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return 1 / COND_M_AIR
    
    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return float('inf')

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)


@dataclass
class CellInsideEndPointAirLayer(CellMaterial):

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        return self.x / self.material.lambda_h

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return 1 / COND_H_AIR

    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m
    
    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return 1 / COND_M_AIR

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return float('inf')


@dataclass
class CellInterior(CellMaterial):

    @property
    def r_h_mns(self) -> float:
        """熱抵抗（ー）, (m2 K)/W"""
        return self.x / self.material.lambda_h / 2

    @property
    def r_h_pls(self) -> float:
        """熱抵抗（＋）, (m2 K)/W"""
        return self.x / self.material.lambda_h / 2

    @property
    def r_m_mns(self) -> float:
        """透湿抵抗（ー）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m / 2

    @property
    def r_m_pls(self) -> float:
        """透湿抵抗（＋）, (m2 Pa)/(kg/s)"""
        return self.x / self.material.lambda_dsh_m / 2

    def r_liq_wp_mns(self, state: State) -> float:
        """水分伝導抵抗（－）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)

    def r_liq_wp_pls(self, state: State) -> float:
        """水分伝導抵抗（＋）, (m2 (J/kg))/(kg/s)"""
        return self.x / self.material.get_lambda_dsh_mu_l(rh=state.rh, t=state.t)
