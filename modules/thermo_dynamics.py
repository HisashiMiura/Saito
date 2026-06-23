import math
from modules.materials import Material


# 絶対温度, K
ATP = 273.15

# 水蒸気の比気体定数, J / (kg K)
RG = 461.5

# 乾き空気の平均分子量(約28.966)に対する水蒸気の分子量(約18.015)の比
RMVA = 0.6217

ROW_CP = 1300.0
"""空気の容積比熱, J / (m3 K)"""

# 水の蒸発潜熱, J / kg
RW = 2.512 * 10**6

CPL = 4200.0
"""水の比熱, J / (kg K)"""

# 水の密度, kg / m3
ROW = 998.0

# 標準大気圧, Pa
P_ATM = 101325.0

# 室外側表面熱伝達率, W/(m2 K)
COND_H_O = 22.4

# 室内側表面熱伝達率, W/(m2 K)
COND_H_I = 9.2

# 通気層表面熱伝達率, W/(m2 K)
COND_H_AIR = 9.2

# 室外側表面湿気伝達率, (kg/s)/(m2 Pa)
COND_M_O = 2.0e-11

# 室内側表面湿気伝達率, (kg/s)/(m2 Pa)
COND_M_I = 3.43e-08

# 通気層表面湿気伝達率, (kg/s)/(m2 Pa)
COND_M_AIR = 3.43e-08

def get_x(p_v: float) -> float:
    """絶対湿度を求める。

    Args:
        p_v: 水蒸気圧, Pa

    Returns:
        絶対湿度, kg/kg(DA)
    """

    return 0.622 * p_v / (P_ATM - p_v)

def get_rho(t: float) -> float:
    """絶対温度tにおける空気の密度を求める。

    Args:
        t: 絶対温度, K
    Returns:
        float: 空気の密度, kg / m3
    
    Notes:
        PV = nRT より、密度ρは以下のように求められる。
        P: 大気圧, Pa (約101325 Pa)
        V: 体積, m3
        n: モル数, mol
        R: 気体定数, J/(mol K) (約8.314 J/(mol K))
        T: 絶対温度, K
        空気の密度は、
        ρ = m / V
        m: 質量, kg
        空気の質量は、
        m = n * M
        M: 空気の平均分子量, kg/mol (約0.028966 kg/mol)
        よって、
        ρ = (n * M) / V
          = PM / RT
        P=101325 Pa, M=0.028966 kg/mol, R=8.314 J/(mol K) を代入すると、
        ρ = (101325 * 0.028966) / (8.314 * t) ≈ 353.0 / t
    """

    return 353.0 / t


def get_wp(rh: float, t: float) -> float:
    """水分化学ポテンシャルを求める。

    Args:
        rh: 相対湿度, %
        t: 絶対温度, K

    Returns:
        水分化学ポテンシャル, J / kg
    """

    return RG * t * math.log(rh * 0.01)


def GOFF(t: float) -> tuple[float, float]:
    """飽和水蒸気圧を求める。
    
    Args:
        TMP: 絶対温度, K
    
    Returns:
        FS: mmHg
        VP: Pa
    """

    ew = - 6096.938 / t + 21.2409642 - t * 2.711193E-2 + t * t * 1.673952E-5 + 2.433502 * math.log(t) 
    
    vp = math.exp(ew)
    
    FS = vp / 133.322 

    return FS, vp


def FUNCX(rh: float, vp: float) -> float:
    """混合比（乾き空気 1 kg あたりの水蒸気重量 [kg]）を求める。

    Args:
        rh: 相対湿度, %
        vp: 飽和水蒸気圧, Pa

    Returns:
        混合比, kg/kg(DA)
    """

    # 133.322 Pa / mmHg
    # 大気圧 101324.72 Pa = 760.0 * 133.322

    return vp * rh * 0.01 * RMVA / (760.0 * 133.322 - vp * rh * 0.01)
    

def get_dpv_dmu(rh: float, t: float) -> float:
    """含水率変化に対する絶対湿度の変化率(DGDU)を計算する。

    Args:
        rh: 相対湿度, %
        tmp: 絶対温度, K
    
    Returns:
        水分化学ポテンシャルに対する水蒸気圧の微分, Pa / (J / kg)
    """

    # --- DGDU の計算 (湿度差による微分) ---
    rh1 = rh + 0.01
    rh2 = rh - 0.01
    
    # 外部定義されている前提の関数
    # get_wp(相対湿度, 温度, パラメータ)
    wp1 = get_wp(rh1, t)
    wp2 = get_wp(rh2, t)
    
    # GOFF関数(飽和水蒸気圧等)から絶対湿度を計算
    # 以前の回答通り、2つの戻り値を想定
    _, vp = GOFF(t)
    x1 = vp * rh1 * 0.01
    x2 = vp * rh2 * 0.01
    
    # ゼロ除算を回避しつつ微分布を計算
    if abs(wp1 - wp2) > 1e-15:
        dgdu = abs((x1 - x2) / (wp1 - wp2))
    else:
        dgdu = 0.0

    return dgdu


def get_dpv_dt(rh: float, t: float) -> float:
    """温度変化に対する絶対湿度の変化率(DGDT)を計算する。

    Args:
        rh: 相対湿度, %
        tmp: 絶対温度, K
    
    Returns:
        絶対温度に対する水蒸気圧の微分, Pa / K

    """

    # --- DGDT の計算 (温度差による微分) ---
    tp1 = t + 0.1
    tp2 = t - 0.1
    
    # 温度を変えて飽和水蒸気圧を取得
    _, vp1 = GOFF(tp1)
    _, vp2 = GOFF(tp2)
    
    x1_t = vp1 * rh * 0.01
    x2_t = vp2 * rh * 0.01
    
    # 温度勾配による絶対湿度の変化率
    # tp1 - tp2 は常に 0.2 なのでゼロ除算の心配はほぼありません
    dgdt = abs((x1_t - x2_t) / (tp1 - tp2))

    return dgdt


def get_rh(mu: float, t: float) -> float:
    """水分化学ポテンシャルから平衡相対湿度 (RH) を計算する。

    Args:
        mu: 水分化学ポテンシャル, J / kg
        t: 絶対温度, K
    
    Returns:
        相対湿度, %

    """
    
    rh = math.exp(mu / RG / t) * 100.0
    
    # 相対湿度は100%を超えることはないので、上限を設定する
    return min(rh, 100.0) 


def DIFF(rh:float, t: float, material: Material):
    """水分伝導率(RML)の計算

    Args:
        rh: 相対湿度, %
        t: 絶対温度, K
        material: 
    Returns: 水分伝導率, (kg/ms) / (J/kg)
    """

    # 含水率, %
    wd = material.get_u(rh=rh) * 100
    
    if 45 < wd < 110:
        d1 = wd * 0.01
        # DW = e^(a + bX + cX^2)
        # 拡散係数, m2/s
        # 含水率差による移動係数
        # 拡散係数は、水分伝導率 (kg/s) / m (kg/kg)　を　水の密度（kg/m3）でわったもの。 
        # m3/s / m kg/kg
        dw = math.exp(-30.91 + 4.2967 * d1 - 0.22017 * (d1**2))
        
        rh1 = rh + 0.005
        rh2 = rh - 0.005
        
        # 微分を差分で近似している計算
        wp1 = get_wp(rh1, t)
        wp2 = get_wp(rh2, t)
        wd1 = material.get_u(rh=rh1) * 100
        wd2 = material.get_u(rh=rh2) * 100
        
        # ゼロ除算のチェック
        if abs(wp1 - wp2) > 1e-12:
            # 拡散係数は、m3/s / m kg/kg
            # これを質量基準に直すために水の密度をかける。
            # 含水率差基準(kg/kg)ではなくて水分化学ポテンシャル差基準(J/kg)になおす。
            # 水の密度(kg/m3) * 拡散係数(m2/s) * 含水率(kg/kg) * 相対湿度差(0.01) / 水分化学ポテンシャル(J/kg)
            # 移動量(kg/s) / (m (J/kg)) 
            # Δwp/Δu = Δwp/Δrh / Δu/Δrh
            # 0.01 は必要かどうか？（要チェック）
            # 998.0: 水の密度, kg/m3
            # dw: 拡散係数, (m3/s) / (m kg/kg) = m2/s
            # wd * 0.01: kg/kg
            # wp: J/kg
            # (kg/ms) / (J/kg)
            rml = 998.0 * dw * abs((wd1 - wd2) * 0.01 / (wp1 - wp2))
        else:
            rml = 0.0
    else:
        rml = 0.0
        
    return rml


def get_dpsi_dmu(mu: float, t: float, gma: float, get_u: callable) -> float:
    """水分化学ポテンシャル変化に対する含水率変化 (DPDU) を計算する。
        (m3/m3)/(J/kg)
        分子：含水率, m3/m3
        分母：水分化学ポテンシャル, J/kg
    
    Args:
        mu: 水分化学ポテンシャル, J/kg
        t: 絶対温度, K
        gma: 材料密度, kg/m3
        get_u: 含水率を計算する関数
    """

    # 差分幅の計算 (d1の1%), J/kg
    dmu = abs(mu * 0.01)
    
    # 微小変化させた含水率, J/kg
    mu1 = mu + dmu
    mu2 = mu - dmu
    
    # 相対湿度, %
    rh1 = get_rh(mu1, t)
    rh2 = get_rh(mu2, t)
    
    # 質量基準の含水率外部定義されている前提の AHGANS, 含水率
    u1 = get_u(rh1)
    u2 = get_u(rh2)
    
    # 体積基準の含水率, m3/m3
    # gma: 材料密度　kg/m3
    # ROW: 水の密度　kg/m3
    psi1 = u1 * gma / ROW
    psi2 = u2 * gma / ROW
    
    # 中央差分による勾配(微分値)の近似
    # 0.5 * (VGT1 - VGT2) / DW
    # dw を2回たしているので2でわっている。
    dpsi_dmu = (psi1 - psi2) / (2 * dmu) 
        
    return dpsi_dmu
