import math


# 絶対温度, K
ATP = 273.15

# 水蒸気の比気体定数, J / (kg K)
RG = 461.5

# 乾き空気の平均分子量(約28.966)に対する水蒸気の分子量(約18.015)の比
RMVA = 0.6217

# 空気の容積比熱, J / (m3 K)
ROW_CP = 1300.0

# 水の蒸発潜熱, J / kg
RW = 2.512 * 10**6

# 水の比熱, J / (kg K)
CPL = 4200.0

# 水の密度, kg / m3
ROW = 998.0


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
    

def get_dgdu(rh: float, t: float) -> float:
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


def get_dgdt(rh: float, t: float) -> float:
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


def WPTRE(wp: float, t: float) -> float:
    """水分化学ポテンシャルから平衡相対湿度 (RH) を計算する。

    Args:
        wp: 水分化学ポテンシャル, J / kg
        t: 絶対温度, K
    
    Returns:
        相対湿度, %

    """
    
    rh = math.exp(wp / RG / t) * 100.0
    
    # 相対湿度は100%を超えることはないので、上限を設定する
    return min(rh, 100.0) 


def AHGANS(rhm, ml0):
    """含水率計算関数 (AHGANS)
    
    Args:
        rhm : 相対湿度 (%)
        ml0 : 材料種別コード

    Returns:
        含水率
    """

    d1 = rhm * 0.01

    # ML0の値に応じて処理を分岐
    # FortranのGO TO (11,10,13,14,15,16,17,18,19,20,21,22,13,13,13,13,13,28,29),ML0 に対応
    if ml0 == 1:
        # ラベル 11: SEASING BOARD
        if rhm <= 89:
            wd = -math.log(1 - 0.01 * rhm) / 0.1612944
        elif 89 < rhm < 95:
            wd = 639.22 - 14.988 * rhm + 0.08943 * rhm * rhm
        else: # rhm >= 95
            wd = 2.446 * rhm - 209.9

    elif ml0 in [3, 13, 14, 15, 16, 17]:
        # ラベル 13: GLASS WOOL
        # ※元のコードでML0が3, 13-17の際にラベル13へ飛ぶ設定
        d2 = math.exp(-1.3016 * (1 - d1**101.667))
        d3 = 0.16 * d2 * d1**7.6448
        if d1 > 0.98:
            d3 = 47.79595 * d1 - 46.79595
        wd = d3 * 100.0

    elif ml0 == 4:
        # ラベル 14: WOOD
        d2 = math.exp(-2.1708 * (1 - d1**90.988))
        d3 = 1.58 * d2 * d1**1.3664
        wd = d3 * 100.0

    elif ml0 == 5:
        # ラベル 15: PLY WOOD
        d2 = math.exp(-1.6471 * (1 - d1**32.051))    #ASHRAE PLYWOOD 1
        d3 = 1.0625 * d2 * d1**1.5373
        wd = d3 * 100.0

    elif ml0 == 6:
        # ラベル 16: PLASTER BOARD
        d2 = math.exp(-1.6445 * (1 - d1**75.811))    #ASHRAE GYPSUM
        d3 = 0.56978 * d2 * d1**0.18552
        wd = d3 * 100.0

    elif ml0 == 7:
        # ラベル 17: ALC
        d2 = math.exp(-0.66866 * (1 - d1**17.316))
        d3 = 0.071101 * d2 * d1**1.0146
        wd = d3 * 100.0

    elif ml0 == 8:
        # ラベル 18: NANSHITSUSENIBAN
        d2 = math.exp(-0.93059 * (1 - d1**4.6145))
        d3 = 0.35093 * d2 * d1**0.48739
        wd = d3 * 100.0

    elif ml0 == 9:
        # ラベル 19: THERMOPLY
        d2 = math.exp(-1.2566 * (1 - d1**5.483))
        d3 = 0.4367 * d2 * d1**0.3513
        wd = d3 * 100.0

    elif ml0 == 10:
        # ラベル 20: SAIDHING
        d2 = math.exp(-1.1225 * (1 - d1**15.878)) #ASHRAE Siding No.38
        d3 = 0.635 * d2 * d1**2.5972
        wd = d3 * 100.0

    elif ml0 == 11:
        # ラベル 21: tutikabe
        d2 = math.exp(-0.69 * (1 - d1**25.06))
        d3 = 0.0645 * d2 * d1**0.658
        wd = d3 * 100.0

    elif ml0 == 12:
        # ラベル 22: 軽量モルタル
        d2 = math.exp(-1.75 * (1 - d1**2.015))
        d3 = 0.203 * d2 * d1**(-0.120)
        wd = d3 * 100.0

    elif ml0 == 18:
        # ラベル 28: 集成材（OMソーラー）
        d2 = math.exp(-0.849 * (1 - d1**9.527))
        d3 = 0.381 * d2 * d1**(0.738)
        wd = d3 * 100.0

    elif ml0 == 19:
        # ラベル 29: 構造用合板（OMソーラー）
        d2 = math.exp(-0.907 * (1 - d1**9.412))
        d3 = 0.382 * d2 * d1**(0.670)
        wd = d3 * 100.0

    # ラベル 10: CONTINUE (処理の終了)
    return wd


def DIFF(n_mat: int, rh:float, k: float):
    """液体伝導率(RML)の計算

    Args:
        n_mat: 材料番号
        rh: 相対湿度, %
        k: 絶対温度, K
    Returns:
    """

    # 含水率
    wd = AHGANS(rhm=rh, ml0=n_mat)
    
    if 45 < wd < 110:
        d1 = wd * 0.01
        # DW = e^(a + bX + cX^2)
        dw = math.exp(-30.91 + 4.2967 * d1 - 0.22017 * (d1**2))
        
        rh1 = rh + 0.005
        rh2 = rh - 0.005
        
        # 微分を差分で近似している計算
        wp1 = get_wp(rh1, k)
        wp2 = get_wp(rh2, k)
        wd1 = AHGANS(rhm=rh1, ml0=n_mat)
        wd2 = AHGANS(rhm=rh2, ml0=n_mat)
        
        # ゼロ除算のチェック
        if abs(wp1 - wp2) > 1e-12:
            rml = 998.0 * dw * abs((wd1 - wd2) * 0.01 / (wp1 - wp2))
        else:
            rml = 0.0
    else:
        rml = 0.0
        
    return rml


def CALDPDU(wpt, tp, gma, ml0, row):
    """
    含水率変化に対するポテンシャル変化率 (DPDU) を計算する。(J/kg)/K
    gma:kg/m3
    TODO: 単位がよくわからない。
    row: kg/m3
    """
    d1 = wpt
    
    # 元のFortranのロジック: 正の値の場合は -100.0 に強制
    if d1 > 0.0:
        d1 = -100.0
    
    # 差分幅の計算 (d1の1%), J/kg
    dw = abs(d1 * 0.01)
    
    # 微小変化させた含水率, J/kg
    w1 = d1 + dw
    w2 = d1 - dw
    
    # 以前定義した WPTRE (平衡相対湿度計算) を呼び出し, %
    rh1 = WPTRE(w1, tp)
    rh2 = WPTRE(w2, tp)
    
    # 外部定義されている前提の AHGANS, 含水率
    wd1 = AHGANS(rhm=rh1, ml0=ml0)
    wd2 = AHGANS(rhm=rh2, ml0=ml0)
    
    # VGTの計算 (0.01は%を小数に戻す係数と推測)
    # kg/m3 / 
    vgt1 = 0.01 * wd1 * gma / row
    vgt2 = 0.01 * wd2 * gma / row
    
    # 中央差分による勾配(微分値)の近似
    # 0.5 * (VGT1 - VGT2) / DW
    if dw != 0:
        dpdu = 0.5 * (vgt1 - vgt2) / dw
    else:
        dpdu = 0.0
        
    return dpdu
