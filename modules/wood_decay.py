import numpy as np

from modules.config import ROTOMOG


# 三浦コメント：このパラメータが何なのか不明
# OSB
WS = -1.0

# 質量減少が起きる場合の相対湿度の閾値か？
# 齋藤先生に要確認
# %
RH_GROWTH = 98.0


def get_mass_loss_is(l_stage_is: list[bool], theta_is: list[float], rh_is: list[float]):

    mass_loss_is = np.zeros_like(a=l_stage_is, dtype=float)

    for i in range(len(l_stage_is)):

        # 両側が発芽しているかフラグが立っていれば b=True にする。
        # この判定はセルの分割の仕方に大きく依存するため、物理的にもっと普遍的な記述がないか？
        if i == 0:
            if l_stage_is[i] and l_stage_is[i+1]:
                b = True
            else:
                b = False
        elif i == len(l_stage_is) - 1:
            if l_stage_is[i-1] and l_stage_is[i]:
                b = True
            else:
                b = False
        else:
            if l_stage_is[i-1] and l_stage_is[i] and l_stage_is[i+1]:
                b = True
            else:
                b = True
        
        if b:
            mass_loss_is[i] = get_mass_reduction(rh=rh_is[i], theta=theta_is[i]) * ROTOMOG
        else:
            mass_loss_is[i] = 0.0

    return mass_loss_is


def update_time_s(theta: float, rh: float, time_s: float) -> float:

    # 温度ごとに定義された相対湿度の閾値。
    rhc = rh_threshold(theta=theta)

    # 相対湿度の閾値を下回る場合は積算時間がリセットされる。
    # もともと0℃以下だと積算時間がリセットという記述があったが、0℃以下の場合の相対湿度の閾値は100%に設定されるため、
    # その記述は不要と考えて削除した。
    if rh < rhc:
    
        return 0.0
    
    else:

        # 時間経過の蓄積
        time_s += 24 * 3600
        
        return time_s


def update_stage(theta: float, rh: float, time_s: float, l_stage: bool):

    # 温度ごとに定義された相対湿度の閾値, %
    rhc = rh_threshold(theta=theta)

    # 温度が0℃以下になると積算時間がリセットされる。
    # 相対湿度の閾値を超えない場合は積算時間がリセットされる。
    # この前に書いてある温度が0℃以下というのは必要ないのではないか？
    if theta <= 0.0 and rh < rhc:
    
        return l_stage
    
    else:

        # Mostafa Nofal, Kumar Kumaran, Biological damage function models for durability assessments of wood and woo-based products in building envelopes
        # Eur. J. Wood Prod. 2011 69:619-631
        fc = (0.1384 * theta + 0.4370 * rh - 42.9450 + WS * (0.034 * theta - 0.021 * rh + 1.721))        
        gc = (-2.2270 * theta - 0.0347 * rh + 0.0244 * theta * rh + WS * (-0.504 * theta + 0.0096 * rh + 0.0047 * theta * rh))

        # ゼロ除算回避
        if fc != 0:
            d1 = - gc / fc
        else:
            d1 = 100.0 # 仮の大きな値

        # 発芽時間の計算 (TIME_INT)
        # d1 がある値より大きくなると、、、、
        # 発芽が始まるまでの時間（相対湿度がある閾値を超えかつある温度を超えた瞬間からの時間）
        if d1 > 0:
            time_int = d1
        elif d1 > -0.59:
            time_int = 0.0
        else:
            time_int = 100.0
        
        # 係数調整
        time_int = max(time_int, 0.5) * 30.0 * 24 * 3600
        
        if time_s > time_int:
            l_stage = True
        
        return l_stage
   

def rh_threshold(theta: float) -> float:
    """発芽する際の相対湿度の閾値を求める。

    Args:
        theta: 温度, ℃

    Returns:
        閾値（相対湿度）, %
    """

    # 閾値を求める関数。最低値は92.5%, 最大値は100.0%
    return np.clip(
        a=100-0.5*theta,
        a_min=92.5,
        a_max=100.0
    )


def get_mass_reduction(rh: float, theta: float):
    """_summary_

    Args:
        rh: 相対湿度, %
        theta: 温度, ℃
    """

    # 相対湿度が閾値以上の場合質量減少が起きる。
    # 温度が0℃より大かつ40℃以下
    if rh >= RH_GROWTH and theta > 0.0 and theta <= 40.0:

        # 腐朽反応速度の計算
        # 齋藤先生の論文の係数か？
        # 反応速度定数, (kg/kg)/s
        # 反応速度定数とは質量減少率（腐朽前の質量に対する腐朽時の質量減少量の比）を時間で除した値。
        k = (2.77 - 3.23 * theta + 0.865 * (theta**2) - 0.0189 * (theta**3)) * 1e-10

        # 1日あたりの質量減少量, (kg/kg)/d
        # 1日あたりにしているが、かえって複雑になるため、時々刻々質量減少量を計算するモデルにした方が良い。
        # 計算量は増えるが余計なif文が減るためかえって計算速度はあがるのではないか？
        return k * 24 * 3600
    
    else:

        # 閾値を超えなければ質量減少は起きない。
        return 0.0


