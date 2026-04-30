import numpy as np


# 三浦コメント：このパラメータが何なのか不明
# OSB
W = -1.0

# 質量減少が起きる場合の相対湿度の閾値か？
# 齋藤先生に要確認
# %
RH_GROWTH = 98.0


def func_1(theta: float, rh: float, time_s: float, l_stage: bool):

    # 温度が0℃以下になると積算時間がリセットされる。
    if theta <= 0.0:
    
        time_s = 0.0
        return time_s, l_stage
    
    else:
    
        # 温度ごとに定義された相対湿度の閾値。
        rhc = rh_threshold(theta=theta)

        # 相対湿度の閾値を超えない場合は積算時間がリセットされる。
        # この前に書いてある温度が0℃以下というのは必要ないのではないか？
        if rh < rhc:

            time_s = 0.0
            return time_s, l_stage

        else:

            # TIN関数の呼び出し (引数wが必要)
            # wという変数があったが固定値だったため関数内に記述した。
            # w が何かは不明。フォートランのコメント「OSB」で-1.0
            gc, fc = tin(theta=theta, rh=rh)

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
            
            # 時間経過の蓄積
            time_s += 24 * 3600

            if time_s > time_int:
                l_stage = True
            
            return time_s, l_stage
   

def func_2():

    

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


def tin(theta: float, rh, w):
    """
    Args:
        theta: 温度, ℃
        rh: 相対湿度, %
    
    TODO: この関数が何なのか、緒言を確認すること。

    """

    fc = (0.1384 * theta + 0.4370 * rh - 42.9450 + W * (0.034 * theta - 0.021 * rh + 1.721))
    
    gc = (-2.2270 * theta - 0.0347 * rh + 0.0244 * theta * rh + w * (-0.504 * theta + 0.0096 * rh + 0.0047 * theta * rh))
    
    return gc, fc


def get_mass_reduction(rh: float, tmp: float):
    """_summary_

    Args:
        rh: 相対湿度, %
        tmp: 温度, ℃
    """

    # 相対湿度が閾値以上の場合質量減少が起きる。
    if rh >= RH_GROWTH:

        # 腐朽反応速度の計算
        # 齋藤先生の論文の係数か？
        # 反応速度定数, (kg/kg)/s
        # 反応速度定数とは質量減少率（腐朽前の質量に対する腐朽時の質量減少量の比）を時間で除した値。
        k = (2.77 - 3.23 * tmp + 0.865 * (tmp**2) - 0.0189 * (tmp**3)) * 1e-10

        # 1日あたりの質量減少量, (kg/kg)/d
        # 1日あたりにしているが、かえって複雑になるため、時々刻々質量減少量を計算するモデルにした方が良い。
        # 計算量は増えるが余計なif文が減るためかえって計算速度はあがるのではないか？
        return k * 24 * 3600
    
    else:

        # 閾値を超えなければ質量減少は起きない。
        return 0.0


