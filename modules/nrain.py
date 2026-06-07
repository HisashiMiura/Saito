from dataclasses import dataclass
import numpy as np
from scipy.stats import rayleigh


# ***********浸水率の設定*********

NRAINS = 14

# 厚壁No
# 座標
# 層No
NRAINPOINT = [
    [1,5,1],
    [1,7,3],
    [2,5,1],
    [2,7,3],
    [3,5,1],
    [3,7,3],
    [4,5,1],
    [4,7,3],
    [5,5,1],
    [5,7,3],
    [6,5,1],
    [6,7,3],
    [7,5,1],
    [7,7,3]
]

# 浸水率%
# 閾値風速m/s
ACOFRAIN = [
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
    [0.5, 0],
]

class NRAIN:

    # 厚壁No.
    wall_no: int

    # 座標
    pos: int

    # 層No
    layer_no: int

    # 浸水率
    ratio: float

    # 閾値風速, m/s
    v: float


def set_default_nrains():

    # 厚壁No
    # 座標
    # 層No
    NRAINPOINT = [
        [1,5,1],
        [1,7,3],
        [2,5,1],
        [2,7,3],
        [3,5,1],
        [3,7,3],
        [4,5,1],
        [4,7,3],
        [5,5,1],
        [5,7,3],
        [6,5,1],
        [6,7,3],
        [7,5,1],
        [7,7,3]
    ]

    # 浸水率%
    # 閾値風速m/s
    ACOFRAIN = [
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
        [0.5, 0],
    ]

    return [
        NRAIN(
            wall_no=n[0],
            pos=n[1],
            layer_no=n[2],
            ratio=a[0],
            v=a[1]
        )
        for (n, a)
        in zip(NRAINPOINT, ACOFRAIN)
    ]

def get_nrains_of_walls(wall_index: int) -> list[NRAIN]:

    nlains_of_wall: list[NRAIN] = []
    for nrain in set_default_nrains():

        if nrain.wall_no == wall_index:
            nlains_of_wall.append(nrain)
    
    return nlains_of_wall


@dataclass
class RainLeakage:

    is_rainpoint: bool

    ratio: float

    wind_threshold: float

    def get_confrain(self, v_mod: float, angle: float) -> float:

        """雨水浸入の割合を計算する。

        Args:
            v_mod: 風速, m/s
            angle: 風向と壁の法線のなす角度, 度

        Returns:
            雨水浸入の割合, -
        """

        if self.is_rainpoint:

            # 風速の分布, m/s, [N]
            x = np.linspace(0.1 * v_mod, 3.0 * v_mod, 30)

            if v_mod > self.wind_threshold:

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
                    self.ratio * x * np.cos(np.radians(angle)) * FE * FD * FL,
                    0.0
                )

                sigma =  v_mod / np.sqrt(np.pi / 2)
                weight = rayleigh.pdf(x=x, scale=sigma)

                return np.sum(weight * d2) / np.sum(weight) * 0.01

            else:

                return 0.0

        else:

            return 0.0
    
