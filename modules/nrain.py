from dataclasses import dataclass

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




