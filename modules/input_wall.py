from dataclasses import dataclass

from .direction import Direction


@dataclass
class InputWall:

    # 垂直方向の長さ, m
    len_vertical: float

    # 水平方向の長さ, m
    len_horizontal: float

    # 方位
    direction: Direction

    # 傾斜角
    angle: float

    # アルベド
    albedo: float

    # 室外側表面の日射吸収率
    absorption: float

    # 室外側表面の長波長放射率
    emissivity: float

    # 壁の下端と上端の高さの差（換気計算に用いられる）, m
    height: float

    # 評価高さ, m
    eva_height: float

    @classmethod
    def read(cls, d: dict):

        # 垂直方向の長さ, m
        len_vertical = d['len_vertical']

        # 水平方向の長さ, m
        len_horizontal = d['len_horizontal']

        # 方位
        direction = Direction(d['direction'])

        # アルベド
        albedo = d.get('albedo', 0.1)

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

        # 室外側表面の日射吸収率
        absorption = d.get('absoption', 0.9)

        # 室外側表面の長波長放射率
        emissivity = d.get('emissivity', 0.9)

        # 壁の下端と上端の高さの差（換気計算に用いられる）, m
        height = d.get('height', len_vertical * np.sin(np.radians(angle)))

        # 評価高さ, m
        eva_height = d['eva_height']

        return InputWall(
            len_vertical=len_vertical,
            len_horizontal=len_horizontal,
            direction=direction,
            angle=angle,
            albedo=albedo,
            absorption=absorption,
            emissivity=emissivity,
            height=height,
            eva_height=eva_height
        )

