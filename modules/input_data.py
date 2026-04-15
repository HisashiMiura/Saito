from dataclasses import dataclass
import numpy as np


import .direction import Direction


@dataclass
class InputRoom:

    # 平均温度, ℃
    average_temperature: float

    # 振幅, ℃
    amplitude_temperature: float

    # 平均湿度
    average_relative_humidity: float

    @classmethod
    def read(cls, d: dict):

        try:
            average_temperature = d.get('average temperature', 22.5)
        except:
            raise ValueError('\'average temperature\' に不正な値が指定されました。')
        
        try:
            amplitude_temperature = d.get('amplitude temperature', 4.5)
        except:
            raise ValueError('\'amplitude temperature\' に不正な値が指定されました。')

        try:
            average_relative_humidity = d.get('average relative humidity', 60.0)
        except:
            raise ValueError('\'average relative humidity\' に不正な値が指定されました。')

        InputRoom(
            average_temperature=average_temperature,
            amplitude_temperature=amplitude_temperature,
            average_relative_humidity=average_relative_humidity
        )


@dataclass
class InputWall:

    # 長辺, m
    len_long: float

    # 短辺, m
    len_short: float

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

        # 長辺, m
        len_long = d['len_long']

        # 短辺, m
        len_short = d['len_short']

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
        height = d.get('height', len_long * np.sin(np.radians(angle)))

        # 評価高さ, m
        eva_height = d['eva_height']

        return InputWall(
            len_long=len_long,
            len_short=len_short,
            direction=direction,
            angle=angle,
            albedo=albedo,
            absorption=absorption,
            emissivity=emissivity,
            height=height,
            eva_height=eva_height
        )

