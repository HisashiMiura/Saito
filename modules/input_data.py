from dataclasses import dataclass
import numpy as np


from modules.direction import Direction


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


