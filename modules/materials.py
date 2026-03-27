import pandas as pd
import os
from dataclasses import dataclass


@dataclass
class Material:

    # 名前
    name: str

    # 熱伝導率, W/(m K) 
    rmd: float

    # 湿気伝導率, (kg/s)/(m Pa)
    rmdd: float

    # 比熱, J/(kg K)
    cgg: float

    # 密度 kg/m3
    gma: float

    # 空隙率, kg/m3
    cggd: float


class Materials:

    def __init__(self):

        self.ms: list[Material] = self.load()
    
    @classmethod
    def load(self):

        file_name = 'material_prop.csv'

        # absolute file path
        path_and_filename = str(os.path.dirname(__file__)) + '/' + file_name

        try:

            if not os.path.isfile(path_and_filename):
                raise FileExistsError('The file does not exist.')
            
            df = pd.read_csv(path_and_filename, skipinitialspace=True, index_col='name')

            # インデックス（name）に重複がないかチェック
            if not df.index.is_unique:
                # 重複している名前を特定してエラーメッセージに入れる
                duplicates = df.index[df.index.duplicated()].unique().tolist()
                raise ValueError(f"CSVの 'name' 列に重複があります: {duplicates}")

        except ValueError as e:
            print(f"エラーが発生しました: {e}")
        
        print("データの読み込みに成功しました。")

        material_list = []

        # 1行ずつ取り出してデータクラスに格納
        for name, row in df.iterrows():
            m = Material(
                name=name,
                rmd=row['thermal_conductivity(W/mK)'],
                rmdd=row['moisture_conductivity(kg/msPa)'],
                cgg=row['specific_heat(J/kgK)'],
                gma=row['density(kg/m3)'],
                cggd=row['porosity']
            )
            material_list.append(m)

        return material_list

    def get_material(self, name: str) -> Material:
        
        ms = []

        for m in self.ms:

            if m.name == name:

                ms.append(m)
        
        if len(ms) == 0:
            raise KeyError('指定の名前の材料が見つかりませんでした。')
        
        elif len(ms) > 2:
            raise Exception()
        
        else:

            return ms[0]
        

