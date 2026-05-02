import pandas as pd
from pandas import Series
import os
from dataclasses import dataclass
import math


@dataclass
class Material:

    # 名前
    name: str

    # 熱伝導率, W/(m K)  ラムダ
    lambda_h: float

    # 湿気伝導率, (kg/s)/(m Pa)　ラムダダッシュ
    lambda_dsh_m: float

    # 比熱, J/(kg K)　シー
    c: float

    # 密度 kg/m3　ロー
    rho: float

    # 空隙率, m3/m3　プサイゼロ
    psi0: float
    
    # 相対湿度（0から1）から含水率（0から1）を求める関数。　
    f_u: callable

    def get_u(self, rh: float):
        """含水率（0.0～1.0）を求める。（質量含水率） ユー JIS A 1476 で建築材料の含水率測定法 が決められている。
           ちなみに、容積含水率の場合はプサイを使う場合が多いらしい。

        Args:
            rh: 相対湿度(0.0～100.0), %

        Returns:
            質量含水率（0.0～1.0）, kg/kg 水分の重量÷材料の乾燥重量
        """

        r = rh * 0.01

        return self.f_u(r)
    
    @classmethod
    def load(cls, name: str, row: Series):

        thermal_conductivity = row['thermal_conductivity(W/mK)']
        moisture_conductivity = row['moisture_conductivity(kg/msPa)']
        specific_heat = row['specific_heat(J/kgK)']
        density = row['density(kg/m3)']
        # エクセルの方もあわせて、列ヘッダ名を porosity(m3/m3) に直すこと。
        porosity = row['porosity']

        n_formula_branches = row['n_formula_branches']
        formula1 = str(row['formula1'])
        formula2 = str(row['formula2'])
        formula3 = str(row['formula3'])
        threshold2 = row['threshold2']
        threshold3 = row['threshold3']

        if n_formula_branches == 1:

            def f_psi_1(r: float):
                return eval(formula1)
            
            f_psi = f_psi_1

        elif n_formula_branches == 2:

            def f_psi_2(r: float):

                if r > threshold2:
                    return eval(formula2)
                else:
                    return eval(formula1)

            f_psi = f_psi_2

        elif n_formula_branches == 3:

            def f_psi_3(r: float):

                if r > threshold3:
                    return eval(formula3)
                elif r > threshold2:
                    return eval(formula2)
                else:
                    return eval(formula1)

            f_psi = f_psi_3

        return Material(
            name=name,
            lambda_h=thermal_conductivity,
            lambda_dsh_m=moisture_conductivity,
            c=specific_heat,
            rho=density,
            psi0=porosity,
            f_u=f_psi
        )




class Materials:

    def __init__(self):

        self.ms: list[Material] = self.load()
    
    @classmethod
    def load(self):

        #file_name = 'material_prop.csv'
        file_name = 'material_prop.xlsx'

        # absolute file path
        path_and_filename = str(os.path.dirname(__file__)) + '/' + file_name

        try:

            if not os.path.isfile(path_and_filename):
                raise FileExistsError('The file does not exist.')
            
            #df = pd.read_csv(path_and_filename, skipinitialspace=True, index_col='name')
            df = pd.read_excel(path_and_filename, sheet_name='material_prop', index_col='name')

            # インデックス（name）に重複がないかチェック
            if not df.index.is_unique:
                # 重複している名前を特定してエラーメッセージに入れる
                duplicates = df.index[df.index.duplicated()].unique().tolist()
                raise ValueError(f"CSVの 'name' 列に重複があります: {duplicates}")

        except ValueError as e:
            print(f"エラーが発生しました: {e}")
        
        print("データの読み込みに成功しました。")

        # 1行ずつ取り出してデータクラスに格納
        return [Material.load(name=name, row=row) for name, row in df.iterrows()]
    


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
    
