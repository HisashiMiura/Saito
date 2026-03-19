# Whamcavity

Heat and Water Transfer Simulation for Building Envelopes based on Water Chemical Potential.

建築外皮における水分化学ポテンシャルに基づく熱・水分移動シミュレーション。

## プロジェクトステータス

| 項目 | 状態 |
|------|------|
| 完成度 | **85%** (Beta移行直前) |
| テスト | 216 passed + 1 xfailed |
| Fortran互換性 | 12件の直接比較テスト合格 |
| 最終評価 | 2026-02-04 第5回討論 |

## 概要

Whamcavityは、多層壁体の1次元熱水分同時移動解析を行うシミュレーションプログラムです。
Hiroaki Saito氏によるFortran 95コード（2021年8月）をPython/NumPyベースに変換したものです。

### 主な機能

- 多層壁体の1次元熱・水分移動解析
- 多階建築外皮の通気層を考慮したモデル
- 風向を考慮した雨水浸入計算
- 木材腐朽（rot）損傷関数

## インストール

```bash
# リポジトリをクローン
git clone <repository-url>
cd whamcavity

# 依存パッケージのインストール
pip install numpy
```

### 依存関係

- Python 3.9+
- NumPy

### テスト実行

```bash
# pytest のインストール（必要な場合）
pip install pytest

# テスト実行
python -m pytest tests/ -v
```

### テスト状況（216 passed + 1 xfailed）

| カテゴリ | テスト数 | 説明 |
|----------|----------|------|
| Fortran直接比較 | 12件 | 同一入力→同一出力を検証 |
| E2Eテスト | 7件 | 完全なシミュレーション実行 |
| 保存則検証 | 4件 | エネルギー・水分保存則 |
| RH影響分析 | 8件 | 高湿度域での挙動検証 |
| 熱力学関数 | 24件 | Goff-Gratch式等 |
| 材料物性 | 18件 | 吸放湿等温線 |
| 木材腐朽 | 16件 | 損傷関数 |
| 日射計算 | 11件 | 太陽位置・相当外気温 |

**重要なテスト:**
- **Fortran直接比較**: Python版とFortran版が同一結果を出力することを確認
- **往復テスト**: `rewpt` ↔ `wptre` が逆関数であることを確認
- **既知値照合**: Goff-Gratch式の飽和蒸気圧（0℃で611Pa、100℃で101325Pa）
- **保存則**: 閉鎖系でのエネルギー・水分保存を確認

## 使用方法

### 対話モード（Fortran互換）

```bash
python -m whamcavity
```

以下のプロンプトが表示されます：
```
形状データファイルを入力して下さい: shape.dat
気象データファイルを入力して下さい: weather.dat
風雨データファイルを入力して下さい: windrain.dat
出力ファイルを入力して下さい（8バイト以内）: output
収束計算あり=1   収束計算なし=0: 1
腐朽緩和係数を入力して下さい: 1.0
水分生成の扱い　無視0　考慮1: 0
```

### バッチモード

```bash
python -m whamcavity shape.dat weather.dat windrain.dat output [options]
```

オプション:
- `--material-file FILE`: 材料物性ファイル（デフォルト: mcoff1.prn）
- `--no-convergence`: 収束計算を無効化
- `--decay-relaxation FLOAT`: 腐朽緩和係数（デフォルト: 1.0）
- `--moisture-generation`: 水分生成を考慮

### ライブラリとして使用

```python
from whamcavity import Simulation, run_simulation

# 簡易実行
run_simulation(
    shape_file="shape.dat",
    weather_file="weather.dat",
    wind_rain_file="windrain.dat",
    output_file="output"
)

# 詳細制御
sim = Simulation()
sim.load_input_files(
    shape_file=Path("shape.dat"),
    weather_file=Path("weather.dat"),
    wind_rain_file=Path("windrain.dat"),
    material_file=Path("mcoff1.prn")
)
sim.use_convergence = True
sim.decay_relaxation = 1.0
sim.initialize("output")
sim.run()
```

## 入力ファイル形式

### 形状データファイル

壁体構成、材料配置、初期条件を定義します。

```
（コメント行）
緯度、経度
35.0 135.0 0.0
（コメント行）
年数 開始月 開始日 終了月 終了日
1 1 1 12 31
...
```

### 気象データファイル

1日25時間分（0時〜24時）の気象データ。

フォーマット: `外気温(℃) 絶対湿度(g/kg) 直達日射 天空日射 夜間放射 ...`

### 風雨データファイル

1日24時間分の風速・降雨量・風向データ。

フォーマット: `風速(m/s) 降雨量(mm/h) 風向(°)`

### 材料物性ファイル (mcoff1.prn)

19種類の材料の物性値。

フォーマット: `番号 熱伝導率 透湿率 比熱 密度`

## モジュール構成

```
whamcavity/
├── __init__.py          # パッケージ初期化、全エクスポート
├── __main__.py          # python -m whamcavity 対応
├── config.py            # 定数・パラメータ定義
├── datatypes.py         # データ構造（dataclass）
├── thermodynamics.py    # 熱力学関数
├── moisture.py          # 水分移動・材料物性
├── solar.py             # 日射計算
├── room.py              # 室内環境計算
├── damage.py            # 木材腐朽モデル
├── io_handlers.py       # ファイル入出力
├── solver.py            # 過緩和法ソルバー
├── simulation.py        # メインシミュレーションクラス
└── main.py              # コマンドラインエントリーポイント
```

### 主要関数

| 関数 | モジュール | 説明 |
|------|-----------|------|
| `goff` | thermodynamics | Goff-Gratch式による飽和蒸気圧計算 |
| `satuwpt` | thermodynamics | 飽和水分ポテンシャル計算 |
| `rewpt` | thermodynamics | 相対湿度→水分ポテンシャル変換 |
| `wptre` | thermodynamics | 水分ポテンシャル→相対湿度変換 |
| `ahgans` | moisture | 19材料の吸放湿等温線 |
| `caldpdu` | moisture | 水分容量係数計算 |
| `diff` | moisture | 液水伝導率計算 |
| `soloct` | solar | 太陽位置計算 |
| `satcal` | solar | 相当外気温計算 |
| `wood_rot` | damage | 木材腐朽損傷関数 |
| `rh_critical` | damage | 腐朽開始臨界湿度 |

## 材料ID一覧

| ID | 材料名 |
|----|--------|
| 1 | シージングボード |
| 2 | 空気層 |
| 3 | グラスウール |
| 4 | 木材 |
| 5 | 合板 |
| 6 | 石膏ボード |
| 7 | ALC |
| 8 | 軟質繊維板 |
| 9 | サーモプライ |
| 10 | サイディング |
| 11 | 土壁 |
| 12 | 軽量モルタル |
| 13-17 | グラスウール（同等） |
| 18 | 集成材 |
| 19 | 構造用合板 |

## 計算パラメータ

| パラメータ | 値 | 説明 |
|-----------|-----|------|
| NXP | 150 | 最大メッシュ点数 |
| KMTLP | 10 | 最大材料層数 |
| NWP | 50 | 最大壁部位数 |
| NRM | 50 | 最大室数 |
| NMTL | 19 | 材料種類数 |
| EPS1 | 0.01 K | 温度収束判定値 |
| EPS2 | 100 Pa | 水分ポテンシャル収束判定値 |
| OMG | 1.2 | 過緩和係数 |
| MXIT | 5000 | 最大反復回数 |

## 出力ファイル

### .DAT ファイル

時刻別の詳細データ（温度、湿度、含水率など）

### .TXT ファイル

日別サマリー（浸水量、平均含水率、腐朽損傷など）

## 理論的背景

### 水分化学ポテンシャル

水分移動の駆動力として水分化学ポテンシャル（水分ポテンシャル）を使用：

```
ψ = R_v × T × ln(RH/100)
```

- ψ: 水分ポテンシャル [J/kg]
- R_v: 水蒸気の気体定数 (461.5 J/(kg·K))
- T: 絶対温度 [K]
- RH: 相対湿度 [%]

### 過緩和法（SOR法）

温度場・水分場の収束計算にGauss-Seidel法と過緩和を適用：

```
φ^(n+1) = φ^(n) + ω × (φ* - φ^(n))
```

- ω: 過緩和係数（1.2）
- φ*: Gauss-Seidel更新値

### 木材腐朽モデル

温度・湿度条件に基づく2段階モデル：

1. **発芽期**: 臨界条件到達までの待機時間
2. **成長期**: 温度依存の反応速度による質量減少

## 既知の制限事項

### RH 98%クランプ

高湿度域（RH > 98%）では、吸放湿等温線モデル（GAB式）が数学的に発散するため、
相対湿度を98%で上限クランプしています。

**影響:**

| 材料 | RH 99%以上での含水率過小評価 |
|------|------------------------------|
| 木材 | 最大17.4%（30%→24.77%） |
| 石膏ボード | 最大21.8%（20%→15.64%） |
| 合板/集成材 | 0%（飽和上限に達済） |

**Fortran版との差異:**
- Fortran版: `IF(RH.GT.99.98) CYCLE` でスキップ
- Python版: `rh = min(rh_percent, 98.0)` でクランプ

**推奨事項:**
- 通常湿度範囲（RH < 95%）での使用を推奨
- 結露・高湿度環境での定量評価には注意が必要
- 定性的傾向分析（腐朽リスク評価等）には適用可能

### 定常熱流束バランス（xfailed）

開放系境界条件下では、定常状態での熱流束バランステストが理論的制限により
完全には成立しません。これは実装の欠陥ではなく、物理モデルの特性です。

## ライセンス

Original Fortran code by Hiroaki Saito (August 2021).

## 参考文献

- 日本建築学会: 建築物の熱・水分性能評価法
- ASHRAE Handbook - Fundamentals
- 伊庭ら: 建築材料の吸放湿特性に関する研究
