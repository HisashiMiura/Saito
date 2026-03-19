import numpy as np
from dataclasses import dataclass

# ////////////////////////////////////
#  Heat and Water transfer based on water chemical potential
#  Building Envelope Simulation for Multi-Story Model 
#  1-dimensional model made by Hiroaki Saito 
#  dxの変更、液水伝導率、通気層、結露量、気象データ、伝達率換算
#  For 漏水あり外壁多層階モデル  風向考慮
#  made by Hiroaki Saito  August 2021
# SI unit
# ///////////////////////////////////

IMPLICIT REAL*8(A-H,O-Z)

CHARACTER MOJI*110,OUTFILE*8,FILENAME0*12

!                    質点数　層数　　　壁数　部屋数

INTEGER,PARAMETER :: NXP=150,KMTLP=10,NWP=50,NRM=50  !,NYP=150,KATLP=70

# サブルーチンで内部保持しているパラメータ（time_s, l_stage）があったため、クラス化した。
wdm = WoodDecayModel(nwp=NWP, nxp=NXP)


INTEGER,PARAMETER :: MXIT=5000,NMTL=19
REAL,PARAMETER :: EPS1=0.01,EPS2=1.0e+2
INTEGER DAY,DAY2,DAY1
COMMON /MS/TEMPO(25),SJD(25),SJS(25),SJN(25),SAT(NWP,25),SJIN(NWP,25),WSJIN(NWP,25),WSJD(NWP,25)

@dataclass
class MS:

    # COMMON / MS で定義されていたグローバル変数
    TEMPO: list # [26]
    SJD: list # [26]
    SJS: list # [26]
    SJN: list # [26]
    SAT: list # [NWP+1, 26]
    SJIN: list # [NWP+1, 26]
    WSJIN: list # [NWP+1, 26]
    WSJD: list # [NWP+1, 26]
    
@dataclass
class SS:

    # COMMON / SS で定義されていたグローバル変数

    W: list #[31]
    DV4: list #[NWP+1, 31]
    BKAKU: float



!//////////COEFFOCIENT//////////////
DIMENSION NX(NWP)
DIMENSION DPDU(NWP,NXP),TMPC(NWP,NXP),WGT(NWP,NXP)
DIMENSION DWGX(NWP,NXP+1),DWLX(NWP,NXP+1),RMDX(NWP,NXP+1),DTGX(NWP,NXP+1),DTLX(NWP,NXP+1),DWX(NWP,NXP+1),DTX(NWP,NXP+1)
DIMENSION ADWGX(NWP,NXP+1),ADWLX(NWP,NXP+1),ADTGX(NWP,NXP+1),ADTLX(NWP,NXP+1),ADWX(NWP,NXP+1),ADTX(NWP,NXP+1)
DIMENSION ALF(NWP,KMTLP,2),ALD(NWP,KMTLP,2),ALDP(NWP,KMTLP,2),ALDT(NWP,KMTLP,2)
DIMENSION ALP_TOTAL(KMTLP),TEMP_CAVITY(KMTLP),HIGHT_CAVITY(KMTLP)
DIMENSION KMTL(KMTLP),WOUTAV(50)
!/////////VARIABLE///////////////
DIMENSION WPU(NWP,NXP),WPS(NWP,NXP),WPW(NWP,NXP),TMP(NWP,NXP),RH(NWP,NXP)
DIMENSION XN(NWP,NXP)
DIMENSION AWPU(NWP,NXP),ATMP(NWP,NXP),BWPU(NWP,NXP),BTMP(NWP,NXP)     !,AWPS(NXP),AWPW(NXP)
DIMENSION HWPU(NWP,NXP),HTMP(NWP,NXP)    !,HWPS(NXP),HWPW(NXP)
! DIMENSION STMP(KMTLP,2),SWPS(KMTLP,2),SWPU(KMTLP,2),SWPW(KMTLP,2)&
DIMENSION QS(NWP,KMTLP,2)  !,SRH(NWP,KMTLP,2),ASWPU(NWP,KMTLP,2)
DIMENSION SCWDAY(NWP,25),MCW(NWP),SCW(NWP),CWV(NWP),CWVDAY(NWP,25)
DIMENSION XOD(25),RHOD(25),SATDV(NWP),QQ(NWP,NXP),XM(NWP,NXP),AXN(NWP,NXP),RHDIS(NWP,NXP,12,31),WGTDIS(NWP,NXP,12,31)
DIMENSION RHAVD(NWP,NXP),WGTAVD(NWP,NXP),DGDUQ(NWP,NXP),DGDTQ(NWP,NXP)
DIMENSION ATMPQ(NWP,NXP),SATA(NWP)
!DIMENSION TDES1(24),TDES2(24),THDES1(24),THDES2(24),TSDES1(24),TSDES2(24),THSDES1(24),THSDES2(24)
DIMENSION TPAVD(NWP,NXP),XNAVD(NWP,NXP),TPDIS(NWP,NXP,12,31),XNDIS(NWP,NXP,12,31)
DIMENSION WLOSS(NWP,NXP),WJW(NWP,NXP)  !DAMAGE FUNC
DIMENSION WJRAIN(NWP,NXP),SWJRAIN(NWP,NXP),COFRAIN(NWP,2),VELOUT(24),RAIN(24),SUMRAIN(NWP,NXP,25)
DIMENSION RN(NWP,NXP),HRN(NWP,NXP),DAYSUMRAIN(NWP,NXP),VELDIR(24)
DIMENSION MDAY(12)
DIMENSION TMPR(NRM),ATMPR(NRM),RHR(NRM),XR(NRM),WPUR(NRM),HTMPR(NRM),XRM(NRM)!,VPR(NRM)
DIMENSION DGDUR(NRM),HWPUR(NRM),DGDTR(NRM),AWPUR(NWP)
DIMENSION WSIN(NWP,25),WSIND(NWP,25)
DIMENSION WSUP(NWP,25),HSUP(NWP,25)

!/////////////CONSTANT/////////////////
DATA RG,RW,CPL,ROW/461.5,2.512D+6,4200,998/! [J/kg],CPL[J/kgK],ROW[kg/m3]
!DATA RG,RW,CPL,ROW/461.5,0.,4200,998/! [J/kg],CPL[J/kgK],ROW[kg/m3]
DATA ATP/273.15/
DATA OMG/1.2/
DATA MDAY/31,28,31,30,31,30,31,31,30,31,30,31/
!RW=0.
!
!**** FILE OPEN ****

  PRINT *,  '　　気象データファイルを入力して下さい'
  OPEN(UNIT=7,FILE="")

!PRINT *,  '　　物性データファイルを入力して下さい'
OPEN(UNIT=11,FILE="mcoff1.prn")

  PRINT *,  '　　風雨データファイルを入力して下さい'
OPEN(UNIT=13,FILE="")

READ(13,210)MOJI
READ(13,210)MOJI
PRINT *,  '　出力ファイルを入力して下さい（8バイト以内）'
READ(5,*)OUTFILE
PRINT *,  '　収束計算あり=1   収束計算なし=0 '
READ(5,*)NREPT
PRINT *,  '腐朽緩和係数を入力して下さい '
READ(5,*)ROTOMG
PRINT *,  '水分生成の扱い　無視0　考慮1 '
READ(5,*)I_HCOFF
PRINT *,  '通気層流量出力ポイント　厚壁No及び層番号1,層番号2 '
READ(5,*)NQOUT1
READ(5,*)NQOUT2
READ(5,*)NQOUT3

# 水分化学ポテンシャルによる建物の温湿度計算

# 緯度
HIDO = 36.7

# 経度
HKEIDO = 137.21

# 方位のズレ
HOI = 0

# 年
LYEAR = 2
 
# 開始月
MON1 = 1
# 開始日
MDAY1 = 1

# 終了月
LMON = 12

# 終了日
LDAY = 31

# 外気温度（初期値）
TMPOC = 15

# 室内温度（初期値）
TMPIC = 15

# 外気湿度（初期値）
RHO = 50

# 室内湿度（初期値）
RHI = 50

# 時間分割（1/DT {h}）
NDVD = 20

# 部屋数
NROOM = 1

@dataclass
class Room:

    # 体積
    volume: float

    # 換気回数
    n_vent: float

    # 計算の有無
    is_calc: bool

    # 平均温度
    t_ave: float

    # 振幅
    t_amp: float

    # 平均湿度
    rh_ave: float

    # 開口数
    n_open: int

    # 空気移動
    flows: list[Flow]

rooms = [
    Room(volume=425.9, n_vent=1, is_calc=False, t_ave=22.5, t_amp=4.5, rh_ave=60.0, n_open=1, flows=[Flow(idx_to=1, idx_from=0, volume=51.5)])
]

@dataclass
class Flow:

    # 自室
    idx_to: int

    # 相手室（外気の場合は0）
    idx_from: int

    # 流入量
    volume: float

# 壁体部位数 NWTYPE 壁TYPE数
NWTYPE = 2

# 部材数
KMTL = [7, 7]

@dataclass
class Layer:

    # 部材No.
    num: int

    # 名称
    name: str

    # 初期温度, ℃
    initial_temp: float

    # 初期湿度, %
    initial_humidity: float

    # 幅, m
    thick: float

    # 分割数
    n_div: int

    # 熱コンダクタンス, W/(m2 K)（外面）
    cond_h_o: float

    # 熱コンダクタンス, W/(m2 K)（内面）
    cond_h_i: float

    # 湿気コンダクタンス, kg/(m s Pa)（外面）
    cond_m_o: float

    # 湿気コンダクタンス, kg/(m s Pa)（内面）
    cond_m_i: float

    # 斎藤先生コメント
    # 伝達率に関しては正確にはコンダクタンスとして扱っており、防湿層がある場合は透湿抵抗値の逆数（ここでは2.00E-11）を入れます。
    # 材料同士が接触する場合のコンダクタンスは、上の表では熱コンダクタンスは50（W/m2K）、湿気コンダクタンスは7.75E-07（kg/msPa）を暫定的に入れています。
    # 塗膜やクロスなどの透湿抵抗も、ここで調整します。

    # 通気層の場合の流量係数α（高さ1mあたりの値）　三浦コメント：これは高さ1mあたりではなく、幅1mあたりの間違いでは？
    alpha: float

    @property
    def dx(self):
        if self.n_div == 1:
            d1 = 1
        else:
            d1 = self.n_div - 1
        
        return self.thick / d1


layers1 = [
    Layer(num=10, name='サイディング', initial_temp=15.0, initial_humidity=50.0, thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
    Layer(num=2, name='通気層', initial_temp=15.0, initial_humidity=50.0, thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
    Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
    Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
    Layer(num=3, name='GW', initial_temp=10.0, initial_humidity=50.0, thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
    Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
    Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
]

layers2 = [
    Layer(num=10, name='サイディング', initial_temp=15.0, initial_humidity=50.0, thick=0.016, n_div=5, cond_h_o=22.4, cond_h_i=9.2, cond_m_o=2.00E-11, cond_m_i=3.43E-08, alpha=0.0),
    Layer(num=2, name='通気層', initial_temp=15.0, initial_humidity=50.0, thick=0.025, n_div=1, cond_h_o=9.2, cond_h_i=9.2, cond_m_o=3.43E-08, cond_m_i=3.43E-08, alpha=0.085),
    Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=9.2, cond_h_i=50.0, cond_m_o=3.43E-08, cond_m_i=7.75E-07, alpha=0.0),
    Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.009, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=7.75E-07, alpha=0.0),
    Layer(num=3, name='GW', initial_temp=10.0, initial_humidity=50.0, thick=0.1, n_div=5, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=7.75E-07, cond_m_i=2.00E-11, alpha=0.0),
    Layer(num=5, name='合板', initial_temp=10.0, initial_humidity=80.0, thick=0.012, n_div=6, cond_h_o=50.0, cond_h_i=50.0, cond_m_o=2.00E-11, cond_m_i=7.75E-07, alpha=0.0),
    Layer(num=6, name='石膏ボード', initial_temp=10.0, initial_humidity=50.0, thick=0.042, n_div=5, cond_h_o=50.0, cond_h_i=9.2, cond_m_o=7.75E-07, cond_m_i=2.60E-10, alpha=0.0)
]

@dataclass
class WallType:

    layers: list[Layer]

    @property
    def n_div_total(self):
        return sum([layer.n_div for layer in self.layers])


walltypes =[
    WallType(layers=layers1),
    WallType(layers=layers2)
]

@dataclass
class Wall:

    # 壁種類
    kwtype: int

    # 方位
    direction: int

    # 隣接室1
    room_sideA: int

    # 隣接室2
    room_sideB: int

    # 長辺
    len_long: float

    # 短辺
    len_short: float

    # 高さ
    height: float

    # ブランク
    blank: float

    # 傾斜角
    angle: float

    # アルベド
    albedo: float

    # 日射吸収率
    absorption: float

    # 放射率
    emissivity: float

    # layers
    layers: list[Layer]

    # 面積
    @property
    def area(self):
        return self.len_long * self.len_short - self.blank
    
    # 相当開口面積（αA）, m2
    @property 
    def alpha_a(self) -> list[float]:
        return [layer.alpha * layer.thick * self.len_short for layer in self.layers]
    

walls = [
    Wall(kwtype=1, direction=4, room_sideA=0, room_sideB=1, len_long=7.00, len_short=0.420, height=7.00, blank=0, angle=90, albedo=0.3, absorption=0.9, emissivity=0.9, layers=walltypes[1]),
    Wall(kwtype=2, direction=4, room_sideA=0, room_sideB=1, len_long=7.00, len_short=0.420, height=7.00, blank=0, angle=90, albedo=0.3, absorption=0.9, emissivity=0.9, layers=walltypes[2]),
    Wall(kwtype=2, direction=4, room_sideA=0, room_sideB=1, len_long=7.00, len_short=0.420, height=7.00, blank=0, angle=90, albedo=0.3, absorption=0.9, emissivity=0.9, layers=walltypes[2]),
    Wall(kwtype=2, direction=4, room_sideA=0, room_sideB=1, len_long=7.00, len_short=0.420, height=7.00, blank=0, angle=90, albedo=0.3, absorption=0.9, emissivity=0.9, layers=walltypes[2]),
    Wall(kwtype=2, direction=4, room_sideA=0, room_sideB=1, len_long=7.00, len_short=0.420, height=7.00, blank=0, angle=90, albedo=0.3, absorption=0.9, emissivity=0.9, layers=walltypes[2]),
    Wall(kwtype=2, direction=4, room_sideA=0, room_sideB=1, len_long=7.00, len_short=0.420, height=7.00, blank=0, angle=90, albedo=0.3, absorption=0.9, emissivity=0.9, layers=walltypes[2]),
    Wall(kwtype=1, direction=4, room_sideA=0, room_sideB=1, len_long=7.00, len_short=0.420, height=7.00, blank=0, angle=90, albedo=0.3, absorption=0.9, emissivity=0.9, layers=walltypes[2])
]


# :::::::::::相当開口面積αAの直列合成:::::::::::

#種類  不要（要修正）
for I in range(len(walltypes)):

    # 層
    for K in range(KMTL(I)):
        D1 = 0.0
        D2 = 0.0
        if walls[I][K].alpha > 1.E-7:
            # 階
            for LW, wall in enumerate(walls):
                IW = wall.kwtype
                D1 = D1 + 1 / (wall.alpha_a[K] * wall.alpha_a[K])
                # 通気層高さの合計
                D2 = D2 + walls[LW].height
            if D1 > 1E-12:
                ALP_TOTAL(K)=1/SQRT(D1) 
                HIGHT_CAVITY(K)=D2

@dataclass
class ThinWall:

    # 方位
    direction: int

    # 隣接室
    room_sideB: int

    # 面積, m2
    area: float

    # 熱貫流率, W/(m2 K)
    u_value: float

    # 遮蔽係数
    shading_factor: float

thin_walls = [
    ThinWall(direction=1, room_sideB=4, area=0.09, u_value=6.5, shading_factor=0.01),
    ThinWall(direction=2, room_sideB=4, area=0.09, u_value=6.5, shading_factor=0.01),
    ThinWall(direction=3, room_sideB=4, area=0.09, u_value=6.5, shading_factor=0.01),
    ThinWall(direction=4, room_sideB=4, area=0.09, u_value=6.5, shading_factor=0.01)
]
# 南壁, 西壁, 北壁, 東壁


# 薄壁部位数
NWIN = 4



# :::::::::結露計算部位入力:::::::

# 結露計算質点
# 出力質点数
NSUWC = 4
# 窓部位
NWCN = [1, 2, 3, 4]

ALDC = 5./3600.

# :::::出力座標::::::::::::::::::::::::: 
# 出力質点
NOUT = 12

NOUTD = [[1, 6], [1, 12], [1, 17], [1, 22], [4, 6], [4, 12], [4, 17], [4, 22], [7, 6], [7, 12], [7, 17], [7, 22]]

# :::::平均値出力材料::::::::::::::::::::
# 平均値出力材料
NOUTAV = 9
NOUTAVD = [[1, 1], [1, 3], [1, 4], [4, 1], [4, 3], [4, 4], [7, 1], [7, 3], [7, 4]]

# :::::腐朽部位出力座標::::::::::::::::::::
# 腐朽部位出力質点
NOUTML = 3
NOUTMLD = [[1, 12, 17], [4, 12, 17], [7, 12, 17]]

# ***********発熱*********
# 室数
NSUP = 1
# 室No. [NSUP]
KSUP = [1]

#発熱量, kW, 発湿量, kg/h, [KSUP, 24]
HSUP = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
WSUP = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

for K in range(NSUP):
    for I in range(24):
        HSUP(KSUP(K),I) = HSUP(KSUP(K),I) * 1000.

# ***********浸水率の設定*********

NRAIN = 14

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

DTM = 1. / NDVD * 3600.
PXCOF = 1. / 133322.


K = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
RMDH = [0.0058, 0.0058, 0.045, 0.0824, 0.16, 0.2192, 1.6, 0.0586, 0.2111, 0.963, 0.522, 0.522, 0.021, 0.0339, 0.028, 0.036, 0.0379, 0.1397, 0.1596]
RMDG = [8.69E-12, 8.69E-12, 0.00E+00, 1.58E-10, 1.86E-12, 4.58E-12, 2.93E-11, 4.58E-12, 2.61E-11, 2.29E-11, 4.17E-10, 2.03E-11, 2.79E-11, 5.00E-14, 3.48E-12, 2.61E-12, 1.30E-10, 1.30E-10, 8.94E-12, 1.36E-12]
GCP = [0, 0, 8.40E+02, 1.62E+03, 1.86E+03, 2.84E+03, 9.00E+02, 1.50E+03, 1.31E+03, 8.79E+02, 1.05E+03, 1.05E+03, 8.78E+02, 1.26E+03, 1.26E+03, 8.37E+02, 8.38E+02, 1.26E+03, 1.30E+03]
GMA = [0, 16, 484, 600, 318, 2300, 306, 715, 1095, 1334, 1091, 38.6, 30, 40, 32, 16, 413, 644]
RMDL = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
GCP = [x * y for (x, y) in zip (GCP, GMA)]


# *********CAL DX,NAFX ***************
for I, walltype in enumerate(walltypes):
    K=0

    for J, layer in enumerate(walltype.layers):

        if J == 1:
            NALX(I, J) = layer.n_div
        else:
            NALX(I, J) = NALX(I, J-1) + layer.n_div

    for J, layer in enumerate(walltype.layers):

        if J == 1:
            NAFX(I, J) = 1
        else:
            NAFX(I, J) = NALX(I, J-1)+1

    NX(I) = walltype.n_div_total

# *******SAT   SINITIAL VALUE *******

for I in range(len(walls)):
    SAT(I, 1) = TMPOC
    SAT(I, 2) = TMPOC

TEMPO(1) = TMPOC
SJD(1) = 0.				
SJS(1) = 0.				
SJN(1) = 0.				
RHOD(1) = walltypes[0][0].initial_humidity
D1 = TMPOC + ATP
RH0 = walltypes[0][0].initial_humidity

FS0, VP = GOFF(D1)

XOD(1) = VP * RH0 * 0.01		
TEMPO(2) = TMPOC
SJD(2) = 0.
SJS(2) = 0.
SJN(2) = 0.
XOD(2) = VP * RH0 * 0.01
RHOD(2) = walltypes[0][0].initial_humidity

# /////////INITIAL CONDITION//////////////

for L in range(len(walls)):
    IW=walls[L].kwtype
    for K in range(KMTL(IW)):
        L1 = NAFX(IW, K)
        L2 = NALX(IW, K)
        L5 = walltypes[IW][K].num
        TP1 = walltypes[IW][K].initial_temp + ATP
        RH0 = walltypes[IW][K].initial_humidity
        WP0 = SATUWPT(TP1)
        WP = REWPT(RH0, TP1, RG)
   
        FS0, VP = GOFF(TP1)
		
        WP1=WP
        WP2=WP0 + WP1
        QQ(L,K)=0.

        for I in range(L1,L2):
            WPW(L,I) = WP2
            WPU(L,I) = WP1
            WPS(L,I) = WP0
            TMP(L,I) = TP1
            ATMPQ(L,I) = TP1 - ATP
            TMPC(L,I) = TP1 - ATP
            RH(L,I) = RH0
            XN(L,I) = VP * RH0 * 0.01
            XM(L,I) = VP * RH0 * 0.01
            AXN(L,I) = VP * RH0 * 0.01
            AWPU(L,I) = WP1
            ATMP(L,I) = TP1
            HWPU(L,I) = WP1
            HTMP(L,I) = TP1
            RHAVD(L,I) = 0.
            WGTAVD(L,I) = 0.
            XNAVD(L,I) = 0.
            TPAVD(L,I) = 0.
            WD = AHGANS(rhm=RH0, ml0=L5)
            WGT(L,I) = WD
            WLOSS(L,I) = 0.

        for I in range(1, 2):
            QS(L,K,I) = 0.

for I in range(NROOM):
    TMPR(I) = TMPIC
    ATMPR(I) = TMPIC + ATP
    HTMPR(I) = TMPIC + ATP
    RHR(I) = RHI
    XR(I) = VP * RHI * 0.01
    XRM(I) = VP * RHI * 0.01
    RH0 = RHI
    TP1 = TMPIC + ATP
    WP1 = REWPT(RH0, TP1, RG)
    WPUR(I) = WP1
    AWPUR(I) = WP1

for L in range(NOUTAV)
    I = NOUTAVD(L,1)
    J = NOUTAVD(L,2)
    IW = walls[I].kwtype
    L1 = NAFX(IW,J)
    WOUTAV(L) = WGT(I,L1)

#DO 444 NYEAR=1,LYEAR
# 三浦コメント：ここから１年のループ計算
for NYEAR in range(LYEAR):

    IF(NYEAR.GT.1) THEN
        READ(7,210)D1
        READ(13,210)D1
        READ(13,210)D1
    END IF
    KDAY=0				
    IF(NYEAR.EQ.1)THEN				
        CALL YOMI(MON1,MDAY1,KDAY)				
        MONF=MON1			
        LMONL=LMON			
        IF(LYEAR.GT.1)LMONL=12			
    ELSE			
        MONF=1			
        LMONL=12
        IF(LYEAR.EQ.NYEAR)LMONL=LMON			
    END IF
    IMQ=0

    # 三浦コメント：ここから１月のループ計算
    DO 999 MON=MONF,LMONL
    
        DAY2=MDAY(MON)
        IF(NYEAR.EQ.LYEAR.AND.MON.EQ.LMON)THEN				
            DAY2=LDAY				
        ELSE				
            DAY2=MDAY(MON)				
        END IF
    
        IF(NYEAR.EQ.1.AND.MON.EQ.MON1)THEN				
            DAY1=MDAY1				
        ELSE				
            DAY1=1				
        END IF
    			
        # 三浦コメント：ここから１日のループ計算
        DO 888 DAY=DAY1,DAY2
    
            PRINT *, '  NYEAR=',NYEAR,'  MON=',MON,' DAY=',DAY,' IT1-3=',IT2,IT3
            KDAY=KDAY+1
            100 FORMAT(F4.1,F3.1,3F3.0,2F3.3,2F4.3)				
            101 FORMAT(24F3.0)				
            102 FORMAT(24I3)
            
            READ(7,*)MOJI

            DO  I=2,25
                READ(7,100)TEMPO(I),XOD(I),SJD(I),SJS(I),SJN(I),D1,D2,D3,D4
                SJD(I)=SJD(I)*1.16                                           !MKS
                SJS(I)=SJS(I)*1.16                                           !MKS
                SJN(I)=SJN(I)*1.16                                           !MKS
                XOD(I)=XOD(I)*0.001				
                D1=TEMPO(I)+ATP				
                RH0=100.				
                FS0, VP = GOFF(D1)				
                XS = FUNCX(FS0,RH0)				
                RHOD(I)=100.0*(XOD(I)/XS)      				
                XOD(I)=VP*RHOD(I)*0.01                                    !CHANGE
            END DO

            !**   外部風速、雨量の読み込み***

            DO I=1,24
                READ(13,*)VELOUT(I),RAIN(I),VELDIR(I)
            END DO
            
            # *******************************

            # 三浦コメント：NWNは壁体の数。何か壁体ごとに処理をしているのか？
            DO LW=1,len(walls)
                AKASYA1=walls[LW].angle
                NORIENT=walls[LW].direction
                RLF1=walls[LW].albedo
                AS1=walls[LW].absorption
                EMI1=walls[LW].emissivity
                IF(NORIENT.LE.4.AND.NORIENT.GE.1)THEN
                    data = {'BKAKU': BKAKU, 'W':  W, 'DV4': DV4}
                    data = SOLOCT(KDAY,HIDO,HKEIDO,HOI,AKASYA1,NORIENT,LW, data)
                    BKAKU, W, DV4 = data['BKAKU'], data['W'], data['DV4']
                    ms = MS(TEMPO=TEMPO, SJD=SJD, SJS=SJS, SJN=SJN, SAT=SAT, SJIN=SJIN, WSJIN=WSJIN, WSJD=WSJD)
                    ss = SS(W=W, DV4=DV4, BKAKU=BKAKU)
                    ms, ss = SATCAL(RLF1,AS1,EMI1,LW, ms, ss)
                    TEMPO, SJD, SJS, SJN, SAT, SJIN, WSJIN, WSJD, W, DV4, BKAKU = ms.TEMPO, ms.SJD, ms.SJS, ms.SJN, ms.SAT, ms.SJIN, ms.WSJIN, ms.WSJD, ss.W, ss.DV4, ss.BUKAKU
                END IF
            END DO

            DO I=1,NWIN
                J=I+len(walls)
                NORIENT=thin_walls[I].direction
                AKASYA1=90.
                data = {'BKAKU': BKAKU, 'W':  W, 'DV4': DV4}
                data = SOLOCT(KDAY,HIDO,HKEIDO,HOI,AKASYA1,NORIENT,J, data)	
                BKAKU, W, DV4 = data['BKAKU'], data['W'], data['DV4']
                ms = MS(TEMPO=TEMPO, SJD=SJD, SJS=SJS, SJN=SJN, SAT=SAT, SJIN=SJIN, WSJIN=WSJIN, WSJD=WSJD)
                ss = SS(W=W, DV4=DV4, BKAKU=BKAKU)
                ms, ss = SATCAL(RLF1,AS1,EMI1,J)
                TEMPO, SJD, SJS, SJN, SAT, SJIN, WSJIN, WSJD, W, DV4, BKAKU = ms.TEMPO, ms.SJD, ms.SJS, ms.SJN, ms.SAT, ms.SJIN, ms.WSJIN, ms.WSJD, ss.W, ss.DV4, ss.BUKAKU
                DO K=2,25
                    WSIN(I,K)=WSJIN(J,K)*thin_walls[I].area * thin_walls[i].shading_factor    !拡散成分
                    WSIND(I,K)=WSJD(J,K)*thin_walls[I].area * thin_walls[i].shading_factor     ! 直達成分
                END DO
            END DO

            DO LW=1,len(walls)
                IW=walls[LW].kwtype
                DO I=1,NX(IW)
                    WGTDIS(LW,I,MON,DAY)=WGT(LW,I)
                END DO
            END DO

            # *************時間ループ*****************
            # 三浦コメント：ここから１時間のループ計算
            DO 777 IM=1,24

                # /////KETURO/////////////
                DO LW=1,NWIN
                    CWVDAY(LW,IM)=0.
                END DO

                # /////浸水率の計算（ワイブル分布を考慮した風速の関数）////////
                IF(NRAIN.GT.0)THEN
                    DO K=1,NRAIN
                        I=NRAINPOINT(K,1)
                        J=NRAINPOINT(K,2)
                        SUMRAIN(I,J,IM)=0.
                        HIGHT_SER=(I-0.5)*walls[LW].height+0.7  !高さを厚壁中央に設定 基礎高0.7m
                        VEL_MOD=VELOUT(IM)*(HIGHT_SER/6.5)**0.25    !風速の高さ補正（基準風速h=6.5m）0.25乗
                        D1=0.
                        D2=0.
                        D3=0.
                        D4=0.

                        IF(VELOUT(IM).GE.0.1)THEN
                            DO I=1,30
                                D1=I*0.1*VEL_MOD 
                                
                                # D2=ACOFRAIN(K,1)*D1*COS((VELDIR(IM)-(180+90*(NORI(NRAINPOINT(K,1),1)-1)))/180*3.1415)*1.5*0.2   !風向の考慮 ASHRAEの暴露係数(SEVERE)を利用
                                D2=ACOFRAIN(K,1)*D1*COS((VELDIR(IM)-(180+90*(walls[NRAINPOINT(K,1)].direction -1)))/180*3.1415)*1.5*0.2   !風向の考慮 ASHRAEの暴露係数(SEVERE)を利用
                                
                                IF(D2.LT.0.)D2=0.
                                D1=3.14*0.5*D1/VEL_MOD/VEL_MOD*EXP(3.14*0.25*(D1/VEL_MOD)**2)      !f(v)
                                D3=D3+D1*D2                                                                ! sig(f(v)*g(v))
                                D4=D4+D1                                                                   ! sig(f(v))
                            END DO
                            COFRAIN(K,1)=D3/D4
                        ELSE 
                            COFRAIN(K,1)=0
                        END IF
                    END DO
                END IF
                
                # //////////////////////////////////////////////////////////
                IMM=IM+1 
                IMQ=IMQ+1
                #   室内側壁及び床面の透過日射量の計算(未完成)
                WSINSUM=0.
                WSINDSUM=0.

                DO  I=1,NWIN
                    WSINSUM=WSINSUM+WSIN(I,IM)     !拡散成分（各面拡散）
                    WSINDSUM=WSINDSUM+WSIND(I,IM)  ! 直達成分（床面入射）
                END DO

                DO LW=1,len(walls)
                    IW=walls[LW].kwtype
                    K=KMTL(IW)
                END DO

                IF(NYEAR.GE.1)THEN
                    IF(IM.EQ.1)THEN
                        I=24
                        ELSE
                        I=IM-1
                    END IF
                    D2=SUMRAIN(NRAINPOINT(1,1),NRAINPOINT(1,2),I)
                    D1=SWJRAIN(NRAINPOINT(1,1),NRAINPOINT(1,2))*3600.
                    D3=RN(NRAINPOINT(1,1),NRAINPOINT(1,2))
                    WRITE(8,*)NYEAR,MON,DAY,IM-1,RAIN(I),VELOUT(I),D1,D2,D3& !   Year Mon Day IM 降雨量(mm/h) 風速(m/s) 瓦浸水量(kg/m2h) 防水層浸水量(kg/m2h) 空気層水分量(kg/m2) 
                            ,(TMPR(K),RHR(K),XR(K),K=1,NROOM),(TMPC(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(XN(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(RH(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT)&
                            ,(WGT(NOUTD(K,1),NOUTD(K,2)),K=1,NOUT),(SATDV(K),K=1,len(walls)),(WOUTAV(K),K=1,NOUTAV),TEMPO(I),XOD(I),QQ(NQOUT1,NQOUT2),QQ(NQOUT1,NQOUT3)
                END IF

                # 三浦コメント：おそらくここから1時間の分割す数に応じたループがはじまる。
                DO 666 NF=1,NDVD
    
                    DO I=1,len(walls)
                        IF(walls[I].direction.GE.1.AND.walls[I].direction.LE.4)THEN
                            D1=(SAT(I,IMM)-SAT(I,IM))*DTM/3600
                            SATDV(I)=SAT(I,IM)+(NF-1)*D1
                            END IF
                        IF(walls[I].direction.EQ.5)THEN
                            D5=(TEMPO(IMM)-TEMPO(IM))*DTM/3600     !  床下温度（低下率0.7）
                            D5=TEMPO(IM)+(NF-1)*D5
                            SATDV(I)=TMPR(walls[I].room_sideB)-0.7*(TMPR(walls[I].room_sideB)-D5)
                            SAT(I,IM)=SATDV(I)
                        END IF
                        SATA(I)=SATDV(I)+ATP
                    END DO
    
                    D5=(TEMPO(IMM)-TEMPO(IM))*DTM/3600				
                    D6=(XOD(IMM)-XOD(IM))*DTM/3600			
                    D7=(RHOD(IMM)-RHOD(IM))*DTM/3600				
                    TPO=TEMPO(IM)+(NF-1)*D5				
                    XO=XOD(IM)+(NF-1)*D6				
                    RHO=RHOD(IM)+(NF-1)*D7	

                    # ***************************
                    #   外気の化学ポテンシャルの計算
                    TMPOC=TPO
                    TP1=TMPOC+ATP
                    WP1 = REWPT(RHO, TP1, RG)
                    WPUO=WP1
                    TMPO=TMPOC+ATP
                    D1, D2 = CALDGDU(RHO,TMPO,RG)
                    DGDUO=1.2*PXCOF*D1
                    DGDTO=1.2*PXCOF*D2

                    IF(NRAIN.GT.0)THEN
                        DO K=1,NRAIN
                            I=NRAINPOINT(K,1)
                            J=NRAINPOINT(K,2)
                            IF(COFRAIN(K,1).GE.0.)THEN
                                SWJRAIN(I,J)=0.01*RAIN(IM)*COFRAIN(K,1)/3600.     !材表面への浸水量kg/m2S 
                            END IF
                        END DO
                    END IF

                    # ////////PUT IN PREBIOUS VALUE///////////
                    DO I=1,NROOM
                        HWPUR(I)=WPUR(I)
                        HTMPR(I)=ATMPR(I)
                        XRM(I)=XR(I)
                    END DO
                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype

                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                HWPU(LW,I)=WPU(LW,I)
                                HTMP(LW,I)=TMP(LW,I)
                                XM(LW,I)=XN(LW,I)
                                HRN(LW,I)=RN(LW,I)
                            END DO
                        END DO
                    END DO

                    !///////////////COMBAIN METHIOD//////////////////////////
                    IT3=0
                    810 AMAX3=0.
                    AMAX4=0.
                    !**** CULCURATION CONDUCTANCE****

                    DO LW=1,len(walls)

                        IW=walls[LW].kwtype

                        DO K=1,KMTL(IW)

                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num

                            DO I=L1,L2

                                IF(RH(LFW,I).GT.99.98)CYCLE

                                IF(RH(LW,I).GT.90.0.AND.L5.EQ.5)THEN
                                    D1=RH(LW,I)
                                    D2=TMP(LW,I)
                                    # Liquid conductivity (kg/msJ/kg)
                                    D3 = DIFF(L5,D1,D2,RG)
                                    RMDL(L5)=D3
                                    ELSE 
                                    RMDL(L5)=0.
                                END IF

                                !*****************************
                                RH1=RH(LW,I)
                                TP1=TMP(LW,I)
                                DGDU, DGDT = CALDGDU(RH1,TP1,RG)

                                D3=walls[LW].area

                                IF(I.EQ.L1)THEN
                                    ALD(LW,K,1)=walls[IW][K].cond_m_o*DGDU*D3
                                    ALDT(LW,K,1)=walls[IW][K].cond_m_o*DGDT*D3
                                    ALDP(LW,K,1)=walls[IW][K].cond_m_o*D3
                                    ALF(LW,K,1)= walls[IW][K].cond_h_o *D3
                                END IF

                                IF(I.EQ.L2)THEN
                                    ALD(LW,K,2)=walls[IW][K].cond_m_i*DGDU*D3
                                    ALDT(LW,K,2)=walls[IW][K].cond_m_i*DGDT*D3
                                    ALDP(LW,K,2)=walls[IW][K].cond_m_i*D3
                                    ALF(LW,K,2)= walls[IW][K].cond_h_i *D3
                                END IF

                                IF(L5.EQ.2)THEN   !空気層
                                    DGDUQ(LW,K)=1.2*DGDU*PXCOF
                                    DGDTQ(LW,K)=1.2*DGDT*PXCOF 
                                    CYCLE
                                END IF
                                !
                                RMDX(LW,I)=RMDH(L5)*D3/ walltypes[IW].layers[K].dx

                                WP=WPU(LW,I)
                                TP=TMP(LW,I)
                                D1=GMA(L5)

                                IF(L5.NE.2)THEN
                                    IF(L5.NE.3)THEN
                                        D2 = CALDPDU(WP,TP,D1,L5,RG,ROW)
                                        DPDU(LW,I)=D2
                                    ELSE
                                        IF(RH(LW,I).LE.99.98)THEN
                                            D2 = CALDPDU(WP,TP,D1,L5,RG,ROW)
                                            DPDU(LW,I)=D2
                                        END IF
                                    END IF

                                END IF

                                D2=RMDG(L5)
                                D1=RH(LW,I)*0.01

                                IF(L5.EQ.5)D2=1.87E-11*D1**2.4019*EXP(-0.78864*(1-D1**1.1471))  !*3.45   !  Moisture conductivity (kg/msPa)  ASHRAE
                                ADWGX(LW,I)=D2*DGDU/ walltypes[IW].layers[K].dx *D3
                                ADWLX(LW,I)=RMDL(L5)*D3/ walltypes[IW].layers[K].dx
                                ADWX(LW,I)=ADWGX(LW,I)+ADWLX(LW,I)
                                ADTGX(LW,I)=D2*DGDT*D3/ walltypes[IW].layers[K].dx
                                ADTLX(LW,I)=0.*RMDL(L5)*DPDU(LW,I)*D3/ walltypes[IW].layers[K].dx
                                ADTX(LW,I)=ADTGX(LW,I)+ADTLX(LW,I)
                            END DO

                            DO I=L1+1,L2
                                DWGX(LW,I)=(ADWGX(LW,I)+ADWGX(LW,I-1))*0.5
                                DWLX(LW,I)=(ADWLX(LW,I)+ADWLX(LW,I-1))*0.5
                                DWX(LW,I)=(ADWX(LW,I)+ADWX(LW,I-1))*0.5
                                DTGX(LW,I)=(ADTGX(LW,I)+ADTGX(LW,I-1))*0.5
                                DTLX(LW,I)=(ADTLX(LW,I)+ADTLX(LW,I-1))*0.5
                                DTX(LW,I)=(ADTX(LW,I)+ADTX(LW,I-1))*0.5
                            END DO
                        END DO
                    END DO

                    DO I=1,NROOM
                        RH1=RHR(I)
                        TP1=ATMPR(I)
                        DGDU, DGDT = CALDGDU(RH1,TP1,RG)
                        DGDUR(I)=1.2*DGDU*PXCOF
                        DGDTR(I)=1.2*DGDT*PXCOF
                    END DO

                    !////////////CAL_TEMP BY OVER RELAXATION METHOD///////////
                    IT1=0                                                                     
                    ITQ=0
                    IT4=IT4+1

                    !***** conbained  temp and ventiration****
    
                    # GOTO 文でここに戻ってくる
                    217 CONTINUE

                    ITQ=ITQ+1
                    IF(ITQ.GT.100)THEN
                        PRINT *, '  no convagence   ventiration'
                        GO TO 218
                    END IF
                    !*****************************************
                    
                    # GOTO 文でここに戻ってくる
                    510 AMAX1=0. 

                    # ******ROOM TEMP*****
                    DO I=1,NROOM
    
                        UHEN=1300.* rooms[I].volume /DTM*HTMPR(I)
                        SAHEN=1300.* rooms[I].volume /DTM

                        DO LW=1,len(walls)     !   厚壁からの取得熱量

                            IW=walls[LW].kwtype

                            IF(walls[LW].room_sideB.EQ.I)THEN
                                UHEN=UHEN+ALF(LW,KMTL(IW),2)*TMP(LW,NX(IW))!+RW*ALDP(LW,KMTL(IW),2)*(XN(LW,NX(IW))-XR(I)) ! SENSIBLE
                                SAHEN=SAHEN+ALF(LW,KMTL(IW),2)
                            END IF

                            IF(walls[LW].room_sideA.EQ.I)THEN
                                UHEN=UHEN+ALF(LW,1,1)*TMP(LW,1)!+RW*ALDP(LW,1,1)*(XN(LW,1)-XR(I)) ! SENSIBLE
                                SAHEN=SAHEN+ALF(LW,1,1)
                            END IF

                        END DO

                        DO LW=1,NWIN     !   薄壁からの取得熱量
                            IF(thin_walls[LW].room_sideB.EQ.I)THEN
                                UHEN=UHEN+ thin_walls[LW].u_vlue *thin_walls[LW].area*TMPO
                                SAHEN=SAHEN+thin_walls[LW].u_vlue*thin_walls[LW].area
                            END IF
                        END DO

                        IF(rooms[I].n_open .NE.0)THEN

                            DO J=1, rooms[I].n_open
                                IF(rooms[I].flows[J].idx_from .EQ.0)THEN
                                    UHEN=UHEN+1300.* rooms[I].flows[rooms[I].flows[J].idx_from].volume * TMPO/3600
                                ELSE 
                                    UHEN=UHEN+1300.* rooms[I].flows[rooms[I].flows[J].idx_from].volume * ATMPR(rooms[I].flows[J].idx_from)/3600      !  換気負荷の積算
                                END IF

                                SAHEN=SAHEN+1300.* rooms[I].flows[rooms[I].flows[J].idx_from].volume )/3600.
                            END DO
                        END IF

                        UHEN=UHEN+HSUP(I,IM)
                        ATMPR(I)=UHEN/SAHEN
                        TMPR(I)=ATMPR(I)-ATP

                        # ****Constant temp********************
                        IF(rooms[I].is_calc == False)THEN
                            D1 = rooms[i].t_ave
                            D2 = rooms[i].t_amp
                            D3 = rooms[i].rh_ave
                            TRR, RR, XRR = ROOM(KDAY,D1,D2,D3)
                            TMPR(I)=TRR
                            ATMPR(I)=TMPR(I)+ATP
                        END IF

                    END DO

                    DO LW=1,len(walls)

                        IW=walls[LW].kwtype

                        DO K=1,KMTL(IW)

                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num

                            !*******  VENTED CAVITY  ******

                            IF(L5.EQ.2)THEN

                                I=L1
                                D1=1300.* walltypes[IW].layers[K].dx *walls[LW].area/DTM
                                D2=ALF(LW,K,1)+ALF(LW,K,2)+RW*(ALDT(LW,K,1)+ALDT(LW,K,2))
                                D4=(ALF(LW,K,1)+RW*ALDT(LW,K,1))*TMP(LW,I-1)+(ALF(LW,K,2)+RW*ALDT(LW,K,2))*TMP(LW,I+1)
                                D5=RW*(ALD(LW,K,1)*(WPU(LW,I-1)-WPU(LW,I))+ALD(LW,K,2)*(WPU(LW,I+1)-WPU(LW,I)))

                                IF(LW.EQ.1)THEN
                                    UHEN=D1*HTMP(LW,I)+D4+D5+1300.*QQ(LW,K)*TMPO
                                ELSE
                                    UHEN=D1*HTMP(LW,I)+D4+D5+1300.*QQ(LW,K)*TMP(LW-1,I)
                                END IF

                                SAHEN=D1+D2+1300.*QQ(LW,K)
                                TMP(LW,I)=UHEN/SAHEN
                                CYCLE
                            END IF

                            # *******************************
                            DO I=L1,L2

                                IF(I.EQ.L1)THEN               !*********端部外

                                    D1=0.5*GCP(L5)* walltypes[IW].layers[K].dx *walls[LW].area/DTM

                                    IF(I.EQ.1)THEN                                               !変更
                                        IF(walls[LW].room_sideA.EQ.0)THEN                                         !             隣室外気
                                            D2=ALF(LW,K,1)+RMDX(LW,I+1)+RW*(ALDT(LW,K,1)+DTGX(LW,I+1))
                                            D3=CPL*(DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
                                            D4=ALF(LW,K,1)*SATA(LW)+RW*ALDT(LW,K,1)*TMPO+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
                                            D5=RW*(ALD(LW,K,1)*(WPUO-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
                                            UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I+1)+QS(LW,K,1)
                                            SAHEN=D1+D2+D3
                                        ELSE                                                              !             隣室室内
                                            D2=ALF(LW,K,1)+RMDX(LW,I+1)+RW*(ALDT(LW,K,1)+DTGX(LW,I+1))
                                            D3=CPL*(DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
                                            D4=ALF(LW,K,1)*ATMPR(walls[LW].room_sideA)+RW*ALDT(LW,K,1)*ATMPR(walls[LW].room_sideA)+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
                                            D5=RW*(ALD(LW,K,1)*(WPUR(walls[LW].room_sideA)-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
                                            UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I+1)+QS(LW,K,1)
                                            SAHEN=D1+D2+D3
                                        END IF
                                    ELSE   
                                        D2=ALF(LW,K,1)+RMDX(LW,I+1)+RW*(ALDT(LW,K,1)+DTGX(LW,I+1))
                                        D3=CPL*(DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
                                        D4=(ALF(LW,K,1)+RW*ALDT(LW,K,1))*TMP(LW,I-1)+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
                                        D5=RW*(ALD(LW,K,1)*(WPU(LW,I-1)-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
                                        UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I+1)+QS(LW,K,1)
                                        SAHEN=D1+D2+D3
                                    END IF

                                    TMP(LW,I)=UHEN/SAHEN

                                END IF

                                IF(I.GT.L1.AND.I.LT.L2)THEN   !*********実質部
                                    D1=GCP(L5)* walltypes[IW].layers[K].dx *walls[LW].area/DTM
                                    D2=CPL*(DWLX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DWLX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I))&
                                        +DTLX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+DTLX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I)))
                                    D3=RMDX(LW,I)+RMDX(LW,I+1)+RW*(DTGX(LW,I)+DTGX(LW,I+1))
                                    D4=(RMDX(LW,I)+RW*DTGX(LW,I))*TMP(LW,I-1)+(RMDX(LW,I+1)+RW*DTGX(LW,I+1))*TMP(LW,I+1)
                                    D5=RW*(DWGX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DWGX(LW,I+1)*(WPU(LW,I+1)-WPU(LW,I)))
                                    D6=D2*(TMP(LW,I-1)+TMP(LW,I+1))
                                    D7=D2*2
                                    UHEN=D1*HTMP(LW,I)+D4+D5+D6
                                    SAHEN=D1+D3+D7
                                    TMP(LW,I)=UHEN/SAHEN
                                END IF

                                IF(I.EQ.L2)THEN               !*********端部内
                                    D1=0.5*GCP(L5)* walltypes[IW].layers[K].dx *walls[LW].area/DTM

                                    IF(I.EQ.NX(IW))THEN
                                        D2=ALF(LW,K,2)+RMDX(LW,I)+RW*(ALDT(LW,K,2)+DTGX(LW,I))
                                        D3=CPL*(DWLX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DTLX(LW,I)*(TMP(LW,I-1)-TMP(LW,I)))
                                        D4=(ALF(LW,K,2)+RW*ALDT(LW,K,2))*ATMPR(walls[LW].room_sideB)+(RMDX(LW,I)+RW*DTGX(LW,I))*TMP(LW,I-1)
                                        D5=RW*(ALD(LW,K,2)*(WPUR(walls[LW].room_sideB)-WPU(LW,I))+DWGX(LW,I)*(WPU(LW,I-1)-WPU(LW,I)))
                                        UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I-1)+QS(LW,K,2)
                                        SAHEN=D1+D2+D3
                                    ELSE
                                        D2=ALF(LW,K,2)+RMDX(LW,I)+RW*(ALDT(LW,K,2)+DTGX(LW,I))
                                        D3=CPL*(DWLX(LW,I)*(WPU(LW,I-1)-WPU(LW,I))+DTLX(LW,I)*(TMP(LW,I-1)-TMP(LW,I)))
                                        D4=(ALF(LW,K,2)+RW*ALDT(LW,K,2))*TMP(LW,I+1)+(RMDX(LW,I)+RW*DTGX(LW,I))*TMP(LW,I-1)
                                        D5=RW*(ALD(LW,K,2)*(WPU(LW,I+1)-WPU(LW,I))+DWGX(LW,I)*(WPU(LW,I-1)-WPU(LW,I)))
                                        UHEN=D1*HTMP(LW,I)+D4+D5+D3*TMP(LW,I-1)+QS(LW,K,2)
                                        SAHEN=D1+D2+D3
                                    END IF

                                    TMP(LW,I)=UHEN/SAHEN

                                END IF

                            END DO

                        END DO   !K

                    END DO   !LW

                    !/////////OVER RELAXATION ////// 
                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                TMP(LW,I)=(TMP(LW,I)-ATMP(LW,I))*OMG+ATMP(LW,I)
                                D8=TMP(LW,I)-ATMP(LW,I)
                                IF(ABS(D8).GT.AMAX1)AMAX1=ABS(D8)
                                ATMP(LW,I)=TMP(LW,I)
                            END DO
                        END DO
                    END DO

                    # /////////JUDGEMENT CONVERGENCE//////

                    IT1=IT1+1
                        
                    IF(IT1.GT.MXIT)THEN
                        WRITE(8,*) ' HASSAN TEMP IT1= ',IT1
                        WRITE(8,*)(TMP(1,I),I=1,NX(1))
                        WRITE(8,*)' 実質部 '
                        201 FORMAT(A)
                        # ****************************************
                        GO TO 555
                    END IF

                    IF(AMAX1.GT.EPS1) GO TO 510

                    # *************換気計算収束判定********************
                    IF(ITQ.GT.1)THEN
                        D1=0
                        D2=0
                        DO LW=1,len(walls)
                            IW=walls[LW].kwtype
                            DO K=1,KMTL(IW)
                                L1=NAFX(IW,K)
                                L2=NALX(IW,K)
                                L5=walltypes[IW][K].num
                                IF(L5.EQ.2.AND.walls[IW][K].alpha.GT.0.)THEN
                                    D1=ABS(TMP(LW,L1)-ATP-ATMPQ(LW,L1))
                                    IF(D1.GT.D2)D2=D1
                                END IF
                                ATMPQ(LW,L1)=TMP(LW,L1)-ATP
                            END DO
                        END DO
                        IF(D2.LT.EPS1)GO TO 218
                    END IF

                    !*************換気量の算出(Q=m3/s)****************
                    DO I=1,NWTYPE
                        DO K=1,KMTL(I)
                            D1=0.
                            L1=NAFX(I,K) !!!
                            L2=NALX(I,K)
                            L5=walltypes[I][K].num

                            DO LW=1,len(walls)
                                IW=walls[LW].kwtype
                                IF(L5.EQ.2.AND.walls[IW][K].alpha.GT.0)THEN
                                    D1=D1+TMP(LW,L1)*walls[LW].height
                                END IF
                            END DO

                            IF(HIGHT_CAVITY(K).GT.0.1)THEN
                                TEMP_CAVITY(K)=D1/HIGHT_CAVITY(K)
                                GMAQ1=353.25/(TEMP_CAVITY(K))
                                GMAO=353.25/(TMPO)
                                PQ1=GMAO-GMAQ1                       
                                D1=0.5*ABS(PQ1)*HIGHT_CAVITY(K)*9.8                !中性帯　高さ中央 0.5, 重力加速度 g(ρ-ρ)
                            END IF

                            DO LW=1,len(walls)
                                IW=walls[LW].kwtype
                                IF(L5.EQ.2.AND.walls[IW][K].alpha.GT.0)QQ(LW,K)=DSQRT(D1)*4.*ALP_TOTAL(K)
                            END DO
                        END DO
                    END DO

                    # 三浦コメント：ここで強制的に前の方に戻される。
                    # 戻された部分からここまでの間に、何かの条件を満たせばGOTO218で飛び、ループを抜けられる。
                    GOTO 217

                    218 CONTINUE

                    # ************END VENT****************************

                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                TP=TMP(LW,I)
                                SW = CALL SATUWPT(TP)
                                WPS(LW,I)=SW
                            END DO
                        END DO
                    END DO

                    # ////////////CAL_WATER POTENTIAL BY OVER RELAXATION METHOD/////////
                    IT2=0
                    CWMAX=2.8E-2

                    ITC=0
                    710 AMAX2=0.D0 

                    # ******ROOM MOISTURE*****
                    # 室内湿度計算
                    DO I=1,NROOM

                        UHEN=DGDUR(I)* rooms[I].volume /DTM*HWPUR(I)
                        SAHEN=DGDUR(I)* rooms[I].volume /DTM

                        DO LW=1,len(walls)     !   厚壁からの放湿
                            IW=walls[LW].kwtype
                            IF(walls[LW].room_sideB.EQ.I)THEN
                                UHEN=UHEN+ALD(LW,KMTL(IW),2)*WPU(LW,NX(IW))+ALDT(LW,KMTL(IW),2)*(TMP(LW,NX(IW))-ATMPR(I))
                                SAHEN=SAHEN+ALD(LW,KMTL(IW),2)
                            END IF
                            IF(walls[LW].room_sideA.EQ.I)THEN
                                UHEN=UHEN+ALD(LW,1,1)*WPU(LW,1)+ALDT(LW,1,1)*(TMP(LW,1)-ATMPR(I))
                                SAHEN=SAHEN+ALD(LW,1,1)
                            END IF
                        END DO

                        DO LW=1,NWIN   !窓ガラスの結露
                            IF(thin_walls[LW].room_sideB.EQ.I)THEN
                                IF(MCW(LW).EQ.1)THEN
                                    T1=ATMPR(I)
                                    WP=WPUR(I)
                                    RH0 = WPTRE(WP,T1,RG)
                                    FS, VP = GOFF(T1)
                                    X1 = FUNCX(FS,RH0)
                                    KMC=MCW(LW)
                                    T1=ATMPR(I)-thin_walls[LW].u_vlue/9.3*(ATMPR(I)-TMPO)
                                    KMC = CWIF1(T1,X1,XSUT,KMC)
                                    MCW(LW)=KMC
                                    UHEN=UHEN+ALDC*(XSUT-X1)*walls[LW].area
                                END IF
                            END IF
                        END DO

                        IF(rooms[I].n_open .NE.0)THEN
                            DO J=1, rooms[I].n_open
                                D1= rooms[I].flows[rooms[I].flows[J].idx_from].volume /3600.
                                IF(rooms[I].flows[J].idx_from .EQ.0)THEN
                                    UHEN=UHEN+DGDUO*WPUO*D1+D1*(DGDTO*TMPO-DGDTR(I)*ATMPR(I)) 
                                ELSE 
                                    UHEN=UHEN+DGDUR(rooms[I].flows[J].idx_from)*WPUR(rooms[I].flows[J].idx_from)*D1+D1*(DGDTR(rooms[I].flows[J].idx_from)*ATMPR(rooms[I].flows[J].idx_from)-DGDTR(I)*ATMPR(I))      !  換気負荷の積算
                                END IF
                                SAHEN=SAHEN+DGDUR(I)*D1
                            END DO
                        END IF

                        UHEN=UHEN+WSUP(I,IM)/3600.
                        WPUR(I)=UHEN/SAHEN
                        D1=(WPUR(I)-AWPUR(I))*OMG+AWPUR(I)

                        D8=D1-AWPUR(I)
                        AWPUR(I)=D1

                        # ****Constant humidity********************
                        IF(rooms[I].is_calc == False)THEN
                            D1 = rooms[i].t_ave
                            D2 = rooms[i].t_amp
                            D3 = rooms[i].rh_ave
                            TRR, RR, XRR = CALL ROOM(KDAY,D1,D2,D3)
                            D1=RR
                            D2=TRR+ATP
                            WP = REWPT(D1, ATP, RG)  !  ???? D2 ?
                            AWPUR(I)=WP
                            WPUR(I)=WP
                            RHR(I)=RR
                            D8=0.
                        END IF

                        IF(ABS(D8).GT.AMAX2)AMAX2=ABS(D8)
                    END DO

                    # *************************************
                    # ************水膜の水分保持量及び吸水量の計算**********
                    IF(NRAIN.GT.0)THEN
                        DO K=1,NRAIN
                            I=NRAINPOINT(K,1)
                            J=NRAINPOINT(K,2)
                            IW=walls[I].kwtype
                            L1=walltypes[IW][NRAINPOINT(k,3)].num
                            L5=walltypes[IW][L1].num

                            IF(J.EQ.NALX(IW,NRAINPOINT(K,3)))THEN     !  水膜との隣接質点の選択
                                D2=WPU(I,J-1)
                                D3=XM(I,J+1)
                            ELSE
                                D2=WPU(I,J+1)
                                IF(J.EQ.1)THEN
                                    D3=XO                               !外装材表面、外気との収支
                                ELSE
                                    D3=XM(I,J-1)
                                END IF
                            END IF

                            IF(SWJRAIN(I,J).GE.0.OR.RN(I,J).GE.0)THEN
                                T1=TMP(I,J)
                                FS, VP = GOFF(T1)
                                RH1=RH(I,J)
                                DGDU, DGDT = CALDGDU(RH1,T1,RG)
                                D4=3.43E-08*(D3-VP)*0.3                                       !****濡れ面率0.3  水膜からの蒸発量

                                IF(L5.GE.10.AND.L5.LE.12)THEN                                 !バックシーラー透水抵抗 2.4e+5 m2sPa/kg by　長村
                                    D1=1/( walltypes[IW].layers[J].dx /3.73E-6+2.4E+5/DGDU)                        !飽和時の水分伝導率 3.73e-6 kg/ms(J/kg)　いぶし瓦 by 伊庭　D論
                                ELSE
                                    D1=1/( walltypes[IW].layers[J].dx /(3.73E-6*0.05)+2.4E+5/DGDU)                 !木製品を想定
                                END IF

                                RN(I,J)=(D1*D2+D4+SWJRAIN(I,J))*3600/NDVD+HRN(I,J)            !******水膜の水分量kg/m2
                                WJRAIN(I,J)=-D1*D2*walls[LW].area

                                IF(RN(I,J).LE.0)THEN
                                    RN(I,J)=0.
                                    WJRAIN(I,J)=SWJRAIN(I,J)*walls[LW].area+HRN(I,J)*walls[LW].area                      !******浸水量WJRAIN kg/s
                                END IF
                            END IF
                        END DO
                    END IF

                    # ****************************************

                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                        L1=NAFX(IW,K)
                        L2=NALX(IW,K)
                        L5=walltypes[IW][K].num
                        # *******  VENTED CAVITY  ******

                        IF(L5.EQ.2)THEN
                            I=L1
                            D4=0.
                            D5=0.
                            D1=ALD(LW,K,1)+ALD(LW,K,2)+DGDUQ(LW,K)*QQ(LW,K)
                            D2=ALD(LW,K,1)*WPU(LW,I-1)+ALD(LW,K,2)*WPU(LW,I+1)     !+DGDUQ(LW,K)*QQ(LW,K)*WPUO
                            D3=ALDT(LW,K,1)*(TMP(LW,I-1)-TMP(LW,I))+ALDT(LW,K,2)*(TMP(LW,I+1)-TMP(LW,I))   !+DGDTQ(LW,K)*QQ(LW,K)*(TMPO-TMP(LW,I))
                            IF(LW.EQ.1)THEN
                                D2=D2+DGDUQ(LW,K)*QQ(LW,K)*WPUO                              !最下層
                                D3=D3+DGDTQ(LW,K)*QQ(LW,K)*(TMPO-TMP(LW,I))
                                UHEN=D1*HTMP(LW,I)+D4+D5+1300.*QQ(LW,K)*TMPO                 !D4,D5 不要?
                            ELSE
                                D2=D2+DGDUQ(LW,K)*QQ(LW,K)*WPU(LW-1,I)                        !2階以上
                                D3=D3+DGDTQ(LW,K)*QQ(LW,K)*(TMP(LW-1,I)-TMP(LW,I))
                            END IF
                            IF(RN(LW,I+1).GT.0.)THEN                                       ! D4 内側の水膜蒸発量 濡れ面率0.3
                                T1=TMP(LW,I+1)
                                FS, VP = GOFF(T1)
                                D4=3.43E-08*(VP-XM(LW,L1))*walls[LW].area*0.3
                            END IF
                            IF(RN(LW,I-1).GT.0.)THEN                                       ! D5 外側の水膜蒸発量 濡れ面率0.3
                                T1=TMP(LW,I-1)
                                FS, VP = GOFF(T1)
                                D5=3.43E-08*(VP-XM(LW,L1))*walls[LW].area*0.3
                            END IF
                            UHEN=D2+D3+D4+D5
                            SAHEN=D1
                            WPU(LW,I)=UHEN/SAHEN
                            IF(WPU(LW,I).GE.0.00)WPU(LW,I)=-1.3E-4
                                CYCLE
                            END IF
                            ! **************

                            DO I=L1,L2
                                D7=WJRAIN(LW,I)

                                IF(I.EQ.L1)THEN               !***********端部外
                                    D1=ROW*DPDU(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area/DTM*0.5

                                    IF(I.EQ.1)THEN
                                        IF(walls[LW].room_sideA.EQ.0)THEN
                                            D2=DWX(LW,I+1)+ALD(LW,K,1)
                                            D3=DWX(LW,I+1)*WPU(LW,I+1)+ALD(LW,K,1)*WPUO
                                            D4=DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))+ALDT(LW,K,1)*(TMPO-TMP(LW,I))
                                            SAHEN=D1+D2
                                            UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area
                                        ELSE
                                            D2=DWX(LW,I+1)+ALD(LW,K,1)
                                            D3=DWX(LW,I+1)*WPU(LW,I+1)+ALD(LW,K,1)*WPUR(walls[LW].room_sideA)
                                            D4=DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))+ALDT(LW,K,1)*(ATMPR(walls[LW].room_sideA)-TMP(LW,I))
                                            SAHEN=D1+D2
                                            UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area
                                        END IF
                                    ELSE
                                        D2=DWX(LW,I+1)+ALD(LW,K,1)
                                        D3=DWX(LW,I+1)*WPU(LW,I+1)+ALD(LW,K,1)*WPU(LW,I-1)
                                        D4=DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))+ALDT(LW,K,1)*(TMP(LW,I-1)-TMP(LW,I))
                                        SAHEN=D1+D2
                                        UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area !
                                    END IF

                                END IF

                                IF(I.GT.L1.AND.I.LT.L2)THEN   !*********実質部
                                    D1=ROW*DPDU(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area/DTM
                                    D2=DWX(LW,I)+DWX(LW,I+1)
                                    D3=DWX(LW,I)*WPU(LW,I-1)+DWX(LW,I+1)*WPU(LW,I+1)
                                    D4=DTX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+DTX(LW,I+1)*(TMP(LW,I+1)-TMP(LW,I))
                                    SAHEN=D1+D2
                                    UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area !WJRAIN(LW,I)!
                                END IF
                                IF(I.EQ.L2)THEN                !************端部内
                                    D1=ROW*DPDU(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area/DTM*0.5

                                    IF(I.EQ.NX(IW))THEN
                                        D2=DWX(LW,I)+ALD(LW,K,2)
                                        D3=DWX(LW,I)*WPU(LW,I-1)+ALD(LW,K,2)*WPUR(walls[LW].room_sideB)
                                        D4=DTX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+ALDT(LW,K,2)*(ATMPR(walls[LW].room_sideB)-TMP(LW,I))
                                        SAHEN=D1+D2
                                        UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area  !WJRAIN(LW,I)
                                    ELSE
                                        D2=DWX(LW,I)+ALD(LW,K,2)
                                        D3=DWX(LW,I)*WPU(LW,I-1)+ALD(LW,K,2)*WPU(LW,I+1)
                                        D4=DTX(LW,I)*(TMP(LW,I-1)-TMP(LW,I))+ALDT(LW,K,2)*(TMP(LW,I+1)-TMP(LW,I))
                                        SAHEN=D1+D2
                                        UHEN=D3+D4+D1*HWPU(LW,I)+D7+WJW(LW,I)* walltypes[IW].layers[K].dx *walls[LW].area  !WJRAIN(LW,I)
                                    END IF
                                END IF

                                WPU(LW,I)=UHEN/SAHEN
                                IF(WPU(LW,I).GE.0.00)WPU(LW,I)=-1.3E-3      !-100.

                            END DO
                        END DO   !K
                    END DO   !LW

                    # /////////OVER RELAXATION ////// 
                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                D1=(WPU(LW,I)-AWPU(LW,I))*OMG+AWPU(LW,I)
                                D8=D1-AWPU(LW,I)
                                IF(ABS(D8).GT.AMAX2)AMAX2=ABS(D8)
                                AWPU(LW,I)=D1
                            END DO
                        END DO
                    END DO

                    DO LW=1,len(walls)

                        IW=walls[LW].kwtype

                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num

                            DO I=L1,L2
                                WPW(LW,I)=WPU(LW,I)+WPS(LW,I)
                                WP0=WPU(LW,I)
                                TP=TMP(LW,I)
                                RH0 = WPTRE(WP0,TP,RG)
                                RH(LW,I)=RH0
                                FS, VP = GOFF(TP)

                                XN(LW,I)=RH0*VP*0.01

                                WD = AHGANS(rhm=RH0, ml0=L5)
                                WGT(LW,I)=WD
                                TMPC(LW,I)=TMP(LW,I)-ATP
                            END DO
                        END DO
                    END DO

                    # /////////JUDGEMENT CONVERGENCE//////
                    IT2=IT2+1

                    IF(IT2.GT.MXIT)THEN
                        AMAX2=0.
                    END IF

                    IF(AMAX2.GT.EPS2) GO TO 710

                    DO I=1,NROOM

                        IF(rooms[I].is_calc == False)CYCLE
                        WP0=WPUR(I)
                        TP=ATMPR(I)
                        RH0 = WPTRE(WP0,TP,RG)
                        RHR(I)=RH0
                        FS, VP = GOFF(TP)
                        XR(I)=RH0*VP*0.01
                        TMPR(I)=ATMPR(I)-ATP

                    END DO

                    # ::::::::窓ガラスの結露判定部、収束判定部に追加
                    IF(ITC.GE.30)THEN
                        PRINT *, ' HASSAN KETURO '
                        PAUSE
                        GO TO 555
                    END IF

                    K=0
                    ITC=ITC+1

                    DO I=1,NROOM

                        DO LW=1,NWIN

                            IF(thin_walls[LW].room_sideB.EQ.I)THEN

                                T1=ATMPR(I)
                                RH0=RHR(I)
                                FS, VP = GOFF(T1)
                                X1 = FUNCX(FS,RH0)
                                KMC=MCW(LW)
                                T1=ATMPR(I)-thin_walls[LW].u_vlue/9.3*(ATMPR(I)-TMPO)
                                KMC = CWIF1(T1,X1,XSUT,KMC)
                                MCW(LW)=KMC

                                IF(MCW(LW).EQ.1)THEN
                                    CW=(X1-XSUT)*DTM*ALDC*walls[LW].area    !X1:ROOM, XSUT:WIN
                                    CWV(LW)=CW
                                    SCW(LW)=CW+SCW(LW)
                                    IF(SCW(LW).GT.CWMAX)SCW(LW)=CWMAX
                                    IF(SCW(LW).LE.0)THEN
                                        MCW(LW)=0
                                        SCW(LW)=0.
                                        CWV(LW)=0.
                                    END IF
                                ELSE
                                    KMC=MCW(LW)
                                    KMC = CWIF1(T1,X1,XSUT,KMC)
                                    MCW(LW)=KMC
                                    IF(MCW(LW).EQ.1)THEN
                                        K=1
                                    END IF
                                END IF

                                SCWDAY(LW,IM)=SCW(LW)

                            END IF
                        END DO
                    END DO
    
                    IF(K.EQ.1)GO TO 710
    
                    # ::::::::結露量入力
                    DO LW=1,NWIN
                        IF(CWV(LW).GE.0.)CWVDAY(LW,IM)=CWVDAY(LW,IM)+CWV(LW)
                    END DO

                    # :::::::::結露終了

                    # ///////JUDGEMENT CONVERGENCE OF HEAT AND MOISTURE/////

                    IF(NREPT.EQ.1)THEN  !////COMBINE

                        DO LW=1,len(walls)
                            IW=walls[LW].kwtype
                            DO K=1,KMTL(IW)

                                L1=NAFX(IW,K)
                                L2=NALX(IW,K)
                                L5=walltypes[IW][K].num

                                DO I=L1,L2
                                    D8=WPU(LW,I)-BWPU(LW,I)
                                    D9=TMP(LW,I)-BTMP(LW,I)
                                    BWPU(LW,I)=WPU(LW,I)
                                    BTMP(LW,I)=TMP(LW,I)
                                    IF(ABS(D8).GT.AMAX3)AMAX3=ABS(D8)
                                    IF(ABS(D9).GT.AMAX4)AMAX4=ABS(D9)
                                END DO
                            END DO
                        END DO

                        # /////////JUDGEMENT CONVERGENCE//////
                        IT3=IT3+1

                        IF(IT3.GT.MXIT)THEN
                            PRINT *, "   NO CONVERGENCE   IM=",IM
                            GO TO 556
    
                        END IF

                        IF(AMAX3.GT.EPS2.OR.AMAX4.GT.EPS1) GO TO 810

                    END IF  !/////COMBINE

                    # 556 CONTINUE

                    # ////////////HEAT FLUX J/m2s  MOISTURE FLUX  g/m2s//////////////////
                    # //////////CAL AVERAGE VALUE//////////////
                    DO LW=1,len(walls)
                        IW=walls[LW].kwtype
                        DO K=1,KMTL(IW)
                            L1=NAFX(IW,K)
                            L2=NALX(IW,K)
                            L5=walltypes[IW][K].num
                            DO I=L1,L2
                                RHAVD(LW,I)=RH(LW,I)+RHAVD(LW,I)
                                WGTAVD(LW,I)=WGT(LW,I)+WGTAVD(LW,I)
                                XNAVD(LW,I)=XN(LW,I)+XNAVD(LW,I)
                                TPAVD(LW,I)=TMPC(LW,I)+TPAVD(LW,I)
                            END DO
                        END DO
                    END DO

                    IF(NRAIN.GT.0)THEN      !浸水量の積算 SUMRAIN (kg/m2h)
                        DO K=1,NRAIN
                            I=NRAINPOINT(K,1)
                            J=NRAINPOINT(K,2)
                            SUMRAIN(I,J,IM)=SUMRAIN(I,J,IM)+WJRAIN(I,J)*3600/NDVD/walls[I].area
                        END DO
                    END IF
                # 666 CONTINUE

                IF(NOUTAV.GT.0)THEN                           !材料の平均含水率の計算
                    DO L=1,NOUTAV
                        LW=NOUTAVD(L,1)
                        J=NOUTAVD(L,2)
                        IW=walls[LW].kwtype
                        L1=NAFX(IW,J)
                        L2=NALX(IW,J)
                        D1=walltypes[IW].layers[J].dx
                        D2=WGT(LW,L1)*0.5
                        D2=D2+WGT(LW,L2)*0.5

                        DO K=L1+1,L2-1
                            D2=D2+WGT(LW,K)
                        END DO

                        WOUTAV(L)=D2/(L2-L1+1)
                    END DO
                END IF

            # 777 CONTINUE

            # *******SAT   IREKAE *******
            DO I=1,len(walls)
                SAT(I,1)=SAT(I,25)
                SJIN(I,1)=SJIN(I,25)
            END DO

            TEMPO(1)=TEMPO(25)				

            SJD(1)=SJD(25)				
            SJS(1)=SJS(25)				
            SJN(1)=SJN(25)				
            XOD(1)=XOD(25)				
            RHOD(1)=RHOD(25)				

            # ***************************************
            DO LW=1,len(walls)
                IW=walls[LW].kwtype
                DO K=1,KMTL(IW)
                    L1=NAFX(IW,K)
                    L2=NALX(IW,K)
                    L5=walltypes[IW][K].num
                    DO I=L1,L2
                        RHDIS(LW,I,MON,DAY)=RHAVD(LW,I)/24./NDVD
                        XNDIS(LW,I,MON,DAY)=XNAVD(LW,I)/24./NDVD
                        TPDIS(LW,I,MON,DAY)= TPAVD(LW,I)/24./NDVD
                        RHAVD(LW,I)=0.
                        WGTAVD(LW,I)=0.
                        XNAVD(LW,I)=0.
                        TPAVD(LW,I)=0.
                    END DO
                END DO 
            END DO

            # //////////雨水浸入量積算値DAYSUMRAIN(kg/m2DAY)/////////
            DO K=1,NRAIN
                I=NRAINPOINT(K,1)
                J=NRAINPOINT(K,2)
                DAYSUMRAIN(I,J)=0.
                D1=0.
                DO L=1,24
                    D1=SUMRAIN(I,J,L)+D1
                END DO
                DAYSUMRAIN(I,J)=D1
            END DO

            # ////////Calculation Mass Loss（日平均値による算出）//////////

            DO LW=1,len(walls)
                IW=walls[LW].kwtype

                DO K=1,KMTL(IW)

                    L1=NAFX(IW,K)
                    L2=NALX(IW,K)
                    L5=walltypes[IW][K].num
                    WLOSSMAX=0.6

                    IF(L5.EQ.4.OR.L5.EQ.5)THEN

                        IF(I_HCOFF.EQ.1)THEN
                            HCOFF=0.319  			  !水分生成量
                        ELSE
                            HCOFF=0.
                        END IF
    
                        DO I=L1,L2
                            D1=TPDIS(LW,I,MON,DAY)
                            D2=RHDIS(LW,I,MON,DAY)
                            DLOSS = wdm.WOOD_ROT(LW, I, D1, D2, ROTOMG)
                            WLOSS(LW,I)=WLOSS(LW,I)+DLOSS
                            IF(WLOSS(LW,I) > WLOSSMAX)DLOSS=0.
                            WJW(LW,I)=HCOFF*DLOSS*GMA(L5)/86400.   ! Time unit:h=24,s=86400
                        END DO 
                    END IF 

                END DO 
            END DO

            # /////////////End Mass Loss//////////////

            WRITE(12,*)NYEAR,MON,DAY,(DAYSUMRAIN(NRAINPOINT(K,1),NRAINPOINT(K,2)),K=1,NRAIN),(WOUTAV(K),K=1,NOUTAV),((WGTDIS(NOUTMLD(L,1),K,MON,DAY),K=NOUTMLD(L,2),NOUTMLD(L,3)),L=1,NOUTML),&
                    ((WLOSS(NOUTMLD(L,1),K),K=NOUTMLD(L,2),NOUTMLD(L,3)),L=1,NOUTML),((WJW(NOUTMLD(L,1),K)*86400,K=NOUTMLD(L,2),NOUTMLD(L,3)),L=1,NOUTML)
        # 888 CONTINUE

    # 999 CONTINUE

    REWIND 7
    REWIND 13

    WRITE(12,*)' '
    WRITE(12,*)'  MOISTURE CONTENT DISTRIBUTION (mass%) '
    WRITE(12,*)'  YEAR MON I '
    DO K=1,12
        J=1
        WRITE(12,*)NYEAR,K,(WGTDIS(1,I,K,J),I=1,NX(2))
    END DO
    WRITE(12,*)' '

# 444 CONTINUE


WRITE(12,*)' '
    WRITE(12,*)'  '
    PRINT *,'  POINT,TMPC(I),RH(I),WGT(I),XN(I)'
    DO I=1,NX(1)
        PRINT *,I,TMPC(1,I),RH(1,I),WGT(1,I),XN(1,I)
    END DO

# *******************************************

PRINT *, MOJI,WPW(1,1),WGT(1,1),XN(1,1),XO,SCWDAY(1,1)
PRINT *,XO,AXN(1,1),XM(1,1),ATMPQ(1,1),HWPUR(1),DGDUR(1)

555 CONTINUE

CLOSE(8)
CLOSE(7)
CLOSE(12)
CLOSE(13)

PAUSE

END


def SATUWPT(TP)

    # J/(kg K)
    CPW = ( 30.36 + 0.009615 * TP + 0.00000118 * TP * TP ) / ( 0.018016 )

    FS, VP = GOFF(TP)

    # J/ kg
    SW=(644243)+CPW*(TP-273.15)-TP*CPW*DLOG(TP/273.15)+461.5*TP*DLOG(VP/101325)

    return SW


def GOFF(TMP)
    """
    
    Args:
        TMP: K
    
    Returns:
        FS: mmHg
        VP: Pa
    """
    EW = -6096.938 / TMP + 21.2409642 - TMP * 2.711193E-2 + TMP * TMP * 1.673952E-5 + 2.433502 * np.log(TMP) 
    
    VP = np.exp(EW)
    
    FS = VP / 133.322 

    return FS, VP


def FUNCX(FS,RH)
    
    FD = FS * RH * 0.01
    X = FD * 0.6217 / (760.0 - FD)
    return X

def REWPT(RH, TP, RG)

    RH0=RH*0.01
    WP=RG*TP*DLOG(RH0)
    return WP


def WPTRE(wp, tp, rg):
    """
    含水率などから平衡相対湿度 (RH) を計算する。
    """
    
    d1 = wp / rg / tp
    d2 = math.exp(d1)
    rh = d2 * 100.0
    
    # 相対湿度は100%を超えることはないので、上限を設定することが多い
    # rh = min(rh, 100.0) 
    
    return rh


def CALDPDU(wpt, tp, gma, ml0, rg, row):
    """
    含水率変化に対するポテンシャル変化率 (DPDU) を計算する。
    """
    d1 = wpt
    
    # 元のFortranのロジック: 正の値の場合は -100.0 に強制
    if d1 > 0.0:
        d1 = -100.0
    
    # 差分幅の計算 (d1の1%)
    dw = abs(d1 * 0.01)
    
    # 微小変化させた含水率
    w1 = d1 + dw
    w2 = d1 - dw
    
    # 以前定義した WPTRE (平衡相対湿度計算) を呼び出し
    rh1 = WPTRE(w1, tp, rg)
    rh2 = WPTRE(w2, tp, rg)
    
    # 外部定義されている前提の AHGANS
    wd1 = ahgans(rhm=rh1, ml0=ml0)
    wd2 = ahgans(rhm=rh2, ml0=ml0)
    
    # VGTの計算 (0.01は%を小数に戻す係数と推測)
    vgt1 = 0.01 * wd1 * gma / row
    vgt2 = 0.01 * wd2 * gma / row
    
    # 中央差分による勾配(微分値)の近似
    # 0.5 * (VGT1 - VGT2) / DW
    if dw != 0:
        dpdu = 0.5 * (vgt1 - vgt2) / dw
    else:
        dpdu = 0.0
        
    return dpdu


def CALDGDU(rh, tmp, rg):
    """
    含水率変化に対する絶対湿度の変化率(DGDU)および
    温度変化に対する絶対湿度の変化率(DGDT)を計算する。
    """
    # --- DGDU の計算 (湿度差による微分) ---
    rh1 = rh + 0.01
    rh2 = rh - 0.01
    
    # 外部定義されている前提の関数
    # rewpt(相対湿度, 温度, パラメータ)
    wp1 = REWPT(rh1, tmp, rg)
    wp2 = REWPT(rh2, tmp, rg)
    
    # GOFF関数(飽和水蒸気圧等)から絶対湿度を計算
    # 以前の回答通り、2つの戻り値を想定
    fs1, vp = GOFF(tmp)
    x1 = vp * rh1 * 0.01
    x2 = vp * rh2 * 0.01
    
    # ゼロ除算を回避しつつ微分布を計算
    if abs(wp1 - wp2) > 1e-15:
        dgdu = abs((x1 - x2) / (wp1 - wp2))
    else:
        dgdu = 0.0

    # --- DGDT の計算 (温度差による微分) ---
    tp1 = tmp + 0.1
    tp2 = tmp - 0.1
    
    # 温度を変えて飽和水蒸気圧を取得
    _, vp1 = GOFF(tp1)
    _, vp2 = GOFF(tp2)
    
    x1_t = vp1 * rh * 0.01
    x2_t = vp2 * rh * 0.01
    
    # 温度勾配による絶対湿度の変化率
    # tp1 - tp2 は常に 0.2 なのでゼロ除算の心配はほぼありません
    dgdt = abs((x1_t - x2_t) / (tp1 - tp2))

    return dgdu, dgdt


def AHGANS(rhm, ml0):
    """
    含水率計算関数 (AHGANS)
    wd  : 含水率 (結果を返す)
    rhm : 相対湿度 (%)
    ml0 : 材料種別コード
    """

    d1 = rhm * 0.01

    # ML0の値に応じて処理を分岐
    # FortranのGO TO (11,10,13,14,15,16,17,18,19,20,21,22,13,13,13,13,13,28,29),ML0 に対応
    if ml0 == 1:
        # ラベル 11: SEASING BOARD
        if rhm <= 89:
            wd = -math.log(1 - 0.01 * rhm) / 0.1612944
        elif 89 < rhm < 95:
            wd = 639.22 - 14.988 * rhm + 0.08943 * rhm * rhm
        else: # rhm >= 95
            wd = 2.446 * rhm - 209.9

    elif ml0 in [3, 13, 14, 15, 16, 17]:
        # ラベル 13: GLASS WOOL
        # ※元のコードでML0が3, 13-17の際にラベル13へ飛ぶ設定
        d2 = math.exp(-1.3016 * (1 - d1**101.667))
        d3 = 0.16 * d2 * d1**7.6448
        if d1 > 0.98:
            d3 = 47.79595 * d1 - 46.79595
        wd = d3 * 100.0

    elif ml0 == 4:
        # ラベル 14: WOOD
        d2 = math.exp(-2.1708 * (1 - d1**90.988))
        d3 = 1.58 * d2 * d1**1.3664
        wd = d3 * 100.0

    elif ml0 == 5:
        # ラベル 15: PLY WOOD
        d2 = math.exp(-1.6471 * (1 - d1**32.051))    #ASHRAE PLYWOOD 1
        d3 = 1.0625 * d2 * d1**1.5373
        wd = d3 * 100.0

    elif ml0 == 6:
        # ラベル 16: PLASTER BOARD
        d2 = math.exp(-1.6445 * (1 - d1**75.811))    #ASHRAE GYPSUM
        d3 = 0.56978 * d2 * d1**0.18552
        wd = d3 * 100.0

    elif ml0 == 7:
        # ラベル 17: ALC
        d2 = math.exp(-0.66866 * (1 - d1**17.316))
        d3 = 0.071101 * d2 * d1**1.0146
        wd = d3 * 100.0

    elif ml0 == 8:
        # ラベル 18: NANSHITSUSENIBAN
        d2 = math.exp(-0.93059 * (1 - d1**4.6145))
        d3 = 0.35093 * d2 * d1**0.48739
        wd = d3 * 100.0

    elif ml0 == 9:
        # ラベル 19: THERMOPLY
        d2 = math.exp(-1.2566 * (1 - d1**5.483))
        d3 = 0.4367 * d2 * d1**0.3513
        wd = d3 * 100.0

    elif ml0 == 10:
        # ラベル 20: SAIDHING
        d2 = math.exp(-1.1225 * (1 - d1**15.878)) #ASHRAE Siding No.38
        d3 = 0.635 * d2 * d1**2.5972
        wd = d3 * 100.0

    elif ml0 == 11:
        # ラベル 21: tutikabe
        d2 = math.exp(-0.69 * (1 - d1**25.06))
        d3 = 0.0645 * d2 * d1**0.658
        wd = d3 * 100.0

    elif ml0 == 12:
        # ラベル 22: 軽量モルタル
        d2 = math.exp(-1.75 * (1 - d1**2.015))
        d3 = 0.203 * d2 * d1**(-0.120)
        wd = d3 * 100.0

    elif ml0 == 18:
        # ラベル 28: 集成材（OMソーラー）
        d2 = math.exp(-0.849 * (1 - d1**9.527))
        d3 = 0.381 * d2 * d1**(0.738)
        wd = d3 * 100.0

    elif ml0 == 19:
        # ラベル 29: 構造用合板（OMソーラー）
        d2 = math.exp(-0.907 * (1 - d1**9.412))
        d3 = 0.382 * d2 * d1**(0.670)
        wd = d3 * 100.0

    # ラベル 10: CONTINUE (処理の終了)
    return wd

#::::::::サブルーチン追加
def CWIF1(T1,X1,X0,mcw)
    
    RH0=100.0
    
    FS0, VP = GOFF(T1)
    
    X0 = FUNCX(FS0,RH0)

    if x1 >= x0 and mcw == 0:
        mcw = 1

    return mcw

# 三浦メモ　このサブルーチンは使われていなかったため、削除した。
# SUBROUTINE DOSU(NYEAR,MON,I,D1)

# 三浦メモ　このサブルーチンは使われていなかったため、削除した。
# SUBROUTINE SIGMA(NYEAR)

!********YOMITOBASI******************				
SUBROUTINE YOMI(MON1,MDAY1,KDAY)
CHARACTER MOJI*50				
DIMENSION MDAY(12)				
DATA MDAY/31,28,31,30,31,30,31,31,30,31,30,31/				
L1=0				
KDAY=0				
DO 60 I = 1,MON1-1				
    L1=MDAY(I)+L1				
60 CONTINUE      				
L2=MDAY1+L1-1				
KDAY=L2				
L3=L2*25				
PRINT *,' YOMITOBASI=',L3				
DO 70 I = 1,L3+1				
    READ(7,*)MOJI
70 CONTINUE

L4=L2*24

DO I=1,L4
    READ(13,*)D1,D2
END DO

PRINT *, '  MOJI ',MOJI,D1,D2				
RETURN				
END				
!********************************				


def SOLOCT(kw, sido, skeido, hoi, akasya, norient, l, data):
    """
    太陽位置と壁面入射角の計算
    kw      : 年通算日 (1-366)
    sido    : 緯度 [度]
    skeido  : 経度 [度]
    hoi     : 方位角 [度]
    akasya  : 傾斜角 [度]
    norient : 方位番号 (1:南, 2:西, 3:北, 4:東 など)
    l       : 地点インデックス
    data    : 共通データを格納した辞書
    """
    
    # 定数 (DATA文相当)
    C = [0.006322, 0.405748, 0.153231, 0.00588, 0.207099, 0.003233, 0.620129]
    SD = [0.000279, 0.122772, 1.498311, 0.165458, 1.261546, 0.005354, 1.1571]
    
    # 単位変換 [度 -> ラジアン]
    hoim = hoi * math.pi / 180.0
    bkaku = akasya * math.pi / 180.0
    hido = sido * math.pi / 180.0
    
    # 共通データへの保存 (他サブルーチンで使用するため)
    data['BKAKU'] = bkaku
    
    # 太陽赤緯(SEKI)と均時差(EKNJ)の計算
    omega = 2.0 * math.pi * kw / 366.0
    seki = C[0] - C[1]*math.cos(omega + C[2]) - C[3]*math.cos(2.0*omega + C[4]) - C[5]*math.cos(3.0*omega + C[6])
    eknj = -SD[0] + SD[1]*math.cos(omega + SD[2]) - SD[2+1]*math.cos(2.0*omega - SD[4]) - SD[5]*math.cos(3.0*omega - SD[6])
    
    # 時刻ループ (I=4時〜20時)
    for i in range(4, 21):
        k = i + 1 # Fortranの添字合わせ
        
        # 時角 TKU の計算 (日本標準時135度基準)
        tku = math.pi * ((i + eknj - 12.0) + (skeido - 135.0) / 15.0) / 12.0
        
        # 太陽高度角に関連する成分
        w_val = math.sin(hido) * math.sin(seki) + math.cos(hido) * math.cos(seki) * math.cos(tku)
        v2_val = math.cos(seki) * math.sin(tku)
        v3_val = -math.sin(seki) * math.cos(hido) + math.cos(seki) * math.sin(hido) * math.cos(tku)
        
        # 負の値をカット (日没後)
        if v3_val < 0: v3_val = 0.0
        if w_val < 0: w_val = 0.0
        
        data['W'][k] = w_val
        hkodo = math.asin(w_val)
        
        # 壁面方位の計算
        j = norient
        e1 = (j - 1) * 90.0
        e2 = e1 * math.pi / 180.0
        hoim2 = hoim + e2
        
        # 太陽方位角成分の計算 (分母が0に近い場合の回避)
        cos_h = math.cos(hkodo)
        if cos_h < 1e-6:
            s1 = 0.0
            s3 = 0.0
        else:
            s1 = v2_val / cos_h
            s3 = v3_val / cos_h
        
        # 壁面入射角の余弦 (DV4) の計算
        s4 = s3 * math.cos(hoim2) + s1 * math.sin(hoim2)
        
        # 最終的な入射角余弦
        res = math.cos(bkaku) * math.sin(hkodo) + math.sin(bkaku) * math.cos(hkodo) * s4
        
        if res <= 0.0:
            res = 0.0
            
        data['DV4'][l, k] = res

        return data


def SATCAL(rlf, as_coef, emi, l_idx, ms: MS, ss: SS):
    """
    SAT(相当外気温度)の計算
    rlf: 地面反射率
    as_coef: 日射吸収率 (FortranのAS)
    emi: 長波長放射率 (FortranのEMI)
    l_idx: 計算対象の層/地点インデックス (FortranのL)
    """
    afo = 22.4  # 表面熱伝達率の逆数に関連する定数
    
    # D1: 形態係数（天空率）の計算
    d1 = (1 + np.cos(SS.BKAKU)) / 2
    
    # I=2から25までループ
    for i in range(2, 26):

        # 直達日射受熱量
        dsrs = MS.SJD[i] * SS.DV4[l_idx, i]
        
        # 天空日射受熱量
        srs = d1 * MS.SJS[i]
        
        # 地面反射日射受熱量の計算用
        dsrh = MS.SJD[i] * SS.W[i]
        hnisya = MS.SJS[i] + dsrh # 水平全天日射量
        
        # 反射日射受熱量 (1-D1は地面に対する形態係数)
        rsrs = (1 - d1) * rlf * hnisya
        
        # 全受熱日射量
        snisya = dsrs + srs + rsrs
        
        # SAT (相当外気温度)
        # SAT = 気温 + (1/AFO) * (吸収日射量 - 夜間放射の影響)
        radiation_loss = d1 * MS.SJN[i] * emi
        net_gain = as_coef * snisya - radiation_loss
        
        MS.SAT[l_idx, i] = MS.TEMPO[i] + (1.0 / afo) * net_gain
        
        # 各成分の保存
        MS.SJIN[l_idx, i] = net_gain
        MS.WSJIN[l_idx, i] = as_coef * (srs + rsrs)
        MS.WSJD[l_idx, i] = as_coef * dsrs
    
    return ms, ss


def DIFF(k, rh, tmp, rg):
    """
    液体伝導率(RML)の計算
    """
    # 外部関数の呼び出しを想定
    # ahgans(rhm, ml0), rewpt(rh, tmp, rg)
    wd = AHGANS(rhm=rh, ml0=k)
    
    if 45 < wd < 110:
        d1 = wd * 0.01
        # DW = e^(a + bX + cX^2)
        dw = math.exp(-30.91 + 4.2967 * d1 - 0.22017 * (d1**2))
        
        rh1 = rh + 0.005
        rh2 = rh - 0.005
        
        # 微分を差分で近似している計算
        wp1 = REWPT(rh1, tmp, rg)
        wp2 = REWPT(rh2, tmp, rg)
        wd1 = AHGANS(rhm=rh1, ml0=k)
        wd2 = AHGANS(rhm=rh2, ml0=k)
        
        # ゼロ除算のチェック
        if abs(wp1 - wp2) > 1e-12:
            rml = 998.0 * dw * abs((wd1 - wd2) * 0.01 / (wp1 - wp2))
        else:
            rml = 0.0
    else:
        rml = 0.0
        
    return rml


def ROOM(kday, trav, trdt, rrav):
    """
    室内の温湿度計算
    kday: 通算日(または経過時間)
    trav: 平均気温
    trdt: 気温の振幅
    rrav: 平均相対湿度
    """
    # D1 = 2 * PI * (KDAY - 212) * 24 / 8760
    # 212日はおそらく8月1日付近を基準（位相）にしている
    d1 = 2 * math.pi * (kday - 212) * 24.0 / 8760.0
    
    # 室温の計算 (余弦曲線での近似)
    tr = trdt * math.cos(d1) + trav
    rr = rrav
    
    t0 = tr + 273.15  # 摂氏 -> ケルビン
    
    # GOFF関数から飽和水蒸気圧等を取得
    fs0, vp = GOFF(t0)
    
    # 水蒸気圧の計算 (VP * 相対湿度[%] * 0.01)
    xr = vp * rr * 0.01
    
    return tr, rr, xr


class WoodDecayModel:
    """
        Biological damage function for wood rot decay 
        File name:hdamage.f95
        Reaction model
        Hiroaki Saito
    """

    def __init__(self, nwp=50, nxp=150):
        # AIコメント：状態を保持するための配列 (FortranのDIMENSION相当)
        self.time_s = np.zeros((nwp + 2, nxp + 2))   # インデックス余裕を持たせる
        self.l_stage = np.zeros((nwp + 2, nxp + 2), dtype=int)

    def WOOD_ROT(self, k, i, tmp, rh, rotomg):
        """
        木材の腐朽被害関数
        k, i: 現在のグリッド等のインデックス
        tmp: 温度
        rh: 相対湿度
        rotomg: 腐朽速度係数
        """
        w = -1.0  # OSB (TIN関数に渡すパラメータと推測)
        rh_growth = 98.0
        dloss = 0.0
        reaction_k = 0.0
        dt_damage = 60 * 60 * 24  # 1日 (秒)

        # --- 発芽期 (Initial response time) ---
        if tmp <= 0.0:
            self.time_s[k, i] = 0.0
        else:
            # 以前定義した RH_CRITICAL を呼び出し
            rhc = RH_CRITICAL(tmp)
            
            if rh < rhc:
                self.time_s[k, i] = 0.0
            else:
                # TIN関数の呼び出し (引数wが必要)
                gc, fc = TIN(tmp, rh, w)
                
                # ゼロ除算回避
                if fc != 0:
                    d1 = -gc / fc
                else:
                    d1 = 100.0 # 仮の大きな値
                
                # 発芽時間の計算 (TIME_INT)
                if d1 > 0:
                    time_int = d1
                elif d1 > -0.59:
                    time_int = 0.0
                else:
                    time_int = 100.0
                
                # 係数調整
                if time_int <= 0.5:
                    time_int = 0.5 * dt_damage * 30.0
                else:
                    time_int = time_int * dt_damage * 30.0
                
                # 時間経過の蓄積
                self.time_s[k, i] += dt_damage
                if self.time_s[k, i] > time_int:
                    self.l_stage[k, i] = 1

        # --- 成長期 (Growth Stage) ---
        if 0 < tmp <= 40:
            # 自身または隣接するノードが発芽しているかチェック
            # (k-1, k+1 の範囲エラーを防ぐため境界チェックが必要)
            if (self.l_stage[k, i] == 1 or 
                self.l_stage[k-1, i] == 1 or 
                self.l_stage[k+1, i] == 1):
                
                if rh >= rh_growth:
                    # 腐朽反応速度の計算
                    reaction_k = (2.77 - 3.23 * tmp + 0.865 * (tmp**2) - 0.0189 * (tmp**3)) * 1e-10 * rotomg
                    dloss = reaction_k * self.dt_damage

        return dloss
	

def RH_CRITICAL(tmp):
    """
    AIのコメント：温度(tmp)に基づき、臨界相対湿度(rhc)を計算する。
    """

    if tmp <= 15.0:
        rhc = -0.5 * tmp + 100.0
    else:
        rhc = 92.5
        
    # 上限を100に制限
    if rhc >= 100.0:
        rhc = 100.0
        
    return rhc


def TIN(t, rh, w):
    """
    AIのコメント：温度(t), 相対湿度(rh), 含水率(?)w に基づき、FCとGCを計算する。
    """
    
    fc = (0.1384 * t + 0.4370 * rh - 42.9450 + w * (0.034 * t - 0.021 * rh + 1.721))
    
    gc = (-2.2270 * t - 0.0347 * rh + 0.0244 * t * rh + w * (-0.504 * t + 0.0096 * rh + 0.0047 * t * rh))
          
    return gc, fc
