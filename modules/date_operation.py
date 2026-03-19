from dataclasses import dataclass


def get_month_days():

    return [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def get_step_d(month: int, day: int) -> int:
    """通算日を求める。
    
    Args:
        month: 月 (1〜12)
        day: 日 (1~28,30,31)

    Returns:
        通算日(1～365)
    
    """

    # 平年の各月の日数リスト（1月〜12月）
    month_days = get_month_days()

    # 指定した月の前月までの日数を合計し、現在の月の日数を足す
    day_of_year = sum(month_days[:month-1]) + day

    return day_of_year


def get_step_d_t(month: int, day: int, hour: int) -> int:
    """通算ステップを求める。

    Args:
        month: 月 (1〜12)
        day: 日 (1~28,30,31)
        hour: 時刻 (0〜23)
    
    Returns:
        通算ステップ (0〜8759)

    """

    # 通算日
    day_of_year = get_step_d(month=month, day=day)

    n = (day_of_year - 1) * 24 + hour

    return n


def get_step_n(month: int, day: int, hour: int, n_hour: int, n_div: int):

    n_d_t = get_step_d(month=month, day=day, hour=hour)

    return n_d_t * n_div + n_hour


def get_step_total(n_div: int):
    """年間の総ステップ数を取得する。

    Args:
        n_div: 1時間の分割数

    Returns:
        年間の総ステップ数
    """

    return n_div * 8760




@dataclass
class Period:

    # 計算年
    LYEAR: int = 2

    # 開始月
    MON1: int = 1

    # 開始日
    MDAY1: int = 1

    # 終了月
    LMON: int = 12

    # 終了日
    LDAY: int = 31

    def get_start_month(self, year: int) -> int:

        # 初年度は指定した月から開始し、2年目からは1月から開始する。
        return self.MON1 if year == 1 else 1
    
    def get_end_month(self, year: int):

        if self.LYEAR == year:
            return self.LMON
        else:
            return 12
    
    def get_start_day(self, year: int) -> int:

        return self.NDAY1 if year == 1 else 1
    
    def get_end_day(self, year: int, month: int) -> int:

        if year == self.LYEAR and month == self.LMON:

            return self.LDAY

        else:

            return get_month_days()[month-1]


