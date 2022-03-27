from statistics import NormalDist
from math import log, exp


class SeismicDisplacement:
    def __init__(self, data) -> None:
        self.ky = data["Yield Coefficiend"]
        self.Ts = data["Fundamental Period"]
        self.Mw = data["Mw"]
        self.Sa_15 = data["Sa(1.5Ts)"]  # g
        self.P1 = data["P1"]  #%
        self.P2 = data["P2"]  #%
        self.P3 = data["P3"]  #%
        self.d_threshold = data["Displacement Threshold"]  # cm
        self.sigma = 0.66

    def get_zero_disp(self):
        a = 0.22 if self.Ts <= 0.05 else -1.1
        D0 = exp(
            a
            + -2.83 * log(self.ky)
            - 0.333 * log(self.ky) ** 2
            + 0.566 * log(self.ky) * log(self.Sa_15)
            + 3.04 * log(self.Sa_15)
            - 0.244 * log(self.Sa_15) ** 2
            + 1.5 * self.Ts
            + 0.278 * (self.Mw - 7)
        )

        return D0

    def get_P0(self):
        inside = (
            -1.76
            - 3.22 * log(self.ky)
            - 0.484 * self.Ts * log(self.ky)
            + 3.52 * log(self.Sa_15)
        )

        P0 = 1 - NormalDist().cdf(inside)

        return P0

    def get_displ(self, probability):
        P0 = self.get_P0()
        D0 = self.get_zero_disp()
        try:
            inverse_dist = NormalDist().inv_cdf(1 - (probability / 100) / (1 - P0))
            disp = exp(log(D0) + self.sigma * inverse_dist)
        except:
            disp = "<1"

        return disp

    def get_P_threshold(self):
        P0 = self.get_P0()
        D0 = self.get_zero_disp()
        P_threshold = (1 - P0) * (
            1 - NormalDist().cdf((log(self.d_threshold) - log(D0)) / self.sigma)
        )
        return P_threshold


data = {
    "Yield Coefficiend": 0.28,
    "Fundamental Period": 0.22,
    "Mw": 6.5,
    "Sa(1.5Ts)": 0.85,
    "P1": 84,
    "P2": 50,
    "P3": 16,
    "Displacement Threshold": 5,
}

SD = SeismicDisplacement(data)
SD.get_displ(84)
SD.get_displ(50)
SD.get_displ(16)
SD.get_P_threshold()
