# CONSTANTS.PY - Encapsulates global constants as an object to be called

class Constants:
    def __init__(self):
        # ==== Define Constants ====
        self.g0 = 9.81 # Gravity (m/s^2)
        self.bar = 100000 # 1 Bar in Pa
        self.Ru = 8314.46 # Universal Gas Constant (J/kmolK)
        self.psi_to_bar = 14.5038 # The number of psi per bar
        self.psi_to_atm = 14.69594878 # The number of psi per atm
        self.P_amb = 1.01325 # Ambient pressure (bar)
        self.C_to_K = 273.15 # Celcius to Kelvin conversion
        self.ft_to_m = 0.3048 # Feet to meters conversion
