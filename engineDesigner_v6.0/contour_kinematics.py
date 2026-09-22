import math
import csv

class Contour:
    def __init__(self, file):
        self.x = []
        self.y = []

        with open(file, 'r') as file:
            reader = csv.reader(file)
            next(reader)
            for row in reader:
                self.x.append(float(row[0]))
                self.y.append(float(row[1]))

        #getting throat ardius with min()
        self.radius = min(self.y)
        self.throat_area = math.pi * (self.radius**2)

    #linear interpolation
    def local_radius(self, x):
        for i in range(len(self.x) - 1):
            if self.x[i] <= x <= self.x[i+1]:
                x1 = self.x[i]
                x2 = self.x[i+1]
                y1 = self.y[i]
                y2 = self.y[i+1]
                return y1 + (y2 - y1) * ((x - x1) / (x2 - x1))
        return self.y[-1]

def gas_state(areaL, areaT, Pstag, Tstag, gamma, R):
    area_ratio = areaL / areaT
    #solving for local chamber mach number with area-mach relation
    #uses bisection root-finding

    low_mach_bound = 0.01
    high_mach_bound = 1.0

    for i in range(30):
        mach_local = (low_mach_bound + high_mach_bound) / 2
        t1 = 2 / (gamma + 1)
        t2 = 1 + ((gamma - 1) / 2) * (mach_local **2)
        exponent = (gamma + 1) / (2 * (gamma -1))
        calc_area_ratio = (1 / mach_local) * ((t1 * t2) ** exponent)

        if calc_area_ratio > area_ratio:
            low_mach_bound = mach_local
        else:
            high_mach_bound = mach_local

    Tstat = Tstag / (1 + ((gamma - 1) / 2) * (mach_local **2))
    Pstat = Pstag / ((1 + ((gamma - 1) / 2) * (mach_local **2)) ** (gamma / (gamma -1)))

    c = math.sqrt(gamma * R * Tstat)
    local_velo = mach_local * c

    return Pstat, Tstat, local_velo, mach_local        



