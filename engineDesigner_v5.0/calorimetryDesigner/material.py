import numpy as np

class MaterialProperties:
    def __init__(self, conductivity: float | list[np.ndarray], yield_strength: float, ultimate_strength: float):
        if (isinstance(conductivity, (float))):
            self._k = [np.array([300]), np.array([conductivity])]
        else:
            assert(len(conductivity) == 2) and (len(conductivity[0]) == len(conductivity[1])), "Temperature varying conductivity must be a list of two arrays of equal length."
            self._k = conductivity
        self.yield_strength = yield_strength
        self.ultimate_strength = ultimate_strength

    def conductivity(self, temperature: float) -> float:
        idx = np.abs(self._k[0] - temperature).argmin()
        if idx != 0 and idx != len(self._k[0]) - 1:
            # Linear interpolation
            t1, t2 = self._k[0][idx], self._k[0][idx + 1]
            k1, k2 = self._k[1][idx], self._k[1][idx + 1]
            return k1 + (k2 - k1) * (temperature - t1) / (t2 - t1)
        else:
            return self._k[1][idx]
    
    def yield_strength_value(self) -> float:
        return self.yield_strength

    def ultimate_strength_value(self) -> float:
        return self.ultimate_strength
    

copper = MaterialProperties(
    # conductivity data: https://www.engineeringtoolbox.com/copper-density-specific-heat-thermal-conductivity-vs-temperature-d_2223.html (not sure how credible this is but it cross checks fine to charts)
    conductivity=[
        np.array([100, 150, 200, 250, 300, 400, 600, 800, 1000, 1200]), # K
        np.array([480, 429, 413, 406, 401, 393, 379, 366,  352,  339]) # W/m-K
    ],
    yield_strength=65e6, # MPa
    ultimate_strength=200e6 # MPa
)

inconel718 = MaterialProperties(
    # conductivity data: https://iopscience.iop.org/article/10.1088/1742-6596/1382/1/012175/pdf
    conductivity=[
        np.array([ 298,   400,   500,   600,   700,   800,  1100,  1200,  1300,  1400]), # K
        np.array([9.94, 11.59, 13.24, 14.91, 16.61, 18.43, 22.72, 23.61, 24.47, 25.32]) # W/m-K
    ],
    yield_strength=1100e6, # MPa
    ultimate_strength=1375e6 # MPa
)


if __name__ == "__main__":
    # Example usage
    temp = 10000  # Example temperature in K
    print(f"Copper conductivity at {temp} K: {copper.conductivity(temp)} W/m-K")