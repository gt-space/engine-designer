from rocketcea.cea_obj import CEA_Obj
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
import re
import math


class Engine:
      
    @property
    def fuel_name(self) -> str:
        return self._fuel_name
    @property
    def ox_name(self) -> str:
        return self._ox_name
    @property
    def thrust(self) -> float:
        return self._thrust
    @property
    def injector_pressure(self) -> float:
        return self._injector_pressure
    @property
    def exit_pressure(self) -> float:
        return self._exit_pressure
    @property
    def mixture_ratio(self) -> float:
        return self._mixture_ratio
    @property
    def fac_contraction_ratio(self) -> float:
        return self._fac_contraction_ratio
    @property
    def characteristic_length(self) -> float:
        return self._characteristic_length
    @property
    def nozzle_diverging_angle(self) -> float:
        return self._nozzle_diverging_angle


    @property
    def exit_velocity(self) -> float:
        return self._exit_velocity
    @property
    def total_mass_flow_rate(self) -> float:
        return self._total_mass_flow_rate
    @property
    def ox_mass_flow_rate(self) -> float:
        return self._ox_mass_flow_rate
    @property
    def fuel_mass_flow_rate(self) -> float:
        return self._fuel_mass_flow_rate
    @property
    def throat_area(self) -> float:
        return self._throat_area
    @property
    def throat_radius(self) -> float:
        return self._throat_radius


    @property
    def num_points(self) -> int:
        return self._num_points
    @property
    def converging_start_index(self) -> int | None:
        return self._converging_start_index   
    @property
    def throat_index(self) -> int | None:
        return self._throat_index


    def __len__(self):
        return self._num_points    

    def __init__(self, num_points:int, fuel_name:str, ox_name:str,
                thrust:float, injector_pressure:float, exit_pressure:float,
                mixture_ratio:float, fac_contraction_ratio:float, characteristic_length:float,
                nozzle_diverging_angle:float):
        
        self._num_points = num_points
        self._fuel_name = fuel_name
        self._ox_name = ox_name
        self._thrust = thrust
        self._injector_pressure = injector_pressure
        self._exit_pressure = exit_pressure
        self._mixture_ratio = mixture_ratio
        self._fac_contraction_ratio = fac_contraction_ratio
        self._characteristic_length = characteristic_length
        self._nozzle_diverging_angle = nozzle_diverging_angle

        # Initialize CEA obj (FIGURE OUT UNITS LATER) 
        self.cea_obj = CEA_Obj(oxName=self._ox_name, fuelName=self._fuel_name, fac_CR=self._fac_contraction_ratio)

        self._get_main_properties()
        self._generate_contour()
        self._get_property_arrays()



    def _get_main_properties(self):
        
        # Calculate important ratios
        total_pressure_ratio = self._injector_pressure / self._exit_pressure
        throat_pressure_ratio = self.cea_obj.get_Throat_PcOvPe(self._injector_pressure, self._mixture_ratio) # injector pressure / throat pressure
        expansion_ratio = self.cea_obj.get_eps_at_PcOvPe(self._injector_pressure, self._mixture_ratio, total_pressure_ratio)

        # Determine exit velocity
        exit_mach_number = self.cea_obj.get_MachNumber(self._injector_pressure, self._mixture_ratio, expansion_ratio)
        exit_sonic_velocity = self.cea_obj.get_SonicVelocities(self._injector_pressure, self._mixture_ratio, expansion_ratio)[2]
        self._nozzle_correction_factor = 0.5 * (1 + math.cos(self._nozzle_diverging_angle * math.pi/180))
        self._exit_velocity = exit_mach_number * exit_sonic_velocity # * self.cstar_eff (NEED TO IMPLEMENT LATER)

        # Determine required mass flow rate
        self._total_mass_flow_rate = self.thrust / (self._nozzle_correction_factor * self._exit_velocity)
        self._ox_mass_flow_rate = (self._total_mass_flow_rate / (self._mixture_ratio + 1)) * self._mixture_ratio
        self._fuel_mass_flow_rate = self._total_mass_flow_rate - self._ox_mass_flow_rate 

        # Get throat conditions
        throat_temperature = self.cea_obj.get_Temperatures(self._injector_pressure, self._mixture_ratio, expansion_ratio, frozen=0, frozenAtThroat=0)[1]
        (throat_molecular_weight, throat_gamma) = self.cea_obj.get_Throat_MolWt_gamma(self._injector_pressure, self._mixture_ratio)
        universal_gas_constant = 8.314 # J/mol-K
        self._throat_area = math.sqrt(throat_temperature * universal_gas_constant / (throat_molecular_weight * throat_gamma)) * self._total_mass_flow_rate / (self._injector_pressure / throat_pressure_ratio)
        self._throat_radius = math.sqrt(self._throat_area / math.pi)



    
    def _generate_contour(self):
        # ******************************************************************
        # HARD CODED FOR TESTING: ADD CONTOUR GENERATION IMPLEMENTATION HERE
        # ******************************************************************
        self.distances = np.linspace(0, 10, self._num_points)
        self.radii = np.ones([len(self.distances)])
        for i in range(len(self.distances)//2, len(self.distances)):
            self.radii[i] = 0.2 + 0.8 * abs(i - 3*len(self.distances)//4) / (len(self.distances)//4)

        # ******************************************************************

        self.distances.setflags(write=False) # Enforce immutability
        self.radii.setflags(write=False)

        self.areas = np.pi * self.radii**2
        self.areas.setflags(write=False)

        # Determine and store converging start index by determining the first index where the radius changes
        radius_diff = np.abs(np.diff(self.radii))
        mask = radius_diff < 1e-6
        self._converging_start_index = np.argmax(mask) if np.any(mask) else None

        # Determine and store throat index by determining the index with the minimum raidus
        self._throat_index = np.argmin(self.radii)

        # Calculate area ratios = areas / throat area 
        self.area_ratios = self.areas / self.areas[self._throat_index]
        self.area_ratios.setflags(write=False)
        
        # Calculate FAC contraction ratio = area at injector / throat area
        self._fac_contraction_ratio = self.areas[0] / self.areas[self._throat_index]



    def _get_property_arrays(self):
        # Initialize data arrays
        self.pressure_ratios = np.full(self._num_points, float('nan'))
        self.mach_numbers = np.full(self._num_points, float('nan'))
        self.pressures = np.full(self._num_points, float('nan'))
        self.temperatures = np.full(self._num_points, float('nan'))
        self.densities = np.full(self._num_points, float('nan'))
        self.enthalpies = np.full(self._num_points, float('nan'))
        self.internal_energies = np.full(self._num_points, float('nan'))
        self.molecular_weights = np.full(self._num_points, float('nan'))
        self.cp_heat_capacities = np.full(self._num_points, float('nan'))
        self.heat_capacity_ratios = np.full(self._num_points, float('nan'))
        self.sonic_velocities = np.full(self._num_points, float('nan'))
        self.viscosities = np.full(self._num_points, float('nan'))
        self.conductivities = np.full(self._num_points, float('nan'))
        self.prandtl_numbers = np.full(self._num_points, float('nan'))

        # TEMP DELETE LATER AFTER CONTOUR GENERATION IS IMPLEMENTED
        self.cea_obj = CEA_Obj(oxName=self._ox_name, fuelName=self._fuel_name, fac_CR=self._fac_contraction_ratio)
        cache = {}

        # Makes a call to get_full_cea_output with the given area ratio and whether it is converging/diverging 
        def _call_cea(side: str, a_ratio: float):
            key = (side, round(float(a_ratio), 4))
            if key in cache:
                return cache[key]
            if side == "sub":
                out = self.cea_obj.get_full_cea_output(
                    self._injector_pressure,
                    self._mixture_ratio,
                    subar=a_ratio,
                    eps=None,
                    short_output=1,
                    show_transport=1,
                    output="siunits",
                )
            else:
                out = self.cea_obj.get_full_cea_output(
                    self._injector_pressure,
                    self._mixture_ratio,
                    eps=a_ratio,
                    frozenAtThroat=1,
                    short_output=1,
                    show_transport=1,
                    output="siunits",
                )
            cache[key] = out
            return out

        # Call CEA for all datapoints and organize in arrays
        converging_outputs = [
            (i, _call_cea("sub", self.area_ratios[i])) for i in range(self._converging_start_index, self._throat_index)
        ]
        diverging_outputs = [
            (i, _call_cea("sup", self.area_ratios[i])) for i in range(self._throat_index, self._num_points)
        ]

        # Parses CEA output to get numerical values for the given row label
        def _grab_row_values(full_txt: str, label: str) -> list[float]:
            _num_pat = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
            for raw in full_txt.splitlines():
                line = raw.strip()
                if line.startswith(label):
                    tail = raw.split(label, 1)[1]
                    # clean up RocketCEA formatting quirks
                    tail = re.sub(r"(\d)\s+0\s+", r"\1 ", tail)
                    tail = re.sub(r"(\d)-(\d)", r"\1e-\2", tail)
                    vals = [float(x) for x in _num_pat.findall(tail)]
                    return vals
            return []
        
        label_map = {
            "Pinj/P ": self.pressure_ratios,
            "MACH NUMBER": self.mach_numbers,
            "P, BAR": self.pressures,
            "T, K": self.temperatures,
            "RHO, KG/CU M": self.densities, 
            "H, KJ/KG": self.enthalpies,
            "U, KJ/KG": self.internal_energies,
            "M, (1/n)": self.molecular_weights,
            "Cp, KJ/(KG)(K)": self.cp_heat_capacities, 
            "GAMMAs": self.heat_capacity_ratios, 
            "SON VEL,M/SEC": self.sonic_velocities, 
            "VISC,MILLIPOISE": self.viscosities, 
            "CONDUCTIVITY  ": self.conductivities,
            "PRANDTL NUMBER": self.prandtl_numbers
        }

        # Assign CEA results to data arrays
        for i, txt in converging_outputs + diverging_outputs:
            for label, array in label_map.items():
                vals = _grab_row_values(txt, label)
                # Expected format is [INJECTOR, COMB END, THROAT, EXIT]. "RHO, KG/CU M" has weird formatting so we get an extraneous value
                if (len(vals) == 4) or (len(vals) == 5 and label == "RHO, KG/CU M"):
                    array[i] = vals[3]
                else:
                    print("Could not parse", label, "at index", i, "(got values", vals, ")")

        # Just apply linear interpolation for chamber properties
        if self._converging_start_index > 0 and converging_outputs:
            seed_txt = converging_outputs[0][1]
            z0 = float(self.distances[0])
            z1 = float(self.distances[self._converging_start_index])
            dz = (z1 - z0) if (z1 != z0) else 1.0
            for i in range(self._converging_start_index):
                s = (float(self.distances[i]) - z0) / dz
                s = max(0.0, min(1.0, s))
                for label, array in label_map.items():
                    vals = _grab_row_values(seed_txt, label)
                    if len(vals) >= 2:
                        INJECTOR, COMB_END = vals[0], vals[1]
                        array[i] = (1.0 - s) * INJECTOR + s * COMB_END

        # Make data arrays readonly
        self.pressure_ratios.setflags(write=False)
        self.mach_numbers.setflags(write=False)
        self.pressures.setflags(write=False)
        self.temperatures.setflags(write=False)
        self.densities.setflags(write=False)
        self.enthalpies.setflags(write=False)
        self.internal_energies.setflags(write=False)
        self.molecular_weights.setflags(write=False)
        self.cp_heat_capacities.setflags(write=False)
        self.heat_capacity_ratios.setflags(write=False)
        self.sonic_velocities.setflags(write=False)
        self.viscosities.setflags(write=False)
        self.conductivities.setflags(write=False)
        self.prandtl_numbers.setflags(write=False)

        


# Testing
if __name__ == "__main__":

    engine = Engine(num_points=50, fuel_name="RP1", ox_name="LOX", thrust=1000, injector_pressure=500, exit_pressure=10,
                    mixture_ratio=2.6, fac_contraction_ratio=1, characteristic_length=1, nozzle_diverging_angle=1)

    plt.plot(engine.distances, engine.radii)
    plt.plot(engine.distances, engine.mach_numbers)
    plt.plot(engine.distances, engine.pressures)
    plt.plot(engine.distances, engine.densities)
    plt.show()
