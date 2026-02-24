import numpy as np
from typing import Union, List, Dict, Optional, Any
import warnings

# --- 1. The Smart Property Wrapper ---
class Prop:
    """
    Represents a single physical property (e.g., Yield Strength).
    Can be a constant float or a temperature-dependent lookup table.
    """
    def __init__(self, name: str, data: Union[float, List[np.ndarray]], units: str = ""):
        self.name = name
        self.units = units
        self._data = data
        self._is_constant = isinstance(data, (float, int))
        
        # Validation
        if not self._is_constant:
            if not isinstance(data, list) or len(data) != 2:
                raise ValueError(f"Property '{name}' must be a float or [Temp_Array, Value_Array]")

    def get(self, T_kelvin: float) -> float:
        """Get value at specific temperature."""
        if self._is_constant:
            return float(self._data)
        
        temps, values = self._data[0], self._data[1]
        # np.interp handles the linear interpolation safely
        return float(np.interp(T_kelvin, temps, values))

    def __repr__(self):
        val_str = f"{self._data}" if self._is_constant else "Temp_Dependent_Array"
        return f"<{self.name}: {val_str} {self.units}>"

# --- 2. The Fatigue Handler ---
class FatigueProfile:
    """
    Handles S-N curves (Stress vs Cycles).
    Ideally, you have different curves for different temperatures.
    This implementation picks the curve closest to the requested Temp.
    """
    def __init__(self, curves: Dict[float, List[np.ndarray]]):
        # Format: { Temperature_K : [Cycles_Array, Stress_Array] }
        self.curves = curves

    def get_limit(self, cycles: float, T_kelvin: float) -> float:
        """Returns max stress for a given number of cycles at Temp T."""
        # 1. Find the closest temperature curve available
        avail_temps = np.array(list(self.curves.keys()))
        idx = (np.abs(avail_temps - T_kelvin)).argmin()
        closest_T = avail_temps[idx]
        
        if abs(closest_T - T_kelvin) > 50:
            warnings.warn(f"Fatigue data unavailable for {T_kelvin}K. Using curve for {closest_T}K.")

        # 2. Interpolate Stress vs Cycles for that curve
        cycle_arr, stress_arr = self.curves[closest_T]
        
        # Log-Log interpolation is usually standard for S-N, but Linear is used here for simplicity
        return float(np.interp(cycles, cycle_arr, stress_arr))

# --- 3. The Material Class ---
class Material:
    def __init__(self, name: str, category: str = "General", default_condition: str = "Standard"):
        self.name = name
        self.category = category
        self.default_condition = default_condition
        
        # Structure: { prop_name: { condition_name: Prop_Object } }
        self.properties: Dict[str, Dict[str, Prop]] = {}
        
        # Structure: { condition_name: FatigueProfile }
        self.fatigue: Dict[str, FatigueProfile] = {}
        
        self.metadata: Dict[str, Any] = {} 

    def add_prop(self, key: str, data, units: str = "", condition: str = None):
        """Add a standard table/value property for a specific condition."""
        cond = condition if condition else self.default_condition
        
        # Initialize the dictionary for this property if it doesn't exist
        if key not in self.properties:
            self.properties[key] = {}
            
        # Create and store the Prop object
        self.properties[key][cond] = Prop(key, data, units)

    def add_custom_prop(self, prop_object: Prop, condition: str = None):
        """
        Use this for your NIST/Polynomial properties.
        """
        cond = condition if condition else self.default_condition
        key = prop_object.name
        
        if key not in self.properties:
            self.properties[key] = {}
            
        self.properties[key][cond] = prop_object

    def add_fatigue(self, curve_data: Dict[float, List[np.ndarray]], condition: str = None):
        """Add S-N curve data for a specific condition."""
        cond = condition if condition else self.default_condition
        self.fatigue[cond] = FatigueProfile(curve_data)

    def add_meta(self, key: str, value):
        """Metadata is usually static across conditions, but you can overwrite if needed."""
        self.metadata[key] = value

    def get(self, prop_name: str, T: float = 298.0, condition: str = None) -> float:
        """
        Get value. If condition is NOT provided, uses the material's default.
        """
        # 1. Determine which condition to look for
        target_cond = condition if condition else self.default_condition
        
        # 2. Check if property exists
        if prop_name not in self.properties:
            raise KeyError(f"Material '{self.name}' has no property '{prop_name}'")
            
        # 3. Check if the specific condition exists for that property
        if target_cond not in self.properties[prop_name]:
            raise KeyError(f"Property '{prop_name}' found, but data for condition '{target_cond}' is missing.")

        # 4. Fetch and Calculate
        return self.properties[prop_name][target_cond].get(T)

    def __repr__(self):
        return f"Material(Name='{self.name}', Default='{self.default_condition}', Props={len(self.properties)})"

# --- 4. The Database Registry ---
class MaterialRegistry:
    def __init__(self):
        self._db: Dict[str, Material] = {}

    def add_material(self, material: Material):
        # Store by a unique key, typically Name + Condition
        key = f"{material.name}_{material.default_condition}".replace(" ", "_").lower()
        self._db[key] = material
        print(f"Registered: {material.name} ({material.default_condition})")

    def get_material(self, name_key: str) -> Material:
        return self._db.get(name_key.lower())

    def list_materials(self):
        return list(self._db.keys())
