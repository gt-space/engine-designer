from dataclasses import dataclass, field
from typing import Any, Dict, Optional
from injector import Injector
from inputs import DesignChoices, EngineInputs, PropellantInputs, FeedSystemInputs
import math

@dataclass
class ImpingingDesign(DesignChoices):
    N_elements: int = 16
    impingement_angle: int = 60 #degrees
    
    LOX_Cd: float = 0.85
    LOX_orifice_L_to_d: int = 6
    LOX_impinging_L_to_d: int = 7
    
    fuel_Cd: float = 0.85
    fuel_orifice_L_to_d: int = 6

@dataclass(frozen=True)
class ImpingingSizingResult:
    LOX_orifices: int = field(metadata={"unit": ""})
    LOX_total_orifice_area: float = field(metadata={"unit": "in^2"})
    LOX_single_orifice_area: float = field(metadata={"unit": "in^2"})
    LOX_orifice_diameter: float = field(metadata={"unit": "in"})
    LOX_orifice_length: float = field(metadata={"unit": "in"})
    LOX_orifice_height: float = field(metadata={"unit": "in"})
    LOX_impinging_length: float = field(metadata={"unit": "in"})
    LOX_orifice_separation_entrance: float = field(metadata={"unit": "in"})
    LOX_orifice_separation_exit: float = field(metadata={"unit": "in"})

    fuel_orifices: int = field(metadata={"unit": ""})
    fuel_total_orifice_area: float = field(metadata={"unit": "in^2"})
    fuel_single_orifice_area: float = field(metadata={"unit": "in^2"})
    fuel_orifice_diameter: float = field(metadata={"unit": "in"})
    fuel_orifice_length: float = field(metadata={"unit": "in"})
    fuel_orifice_length_unbounded: float = field(metadata={"unit": "in"})
    fuel_impinging_length: float = field(metadata={"unit": "in"})

class ImpingingInjector(Injector):
    def __init__(
        self,
        engine: EngineInputs,
        prop: PropellantInputs,
        feed: FeedSystemInputs,
        design: ImpingingDesign,
    ) -> None:
        super().__init__(engine, prop, feed, design)
        self.imp = design
        self._result: Optional[ImpingingSizingResult] = None

    def size(self) -> ImpingingSizingResult:

        LOX_orifices = self.imp.N_elements * 2

        #A = mdot / (Cd * sqrt(2 * rho * dP))
        #(39.3701**2) * (1*10**(-4)): convert m^2 to in^2
        LOX_total_orifice_area = (39.3701**2) * (1*10**(-4)) * self.prop.LOX_m_dot / (self.imp.LOX_Cd * math.sqrt(2 * self.prop.LOX_density * self.feed.LOX_pressure_drop))
        LOX_single_orifice_area = LOX_total_orifice_area / LOX_orifices
        LOX_orifice_diameter = 2 * math.sqrt((LOX_single_orifice_area) / math.pi)
        LOX_orifice_length = self.imp.LOX_orifice_L_to_d * LOX_orifice_diameter
        LOX_orifice_height = LOX_orifice_length * math.cos(math.radians(self.imp.impingement_angle / 2))
        LOX_impinging_length = self.imp.LOX_impinging_L_to_d * LOX_orifice_diameter

        #2 * (Lio + Lo) * sin(theta/2)
        LOX_orifice_separation_entrance = 2 * (LOX_impinging_length + LOX_orifice_length) * math.sin(math.radians(self.imp.impingement_angle / 2))
        #2 * Lio * sin(theta/2)
        LOX_orifice_separation_exit = 2 * LOX_impinging_length * math.sin(math.radians(self.imp.impingement_angle / 2))

        fuel_orifices = self.imp.N_elements

        #A = mdot / (Cd * sqrt(2 * rho * dP))
        #(39.3701**2) * (1*10**(-4)): convert m^2 to in^2
        fuel_total_orifice_area = (39.3701**2) * (1*10**(-4)) * self.prop.fuel_m_dot / (self.imp.fuel_Cd * math.sqrt(2 * self.prop.fuel_density * self.feed.fuel_pressure_drop))
        fuel_single_orifice_area = fuel_total_orifice_area / fuel_orifices
        fuel_orifice_diameter = 2 * math.sqrt((fuel_single_orifice_area) / math.pi)
        fuel_orifice_length = LOX_orifice_height
        fuel_orifice_length_unbounded = self.imp.fuel_orifice_L_to_d * fuel_orifice_diameter
        fuel_impinging_length = LOX_impinging_length * math.cos(math.radians(self.imp.impingement_angle / 2))
        
        result = ImpingingSizingResult(LOX_orifices=LOX_orifices, LOX_total_orifice_area=LOX_total_orifice_area, 
                                        LOX_single_orifice_area=LOX_single_orifice_area, 
                                        LOX_orifice_diameter=LOX_orifice_diameter, LOX_orifice_length=LOX_orifice_length, 
                                        LOX_orifice_height=LOX_orifice_height, LOX_impinging_length=LOX_impinging_length,
                                        LOX_orifice_separation_entrance=LOX_orifice_separation_entrance, 
                                        LOX_orifice_separation_exit=LOX_orifice_separation_exit,
                                        fuel_orifices=fuel_orifices, fuel_total_orifice_area=fuel_total_orifice_area, 
                                        fuel_single_orifice_area=fuel_single_orifice_area, fuel_orifice_diameter=fuel_orifice_diameter, 
                                        fuel_orifice_length=fuel_orifice_length, fuel_orifice_length_unbounded=fuel_orifice_length_unbounded, 
                                        fuel_impinging_length=fuel_impinging_length)
        
        self._result = result
        return result

    def check_constraints(self) -> None:
        print("[ImpingingInjector.check_constraints] Placeholder checks passed.")