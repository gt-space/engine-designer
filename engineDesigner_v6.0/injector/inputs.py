from dataclasses import dataclass
from typing import Optional


@dataclass
class EngineInputs:
    Pc: float = 21.6 #bar
    Dc: float = 0.105 #m
    MR: float = 2 # LOX/fuel


@dataclass
class PropellantInputs:
    LOX_density: float = 1.141 #kg/L
    LOX_temperature: float = 90 #K
    LOX_m_dot: float = 3.24 #kg/s
    LOX_volume_flow: float = 2.8396 #L/s

    fuel_density: float = 0.8 #kg/L
    fuel_temperature: float = 293 #K
    fuel_m_dot: float = 1.62 #kg/s
    fuel_volume_flow: float = 2.0250 #L/s


@dataclass
class FeedSystemInputs:
    LOX_inlet_pressure: float = 25.92 #bar
    LOX_pressure_drop: float = 4.32 #bar

    fuel_inlet_pressure: float = 25.92 #bar
    fuel_pressure_drop: float = 4.32 #bar


@dataclass
class DesignChoices:
    pass