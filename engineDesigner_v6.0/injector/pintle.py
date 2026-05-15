from dataclasses import dataclass
from typing import Any, Dict, Literal

from injector import Injector
from inputs import DesignChoices, EngineInputs, PropellantInputs, FeedSystemInputs


@dataclass
class PintleDesign(DesignChoices):
    pintle_type: Literal["lox_center", "fuel_center"] = "lox_center"


class PintleInjector(Injector):
    def __init__(
        self,
        engine: EngineInputs,
        prop: PropellantInputs,
        feed: FeedSystemInputs,
        design: PintleDesign,
    ) -> None:
        super().__init__(engine, prop, feed, design)
        self.pint = design

    def size(self) -> Dict[str, Any]:
        print("[PintleInjector.size] Placeholder sizing.")
        return {
            "pintle_type": self.pint.pintle_type,
            "example_gap_m": 0.0005,
        }

    def check_constraints(self) -> None:
        print("[PintleInjector.check_constraints] Placeholder checks passed.")