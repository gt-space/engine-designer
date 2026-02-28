from dataclasses import dataclass
from typing import Any, Dict, Literal

from injector import Injector
from inputs import DesignChoices, EngineInputs, PropellantInputs, FeedSystemInputs


@dataclass
class CoaxDesign(DesignChoices):
    coax_type: Literal["shear", "swirl"] = "shear"


class CoaxInjector(Injector):
    def __init__(
        self,
        engine: EngineInputs,
        prop: PropellantInputs,
        feed: FeedSystemInputs,
        design: CoaxDesign,
    ) -> None:
        super().__init__(engine, prop, feed, design)
        self.coax = design

    def size(self) -> Dict[str, Any]:
        print("[CoaxInjector.size] Placeholder sizing.")
        return {
            "coax_type": self.coax.coax_type,
            "example_inner_d_m": 0.002,
        }

    def check_constraints(self) -> None:
        print("[CoaxInjector.check_constraints] Placeholder checks passed.")