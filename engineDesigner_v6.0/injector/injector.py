from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any, Dict, Tuple
from dataclasses import asdict, fields

from inputs import EngineInputs, PropellantInputs, FeedSystemInputs, DesignChoices


class Injector(ABC):
    """Abstract base injector class."""

    def __init__(
        self,
        engine: EngineInputs,
        prop: PropellantInputs,
        feed: FeedSystemInputs,
        design: DesignChoices,
    ) -> None:
        self.engine = engine
        self.prop = prop
        self.feed = feed
        self.design = design

    def run(self) -> Dict[str, Any]:
        # self.validate_inputs()
        geom = self.size()
        # self.check_constraints()
        return self.report(geom)

    def validate_inputs(self) -> None:
        #Catch errors or missing data BEFORE running size()
        print("[Injector.validate_inputs] This is the base validation (placeholder).")

    @abstractmethod
    def size(self) -> Dict[str, Any]:
        #Main method to size the injector
        raise NotImplementedError

    @abstractmethod
    def check_constraints(self) -> None:
        #Validate the results AFTER running size(). Sanity check and manufacturable
        raise NotImplementedError


    def report(self, geometry) -> Dict[str, Any]:
        formatted = {}
        inserted_break = False

        for f in fields(geometry):
            value = getattr(geometry, f.name)
            
            if f.name.startswith("fuel_") and not inserted_break:
                print()
                inserted_break = True
            
            if isinstance(value, float):
                value = round(value, 4)

            unit = f.metadata.get("unit", "")
            formatted_value = f"{value} {unit}".strip()
            formatted[f.name] = formatted_value
            print(f"{f.name}: {formatted_value}")

        return formatted