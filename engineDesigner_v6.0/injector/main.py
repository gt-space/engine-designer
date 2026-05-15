from inputs import EngineInputs, PropellantInputs, FeedSystemInputs
from impinging import ImpingingDesign, ImpingingInjector
# from .coax import CoaxDesign, CoaxInjector
# from .pintle import PintleDesign, PintleInjector

# python3 engineDesigner_v6.0/injector/main.py
def main() -> None:
    engine = EngineInputs
    prop = PropellantInputs
    feed = FeedSystemInputs
    impDesign = ImpingingDesign
    
    print("\n" + "=" * 60, "\n")
    impinging = ImpingingInjector(engine, prop, feed, impDesign)
    impinging.run()
    print("\n" + "=" * 60, "\n")


if __name__ == "__main__":
    main()