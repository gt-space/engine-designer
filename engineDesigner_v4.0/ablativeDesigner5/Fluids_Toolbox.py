from pyfluids import Fluid, FluidsList

# https://github.com/portyanikhin/PyFluids
p = 140
p_Pa = 6894.76 * p
nitrogen = Fluid(FluidsList.Nitrogen).dew_point_at_pressure(p_Pa)
print(str(nitrogen.temperature+273.15) + "deg K ")
print(nitrogen.units_system)