from CoolProp.CoolProp import PropsSI

#NEED TO UPDATE. coolprop doesn't have RP-1, so n-Dodecane is used as substitute
fluid = "n-Dodecane"

def get_coolant_h(T, P):
    return PropsSI("Hmass", "T", T, "P", P, fluid)


def get_coolant_props(P, h):
    T = PropsSI("T", "P", P, "Hmass", h, fluid)
    rho = PropsSI("Dmass", "P", P, "Hmass", h, fluid)
    mu = PropsSI("V", "P", P, "Hmass", h, fluid)
    k = PropsSI("L", "P", P, "Hmass", h, fluid)
    cp = PropsSI("Cpmass", "P", P, "Hmass", h, fluid)

    Pr = cp * mu / k

    return {
            "T": T,
            "rho": rho,
            "mu": mu,
            "k": k,
            "cp": cp,
            "Pr": Pr,
        }


