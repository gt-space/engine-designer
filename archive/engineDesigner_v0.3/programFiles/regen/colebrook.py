# F = COLEBROOK(R,K) fast, accurate and robust computation of the
#     Darcy-Weisbach friction factor F according to the Colebrook equation:
#                             -                       -
#      1                     |    K        2.51        |
#  ---------  =  -2 * Log_10 |  ----- + -------------  |
#   sqrt(F)                  |   3.7     R * sqrt(F)   |
#                             -                       -

def colebrook(Re, d_hyd, k):
    import numpy
    from scipy.optimize import root

    def f(x):
        return -2*numpy.log10(k/(3.71*d_hyd) + 2.51/(Re*numpy.sqrt(x))) - 1.0/numpy.sqrt(x)

    fric = root(f, 0.04).x[0]
    return fric
