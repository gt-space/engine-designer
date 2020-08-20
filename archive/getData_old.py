# Originally I had planned to use the many functions that RocketCEA has to offer
# to gather data, as they give values to more decimal places and would not require
# complex string parsing. However, these functions do not offer querying of all
# the data in the main output string. As a result, I sill had to intermingle string
# parsing to get all of the data and accont for unit conversions. The added
# complexity made me decide to just go full stirng parsing so I could make the code
# cleaner as a loop.

# This old code does have the benefit of having more decial places for some values,
# but the worst the main output gets is 4 sig figs, which itself is at worst 0.05% off.
# I have assumed we are willing to tolerate that, so this code has been archived.

# NOTE: THE CODE IS ALSO INCOMPLETE

# ==== Calculate Properties ====
# The array generated shall have the following header format:
# R, Z, cRat, pip, aeat, mach, cf, ivac, isp, p, t, rho, h, u, mw, cp, gam, son,
#   vis, cond, condfz, pran, pranfz
# A description of each of these parameters can be found at:
# https://cearun.grc.nasa.gov/cgi-bin/CEARUN/plotParameters.cgi

# Helper function for parsing the large output string. It might be possible to
# optimize all the loops, but for now it seems to work well enough.
def getValue(str, i=1):
    # i is the col you're trying to get from the row of string (first item is 1)
    ind = output.find(str) + len(str)
    # l is the length of the value we're finding
    l = 1
    # Some numbers have exponents, so that must be taken into account
    exp = 0
    while i > 1:
        # When we're not yet at the point
        while output[ind] == " ":
            ind += 1
        # passing some earlier numbers
        while output[ind] != " ":
            ind += 1
        i -= 1
    # figure out and return value if you're at its index:
    while output[ind] == " ":
        ind += 1
    while output[ind + l] != " ":
        l += 1
    return float(output[ind:ind+l])

# Find the Chamber Barrel Properties. Idealized as linear from inj to comb end
# Chamber barrel length
numPTS = len(chBarrel)
# full CEA output (used when other CEA methods are not available to gather specific data)
output = ispObj.get_full_cea_output(P_inj * 14.5038, MR, short_output=1, output="siunits")
print(output)
# Define chamber props while adding pip column to the data
pip_combEnd = getValue("Pinj/P ", 2)
chamberProps = np.column_stack((chBarrel, np.linspace(1, pip_combEnd, numPTS)))
# Add aeat col
chamberProps = np.column_stack((chamberProps, np.ones(numPTS) * 6))
# Add mach col
chamberProps = np.column_stack((chamberProps, np.linspace(0, ispObj.get_Chamber_MachNumber(P_inj_psi, MR), numPTS)))
# Add cf col
CF = getValue("CF")
chamberProps = np.column_stack((chamberProps, np.linspace(0, CF, numPTS)))
# Add ivac col
Ivac = getValue("Ivac, M/SEC")
chamberProps = np.column_stack((chamberProps, np.linspace(0, Ivac, numPTS)))
# Add isp col
isp = C_star * CF / 0.975
chamberProps = np.column_stack((chamberProps, np.linspace(0, isp, numPTS)))
# Add p col
P_end = (P_inj / pip_combEnd) #Pressure at combustion end (bar)
chamberProps = np.column_stack((chamberProps, np.linspace(P_inj, P_end, numPTS)))
# Add t col
t_inj = getValue("T, K")
t_end = ispObj.get_Temperatures(P_inj_psi, MR)[0] * 0.555556 # Convert R to K
chamberProps = np.column_stack((chamberProps, np.linspace(t_inj, t_end, numPTS)))
# Add rho col
rho_inj = getValue("RHO, KG/CU M")
rho_end = ispObj.get_Chamber_Density(P_inj_psi, MR) * 16.0185  # Convert lb/ft^3 to kg/m^3
chamberProps = np.column_stack((chamberProps, np.linspace(rho_inj, rho_end, numPTS)))
# Add h col
h_inj = getValue("H, KJ/KG")
h_end = ispObj.get_Chamber_H(P_inj_psi, MR) * 2.326 # Convert btu/lb to kJ/kg
chamberProps = np.column_stack((chamberProps, np.linspace(h_inj, h_end, numPTS)))
# Add u col
u_inj = getValue("U, KJ/KG")
u_end = getValue("U, KJ/KG", 2)
chamberProps = np.column_stack((chamberProps, np.linspace(u_inj, u_end, numPTS)))
# Add mw col
mw_inj = getValue("M, (1/n)")
mw_end = ispObj.get_Chamber_MolWt_gamma(P_inj_psi, MR)[0]
chamberProps = np.column_stack((chamberProps, np.linspace(mw_inj, mw_end, numPTS)))
# Add cp col
cp_inj = getValue("Cp, KJ/(KG)(K)")
cp_end = ispObj.get_Chamber_Cp(P_inj_psi, MR) * 4.1868 # btu/(lb R) to KJ/(kg K)
print(cp_end)
chamberProps = np.column_stack((chamberProps, np.linspace(cp_inj, cp_end, numPTS)))
# Add gam col
# Add son col
print(chamberProps)
