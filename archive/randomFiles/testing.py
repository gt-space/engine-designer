# this was a preliminary experimental script. Probably not of much use

pRatStations = [1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.7,1.8,1.9]

# pip is the ratio between injector pressure and local pressure
pip = np.array([1.00000])
pip = np.append(pip, [pip_combEnd, pip_t])
pip = np.append(pip, pRatStations)
A_t = ispObj.get_eps_at_PcOvPe(P_c_psi, MR, PcOvPe=pip_t)

aeat = np.array([0, conRat, 1.0000])
M_combEnd = ispObj.get_Chamber_MachNumber(P_c_psi, MR)
mach = np.array([0, M_combEnd, 1.0000])

resolution = 10 #How many points along the engine?
# Main combustion properties:
for item in pRatStations:
    s = ispObj.get_full_cea_output(P_c, MR, PcOvPe=item, subar=[3, 2, 1], short_output=1, pc_units='bar')
    aRat = ispObj.get_eps_at_PcOvPe(P_c_psi, MR, PcOvPe=item)
    # print(aRat)
    aeat = np.append(aeat, [aRat])
    M = ispObj.get_MachNumber(P_c_psi, MR, aRat)
    # print(M)
    mach = np.append(mach, [M])
    C_f = ispObj.get_PambCf(14.7, P_end, MR, aRat)
    # cf = np.array([0, , ])
    x = ispObj.get_PcOvPe(P_c_psi, MR, eps=aRat)

    # print(s)
    # print(x)
    # print(C_f)
    # print(M)

    # print(pip_combEnd)
    # conList.append(s)
    # print(s.find("Pinf"))

# print(len())
# print(mach)
# subar=[3, 2, 1.5]
s = ispObj.get_full_cea_output(P_c, MR, subar=[5], short_output=1, pc_units='bar')
# print(s)
