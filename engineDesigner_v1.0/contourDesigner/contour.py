# CONTOUR.PY – Does the geoemetry to construct a 200 point contour describing the
# engine given various input parameters

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

def get_contour(R_t, con_rat, exp_rat, L_star, adv_data):
    # Generate Nozzle Contour given the following input parameters:
    # R_t - Throat Radius (m)
    # con_rat - Contraction Ratio for chamber
    # exp_rat - Expansion ratio for exit
    # L_star - L* Value (derived from Propellants)
    # adv_data - Dictionary with additional sizing parameters

    # Unpack advanced parameters from adv_data dictionary
    div_ang = adv_data["div_ang"] * math.pi / 180 # Nozzle exit half angle (conical) [degrees to rad]
    con_ang = adv_data["con_ang"] * math.pi / 180 # Nozzle contraction angle [degrees to rad]

    throatLeadInRadius = adv_data["lead_in_factor"] * R_t # Throat lead in radius after contraction angle
    throatLeadOutRadius = adv_data["lead_out_factor"] * R_t # Nozzle lead in radius

    rad_rat = adv_data["rad_rat"] # Radius Ratio conLeadInRadius/conLdRadMax for chamber convergence section lead in radius

    solve_bell = adv_data["solve_bell"] # Bool for solving bell vs conical diverging section
    theta_i = adv_data["theta_i"] * math.pi / 180 # Angle leaving throat [deg to rad] (Should be between 20 and 50)
    theta_e = adv_data["theta_e"] * math.pi / 180 # Exit angle [deg to rad] (Sutton recommends <= 10)
    percent_of_conical = adv_data["percent_of_conical"] # Percent length compared to conical alternative (Should be ~80%)

    R_tCurve = (throatLeadInRadius + throatLeadOutRadius)/2 # Throat radius of curvature (for Bartz)
    throatArea = math.pi * (R_t**2) # Cross-sectional area of the throat [m^2]
    chamberArea = con_rat * throatArea
    chamberRad = math.sqrt(con_rat) * R_t


    # ==== Nozzle Critical Points ====
    #     - 0: Chamber Start, Injector Plane
    #     - 1: Converging radius start
    #     - 2: Converging radius end, angular section start (imagine a frustum)
    #     - 3: Angular section end, throat lead in start
    #     - 4: Throat - End of throat lead in (min radius)
    #     - 5: Start of major nozzle contour
    #     - 6: End of Nozzle

    # ==== Size Angular Section Lead in ====
    throatLeadRDist = throatLeadInRadius-(math.cos(con_ang) * throatLeadInRadius) + R_t #R distance occupied by throat lead in
    convLdRDelta = chamberRad - throatLeadRDist #Delta R for converging lead in
    conLdRadMax = convLdRDelta/(1 - math.cos(con_ang)) #Max lead in rad for meet @con_ang
    conLeadInRadius = rad_rat * conLdRadMax
    a1r = chamberRad - conLeadInRadius #Radial center point for lead out arc
    r0 = chamberRad

    # Chamber lead out start
    z1 = 0
    r1 = chamberRad
    # Chamber lead out end, Angular contraction start
    z2 = z1 + conLeadInRadius * math.sin(con_ang) #axial coordinate of end of lead in rad
    r2 = chamberRad - conLeadInRadius * (1 - math.cos(con_ang)) #Set radial height of ending coord
    a2r = throatLeadInRadius + R_t #Radial center point for lead in
    # Angular contraction end,throat lead in start
    r3 = throatLeadRDist
    z3 = (r3 - r2)/(-math.tan(con_ang)) + z2
    a3r = throatLeadOutRadius + R_t
    # Throat
    z4 = z3 + throatLeadInRadius * math.sin(con_ang)
    r4 = R_t
    # Throat lead out end
    z5 = z4 + (throatLeadOutRadius * math.sin(div_ang))
    r5 = throatLeadOutRadius - math.sqrt(throatLeadOutRadius**2-(throatLeadOutRadius * math.sin(div_ang))**2) + R_t
    # Nozzle Sizer
    # Get conic length
    r6 = math.sqrt(exp_rat) * R_t
    z6 = z5 + (r6 - r5)/(math.tan(div_ang))

    # Get bell length if desired
    if solve_bell:
        z6 = z5 + ((r6 - r5)/(math.tan(div_ang)) * (percent_of_conical/100)) # Scale length
        z5 = z4 + (throatLeadOutRadius * math.sin(theta_i)) # Change lead out end axial location
        r5 = throatLeadOutRadius - math.sqrt(throatLeadOutRadius**2-(throatLeadOutRadius * math.sin(theta_i))**2) + R_t # Change lead out end radial loacation

    # ==== Generate Nozzle Contour Vector ====
    nozzleLength = z6 #Total Length of Nozzle
    numPTS = 100 # Number of points defining each barrel and nozzle
    nozzleZ = np.linspace(0, nozzleLength, num=numPTS)
    nozzleR = np.zeros(numPTS)
    for i in range(len(nozzleZ)):
        if z1 <= nozzleZ[i] and nozzleZ[i] < z2:
            #Chamber Lead Out Rad
            nozzleR[i] = (chamberRad - conLeadInRadius) + math.sqrt(conLeadInRadius**2 - nozzleZ[i]**2)
        elif z2 <= nozzleZ[i] and nozzleZ[i] < z3:
            nozzleR[i] = r2-(nozzleZ[i]-z2) * math.tan(con_ang)
        elif z3 <= nozzleZ[i] and nozzleZ[i] < z4:
            nozzleR[i] = R_t + throatLeadInRadius - math.sqrt(throatLeadInRadius**2 - (z4 - nozzleZ[i])**2)
        elif z4 <= nozzleZ[i] and nozzleZ[i] < z5:
            nozzleR[i] = R_t + throatLeadOutRadius - math.sqrt(throatLeadOutRadius**2 - (nozzleZ[i] - z4)**2)
        elif z5 <= nozzleZ[i] and nozzleZ[i] <= z6:
            if solve_bell:
                # Function defining parabola:
                def f(r):
                    a = (1/math.tan(theta_e) - 1/math.tan(theta_i))/(2*(r6-r5))
                    b = 1/math.tan(theta_i)
                    return a * (r**2) + b * r - (nozzleZ[i] - z5)
                # Solve the root to get radius at that point
                nozzleR[i] = r5 + root(f, 0.04).x[0]
            else:
                nozzleR[i] = r5 + (nozzleZ[i] - z5) * math.tan(div_ang)
        else:
            nozzleR[i] = 0

    # ==== Size Combustion Chamber ====
    chamberVolume = throatArea * L_star
    # Solve for convergent volume:
    convergingVolume = 0
    nozzleZ = np.append(nozzleZ, z4) # Add the throat position to the Z points
    nozzleZ = np.sort(nozzleZ) # Put it in its place
    throatInd = np.where(nozzleZ == z4)[0][0] #Find throat index
    # print(throatInd)
    nozzleR = np.insert(nozzleR, throatInd, r4)
    for i in range(1,throatInd+1):
        deltaVol = math.pi * (nozzleZ[i] - nozzleZ[i-1]) * ((nozzleR[i] + nozzleR[i - 1])/2)**2
        convergingVolume = convergingVolume + deltaVol
    cylVol = chamberVolume - convergingVolume
    cylLength = cylVol/(math.pi * (chamberRad**2))
    chBarrelZ = np.linspace(0, cylLength, numPTS)
    chBarrelR = np.ones(numPTS) * chamberRad

    # ==== Shift Points and Generate Outputs ====
    zCrit = [z1,z2,z3,z4,z5,z6]
    zCrit = [x+cylLength for x in zCrit]
    zCrit = [0] + zCrit #Critical Z points per defs 1-6 above
    rCrit = [r0,r1,r2,r3,r4,r5,r6]
    chBarrel = np.column_stack((chBarrelR, chBarrelZ)) #Points for chamber barrel
    nozzleZ = nozzleZ + cylLength
    nozzleContour = np.column_stack((nozzleR, nozzleZ)) #Points for nozzle section
    nozzleContour = np.delete(nozzleContour, 0, 0) #Drop first point so it isn't repeated
    engineContour = np.concatenate((chBarrel,nozzleContour), axis=0) #Combine into one index

    # Plot the Contour if you want:
    # plt.plot(engineContour[:, 1], engineContour[:, 0])
    # plt.show()

    return (engineContour, chBarrel, nozzleContour, R_tCurve, throatInd, conLeadInRadius)
