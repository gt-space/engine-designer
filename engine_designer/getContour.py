import math
import numpy as np
import matplotlib.pyplot as plt

def getContour(R_t, LStar, conRat, conAng, divAng, radRat, expRat):
    # Generate Nozzle Contour given the following input parameters:
    # R_t - Throat Radius (m)
    # LStar - L* Value (derived from Propellants)
    # conRat - Contraction Ratio for chamber
    # conAng - Convergence Angle for chamber
    # divAng - Divergence Angle for exit
    # radRat - Radius Ratio conLeadInRadius/conLdRadMax for chamber convergence section lead in radius,
    # where R_cyl is the cylindrical chamber radius
    # expRat - Expansion ratio for exit

    # ==== Size Chamber Parameters ====
    exitHalfAngle = divAng * math.pi / 180 #Nozzle Exit half angle, degrees to rad
    throatLeadInRadius = 1.5 * R_t #Throat lead in radius after contraction angle
    throatLeadOutRadius = 0.382 * R_t #Nozzle lead in radius
    R_tCurve = (throatLeadInRadius + throatLeadOutRadius)/2 #Throat radius of curvature (for Bartz)
    throatArea = math.pi * (R_t**2) #Area of the throat
    chamberArea = conRat * throatArea
    chamberRad = math.sqrt(conRat)*R_t

    # ==== Nozzle Critical Points ====
    # (z  h)^2 + (r  k)^2 = r^2
    # Critical Pts:
    #     - 0: Chamber Start, Injector Plane
    #     - 1: Converging radius start
    #     - 2: Converging radius end, angular section start
    #     - 3: Angular section end, nozzle lead in start
    #     - 4: Throat - Start of nozzle lead in (minor contour)
    #     - 5: Start of major nozzle contour
    #     - 6: End of Nozzle

    # ==== Size Angular Section Lead in ====
    throatLeadRDist = throatLeadInRadius-(math.cos(conAng) * throatLeadInRadius) + R_t #R distance occupied by throat lead in
    convLdRDelta = chamberRad - throatLeadRDist #Delta R for converging lead in
    conLdRadMax = convLdRDelta/(1 - math.cos(conAng)) #Max lead in rad for meet @conAng
    conLeadInRadius = radRat * conLdRadMax
    a1r = chamberRad - conLeadInRadius #Radial center point for lead out arc
    r0 = chamberRad
    # Chamber lead out start
    z1 = 0
    r1 = chamberRad
    # Chamber lead out end, Angular contraction start
    z2 = z1 + conLeadInRadius * math.sin(conAng) #axial coordinate of end of lead in rad
    r2 = chamberRad - conLeadInRadius * (1 - math.cos(conAng)) #Set radial height of ending coord
    a2r = throatLeadInRadius + R_t #Radial center point for lead in
    # Angular contraction end,throat lead in start
    r3 = throatLeadRDist
    z3 = (r3 - r2)/(-math.tan(conAng)) + z2
    a3r = throatLeadOutRadius + R_t
    # Throat
    z4 = z3 + throatLeadInRadius * math.sin(conAng)
    r4 = R_t
    # Throat lead out start
    z5 = z4 + (throatLeadOutRadius * math.sin(exitHalfAngle))
    r5 = throatLeadOutRadius - math.sqrt(throatLeadOutRadius**2-(throatLeadOutRadius * math.sin(exitHalfAngle))**2) + R_t
    # Nozzle Sizer
    # Currently conic only - mabye parabolic in the future.
    r6 = math.sqrt(expRat) * R_t
    z6 = z5 + (r6 - r5)/(math.tan(exitHalfAngle))

    # ==== Generate Nozzle Contour Vector ====
    nozzleLength = z6 #Total Length of Nozzle
    numPTS = 100
    nozzleZ = np.linspace(0, nozzleLength, num=numPTS)
    nozzleR = np.zeros(numPTS)
    for i in range(len(nozzleZ)):
        if z1 <= nozzleZ[i] and nozzleZ[i] < z2:
            #Chamber Lead Out Rad
            nozzleR[i] = (chamberRad - conLeadInRadius) + math.sqrt(conLeadInRadius**2 - nozzleZ[i]**2)
        elif z2 <= nozzleZ[i] and nozzleZ[i] < z3:
            nozzleR[i] = r2-(nozzleZ[i]-z2) * math.tan(conAng)
        elif z3 <= nozzleZ[i] and nozzleZ[i] < z4:
            nozzleR[i] = R_t + throatLeadInRadius - math.sqrt(throatLeadInRadius**2 - (z4 - nozzleZ[i])**2)
        elif z4 <= nozzleZ[i] and nozzleZ[i] < z5:
            nozzleR[i] = R_t + throatLeadOutRadius - math.sqrt(throatLeadOutRadius**2 - (nozzleZ[i] - z4)**2)
        elif z5 <= nozzleZ[i] and nozzleZ[i] <= z6:
            nozzleR[i] = r5+(nozzleZ[i] - z5) * math.tan(exitHalfAngle)
        else:
            nozzleR[i] = 0

    # ==== Size Combustion Chamber ====
    chamberVolume = throatArea * LStar
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

    return (engineContour, chBarrel, nozzleContour, R_tCurve, throatInd)
