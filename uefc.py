#!/usr/bin/env

import numpy as np
import matplotlib.pyplot as plt

# DEFINE CONSTANTS #
PI = 3.14159
G = 9.81 # m/s^2

# PLANE VANILLA PARAMETERS #
print "*** PLANE VANILLA PARAMS ***"
B = 0.76*2 # m
RHO_FOAM = 32.0 # kg / m^3
RHO_AIR = 1.225
E_FOAM = 19.3 * (10**6) # Pa
TIP_CHORD = 0.10 # m
ROOT_CHORD = 0.20 # m
TAPER_RATIO = float(TIP_CHORD) / ROOT_CHORD
e = 0.95 # efficiency
EPSILON = 0.03
S_REF = 2.0 * (B / 2.0) * (TIP_CHORD+ROOT_CHORD) / 2.0  # m^2
print "S_REF: ", S_REF

AR = B**2 / S_REF
print "AR: ", AR

T_MAX = 0.7 # N
CDA0 = 0.004 # m^2

# PROFILE DRAG COEFFICIENT STUFF 
c_d_0 = 0.020
c_d_1 = -0.004
c_d_2 = 0.020
c_d_8 = 1.0
c_l_0 = 0.8

# C_L is for Plane Vanilla
C_L = 0.8

WEIGHT_WING = 0.7348 # N 
WEIGHT_FUSE = 2.7 # N
print "W_FUSE: ", WEIGHT_FUSE
print "W_WING (manual): ", WEIGHT_WING

TAU = 0.11 # airfoil thickness ratio

# END PLANE VANILLA PARAMETERS #


def planeVanillaChord(y):
    """
    Chord length as a function of y.
    c(y) = (root_chord - tip_chord) * (1 - 2y/span) + tip_chord
    """
    return (ROOT_CHORD, TIP_CHORD) * (1- (2*y / B)) + TIP_CHORD


def calculateProfileDragCoefficient(c_l, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8):
    c_d = c_d_0 + c_d_1 * (c_l - c_l_0) + c_d_2 * ((c_l-c_l_0)**2) + c_d_8 * ((c_l-c_l_0)**8)
    return c_d

def calculateCoeffDrag(CDA0, S, profDragCoeff, C_L, AR, e):
    t1 = float(CDA0) / S
    t2 = profDragCoeff
    t3 = float(C_L**2) / (PI * AR * e)
    return (t1 + t2 + t3)
    
def calculateCoeffLift(W_fuse, W_wing, W_pay, rho_air, velocity):
    """
    Lift = Weight used to derive
    """
    C_L = (W_fuse + W_wing + W_pay) / (0.5 * rho_air * velocity**2)
    return C_L

def calculateWingWeight(tip_chord, root_chord, span, rho, g, tau):
    delta_c = root_chord - tip_chord

    coeff = 0.2 * span * tau * rho / delta_c
    diff = root_chord**3 - tip_chord**3

    return coeff * diff * g


calcWingWeight = calculateWingWeight(TIP_CHORD, ROOT_CHORD, B, RHO_FOAM, G, TAU)
print "Calculated Wing Weight: ", calcWingWeight


def calculatePayloadWeight(AR, S_ref, CDA0, C_L, c_d, T_MAX, W_fuse, W_wing, e):
    """
    Equation (8) from the notes.
    """
    denom1 = float(CDA0) / (S_ref * C_L)
    denom2 = float(c_d) / C_L
    denom3 = float(C_L) / (PI * AR * e)

    term1 = (float(T_MAX) / (denom1 + denom2 + denom3))
    print "Term 1: ", term1
    W_pay = term1 - W_fuse - W_wing
    return W_pay


def calculateDeltaBRatio(W_fuse, W_pay, E, tau, epsilon, lambda_, aspect_ratio, S_ref, N=1):
    """
    Calculates the sigma / b ratio, using eqn 16 in notes.
    """
    part1 = float(W_fuse + W_pay) / (E * tau * (tau**2 + epsilon**2))
    part2 = (1 + lambda_)**3 * (1 + 2*lambda_)
    part3 = aspect_ratio**3 / S_ref
    return 0.018 * part1 * part2 * part3


def calculatePayloadWeightUpperBoundGivenDeltaBMax(W_fuse, E, tau, epsilon, lambda_, aspect_ratio, \
                                                   S_ref, deltaBMaxRatio, N=1):
    """
    Given a delta/b maximum ratio, determine the maximum value that the W_pay can have.
    """
    num = deltaBMaxRatio * E * tau * (tau**2 + epsilon**2) * S_ref
    den = 0.018 * N * (1 + lambda_)**3 * (1 + 2*lambda_) * aspect_ratio**3
    res = float(num) / den - W_fuse
    return res


def calculate_Maximum_Load_Factor_Given_DeltaBMax(W_fuse, E, tau, epsilon, lambda_, aspect_ratio, \
                                                  S_ref, deltaBMaxRatio, W_pay):

    """
    Given the maximum delta/b ratio constraint, calculate the maxmimum load factor that can be achieved.
    """
    num = deltaBMaxRatio * E * tau * (tau**2 + epsilon**2) * S_ref
    den = 0.018 * (W_fuse + W_pay) * (1 + lambda_)**3 * (1 + 2*lambda_) * aspect_ratio**3
    return num / den


def calculateMaxVelocity(T_max, C_d, rho_air, S_ref):
    """
    Tmax = D
    Solves for the velocity at maximum thrust.
    """
    v_max = ((2 * T_max) / (C_d * rho_air * S_ref)) ** 0.5
    return v_max


def calculateLoadFactorGivenVelocity(velocity, radius, g):
    """
    Rearranged eqn (20) to solve for N
    """
    term1 = (velocity**4) / (g**2 * radius**2) + 1
    return term1**0.5


def calculateBankedLoadFactorN(W_fuse, W_wing, W_pay, g, R_turn, S_ref, C_l, rho_air=1.225):
    """
    Calculats the load factor N during a banked turn of radius R_turn.
    """
    first = 1 - ((W_fuse + W_pay + W_wing) / (0.5 * rho_air * g * R_turn * S_ref * C_l)) ** 2 
    if first < 0:
        print "Error: trying to raise negative number to a fractional power. Returning a load factor of 1.01"
        return 1.01
    return (first ** (-0.5))


def calculateBankedVelocityGivenLoadFactor(G, R, load_factor):
    """
    use eqn (20)
    """
    if load_factor < 1:
        print "Error: load factor of %f is less than 1. Using 1.01." % load_factor
        load_factor = 1.01
    v = (G * R * (load_factor**2 - 1)**0.5) ** 0.5
    return v


def calculateRevolutionTime(W_fuse, rho_foam, g, S_ref, AR, CDA0, c_d, T_MAX, \
                            tip_chord, root_chord, span, tau, C_L=0.8, R_turn=12.5, \
                            N=None, W_wing=None, rho_air=1.225):
    """
    Assumptions made for bending constrained trev:
    CL = 0.8
    Wpay = 0
    (delta/b)max = 0.1
    """

    # calculate the wing weight so that we can get 
    if W_wing == None:
        print "[REV TIME] No wing weight specified, calculating..."
        W_wing = calculateWingWeight(tip_chord, root_chord, span, rho_foam, g, tau)


    # now get Wpay
    #Wpay = calculatePayloadWeight(AR, S_ref, CDA0, C_L, c_d, T_MAX, W_fuse, W_wing, e)
    W_pay = 0

    # calculate N using eqn (22) if it hasn't been specified
    if N==None:
        print "[REV TIME] No load factor specified, calculating..."
        N = calculateBankedLoadFactorN(W_fuse, W_wing, W_pay, g, R_turn, S_ref, C_L)

    # now calculate T_rev^2 using eqn (19) and (20)
    t_rev_squared = float(4 * PI**2 * R_turn) / (g * (N**2 - 1)**0.5)

    # take the square root to get the final answer for T_rev
    return t_rev_squared**0.5



def planeVanillaAnalysis():

    print "*** WPAY ANALYSIS *** "
    # Note: "make the assumption that c_l = C_L"
    c_d = calculateProfileDragCoefficient(C_L, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8)
    print "c_d: ", c_d 

    C_D = calculateCoeffDrag(CDA0, S_REF, c_d, C_L, AR, e)
    print "C_D: ", C_D

    W_pay_calc = calculatePayloadWeight(AR, S_REF, CDA0, C_L, c_d, T_MAX, WEIGHT_FUSE, WEIGHT_WING, e)
    print "Payload weight: ", W_pay_calc, "N"
    sigma_b_ratio = calculateDeltaBRatio(WEIGHT_FUSE, W_pay_calc, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, S_REF)
    print "Sigma B Ratio: ", sigma_b_ratio

    bending_constrained_wpay = calculatePayloadWeightUpperBoundGivenDeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, 0.228, 0.1)
    print "Bending constrained Wpay: ", bending_constrained_wpay


    print "\n *** TREV ANALYSIS ***"

    bending_constrained_N = calculate_Maximum_Load_Factor_Given_DeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, \
                                                                          S_REF, 0.1, 0)

    print "Bending constrained Load Factor: ", bending_constrained_N


    # calculate revolution time
    min_Nconstrained_rev_time = calculateRevolutionTime(WEIGHT_FUSE, RHO_FOAM, G, S_REF, AR, CDA0, c_d, T_MAX, \
                                0.1, 0.2, B, TAU, C_L=0.8, R_turn=12.5, N=bending_constrained_N, \
                                W_wing=WEIGHT_WING)

    print "Min. Rev. Time Bending-Constrained: ", min_Nconstrained_rev_time

    N = calculateBankedLoadFactorN(WEIGHT_FUSE, WEIGHT_WING, 0, G, 12.5, S_REF, C_L, rho_air=1.225)
    print "Load factor: ", N

    V = calculateBankedVelocityGivenLoadFactor(G, 12.5, N)

    print "Velocity: ", V

    rev_time = 2 * PI * 12.5 / V
    print "Rev. Time: ", rev_time



def calculateMinRevTimeOptimization(AR, S_REF, C_L, MAX_DELTA_B=0.1):
    
    print "\n *** OPTIMIZATION *** \n "
    print "S_REF: %f  AR: %f  CL: %f" % (S_REF, AR, C_L)
    print "*** Derived Parameters: *** "
    # B and C can be determined by S_REF and AR
    B = (float(S_REF) * AR) ** 0.5
    C = (float(S_REF) / AR) ** 0.5
    print "B: %f  C: %f" % (B, C)

    TIP_CHORD = float(2 * C) / 3
    ROOT_CHORD = float(4 * C) / 3  

    # Now can calculate W_wing
    WING_WEIGHT = calculateWingWeight(TIP_CHORD, ROOT_CHORD, B, RHO_FOAM, G, TAU)
    
    # calculate drag coefficient
    #C_L = 0.8
    c_d = calculateProfileDragCoefficient(C_L, c_l_0, c_d_0, c_d_1, c_d_2, c_d_8)
    C_D = calculateCoeffDrag(CDA0, S_REF, c_d, C_L, AR, e)
    
    # Calculate the maximum velocity given T_max
    V_MAX_THRUST = calculateMaxVelocity(T_MAX, C_D, RHO_AIR, S_REF)

    print "Calculated Wing Weight: ", WING_WEIGHT
    print "c_d: ", c_d 
    print "C_D: ", C_D
    print "TIP_CHORD: %f  ROOT_CHORD: %f" % (TIP_CHORD, ROOT_CHORD)
    print "V_max due to thrust constraint: ", V_MAX_THRUST

    # calculate lift coefficient using maximum velocity constrained by thrust
    #C_L = calculateCoeffLift(WEIGHT_FUSE, WING_WEIGHT, 0, RHO_AIR, V_MAX_THRUST)
    #C_L = 0.8

    # Using wing weight, can calculate the maxmimum load factor (not constrained by thrust)
    N_MAX_BENDING = calculate_Maximum_Load_Factor_Given_DeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, \
                                                                          S_REF, MAX_DELTA_B, 0) # Wpay is zero

    #V_MAX_BENDING = calculateBankedVelocityGivenLoadFactor(G, 12.5, bending_constrained_N)

    print "Bending Constrained Load Factor: ", N_MAX_BENDING
    #print "Bending Constrained Velocity: ", V_MAX_BENDING

    # N can be at most N_MAX_BENDING
    # V can be at most V_MAX_THRUST
    N = calculateBankedLoadFactorN(WEIGHT_FUSE, WING_WEIGHT, 0, G, 12.5, S_REF, C_L, rho_air=1.225)

    N_IS_BENDING_CONSTRAINED = False
    if N >= N_MAX_BENDING:
        N_IS_BENDING_CONSTRAINED = True
        print "Load factor of %f for AR:%f and SREF:%f exceeds the maximum load factor contrained by bending." % (N, AR, S_REF)
        N = N_MAX_BENDING

    V = calculateBankedVelocityGivenLoadFactor(G, 12.5, N)
    V_IS_THRUST_CONTRAINED = False
    if V > V_MAX_THRUST:
        V_IS_THRUST_CONTRAINED = True
        print "Velocity of %f for AR:%f and SREF:%f exceeds the maximum velocity at full thrust." % (V, AR, S_REF)
        V = V_MAX_THRUST

    rev_time = 2 * PI * 12.5 / V

    print "Load factor: ", N
    print "Velocity: ", V
    print "Rev. Time: ", rev_time

    return rev_time, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, V_IS_THRUST_CONTRAINED


def testCalculateMinRevTimeOptimization():
    AR = 10.33
    S = 0.228
    C_L = 0.8
    rev_time, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, V_MAX_THRUST, \
    N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, V_IS_THRUST_CONTRAINED = calculateMinRevTimeOptimization(AR, S, C_L)

    print "[TEST] Rev. Time: ", rev_time
    print "[TEST] Velocity: ", V
    print "[TEST] TIP_CHORD: ", TIP_CHORD
    print "[TEST] ROOT_CHORD: ", ROOT_CHORD


def optimizeRevTime():

    aspectRatios = [0.5 * i for i in range(2, 20)]
    areas = [0.005 * i for i in range(1,80)]
    C_Ls = [0.3 + 0.05 * i for i in range(1, 50)]

    optimalRevTime = 1000
    optimalS_REF = None
    optimalAR = None
    optimalCL = None
    loadFactorAtOptimal = None
    velocityAtOptimal = None
    wingWeightAtOptimal = None
    spanAtOptimal = None
    averageChordAtOptimal = None
    rootChordAtOptimal = None
    tipChordAtOptimal = None
    maxVelocityAtThrust = None
    bendingConstrainedLoadFactor = None
    loadFactorIsBendingConstrained = None
    velocityIsThrustConstrained = None

    for aspect_ratio in aspectRatios:
        for area in areas:
            for C_L in [0.65]:
            #for C_L in C_Ls:

                minimumRevTime, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, \
                V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, \
                 V_IS_THRUST_CONTRAINED= calculateMinRevTimeOptimization(aspect_ratio, area, C_L, MAX_DELTA_B=0.1)

                if minimumRevTime < optimalRevTime:
                    optimalRevTime = minimumRevTime
                    optimalS_REF = area
                    optimalAR = aspect_ratio
                    optimalCL = C_L

                    loadFactorAtOptimal = N
                    velocityAtOptimal = V
                    wingWeightAtOptimal = WING_WEIGHT
                    spanAtOptimal = B
                    averageChordAtOptimal = C
                    rootChordAtOptimal = ROOT_CHORD
                    tipChordAtOptimal = TIP_CHORD
                    maxVelocityAtThrust = V_MAX_THRUST
                    bendingConstrainedLoadFactor = N_MAX_BENDING
                    loadFactorIsBendingConstrained = N_IS_BENDING_CONSTRAINED
                    velocityIsThrustConstrained = V_IS_THRUST_CONTRAINED

    print "\n *** FINAL RESULTS ***"
    print "Optimal Rev. Time: ", optimalRevTime, "sec"
    print "Optimal S_ref: ", optimalS_REF, "m^2"
    print "Optimal AR: ", optimalAR,
    print "Optimal CL: ", optimalCL
    print ""
    print "Derived Params at Optimal:"
    print "Load Factor at Optimal: ", loadFactorAtOptimal
    print "Velocity at Optimal: ", velocityAtOptimal
    print "Wing Weight at Optimal: ", wingWeightAtOptimal
    print "Span at Optimal: ", spanAtOptimal
    print "Avg. Chord at Optimal: ", averageChordAtOptimal
    print "Root Chord at Opt:", rootChordAtOptimal
    print "Tip Chord at Opt:", tipChordAtOptimal
    print "Velocity at Max Thrust:", maxVelocityAtThrust
    print "Bending Constrained Load Factor:", bendingConstrainedLoadFactor
    print "Load Factor is Bending Constrained? : ", loadFactorIsBendingConstrained
    print "Velocity is Thrust Constrained? : ", velocityIsThrustConstrained


def plotResults(AR, S):

    optimalRevTime = 1000
    optimalS_REF = None
    optimalAR = None
    optimalCL = None
    loadFactorAtOptimal = None
    velocityAtOptimal = None
    wingWeightAtOptimal = None
    spanAtOptimal = None
    averageChordAtOptimal = None
    rootChordAtOptimal = None
    tipChordAtOptimal = None
    maxVelocityAtThrust = None
    bendingConstrainedLoadFactor = None
    loadFactorIsBendingConstrained = None
    velocityIsThrustConstrained = None

    xlist = np.linspace(1.0, 15.0, 20) # aspect ratios
    ylist = np.linspace(0.005, 0.3, 80) # S refs
    X, Y = np.meshgrid(xlist, ylist)
    Z = np.zeros(shape=(len(ylist), len(xlist)))

    #print type(xlist)
    for x in range(len(xlist)):
        for y in range(len(ylist)):
            #print "indices:", x, y

            minimumRevTime, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, \
                V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, \
                 V_IS_THRUST_CONTRAINED = calculateMinRevTimeOptimization(xlist[x], ylist[y], 0.65, MAX_DELTA_B=0.1)

            if minimumRevTime < optimalRevTime:
                optimalRevTime = minimumRevTime
                optimalS_REF = ylist[y]
                optimalAR = xlist[x]
                optimalCL = 0.65

                loadFactorAtOptimal = N
                velocityAtOptimal = V
                wingWeightAtOptimal = WING_WEIGHT
                spanAtOptimal = B
                averageChordAtOptimal = C
                rootChordAtOptimal = ROOT_CHORD
                tipChordAtOptimal = TIP_CHORD
                maxVelocityAtThrust = V_MAX_THRUST
                bendingConstrainedLoadFactor = N_MAX_BENDING
                loadFactorIsBendingConstrained = N_IS_BENDING_CONSTRAINED
                velocityIsThrustConstrained = V_IS_THRUST_CONTRAINED

            Z[y][x] = minimumRevTime


    print "\n *** FINAL RESULTS ***"
    print "Optimal Rev. Time: ", optimalRevTime, "sec"
    print "Optimal S_ref: ", optimalS_REF, "m^2"
    print "Optimal AR: ", optimalAR,
    print "Optimal CL: ", optimalCL
    print ""
    print "Derived Params at Optimal:"
    print "Load Factor at Optimal: ", loadFactorAtOptimal
    print "Velocity at Optimal: ", velocityAtOptimal
    print "Wing Weight at Optimal: ", wingWeightAtOptimal
    print "Span at Optimal: ", spanAtOptimal
    print "Avg. Chord at Optimal: ", averageChordAtOptimal
    print "Root Chord at Opt:", rootChordAtOptimal
    print "Tip Chord at Opt:", tipChordAtOptimal
    print "Velocity at Max Thrust:", maxVelocityAtThrust
    print "Bending Constrained Load Factor:", bendingConstrainedLoadFactor
    print "Load Factor is Bending Constrained? : ", loadFactorIsBendingConstrained
    print "Velocity is Thrust Constrained? : ", velocityIsThrustConstrained

    plt.figure()
    contour = plt.contourf(X, Y, Z)
    #plt.clabel(contour, colors = 'k', fmt = '%2.1f', fontsize=12)
    #c = ('#ff0000', '#ffff00', '#0000FF', '0.6', 'c', 'm')
    #contour_filled = plt.contourf(X, Y, Z, colors=c)
    contour_filled = plt.contour(X, Y, Z)
    plt.colorbar(contour)
    plt.title('Filled Contours Plot')
    plt.xlabel('Aspect Ratio')
    plt.ylabel('S_REF')
    plt.show()


if __name__ == '__main__':
    #planeVanillaAnalysis()

    #optimizeRevTime()
    #testCalculateMinRevTimeOptimization()

    plotResults(None, None)