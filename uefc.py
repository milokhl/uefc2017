### By Milo Knowles ###
### March 15, 2017 ###

#!/usr/bin/env

import numpy as np
import math
import matplotlib.pyplot as plt

# # DEFINE CONSTANTS #
# PI = 3.14159
# G = 9.81 # m/s^2

# # PLANE VANILLA PARAMETERS #
# print "*** PLANE VANILLA PARAMS ***"
# B = 0.76*2 # m
# RHO_FOAM = 32.0 # kg / m^3
# RHO_AIR = 1.225
# E_FOAM = 19.3 * (10**6) # Pa
# TIP_CHORD = 0.10 # m
# ROOT_CHORD = 0.20 # m
# TAPER_RATIO = float(TIP_CHORD) / ROOT_CHORD
# e = 0.95 # efficiency
# EPSILON = 0.03
# S_REF = 2.0 * (B / 2.0) * (TIP_CHORD+ROOT_CHORD) / 2.0  # m^2
# print "S_REF: ", S_REF

# AR = B**2 / S_REF
# print "AR: ", AR

# T_MAX = 0.7 # N
# CDA0 = 0.004 # m^2

# # PROFILE DRAG COEFFICIENT STUFF 
# c_d_0 = 0.020
# c_d_1 = -0.004
# c_d_2 = 0.020
# c_d_8 = 1.0
# c_l_0 = 0.8

# # C_L is for Plane Vanilla
# C_L = 0.8

# WEIGHT_WING = 0.7348 # N 
# WEIGHT_FUSE = 2.7 # N
# print "W_FUSE: ", WEIGHT_FUSE
# print "W_WING (manual): ", WEIGHT_WING

# TAU = 0.11 # airfoil thickness ratio
# # END PLANE VANILLA PARAMETERS #

def calculateMaxVelocityQuadratic(C_d, rho_air, S_ref, T0=1.0, T1=-0.08, T2=-0.0028):
    """
    Solves a quadratic eqn to find the maximum velocity.
    """
    a = T2 - 0.5*rho_air*C_d*S_ref
    b = T1
    c = T0

    if (b**2 - 4*a*c) < 0: # cannot have negative inside root
        print "Error: negative argument inside of root in calculateMaxVelocityQuadratic"

    ans1 = (-a + math.sqrt(b**2 - 4*a*c)) / (2*a)
    ans2 = (-a - math.sqrt(b**2 - 4*a*c)) / (2*a)

    if ans1<0 and ans2<0:
        print "Error: both solutions to quadratic are negative"

    return max(ans1, ans2)

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
    PI = 3.14159
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


def calculatePayloadWeight(AR, S_ref, CDA0, C_L, c_d, T_MAX, W_fuse, W_wing, e):
    """
    Equation (8) from the notes.
    Wpay = T/((CDA0/S)/Cl + cd_prof/Cl + Cl/(math.pi*AR*e)) - Wfuse - Wwing
    """
    PI = 3.14159
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
                                                   S_ref, deltaBMaxRatio=0.1, N=1):
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


def calculateBankedLoadFactorN(W_fuse, W_wing, W_pay, g, R_turn, S_ref, C_l, rho_air=1.225, handle_ex=True):
    """
    Calculats the load factor N during a banked turn of radius R_turn.
    """
    first = 1 - ((W_fuse + W_pay + W_wing) / (0.5 * rho_air * g * R_turn * S_ref * C_l)) ** 2 

    if first < 0 and handle_ex==True: # only handle this exception is desired
        print "Error: trying to raise negative number to a fractional power. Returning a load factor of 1.01"
        return 1.01

    return (first ** (-0.5))


def calculateBankedVelocityGivenLoadFactor(G, R, load_factor, handle_ex=True):
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


    calcWingWeight = calculateWingWeight(TIP_CHORD, ROOT_CHORD, B, RHO_FOAM, G, TAU)
    print "Calculated Wing Weight: ", calcWingWeight

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



def calculateMinRevTimeOptimization(AR, S_REF, C_L, WEIGHT_FUSE=1.2, MAX_DELTA_B=0.1, verbose=False):
    """
    verbose: should the function output dense print statements?

    The subroutine called by plotResults()
    """
    if verbose:
        print "\n *** OPTIMIZATION *** \n "

    # B and C can be determined by S_REF and AR
    B = (float(S_REF) * AR) ** 0.5
    C = (float(S_REF) / AR) ** 0.5
    
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

    if verbose:
        print "Calculated Wing Weight: ", WING_WEIGHT
        print "c_d: ", c_d 
        print "C_D: ", C_D
        print "B: %f  C: %f" % (B, C)
        print "TIP_CHORD: %f  ROOT_CHORD: %f" % (TIP_CHORD, ROOT_CHORD)
        print "V_max due to thrust constraint: ", V_MAX_THRUST

    # Using wing weight, can calculate the maxmimum load factor (not constrained by thrust)
    N_MAX_BENDING = calculate_Maximum_Load_Factor_Given_DeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, \
                                                                          S_REF, MAX_DELTA_B, 0) # Wpay is zero

    N = calculateBankedLoadFactorN(WEIGHT_FUSE, WING_WEIGHT, 0, G, 12.5, S_REF, C_L, rho_air=1.225)
    N_IS_BENDING_CONSTRAINED = False
    if N >= N_MAX_BENDING:
        N_IS_BENDING_CONSTRAINED = True
        if verbose:
            print "Load factor of %f for AR:%f and SREF:%f exceeds the maximum load factor contrained by bending." % (N, AR, S_REF)
        N = N_MAX_BENDING

    V = calculateBankedVelocityGivenLoadFactor(G, 12.5, N, handle_ex=True)
    V_IS_THRUST_CONTRAINED = False
    if V > V_MAX_THRUST:
        V_IS_THRUST_CONTRAINED = True
        if verbose:
            print "Velocity of %f for AR:%f and SREF:%f exceeds the maximum velocity at full thrust." % (V, AR, S_REF)
        V = V_MAX_THRUST

    rev_time = 2 * PI * 12.5 / V

    if verbose:
        print "Bending Constrained Load Factor: ", N_MAX_BENDING
        print "Load factor: ", N
        print "Velocity: ", V
        print "Rev. Time: ", rev_time

    # this is always printed
    print "S_REF: %f  AR: %f  CL: %f  Time: %f" % (S_REF, AR, C_L, rev_time)

    return rev_time, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, \
    V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, V_IS_THRUST_CONTRAINED, C_D


def calculateMinRevTimeOptimization2(AR, S_REF, C_L, B, C, TIP_CHORD, ROOT_CHORD, WING_WEIGHT, c_d, C_D, T_MAX, \
                                    RHO_AIR, E_FOAM, TAU, EPSILON, TAPER_RATIO, WEIGHT_FUSE=2.7, MAX_DELTA_B=0.1, \
                                    verbose=False):
    """
    verbose: should the function output dense print statements?
    No parameters are derived inside of this function, so it takes in more arguments.
    """
    G = 9.81 # m/s^2
    PI = 3.14159

    # Calculate the maximum velocity given T_max
    # C_d, rho_air, S_ref, T0=1.0, T1=-0.08, T2=-0.0028
    V_MAX_THRUST = calculateMaxVelocityQuadratic(C_D, RHO_AIR, S_REF)

    # Using wing weight, can calculate the maxmimum load factor (not constrained by thrust)
    N_MAX_BENDING = calculate_Maximum_Load_Factor_Given_DeltaBMax(WEIGHT_FUSE, E_FOAM, TAU, EPSILON, TAPER_RATIO, AR, \
                                                                          S_REF, MAX_DELTA_B, 0) # Wpay is zero

    # check if the banked load factor exceeds the bending constraint
    N = calculateBankedLoadFactorN(WEIGHT_FUSE, WING_WEIGHT, 0, G, 12.5, S_REF, C_L, rho_air=1.225)
    N_IS_BENDING_CONSTRAINED = False
    if N >= N_MAX_BENDING:
        N_IS_BENDING_CONSTRAINED = True
        if verbose:
            print "Load factor of %f for AR:%f and SREF:%f exceeds the maximum load factor contrained by bending." % (N, AR, S_REF)
        N = N_MAX_BENDING

    # check if the banked velocity exceeds the maximum velocity give available thrust
    V = calculateBankedVelocityGivenLoadFactor(G, 12.5, N, handle_ex=True)
    V_IS_THRUST_CONTRAINED = False
    if V > V_MAX_THRUST:
        V_IS_THRUST_CONTRAINED = True
        if verbose:
            print "Velocity of %f for AR:%f and SREF:%f exceeds the maximum velocity at full thrust." % (V, AR, S_REF)
        V = V_MAX_THRUST

    # calculate the rev. time using velocity
    if V == 0:
        print "[ERROR] Printing zero division error."
        V = 0.01

    rev_time = 2 * PI * 12.5 / V

    if verbose:
        print "Bending Constrained Load Factor: ", N_MAX_BENDING
        print "Load factor: ", N
        print "Velocity: ", V
        print "Revolution Time: ", rev_time

    # this is always printed
    print "S_REF: %f  AR: %f  CL: %f  Time: %f" % (S_REF, AR, C_L, rev_time)

    return rev_time, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, \
    V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, V_IS_THRUST_CONTRAINED, C_D


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
                 V_IS_THRUST_CONTRAINED = calculateMinRevTimeOptimization(aspect_ratio, area, C_L, MAX_DELTA_B=0.1)

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


def plotTvsDeltaB():
    T_vals = []
    deltab_vals = np.linspace(0.03, 0.3, 20)

    xlist = np.linspace(1.0, 12.0, 24)
    ylist = np.linspace(0.005, 0.2, 50)

    for deltab in deltab_vals:
        t = plotResults(AR=xlist, S=ylist, MAX_DELTA_B=deltab, display=False)
        T_vals.append(t)

    plt.plot(deltab_vals, T_vals)
    plt.title('T_rev vs. delta/b Ratio')
    plt.xlabel('delta/b')
    plt.ylabel('T_rev (sec)')
    plt.show()


def plotResults(AR=None, S=None, MAX_DELTA_B=0.1, filled=True, display=True):
    """
    AR: a list of AR values to try
    S: a list of S values to try
    MAX_DELTA_B: the maximum bending constraint (delta/B ratio)
    filled: Whether or not the contour plot should be filled.
    Creates a search grid using the possible values of AR and S, and optimizes over that grid.
    """

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

    if AR == None:
        xlist = np.linspace(1.0, 12.0, 48) # aspect ratios
    else:
        xlist = AR

    if S == None:
        ylist = np.linspace(0.005, 0.5, 150) # S refs
    else:
        ylist = S

    C_L_List = np.linspace(0.4, 1.2, 15)

    X, Y = np.meshgrid(xlist, ylist)
    Z = np.zeros(shape=(len(ylist), len(xlist)))

    for x in range(len(xlist)):
        for y in range(len(ylist)):
            for c in range(len(C_L_List)):

                try:
                    # calculate the min T_rev with given input params
                    minimumRevTime, N, V, WING_WEIGHT, B, C, TIP_CHORD, ROOT_CHORD, \
                        V_MAX_THRUST, N_MAX_BENDING, N_IS_BENDING_CONSTRAINED, \
                         V_IS_THRUST_CONTRAINED, DRAG_COEFF_AT_OPTIMAL = calculateMinRevTimeOptimization(xlist[x], ylist[y], C_L_List[c], MAX_DELTA_B=MAX_DELTA_B)

                    # if the calculate T_rev improves our best T_rev, store the result and the params
                    if minimumRevTime < optimalRevTime:
                        optimalRevTime = minimumRevTime
                        optimalS_REF = ylist[y]
                        optimalAR = xlist[x]
                        optimalCL = C_L_List[c]

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
                        dragCoeffAtOptimal = DRAG_COEFF_AT_OPTIMAL

                    Z[y][x] = minimumRevTime

                except: # if a value fails to calculate (i.e negative argument inside sqrt), set it to NaN
                    Z[y][x] = float('nan')

    if display:
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
        print "Drag Coeff at Opt:", dragCoeffAtOptimal
        print "Velocity at Max Thrust:", maxVelocityAtThrust
        print "Bending Constrained Load Factor:", bendingConstrainedLoadFactor
        print "Load Factor is Bending Constrained? : ", loadFactorIsBendingConstrained
        print "Velocity is Thrust Constrained? : ", velocityIsThrustConstrained

        plt.figure()

        if filled:
            contour = plt.contourf(X, Y, Z)
        else:
            contour = plt.contour(X, Y, Z)
            

        plt.colorbar(contour, label="Rev. Time (sec)")
        plt.title('Revolution Time vs. AR and S_ref (Delta/B=%.2f)' % MAX_DELTA_B)
        plt.xlabel('Aspect Ratio')
        plt.ylabel('Reference Area (m^2)')
        plt.show()

    return optimalRevTime




if __name__ == '__main__':
    """
    Uncomment the functions below to see their results.
    """

    #planeVanillaAnalysis()

    #optimizeRevTime()
    #testCalculateMinRevTimeOptimization()

    #plotResults(MAX_DELTA_B = 0.1)

    plotTvsDeltaB()